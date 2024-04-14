import argparse
import subprocess
from cyvcf2 import VCF, Writer
import os
import sys
import hashlib
from pyfaidx import Fasta

dump = open('/dev/null', 'w')


def read_cluster_vcf(vcf_path, max_af=.01, flank_size=50000):
    # record SV info
    svs = {}
    vcf = VCF(vcf_path)
    for variant in vcf:
        if len(variant.REF) > 30 or len(variant.ALT[0]) > 30:
            svinfo = {}
            svinfo['ref'] = variant.REF
            svinfo['alt'] = variant.ALT[0]
            svinfo['seqn'] = variant.CHROM
            svinfo['start'] = variant.start + 1
            svinfo['end'] = variant.end + 1
            seq = '{}_{}'.format(variant.REF, variant.ALT[0])
            seq = hashlib.sha1(seq.encode())
            svinfo['svid'] = '{}_{}_{}'.format(variant.CHROM,
                                               variant.start,
                                               seq.hexdigest())
            svtype = variant.INFO.get('SVTYPE')
            if svtype is not None and svtype == 'BND':
                continue
            af = variant.INFO.get('AF')
            if af is None or af <= max_af:
                # add variants and ref sequence to input files
                # for the pangenome
                svs[svinfo['svid']] = svinfo
    vcf.close()
    # sort SVs per chromosome
    sv_pos_chr = {}
    for svid in svs:
        svinfo = svs[svid]
        if svinfo['seqn'] not in sv_pos_chr:
            sv_pos_chr[svinfo['seqn']] = []
        sv_pos_chr[svinfo['seqn']].append([svid,
                                           svinfo['start'],
                                           svinfo['end']])
    sv_cls = []
    cur_cl = []
    for seqn in sv_pos_chr:
        pos_sorted = sorted(sv_pos_chr[seqn], key=lambda k: (k[1], k[2]))
        # create SV clusters
        for svpos in pos_sorted:
            if len(cur_cl) == 0:
                # current cluster is empty, so just init with the next SV
                cur_cl.append(svpos[0])
            else:
                # compare positions with last SV in cluster
                cl_sv = cur_cl[-1]
                if svs[cl_sv]['end'] + flank_size > svpos[1] - flank_size:
                    # their regions would overlap, add to the current cluster
                    cur_cl.append(svpos[0])
                else:
                    # they are far enough from each other, save current cluster
                    # and start a new one
                    sv_cls.append(cur_cl)
                    cur_cl = [svpos[0]]
        if len(cur_cl) > 0:
            sv_cls.append(cur_cl)
            cur_cl = []
    vcf.close()
    return ({'svs': svs, 'cls': sv_cls})


def add_sv_to_pan(svs, ref, svcl,
                  vcf_outf, ref_outf, reg_outf, chr_lens, flank_size=50000):
    # broad region of interest (for the reference sequence)
    cl_start = min([svs[svid]['start'] for svid in svcl])
    region_start = max(1, cl_start - flank_size)
    cl_end = max([svs[svid]['end'] for svid in svcl])
    region_end = cl_end + flank_size
    seqn = svs[svcl[0]]['seqn']
    if seqn in chr_lens:
        region_end = min(region_end, chr_lens[seqn] - 1)
    # extract reference sequence for the FASTA
    ref_cont = ref[seqn][region_start:region_end]
    cont_name = 'ref_' + '_'.join([svs[svid]['svid'] for svid in svcl])
    # append to fasta file
    ref_outf.write('>{}\n{}\n'.format(cont_name, str(ref_cont)))
    # add each SV to VCF
    rec_to_format = "{seqn}\t{pos}\t{svid}\t{ref}\t{alt}\t.\t.\t.\n"
    for svid in svcl:
        cont_pos = svs[svid]['start'] - region_start
        vcf_outf.write(rec_to_format.format(seqn=cont_name,
                                            pos=cont_pos,
                                            svid=svid,
                                            ref=svs[svid]['ref'],
                                            alt=svs[svid]['alt']))
    # regions for read extraction use smaller flanks to make sure they align
    flank_size_n = int(flank_size / 5)
    eregion_start = max(1, cl_start - flank_size_n)
    eregion_end = cl_end + flank_size_n
    if seqn in chr_lens:
        eregion_end = min(eregion_end, chr_lens[seqn] - 1)
    # write bed
    reg_outf.write('{}\t{}\t{}\n'.format(seqn,
                                         eregion_start,
                                         eregion_end))


def genotype_svs(gt_vcf, ref_fa, reads_fq):
    # #
    print('      Build pangenome...')
    # pre-process VCF for pangenome construction
    bgzip_args = ["bgzip", '-f', gt_vcf]
    gt_vcf_gz = gt_vcf + '.gz'
    subprocess.run(bgzip_args, check=True)
    tabix_args = ["tabix", gt_vcf_gz]
    subprocess.run(tabix_args, check=True)
    # make graph with SV
    construct_args = ["vg", "construct", "-a", "-m", "1024", "-S",
                      "-t", str(args.t), "-r", ref_fa, "-v", gt_vcf_gz]
    vg_output_path = os.path.join(args.d, "graph.vg")
    with open(vg_output_path, 'w') as file:
        subprocess.run(construct_args, check=True, stdout=file, stderr=dump)
    # #
    print('      Align extracted reads to pangenome...')
    # align reads to pangenome
    convert_args = ["vg", "convert", "-f", vg_output_path]
    gfa_output_path = os.path.join(args.d, "graph.gfa")
    with open(gfa_output_path, 'w') as file:
        subprocess.run(convert_args, check=True, stdout=file, stderr=dump)
    map_args = ["minigraph", "-t", str(args.t),
                "-c", gfa_output_path, reads_fq]
    gaf_output_path = os.path.join(args.d, "reads.gaf")
    with open(gaf_output_path, 'w') as file:
        subprocess.run(map_args, check=True, stdout=file, stderr=dump)
    # #
    print('      Genotype SV calls...')
    pack_output_path = os.path.join(args.d, "reads.pack")
    pack_args = ["vg", "pack", "-t", str(args.t), "-e",
                 "-x", vg_output_path, "-o", pack_output_path,
                 '-a', gaf_output_path]
    try:
        subprocess.run(pack_args, check=True, stdout=sys.stdout)
    except subprocess.CalledProcessError as e:
        print('cmd: {}\nerror code: {}.\noutput: {}\nstderr:'
              ' {}'.format(e.cmd,
                           e.returncode,
                           e.output,
                           e.stderr))
    call_args = ["vg", "call", "-t", str(args.t),
                 "-k", pack_output_path, '-v', gt_vcf_gz, vg_output_path]
    call_output_path = os.path.join(args.d, "called.vcf")
    with open(call_output_path, 'w') as file:
        subprocess.run(call_args, check=True, stdout=file, stderr=dump)
    # delete temporary files
    for ftrm in [gt_vcf_gz, gt_vcf_gz + '.tbi', reads_fq,
                 ref_fa, ref_fa + '.fai',
                 vg_output_path, gfa_output_path, gaf_output_path,
                 pack_output_path]:
        os.remove(ftrm)
    # return path to calls
    return (call_output_path)


parser = argparse.ArgumentParser()
parser.add_argument('-b', help='BAM file (indexed)', required=True)
parser.add_argument('-f', help='reference FASTA file (indexed)', required=True)
parser.add_argument('-v', help='variants in VCF (can be bgzipped)',
                    required=True)
parser.add_argument('-d', help='output directory', default='temp_valsv')
parser.add_argument('-F', default=50000, type=int,
                    help='size of flanking regions added. Default 50000')
parser.add_argument('-m', default=0.01, type=float,
                    help='maximum frequency (AF) to annotate. Default 0.01')
parser.add_argument('-B', default=200, type=int,
                    help='batch size. Default 200')
parser.add_argument('-o', default='out.vcf',
                    help='output (annotated) VCF (will be '
                    'bgzipped if ending in .gz)')
parser.add_argument('-t', default=2, type=int,
                    help='number of threads used by the tools '
                    '(vg and minigraph))')
args = parser.parse_args()

# python3 validate-svs.py -b ../../test/reads.bam -f ../../test/hg38.fa -v ../../test/test.vcf.gz
# args = parser.parse_args('-b ../../test/reads.bam -f ../../test/hg38.fa -v ../../test/test.vcf.gz'.split())

# create temp directory if it doesn't exist
if not os.path.exists(args.d):
    os.makedirs(args.d)

# load reference genome index
ref = Fasta(args.f)
# extract chromosome lengths from the .fai (2nd column)
chr_lens = {}
with open(args.f + '.fai', 'rt') as inf:
    for line in inf:
        line = line.rstrip().split('\t')
        chr_lens[line[0]] = int(line[1])


# Cluster SV regions to avoid multimap issues if they overlap
print('Read input VCF and cluster SVs...')
in_svs = read_cluster_vcf(args.v, max_af=args.m, flank_size=args.F)

# if there are a lot of regions to process, create batches
batches = []
cur_batch = []
for cl in in_svs['cls']:
    if len(cur_batch) < args.B:
        cur_batch.append(cl)
    else:
        batches.append(cur_batch)
        cur_batch = []
if len(cur_batch) > 0:
    batches.append(cur_batch)

print('{} SV(s) in {} cluster(s) will be processed in {} '
      'batch(es).'.format(len(in_svs['svs']),
                          len(in_svs['cls']),
                          len(batches)))
val_score = {}  # to record validation scores for all SVs
for bii, batch in enumerate(batches):
    print('   Batch {}/{}'.format(bii + 1, len(batches)))
    print('      Prepare files for pangenome construction...')
    # prepare 3 files: VCF + FASTA for the pangenome, + bed with regions
    gt_vcf = os.path.join(args.d, "for_construct.vcf")
    reg_bed = os.path.join(args.d, "regions.bed")
    ref_fa = os.path.join(args.d, "ref.fa")
    # open connections to files
    ref_outf = open(ref_fa, 'wt')
    vcf_outf = open(gt_vcf, 'wt')
    reg_outf = open(reg_bed, 'wt')
    # init VCF file
    vcf_outf.write("##fileformat=VCFv4.2\n")
    vcf_outf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    for svcl in batch:
        # add variants and ref sequence to input files for the pangenome
        add_sv_to_pan(in_svs['svs'], ref, svcl,
                      vcf_outf, ref_outf, reg_outf,
                      chr_lens=chr_lens,
                      flank_size=args.F)
    # close connections
    vcf_outf.close()
    ref_outf.close()
    reg_outf.close()
    print('      Extract reads...')
    # extract reads
    reads_bam = os.path.join(args.d, "reads.bam")
    extract_cmd = ["samtools", "view", "-h", '-@', str(args.t - 1),
                   '--region-file', reg_bed,
                   '-o', reads_bam, args.b]
    subprocess.run(extract_cmd, check=True)
    os.remove(reg_bed)
    reads_fq = os.path.join(args.d, "reads.fq.gz")
    extract_cmd = ['samtools', "fastq", '-@', str(args.t - 1),
                   '-0', reads_fq, reads_bam]
    subprocess.run(extract_cmd, check=True, stderr=dump)
    os.remove(reads_bam)
    # align and genotype
    calls_vcf = genotype_svs(gt_vcf, ref_fa, reads_fq)
    # #
    print('      Collect genotyping information...')
    called_vcf = VCF(calls_vcf)
    for variant in called_vcf:
        score = {}
        ad = variant.format('AD')[0]
        dp = ad[0] + ad[1]
        score['ad'] = '{}|{}'.format(ad[0], ad[1])
        if dp == 0:
            score['prop'] = 0
        else:
            score['prop'] = float(ad[1]) / dp
        val_score[variant.ID] = score
    called_vcf.close()
    # delete temporary files
    os.remove(calls_vcf)

# #
print('Annotate input VCF...')
# annotate input VCF
vcf = VCF(args.v)
vcf.add_info_to_header({'ID': 'RS_PROP',
                        'Description': 'Proportion of supporting reads '
                        '(from vg genotyping)',
                        'Type': 'Float', 'Number': '1'})
vcf.add_info_to_header({'ID': 'RS_AD',
                        'Description': 'Number of reads supporting the '
                        'ref|alt (from vg genotyping)',
                        'Type': 'String', 'Number': '1'})
vcf_o = Writer(args.o, vcf)
for variant in vcf:
    if len(variant.REF) > 30 or len(variant.ALT[0]) > 30:
        seq = '{}_{}'.format(variant.REF, variant.ALT[0])
        seq = hashlib.sha1(seq.encode())
        svid = '{}_{}_{}'.format(variant.CHROM,
                                 variant.start,
                                 seq.hexdigest())
        if svid in val_score:
            variant.INFO["RS_PROP"] = val_score[svid]['prop']
            variant.INFO["RS_AD"] = val_score[svid]['ad']
        else:
            variant.INFO["RS_PROP"] = 2
            variant.INFO["RS_AD"] = 'NA'
    vcf_o.write_record(variant)
vcf_o.close()
vcf.close()

dump.close()
