import argparse
from subprocess import run
from cyvcf2 import VCF, Writer
import os
import sys
import hashlib
from pyfaidx import Fasta

dump = open('/dev/null', 'w')
flank_size = 50000


def add_sv_to_pan(svs, ref, bam_path, svcl,
                  vcf_outf, ref_outf, reads_outf, chr_lens):
    # region of interest
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
    # extract reads using smaller flanks to make sure they align
    flank_size_n = int(flank_size / 5)
    eregion_start = max(1, cl_start - flank_size_n)
    eregion_end = cl_end + flank_size_n
    if seqn in chr_lens:
        eregion_end = min(eregion_end, chr_lens[seqn] - 1)
    eregion_coord = '{}:{}-{}'.format(seqn,
                                      eregion_start,
                                      eregion_end)
    # extract reads
    extract_args = ["samtools", "view", "-h", bam_path, eregion_coord]
    sam_view = run(extract_args, check=True, capture_output=True)
    extract_args = ["samtools", "fasta"]
    run(extract_args, check=True, input=sam_view.stdout,
        stdout=reads_outf, stderr=dump)
    # add to VCF
    rec_to_format = "{seqn}\t{pos}\t{svid}\t{ref}\t{alt}\t.\t.\t.\n"
    for svid in svcl:
        cont_pos = svs[svid]['start'] - region_start
        vcf_outf.write(rec_to_format.format(seqn=cont_name,
                                            pos=cont_pos,
                                            svid=svid,
                                            ref=svs[svid]['ref'],
                                            alt=svs[svid]['alt']))
    return (False)


parser = argparse.ArgumentParser()
parser.add_argument('-b', help='BAM file (indexed)', required=True)
parser.add_argument('-f', help='reference FASTA file (indexed)', required=True)
parser.add_argument('-v', help='variants in VCF (can be bgzipped)',
                    required=True)
parser.add_argument('-d', help='output directory', default='temp_valsv')
parser.add_argument('-F', default=0.01, type=float,
                    help='maximum frequency (AF) to annotate. Default 0.01')
parser.add_argument('-o', default='out.vcf',
                    help='output (annotated) VCF (will be '
                    'bgzipped if ending in .gz)')
parser.add_argument('-t', default=2,
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
# record SV info
svs = {}
vcf = VCF(args.v)
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
        if 'AF' not in variant.INFO or variant.INFO['AF'] <= args.F:
            # add variants and ref sequence to input files for the pangenome
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

print('Prepare files for pangenome construction...')
# prepare 3 files: VCF + FASTA for the pangenome, + local reads
gt_vcf = os.path.join(args.d, "for_construct.vcf")
reads_fa = os.path.join(args.d, "reads.fasta")
ref_fa = os.path.join(args.d, "ref.fa")
# open connections to files
ref_outf = open(ref_fa, 'wt')
vcf_outf = open(gt_vcf, 'wt')
reads_outf = open(reads_fa, 'wt')
# init VCF file
vcf_outf.write("##fileformat=VCFv4.2\n")
vcf_outf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
for svcl in sv_cls:
    # add variants and ref sequence to input files for the pangenome
    add_sv_to_pan(svs, ref, args.b, svcl,
                  vcf_outf, ref_outf, reads_outf, chr_lens)
vcf.close()
# close connections
vcf_outf.close()
reads_outf.close()
ref_outf.close()


# #######
print('Build pangenome with {} SV calls in {} clusters to '
      'genotype...'.format(len(svs), len(sv_cls)))
# #######

# pre-process VCF for pangenome construction
bgzip_args = ["bgzip", '-f', gt_vcf]
gt_vcf_gz = gt_vcf + '.gz'
run(bgzip_args, check=True)
tabix_args = ["tabix", gt_vcf_gz]
run(tabix_args, check=True)
# make graph with SV
construct_args = ["vg", "construct", "-a", "-m", "1024", "-S",
                  "-t", str(args.t), "-r", ref_fa, "-v", gt_vcf_gz]
vg_output_path = os.path.join(args.d, "graph.vg")
with open(vg_output_path, 'w') as file:
    run(construct_args, check=True, stdout=file, stderr=dump)


# #######
print('Align extracted reads to pangenome...')
# #######

# align reads to pangenome
convert_args = ["vg", "convert", "-f", vg_output_path]
gfa_output_path = os.path.join(args.d, "graph.gfa")
with open(gfa_output_path, 'w') as file:
    run(convert_args, check=True, stdout=file, stderr=dump)
map_args = ["minigraph", "-t", str(args.t),
            "-c", gfa_output_path, reads_fa]
gaf_output_path = os.path.join(args.d, "reads.gaf")
with open(gaf_output_path, 'w') as file:
    run(map_args, check=True, stdout=file, stderr=dump)


# #######
print('Genotype SV calls...')
# #######

pack_output_path = os.path.join(args.d, "reads.pack")
pack_args = ["vg", "pack", "-t", str(args.t), "-e",
             "-x", vg_output_path, "-o", pack_output_path,
             '-a', gaf_output_path]
run(pack_args, check=True, stdout=sys.stdout, stderr=dump)
call_args = ["vg", "call", "-t", str(args.t),
             "-k", pack_output_path, '-v', gt_vcf_gz, vg_output_path]
call_output_path = os.path.join(args.d, "called.vcf")
with open(call_output_path, 'w') as file:
    run(call_args, check=True, stdout=file, stderr=dump)


# #######
print('Annotate input VCF...')
# #######

# parse genotyped SVs
val_score = {}
called_vcf = VCF(call_output_path)
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

# delete temporary files
for ftrm in [gt_vcf_gz, gt_vcf_gz + '.tbi', reads_fa, ref_fa, ref_fa + '.fai',
             vg_output_path, gfa_output_path, gaf_output_path,
             pack_output_path, call_output_path]:
    os.remove(ftrm)
