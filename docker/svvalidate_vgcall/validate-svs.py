import argparse
from subprocess import run
from cyvcf2 import VCF, Writer
import os
import sys
import hashlib


# function to write a VCF for one SV
# that VCF file is later used to build a pangenome
def write_single_sv_vcf(sv_info, vcf_path):
    outf = open(vcf_path, 'wt')
    outf.write("##fileformat=VCFv4.2\n")
    outf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    rec_to_format = "{seqn}\t{pos}\t.\t{ref}\t{alt}\t.\t.\t.\n"
    outf.write(rec_to_format.format(seqn=sv_info['seqn'],
                                    pos=sv_info['start'],
                                    ref=sv_info['ref'],
                                    alt=sv_info['alt']))
    outf.close()
    bgzip_args = ["bgzip", "-c", vcf_path]
    vcf_path_gz = vcf_path + '.gz'
    with open(vcf_path_gz, 'w') as file:
        run(bgzip_args, check=True, stdout=file,
            stderr=sys.stderr, universal_newlines=True)
    tabix_args = ["tabix", vcf_path_gz]
    run(tabix_args, check=True, stdout=sys.stdout,
        stderr=sys.stderr, universal_newlines=True)
    return (vcf_path_gz)

# function to evaluate a SV
# sv_info is a dict with a 'svid', position and ref/alt sequences information
def evaluate_sv(sv_info, ref_fa_path, bam_path, output_dir,
                debug_mode=False, nb_cores=2, chr_lens={}):
    dump = open('/dev/null', 'w')
    # make VCF with just the one SV
    vcf_path = os.path.join(output_dir, sv_info['svid'] + ".vcf")
    vcf_path_gz = write_single_sv_vcf(sv_info, vcf_path)
    # decide on the region to consider
    # for now take the SV position and add some flanking regions
    # make sure to not go over the chromosome boundaruies though
    flank_size = 50000
    region_start = max(0, sv_info['start'] - flank_size)
    region_end = sv_info['end'] + flank_size
    if sv_info['seqn'] in chr_lens:
        region_end = min(region_end, chr_lens[sv_info['seqn']])
    region_coord_vg = '{}:{}-{}'.format(sv_info['seqn'],
                                        region_start,
                                        region_end)
    # also use smaller flanks to extract the reads to make sure they align
    flank_size = 10000
    region_start = max(0, sv_info['start'] - flank_size)
    region_end = sv_info['end'] + flank_size
    if sv_info['seqn'] in chr_lens:
        region_end = min(region_end, chr_lens[sv_info['seqn']])
    region_coord = '{}:{}-{}'.format(sv_info['seqn'],
                                     region_start,
                                     region_end)
    # make graph with SV
    construct_args = ["vg", "construct", "-a", "-m", "1024", "-S",
                      "-r", ref_fa_path, "-v", vcf_path_gz,
                      "-R", region_coord_vg]
    vg_output_path = os.path.join(output_dir, sv_info['svid'] + ".vg")
    with open(vg_output_path, 'w') as file:
        run(construct_args, check=True, stdout=file,
            stderr=sys.stderr, universal_newlines=True)
    # extract reads
    extract_args = ["samtools", "view", "-h", bam_path, region_coord]
    sam_output_path = os.path.join(output_dir, sv_info['svid'] + ".sam")
    with open(sam_output_path, 'w') as file:
        run(extract_args, check=True, stdout=file,
            stderr=sys.stderr, universal_newlines=True)
    extract_args = ["samtools", "fasta", sam_output_path]
    fa_output_path = os.path.join(output_dir, sv_info['svid'] + ".fasta")
    with open(fa_output_path, 'w') as file:
        run(extract_args, check=True, stdout=file,
            stderr=dump, universal_newlines=True)
    # align reads to pangenome
    convert_args = ["vg", "convert", "-f", vg_output_path]
    gfa_output_path = os.path.join(output_dir, sv_info['svid'] + ".gfa")
    with open(gfa_output_path, 'w') as file:
        run(convert_args, check=True, stdout=file,
            stderr=sys.stderr, universal_newlines=True)
    map_args = ["minigraph", "-t", str(nb_cores),
                "-c", gfa_output_path, fa_output_path]
    gaf_output_path = os.path.join(output_dir, sv_info['svid'] + ".gaf")
    with open(gaf_output_path, 'w') as file:
        run(map_args, check=True, stdout=file,
            stderr=dump, universal_newlines=True)
    # genotype SV
    pack_output_path = os.path.join(output_dir, sv_info['svid'] + ".pack")
    pack_args = ["vg", "pack", "-t", str(nb_cores), "-e",
                 "-x", vg_output_path, "-o", pack_output_path,
                 '-a', gaf_output_path]
    run(pack_args, check=True, stdout=sys.stdout,
        stderr=sys.stderr, universal_newlines=True)
    call_args = ["vg", "call", "-t", str(nb_cores),
                 "-k", pack_output_path, '-v', vcf_path, vg_output_path]
    call_output_path = os.path.join(output_dir,
                                    sv_info['svid'] + ".called.vcf")
    with open(call_output_path, 'w') as file:
        run(call_args, check=True, stdout=file,
            stderr=sys.stderr, universal_newlines=True)
    # update SV information or return a score
    score = {'prop': -1, 'ad': ''}
    for variant in VCF(call_output_path):
        ad = variant.format('AD')[0]
        dp = ad[0] + ad[1]
        score['ad'] = '{}|{}'.format(ad[0], ad[1])
        if dp == 0:
            score['prop'] = 0
        else:
            score['prop'] = float(ad[1]) / dp
    # remove intermediate files
    if not debug_mode:
        for ff in [sam_output_path, vg_output_path,
                   fa_output_path, gfa_output_path,
                   gaf_output_path, pack_output_path, call_output_path,
                   vcf_path, vcf_path_gz, vcf_path_gz + '.tbi']:
            os.remove(ff)
    dump.close()
    return (score)


parser = argparse.ArgumentParser()
parser.add_argument('-b', help='BAM file (indexed)', required=True)
parser.add_argument('-f', help='reference FASTA file (indexed)', required=True)
parser.add_argument('-v', help='variants in VCF (can be bgzipped)',
                    required=True)
parser.add_argument('-d', help='output directory', default='temp_valsv')
parser.add_argument('-o', default='out.vcf',
                    help='output (annotated) VCF (will be bgzipped if ending in .gz)')
parser.add_argument('-t', default=2,
                    help='number of threads used by the tools (vg and minigraph))')
args = parser.parse_args()

DEBUG_MODE = True

# index reference fasta if needed
if not os.path.isfile(args.f + '.fai'):
    index_fasta(args.f)
# extract chromosome lengths from the .fai (2nd column)
chr_lens = {}
with open(args.f + '.fai', 'rt') as inf:
    for line in inf:
        line = line.rstrip().split('\t')
        chr_lens[line[0]] = int(line[1])

vcf = VCF(args.v)
vcf.add_info_to_header({'ID': 'RS_PROP',
                        'Description': 'Proportion of supporting reads (from vg genotyping)',
                        'Type': 'Float', 'Number': '1'})
vcf.add_info_to_header({'ID': 'RS_AD',
                        'Description': 'Number of reads supporting the ref|alt (from vg genotyping)',
                        'Type': 'String', 'Number': '1'})
vcf_o = Writer(args.o, vcf)

# Read VCF and evaluate each SV
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
        score = evaluate_sv(svinfo, args.f, args.b, args.d,
                            DEBUG_MODE, nb_cores=args.t,
                            chr_lens=chr_lens)
        variant.INFO["RS_PROP"] = score['prop']
        variant.INFO["RS_AD"] = score['ad']
    vcf_o.write_record(variant)

vcf_o.close()
vcf.close()
