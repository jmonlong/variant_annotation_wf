from pyfaidx import Fasta
import gzip
import random
import argparse

parser = argparse.ArgumentParser(description="simulate reads")
parser.add_argument('-v', help='VCF file with variants')
parser.add_argument('-f', help='reference Fasta file')
parser.add_argument('-r', help='reference region to extract')
parser.add_argument('-l', type=int, help='read length', default=1000)
parser.add_argument('-n', type=int, help='read number', default=1000)
args = parser.parse_args()

# ref_fa = '../../../../genomedata/hg38.fa'
# region = 'chr7:97357200-104567500'
# var_vcf = 'svs.vcf.gz'

# load genomic region from reference fasta
ref_fa = Fasta(args.f)
region_s = args.r.split(':')
region_r = [int(x) for x in region_s[1].split('-')]
# into two haplotypes
h1 = str(ref_fa[region_s[0]][region_r[0]:region_r[1]])
h1_offset = region_r[0]
h2 = str(ref_fa[region_s[0]][region_r[0]:region_r[1]])
h2_offset = region_r[0]

# read sorted VCF and edit haplotypes
in_vcf = gzip.open(args.v, 'rt')
for line in in_vcf:
    if line[0] == '#':
        continue
    line = line.rstrip().split('\t')
    # find GT
    fmt_v = line[8].split(':')
    samp_v = line[9].split(':')
    gt = samp_v[fmt_v.index('GT')]
    if '/' in gt:
        gt = gt.split('/')
    else:
        gt = gt.split('|')
    gt = [int(x) for x in gt]
    # prepare alleles
    ref = line[3]
    alts = line[4].split(',')
    ref_r = [int(line[1]) - h1_offset]
    ref_r += [ref_r[0] + len(ref)]
    # edit haplotypes
    if gt[0] != '0':
        h1 = h1[:ref_r[0]] + alts[gt[0]-1] + h1[ref_r[1]:]
        h1_offset += len(ref) - len(alts[gt[0]-1])
    if gt[1] != '0':
        h2 = h2[:ref_r[0]] + alts[gt[1]-1] + h1[ref_r[1]:]
        h2_offset += len(ref) - len(alts[gt[1]-1])
in_vcf.close()

# simulate reads from haplotypes
fastq_f = open('reads.fastq', 'wt')
readid = 0
for rr in range(args.n):
    # pick an haplotype randomly
    if random.random() > .5:
        pos = random.randint(0, len(h1) - args.l)
        read = h1[pos:(pos + args.l)]
    else:
        pos = random.randint(0, len(h2) - args.l)
        read = h2[pos:(pos + args.l)]
    fastq_f.write('@r{}\n{}\n+\n{}\n'.format(readid, read, '~'*len(read)))
    readid += 1
fastq_f.close()
