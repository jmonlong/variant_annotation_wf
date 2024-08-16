import argparse
from cyvcf2 import VCF


parser = argparse.ArgumentParser()
parser.add_argument('-v', help='variants in VCF (can be bgzipped)',
                    required=True)
parser.add_argument('-o', default='vcf_stats.tsv',
                    help='output TSV with stats')
args = parser.parse_args()

# args = parser.parse_args(['-v', 'SG190182.hapdiff.svs.vcf.gz', '-o', 'test.tsv'])

vartypes = {}
nhet = {}
nhom = {}
nhet_phased = {}
nhet_notphased = {}
ps_coord = {}

# Read VCF
vcf = VCF(args.v)
for variant in vcf:
    # filter homozygous ref variants
    if variant.num_het == 0 and variant.num_hom_alt == 0:
        continue
    # update the boundary of phase block is PS is present
    if 'PS' in variant.FORMAT:
        ps = variant.format('PS')[0][0]
        if ps in ps_coord:
            ps_coord[ps]['start'] = min(variant.start, ps_coord[ps]['start'])
            ps_coord[ps]['end'] = max(variant.end, ps_coord[ps]['end'])
        else:
            ps_coord[ps] = {'start': variant.start, 'end': variant.end}
    # check the genotype for each allele
    for ali in range(len(variant.ALT)):
        # init with the SVTYPE info, if available
        var_type = variant.INFO.get('SVTYPE')
        # even if SVTYPE is here, check that it's not clearly a SNV/indel
        if var_type is None or \
           (len(variant.REF) <= 30 and len(variant.ALT[ali]) <= 30):
            # guess from allele sizes
            if len(variant.REF) == 1 and len(variant.ALT[ali]) == 1:
                var_type = 'SNV'
            elif len(variant.REF) > 30 or len(variant.ALT[ali]) > 30:
                if len(variant.REF) > len(variant.ALT[ali]):
                    var_type = 'INS'
                else:
                    var_type = 'DEL'
            else:
                var_type = 'indel'
        # potentially init new vartype
        if var_type not in vartypes:
            vartypes[var_type] = True
            nhet[var_type] = 0
            nhom[var_type] = 0
            nhet_phased[var_type] = 0
            nhet_notphased[var_type] = 0
        # count alleles in genotype
        gt = variant.genotypes[0]
        nal_i = 0
        if gt[0] == ali + 1:
            nal_i += 1
        if gt[1] == ali + 1:
            nal_i += 1
        # update het/hom counts
        if nal_i == 1:
            nhet[var_type] += 1
            if nal_i > 0 and gt[2]:
                nhet_phased[var_type] += 1
            if nal_i > 0 and not gt[2]:
                nhet_notphased[var_type] += 1
        elif nal_i == 2:
            nhom[var_type] += 1
vcf.close()

# compute phase block stats
ps_sizes = []
for ps in ps_coord:
    ps_sizes.append(ps_coord[ps]['end'] - ps_coord[ps]['start'])
ps_sizes = sorted(ps_sizes)
half_tot_size = sum(ps_sizes)/2
ps_n50 = 0
ps_csum = 0
for pss in ps_sizes:
    if ps_csum > half_tot_size:
        break
    ps_n50 = pss
    ps_csum += pss

# write output TSV
with open(args.o, 'wt') as outf:
    outf.write('metric\tvalue\n')
    # write phase block metric
    if len(ps_sizes) > 0:
        outf.write('{}\t{}\n'.format('phase_block_n50', ps_n50))
    # write variant counts
    for vtype in vartypes:
        if vtype in nhet:
            outf.write('n_{}_het\t{}\n'.format(vtype, nhet[vtype]))
        if vtype in nhet_phased:
            outf.write('n_{}_het_phased\t{}\n'.format(vtype,
                                                      nhet_phased[vtype]))
        if vtype in nhet_notphased:
            outf.write('n_{}_het_notphased\t{}\n'.format(vtype,
                                                         nhet_notphased[vtype]))
        if vtype in nhom:
            outf.write('n_{}_hom\t{}\n'.format(vtype, nhom[vtype]))
