import fileinput

hd = ''
ofmt = '\t'.join(['{}'] * 7)
print(ofmt.format('seqnames', 'pos', 'ref', 'alt',
                  'sample', 'haplotype', 'phased'))

for line in fileinput.input():
    line = line.rstrip().split('\t')
    # skip headers
    if line[0][0] == '#':
        if line[0] == '#CHROM':
            hd = line
        continue
    # skip already if variant looks short
    if len(line[3]) < 10 and len(line[4]) < 10:
        continue
    # process each variant
    alts = line[4].split(',')
    gt_idx = line[8].split(':').index('GT')
    # check genotype of each sample
    for s_ii in range(9, len(line)):
        s_phased = True
        s_gt_raw = line[s_ii].split(':')[gt_idx]
        s_gt = s_gt_raw.split('|')
        if '/' in s_gt_raw:
            s_gt = s_gt_raw.split('/')
            s_phased = False
        if len(s_gt) == 2 and s_gt[0] == s_gt[1]:
            # homozygous calls are phased in the sense that
            # we know they are in both haplotypes
            s_phased = True
        for al_ii, s_al in enumerate(s_gt):
            if s_al != '0' and s_al != '.':
                if len(line[3]) < 10 and len(alts[int(s_al)-1]) < 10:
                    continue
                # non-ref allele, write record
                if s_phased:
                    rec = ofmt.format(line[0],
                                      line[1],
                                      len(line[3]),
                                      len(alts[int(s_al)-1]),
                                      hd[s_ii],
                                      al_ii,
                                      s_phased)
                    print(rec)
                else:
                    for al_ii, dump in enumerate(s_gt):
                        rec = ofmt.format(line[0],
                                          line[1],
                                          len(line[3]),
                                          len(alts[int(s_al)-1]),
                                          hd[s_ii],
                                          al_ii,
                                          s_phased)
                        print(rec)
