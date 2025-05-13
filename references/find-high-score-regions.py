# Reads input TSV stream and outputs regions where bases have a score greater than a threshold.
import fileinput


# which column to use for the score? 2: raw score; 3: percentile
score_col_idx = 3
# minimum score of bases in the region
min_score = .99
# extend up to X bases below the minimum score
extend_bp = 10
min_extend_score = .9

# information about the region in progress
cur_seqn = ''
cur_pos = -1
cur_start = -1
cur_score_sum = 0
cur_extend_score_sum = 0

for line in fileinput.input():
    line = line.rstrip().split('\t')
    pos = int(line[1])
    score = float(line[score_col_idx])
    if score >= min_score:
        new_region = True
        if line[0] == cur_seqn and pos < cur_pos + 2 + extend_bp:
            # we could extend, let's check the score with extension
            extend_score = (cur_score_sum + score + cur_extend_score_sum) / (pos - cur_start + 1)
            if extend_score >= min_extend_score:
                # extend current region
                cur_pos = pos
                cur_score_sum += score + cur_extend_score_sum
                cur_extend_score_sum = 0
                new_region = False
            else:
                # score would drop too much, don't extent
                # write in-progress region
                mean_score = cur_score_sum / (cur_pos - cur_start + 1)
                print('{}\t{}\t{}\t{}'.format(cur_seqn, cur_start,
                                              cur_pos, round(mean_score, 3)))
        if new_region:
            # first base of a high-score region, init
            cur_seqn = line[0]
            cur_pos = pos
            cur_start = cur_pos
            cur_score_sum = score
            cur_extend_score_sum = 0
    else:
        # low-score base
        if cur_start > 0 and pos > cur_pos + extend_bp:
            # Write in-progress region
            mean_score = cur_score_sum / (cur_pos - cur_start + 1)
            print('{}\t{}\t{}\t{}'.format(cur_seqn, cur_start,
                                          cur_pos, round(mean_score, 3)))
            cur_seqn = ''
            cur_start = -1
            cur_pos = -1
            cur_score_sum = 0
        else:
            # maybe we're trying to extend, so save the scores in case
            cur_extend_score_sum += score

# write last region if one was in progress
if cur_start > 0:
    mean_score = cur_score_sum / (cur_pos - cur_start + 1)
    print('{}\t{}\t{}\t{}'.format(cur_seqn, cur_start,
                                  cur_pos, round(mean_score, 3)))
