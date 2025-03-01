import argparse
import random
import gzip
from heapq import merge
import math

class Error(Exception):
    """Base class for other exceptions"""
    pass

class DifferingLengths(Error):
    """Raised when two file lengths differ"""
    pass

class MissingColumn(Error):
    """Raised when legend is missing a column"""
    pass

class MissingProbs(Error):
    """Raised when legend is missing probs column"""
    pass


def get_bin(bins, val):
    for i in range(len(bins)):
        if val >= bins[i][0] and val <= bins[i][1]:
            return i
    return i+1

def get_args():
    parser = argparse.ArgumentParser()

    parser.add_argument('-m',
                        dest='sparse_matrix',
                        required=True,
                        help='Input sparse matrix path')

    parser.add_argument('-b',
                        dest='exp_bins',
                        help='Input expected bin sizes')

    parser.add_argument('--functional_bins',
                        dest='exp_fun_bins',
                        help='Input expected bin sizes for functional variants')

    parser.add_argument('--synonymous_bins',
                        dest='exp_syn_bins',
                        help='Input expected bin sizes for synonymous variants')

    parser.add_argument('-l',
                        dest='input_legend',
                        help='Input variant site legend')

    parser.add_argument('-L',
                        dest='output_legend',
                        help='Output variant site legend')

    parser.add_argument('-H',
                        dest='output_hap',
                        required=True,
                        help='Output compress hap file')

    parser.add_argument('--f_only',
                        dest='fun_bins_only',
                        help='Input expected bin sizes for only functional variants')

    parser.add_argument('--s_only',
                        dest='syn_bins_only',
                        help='Input expected bin sizes for synonymous variants only')

    parser.add_argument('-z',
                        action='store_true',
                        help='Rows of zeros and pruned rows are removed from input haps file')

    parser.add_argument('-prob',
                        action='store_true',
                        help='Rows are pruned allele by allele given a probability of removal')

    parser.add_argument('--small_sample',
                        action='store_true',
                        help='Override error to allow for simulation of small sample size')

    parser.add_argument('--keep_protected',
                        action='store_true',
                        help='Rows in the legend marked with a 1 in the protected column will be accounted for but not pruned')

    parser.add_argument('--stop_threshold',
                        dest='stop_threshold',
                        default='20',
                        help='Percentage threshold for the pruning process 0-100. Provides a stop to prevent us from going the given % below the expected count for any given bin during pruning. Default value of 20.')

    parser.add_argument('--activation_threshold',
                        dest='activation_threshold',
                        default='10',
                        help='Percentage threshold for activation of the pruning process. Requires that the actual count for a bin must be more than the given percentage different from the expected count to activate pruning on the bin.')

    args = parser.parse_args()

    return args




def prune_bins(bin_assignments, bins, extra_rows, matrix, activation_threshold, stopping_threshold, legend):
    '''
    @param bin_assignments: A map of bin_id -> list of rows assigned to the bin
    @param bins: the provided input bins with lower bound, upper bound, and expected count
    @param extra_rows: list of rows that were removed and never re-added
    @param matrix: sparse matrix of haps
    '''
    reserve_pool = bin_assignments[len(bin_assignments)-1]
    # Loop through the bins from largest to smallest
    for bin_id in range(len(bin_assignments)-1)[::-1]:
        # If there are any rows with too many 1s for the largest bin, they get put into an extra bin,
        # and we do not prune these alleles away
        if bin_id == len(bins):
            continue

        # How many variants to we need and how many do we have
        need = bins[bin_id][2]
        have = len(bin_assignments[bin_id])

        activation = min([(float(activation_threshold) / 100) * need, 10])
        stop = (float(stopping_threshold) / 100) * need

        # If we have more rows in the bin than what we need, remove some rows
        if have - need > activation:
            # Calculate the probability to remove any given row in the bin
            prob_remove = 1 - float(need) / float(have) if have != 0 else 0
            row_ids_to_rem = []
            for i in range(have):
                # If we hit our stopping threshold to stop the pruning process, then stop the pruning process
                if have - len(row_ids_to_rem) <= need - stop:
                    print (f"Stopping threshold reaching for bin {bin_id}")
                    break
                flip = random.uniform(0, 1)
                # If the row is 'chosen' for removal, remove it and add the row to the list of rows that may be used
                # to make up for not having enough rows in later bins
                if flip < prob_remove:
                    row_id = bin_assignments[bin_id][i]
                    if legend[row_id]['protected'] == '1':
                        print(f"WARNING: Attempting to prune a row that is protected. This should not happen and is a bug. RowId = {row_id}")
                    # Add the ith row in the bin to the list of row ids to remove and the list of available
                    # rows to pull from if needed later on
                    row_ids_to_rem.append(row_id)
                    extra_rows.append(row_id)
            for row_id in row_ids_to_rem:
                bin_assignments[bin_id].remove(row_id)

        # If we don't have enough rows in the current bin, pull from the list of excess rows
        elif have < round(need - activation):
            pullFromCommons = False
            if len(extra_rows) < round(abs(need - have)):
                if len(reserve_pool) < round(abs(need-have)):
                    raise Exception(f'ERROR: Current bin has {have} variants, but the model needs {need} variants '
                                    f'and only {len(extra_rows)} excess rows are available to use and prune down '
                                    f'from larger bins and only {len(reserve_pool)} common variants are available '
                                    f'to prune down.')
                else:
                    print(f"WARNING: No more extra pruned rows to pull from. Common variants may be pruned down for bin {bin_id+1} if necessary")
                    pullFromCommons = True

            # Calculate the probability to use any given row from the available list
            prob_add_from_extras = float(need - have) / float(len(extra_rows)) if len(extra_rows) > 0 else 1
            prob_add_from_reserve = float(need - have - len(extra_rows)) / float(len(reserve_pool))
            row_ids_to_add_from_extras = []
            for i in range(len(extra_rows)):
                flip = random.uniform(0, 1)
                # If the current row is 'chosen' for use, we will prune it down to the desired number of variants
                # and add it back in to the current bin and remove it from the list of available rows to pull from
                if flip < prob_add_from_extras:
                    row_id = extra_rows[i]
                    row_ids_to_add_from_extras.append(row_id)
                    num_to_keep = random.randint(bins[bin_id][0], bins[bin_id][1])
                    num_to_rem = matrix.row_num(row_id) - num_to_keep
                    matrix.prune_row(row_id, num_to_rem)
                    if matrix.row_num(row_id) != num_to_keep:
                        print(f"WARNING: Requested to prune row {row_id} down to {num_to_keep} variants, but after the request the row still has {matrix.row_num} variants. This should not happen and the issue likely needs to be raised with a developer.")
                    bin_assignments[bin_id].append(row_id)

            for row_id in row_ids_to_add_from_extras:
                extra_rows.remove(row_id)

            if pullFromCommons:
                row_ids_to_add_from_reserve = []
                for i in range(len(reserve_pool)):
                    flip = random.uniform(0, 1)
                    if flip < prob_add_from_reserve:
                        row_id = reserve_pool[i]
                        row_ids_to_add_from_reserve.append(row_id)
                        num_to_keep = random.randint(bins[bin_id][0], bins[bin_id][1])
                        num_to_rem = matrix.row_num(row_id) - num_to_keep
                        matrix.prune_row(row_id, num_to_rem)
                        if matrix.row_num(row_id) != num_to_keep:
                            print(f"WARNING: Requested to prune row {row_id} down to {num_to_keep} variants, but after the request the row still has {matrix.row_num} variants. This should not happen and the issue likely needs to be raised with a developer.")
                        bin_assignments[bin_id].append(row_id)
                for row_id in row_ids_to_add_from_reserve:
                    reserve_pool.remove(row_id)

def print_bin(bin_assignments, bins):
    print(f"{'Bin':<12}{'Expected':<18}{'Actual':<10}")
    for bin_id in range(len(bin_assignments)):
        if bin_id < len(bins):
            col1 = f"[{str(bins[bin_id][0])},{str(bins[bin_id][1])}]"
            col2 = f"{bins[bin_id][2]:.10f}"
            col3 = f"{str(len(bin_assignments[bin_id]))}"
        else:
            col1 = f"[{str(bins[bin_id-1][1]+1)},\u221E]"
            col2 = "N/A"
            col3 = f"{str(len(bin_assignments[bin_id]))}"
        print(f"{col1:<12}{col2:<18}{col3:<10}")


def read_legend(legend_file_name):
    header = None
    legend = []
    with open(legend_file_name, "r") as f:
        for l in f:
            A = l.rstrip().split()

            if header == None:
                header = A
            else:
                legend.append(dict(zip(header,A)))

    return header, legend

def read_expected(expected_file_name):

    bins = []

    header = None
    with open(expected_file_name) as f:
        for l in f:
            A = l.rstrip().split()

            if header == None:
                header = A
            else:
                bins.append([int(A[0]), int(A[1]), float(A[2])])

    return bins



def get_split(args):
    func_split = False
    fun_only = False
    syn_only = False

    if args.exp_bins is None and not args.prob:
        if args.exp_fun_bins is not None \
            and args.exp_syn_bins is not None:
            func_split = True
        elif args.fun_bins_only is not None:
            fun_only = True
        elif args.syn_bins_only is not None:
            syn_only = True
        else:
            raise Exception('If variants are split by functional/synonymous ' + \
                     'files must be provided for --functional_bins ' + \
                     'and --synonymous_bins')
    return func_split, fun_only, syn_only

def verify_legend(legend, legend_header, M, split, probs):
    if split and 'fun' not in legend_header and not probs:
        raise MissingColumn('If variants are split by functional/synonymous ' + \
                 'the legend file must have a column named "fun" ' + \
                 'that specifies "fun" or "syn" for each site')

    if M.num_rows() != len(legend):
        raise DifferingLengths(f"Lengths of legend {len(legend)} and hap {M.num_rows()} files do not match")

    if probs and 'prob' not in legend_header:
        raise MissingProbs('The legend file needs to have a "prob" column ' + \
                'to indicate the pruning probability of a given row ')




def assign_bins(matrix, bins, legend, is_func_split, is_fun_only, is_syn_only):
    bin_assignments = {}

    if is_func_split:
        bin_assignments['fun'] = {bin_id: [] for bin_id in range(len(bins['fun']) + 1)}
        bin_assignments['syn'] = {bin_id: [] for bin_id in range(len(bins['syn']) + 1)}
    else:
        bin_assignments = {bin_id: [] for bin_id in range(len(bins) + 1)}

    row_i = 0
    for row in range(matrix.num_rows()):
        row_num = matrix.row_num(row)

        if row_num > 0:
            if is_func_split:
                bin_id = get_bin(bins[legend[row_i]['fun']], row_num)
            elif is_syn_only:
                if legend[row_i]['fun'] != 'syn':
                    row_i += 1
                    continue
                bin_id = get_bin(bins, row_num)
            elif is_fun_only:
                if legend[row_i]['fun'] != 'fun':
                    row_i += 1
                    continue
                bin_id = get_bin(bins, row_num)
            else:
                bin_id = get_bin(bins, row_num)

            #Depending on split status, either append to bin_assignments or to just the annotated dictionary
            target_map = bin_assignments

            if is_func_split:
                target_map = bin_assignments[legend[row_i]['fun']]
            if bin_id not in target_map:
                target_map[bin_id] = []


            target_map[bin_id].append(row_i)

        row_i += 1

    return bin_assignments


def write_legend(all_kept_rows, input_legend, output_legend):
    header = None
    legend_file_index = 0
    kept_row_index = 0
    f = open(output_legend, 'w')
    r = open(f'{output_legend}-pruned-variants', 'w')
    with open(input_legend, "r") as in_file:
        for l in in_file.readlines():
            if header == None:
                f.write(l)
                header = True
                continue

            if kept_row_index < len(all_kept_rows) and legend_file_index == all_kept_rows[kept_row_index]:
                f.write(l)
                kept_row_index+=1
            else:
                r.write(l)
            legend_file_index+=1


def write_hap(all_kept_rows, output_file, M):
    step = int(len(all_kept_rows)/10)
    i = 0

    with gzip.open(output_file, 'wb') as f:
        for row_i in all_kept_rows:
            row = []
            for col_i in range(M.row_num(row_i)):
                row.append(M.get(row_i, col_i))

            O = ['0'] * M.num_cols()

            for col_i in row:
                O[col_i] = '1'

            s = ' '.join(O) + '\n'
            f.write(s.encode())

            if (i % step == 0):
                print('.', end='', flush=True)


            i+=1
    print()



def print_frequency_distribution(bins, bin_assignments, func_split, fun_only, syn_only):
    if func_split:
        print('Functional')
        print_bin(bin_assignments['fun'], bins['fun'])
        print('\nSynonymous')
        print_bin(bin_assignments['syn'], bins['syn'])
    elif fun_only:
        print('Functional')
        print_bin(bin_assignments, bins)
    elif syn_only:
        print('Synonymous')
        print_bin(bin_assignments, bins)
    else:
        print_bin(bin_assignments, bins)


def get_all_kept_rows(bin_assignments, R, func_split, fun_only, syn_only, legend):
    all_kept_rows = []

    if func_split:
        for bin_id in range(len(bin_assignments['fun'])):
            all_kept_rows += bin_assignments['fun'][bin_id]
        for bin_id in range(len(bin_assignments['syn'])):
            all_kept_rows += bin_assignments['syn'][bin_id]

    elif fun_only:
        for bin_id in bin_assignments:
            all_kept_rows += bin_assignments[bin_id]
        for row_id in range(len(legend)):
            if legend[row_id]['fun'] == 'syn':
                all_kept_rows.append(row_id)
    elif syn_only:
        for bin_id in bin_assignments:
            all_kept_rows += bin_assignments[bin_id]
        for row_id in range(len(legend)):
            if legend[row_id]['fun'] == 'fun':
                all_kept_rows.append(row_id)
    else:
        for bin_id in range(len(bin_assignments)):
            all_kept_rows += bin_assignments[bin_id]

    all_kept_rows.sort()

    if isinstance(R, dict):
        R = [item for sublist in R.values() for item in sublist]
    R = sorted(R)

    all_kept_rows = list(dict.fromkeys(all_kept_rows))
    return all_kept_rows


def get_expected_bins(args, func_split, fun_only, syn_only):
    bins = None
    if func_split:
        bins = {}
        bins['fun'] = read_expected(args.exp_fun_bins)
        bins['syn'] = read_expected(args.exp_syn_bins)
    elif fun_only:
        bins = read_expected(args.fun_bins_only)
    elif syn_only:
        bins = read_expected(args.syn_bins_only)
    else:
        bins = read_expected(args.exp_bins)
    return bins

def adjust_for_protected_variants(bins, bin_assignments, legend):
    ret = {bin_id : [] for bin_id in range(len(bins))}
    for bin_id in range(len(bins)):
        for row_id in bin_assignments[bin_id]:
            if legend[row_id]['protected'] == '1':
                bin_assignments[bin_id].remove(row_id)
                bins[bin_id][2] -= 1
                ret[bin_id].append(row_id)
    return ret

def add_protected_rows_back(bins, bin_assignments, protected_var_counts_per_bin):
    for bin_id in protected_var_counts_per_bin:
        if bin_id == len(bins): continue
        bins[bin_id][2] += len(protected_var_counts_per_bin[bin_id])
        bin_assignments[bin_id] += protected_var_counts_per_bin[bin_id]
        bin_assignments[bin_id] = sorted(bin_assignments[bin_id])


