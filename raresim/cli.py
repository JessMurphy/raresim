from raresim.common.sparse import SparseMatrixReader, SparseMatrixWriter
from raresim.engine.runner import DefaultRunner
from raresim.engine.config import RunConfig
from raresim.calculate.expected_vars import calc
import argparse
import random
import os
import gzip


def parseCommand():
    parser = argparse.ArgumentParser()

    subparsers = parser.add_subparsers(dest='command', required=True, help='Available commands: sim, convert, extract')

    sim_parser = subparsers.add_parser('sim')
    convert_parser = subparsers.add_parser('convert')
    extract_parser = subparsers.add_parser('extract')
    calc_parser = subparsers.add_parser('calc')

    calc_parser.add_argument('--mac',
                        dest='mac',
                        required=True,
                        help='Provided mac bins with proportions')
    calc_parser.add_argument('-o',
                        dest='output',
                        required=True,
                        help='Output file to be written')
    calc_parser.add_argument('-N',
                        dest='n',
                        required=True,
                        help='Number of expected variants')
    calc_parser.add_argument('--pop',
                        dest='pop',
                        help='Population to use default values for if not providing alpha, beta, omega, phi, b, or target values')
    calc_parser.add_argument('--alpha',
                        dest='alpha',
                        help='Provided alpha value')
    calc_parser.add_argument('--beta',
                        dest='beta',
                        help='Provided beta value')
    calc_parser.add_argument('--omega',
                        dest='omega',
                        help='Provided omega value')
    calc_parser.add_argument('--phi',
                        dest='phi',
                        help='Provided phi value')
    calc_parser.add_argument('-b',
                        dest='b',
                        help='Provided b value')
    calc_parser.add_argument('--nvar_target_data',
                        dest='nvar_target_data',
                        help='Provided target values for nvars')
    calc_parser.add_argument('--afs_target_data',
                        dest='afs_target_data',
                        help='Provided target values for afs')
    calc_parser.add_argument('--reg_size',
                        dest='reg_size',
                        help='Region size in kilobases')
    calc_parser.add_argument('-w',
                        dest='w',
                        default='1.0',
                        help='Weight value to multiple total number of variants by in non-stratified runs. Valid range of values is [0,2] with a default of 1')
    calc_parser.add_argument('--w_fun',
                        dest='w_fun',
                        default='1.0',
                        help='Weight value to multiple total number of functional variants by in stratified runs. Valid range of values is [0,2] with a default of 1')
    calc_parser.add_argument('--w_syn',
                        dest='w_syn',
                        default='1.0',
                        help='Weight value to multiple total number of synonymous variants by in stratified runs. Valid range of values is [0,2] with a default of 1')

    extract_parser.add_argument('-i',
                                dest='input_file',
                                required=True,
                                help='Input cases path')
    extract_parser.add_argument('-o',
                                dest='output_file',
                                required=True,
                                help='Output cases path')
    extract_parser.add_argument('-s', '--seed',
                                dest='seed',
                                type=int,
                                help='Seed for random sample')
    extract_parser.add_argument('-n',
                                dest='num',
                                type=int,
                                required=True,
                                help='Number of haplotypes to extract')

    sim_parser.add_argument('-m',
                        dest='sparse_matrix',
                        required=True,
                        help='Input sparse matrix path, can be a .haps, .sm, or .gz file')

    sim_parser.add_argument('-b',
                        dest='exp_bins',
                        help='Input expected bin sizes')

    sim_parser.add_argument('--functional_bins',
                        dest='exp_fun_bins',
                        help='Input expected bin sizes for functional variants')

    sim_parser.add_argument('--synonymous_bins',
                        dest='exp_syn_bins',
                        help='Input expected bin sizes for synonymous variants')

    sim_parser.add_argument('-l',
                        dest='input_legend',
                        required=True,
                        help='Input variant site legend')

    sim_parser.add_argument('-L',
                        dest='output_legend',
                        help='Output variant site legend')

    sim_parser.add_argument('-H',
                        dest='output_hap',
                        required=True,
                        help='Output compress hap file')

    sim_parser.add_argument('--f_only',
                        dest='fun_bins_only',
                        help='Input expected bin sizes for only functional variants')

    sim_parser.add_argument('--s_only',
                        dest='syn_bins_only',
                        help='Input expected bin sizes for synonymous variants only')

    sim_parser.add_argument('-z',
                        action='store_true',
                        help='Rows of zeros and pruned rows are removed')

    sim_parser.add_argument('-prob',
                        action='store_true',
                        help='Rows are pruned allele by allele given a probability of removal')

    sim_parser.add_argument('--small_sample',
                        action='store_true',
                        help='Override error to allow for simulation of small sample size')

    sim_parser.add_argument('--keep_protected',
                        action='store_true',
                        help='Rows in the legend marked with a 1 in the protected column will be accounted'
                             ' for but not pruned')

    sim_parser.add_argument('--stop_threshold',
                        dest='stop_threshold',
                        default='20',
                        help='Percentage threshold for the pruning process 0-100. Provides a stop to prevent us from going the given % below the expected count for any given bin during pruning. Default value of 20.')

    sim_parser.add_argument('--activation_threshold',
                        dest='activation_threshold',
                        default='10',
                        help='Percentage threshold for activation of the pruning process. Requires that the actual count for a bin must be more than the given percentage different from the expected count to activate pruning on the bin.')

    sim_parser.add_argument('--verbose',
                        action='store_true',
                        help='Rows in the legend marked with a 1 in the protected column will be accounted for but not pruned')

    convert_parser.add_argument('-i',
                                dest='input_file',
                                required=True,
                                help='Input sparse matrix path, can be a .haps, .sm, or .gz file')

    convert_parser.add_argument('-o',
                                dest='output_file',
                                required=True,
                                help='Output haplotype file path')

    args = parser.parse_args()

    return args

def extract(args):
    random.seed(args.seed)
    with gzip.open(args.input_file, 'rt') as f:
        line = f.readline()
        columns = line.split()
    size = len(columns)
    columnsToExtract = random.sample(range(0, size), args.num)
    otherColumns = [i for i in range(size) if i not in columnsToExtract]
    columnsToExtract.sort()
    base, ext = os.path.splitext(args.output_file)
    output_file_name = base
    with gzip.open(f'{output_file_name}-sample.gz', 'wb') as s:
        with gzip.open(f'{output_file_name}-remainder.gz', 'wb') as r:
            with gzip.open(args.input_file, 'rt') as input_haps:
                for l in input_haps.readlines():
                    cols = l.split()
                    sampleLine = [cols[i] for i in columnsToExtract]
                    remainderLine = [cols[i] for i in otherColumns]
                    s.write((" ".join(sampleLine) + "\n").encode())
                    r.write((" ".join(remainderLine) + "\n").encode())


def main():
    command = parseCommand()
    if command.command == 'sim':
        runConfig = RunConfig(command)
        runner: DefaultRunner = DefaultRunner(runConfig)
        runner.run()

    elif command.command == 'convert':
        args = command
        reader = SparseMatrixReader()
        writer = SparseMatrixWriter()
        matrix = reader.loadSparseMatrix(args.input_file)
        writer.writeToHapsFile(matrix, args.output_file, "sm")

    elif command.command == 'extract':
        args = command
        extract(args)

    elif command.command == 'calc':
        args = command
        calc(args)


if __name__ == '__main__':
    main()
