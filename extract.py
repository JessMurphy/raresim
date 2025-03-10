import random
import argparse
import gzip
import os


def get_args():
    parser = argparse.ArgumentParser()

    parser.add_argument('-i',
                        dest='input_file',
                        required=True,
                        help='Input Haplotype file')

    parser.add_argument('-o',
                        dest='output_file',
                        required=True,
                        help='Ouput cases path')
    
    parser.add_argument('--seed',
                        dest='seed',
                        type=int,
                        help='Seed for random sample')

    parser.add_argument('-n',
                        dest='num',
                        type=int,
                        help='Number of haplotypes to extract')

    args = parser.parse_args()

    return args

def main():
    args = get_args()
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


if __name__ == '__main__':
    main()
