import argparse

def get_args():
    parser = argparse.ArgumentParser()

    parser.add_argument('-N',
                        dest='n',
                        required=True,
                        help='Number of expected variants')
    parser.add_argument('--mac',
                        dest='macs',
                        required=True,
                        help='Provided mac bins with proportions')
    parser.add_argument('-o',
                        dest='output',
                        help='Output file to be written')

    args = parser.parse_args()

    return args

def expected_variants(macs, n):
    ret = []
    for tup in macs:
        ret.append((tup[0], tup[1], tup[2]*n))
    return ret

def write_expected_variants(out_file, n, macs):
    with open(out_file, 'w') as output:
        output.writelines("Lower\tUpper\tExpected_var\n")
        for tup in macs:
            line = line.split(',')
            output.writelines(f"{tup[0]}\t{tup[1]}\t{tup[2]*n}\n")

def main():
    args = get_args()
    n = int(args.n)
    macs = []
    with open(args.macs) as macs_file:
        lines = macs_file.readlines()
        for line in lines:
            line = line.split(',')
            macs.append((int(line[0]), int(line[1]), float(line[2])))

    write_expected_variants(args.output, n, macs)


if __name__ == '__main__':
    main()
        
