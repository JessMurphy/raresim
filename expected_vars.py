import argparse
from nvariants import nvariants
from afs import afs
from expected_variants import expected_variants, write_expected_variants
import os

DEFAULT_PARAMS = {
    'AFR': {"phi":0.1576, "omega":0.6247, "alpha": 1.5883, "beta": -0.3083, "b": 0.2872},
    'EAS': {"phi":0.1191, "omega":0.6369, "alpha": 1.6656, "beta": -0.2951, "b": 0.3137},
    'NFE': {"phi":0.1073, "omega":0.6539, "alpha": 1.9470, "beta": -0.1180, "b": 0.6676},
    'SAS': {"phi":0.1249, "omega":0.6495, "alpha": 1.6977, "beta": -0.2273, "b": 0.3564}
}

def get_args():
    parser = argparse.ArgumentParser()

    parser.add_argument('--mac',
                        dest='mac',
                        required=True,
                        help='Provided mac bins with proportions')
    parser.add_argument('-o',
                        dest='output',
                        required=True,
                        help='Output file to be written')
    parser.add_argument('-N',
                        dest='n',
                        required=True,
                        help='Number of expected variants')
    parser.add_argument('--pop',
                        dest='pop',
                        help='Population to use default values for if not providing alpha, beta, omega, phi, b, or target values')
    parser.add_argument('--alpha',
                        dest='alpha',
                        help='Provided alpha value')
    parser.add_argument('--beta',
                        dest='beta',
                        help='Provided beta value')
    parser.add_argument('--omega',
                        dest='omega',
                        help='Provided omega value')
    parser.add_argument('--phi',
                        dest='phi',
                        help='Provided phi value')
    parser.add_argument('-b',
                        dest='b',
                        help='Provided b value')
    parser.add_argument('--nvar_target_data',
                        dest='nvar_target_data',
                        help='Provided target values for nvars')
    parser.add_argument('--afs_target_data',
                        dest='afs_target_data',
                        help='Provided target values for afs')
    parser.add_argument('--reg_size',
                        dest='reg_size',
                        help='Region size in kilobases')

    args = parser.parse_args()

    return args

def read_mac_bins(macs_file):
    root, extension = os.path.splitext(macs_file)
    macs = []
    if extension == '.csv':
        with open(macs_file) as macs_file:
            lines = macs_file.readlines()
            header = lines[0].split(',')
            if header[0].strip() != 'Lower' or header[1].strip() != 'Upper':
                raise Exception("Mac bins file needs to have column names Lower and Upper")

            for line in lines[1:]:
                l = line.strip().split(',')
                macs.append((int(l[0]), int(l[1])))
    elif extension == '.txt':
        with open(macs_file) as macs_file:
            lines = macs_file.readlines()
            header = lines[0].split('\t')
            if header[0].strip() != 'Lower' or header[1].strip() != 'Upper':
                raise Exception("Mac bins file needs to have column names Lower and Upper")

            for line in lines[1:]:
                l = line.strip().split('\t')
                macs.append((int(l[0]), int(l[1])))

    return macs

def main():
    args = get_args()

    # Validate inputs
    if args.pop:
        if args.pop.strip() not in DEFAULT_PARAMS.keys():
            raise Exception(f"{args.pop} is not a valid population. Valid populations are: {','.join(DEFAULT_PARAMS.keys())}")
    elif (args.nvar_target_data is None or args.afs_target_data is None) and (args.alpha is None or args.beta is None or args.omega is None or args.phi is None or args.b is None):
        raise Exception('Error: either a default population should be specified or alpha, beta, omega, phi, and b parameters provided or target values for nvars and afs should be provided')

    # Get alpha, beta, omega, phi, b
    # TODO: Add support for target values
    if args.pop:
        pop = args.pop.strip()
        alpha = DEFAULT_PARAMS[pop]['alpha']
        beta = DEFAULT_PARAMS[pop]['beta']
        omega = DEFAULT_PARAMS[pop]['omega']
        phi = DEFAULT_PARAMS[pop]['phi']
        b = DEFAULT_PARAMS[pop]['b']
    else:
        alpha = float(args.alpha)
        beta = float(args.beta)
        omega = float(args.omega)
        phi = float(args.phi)
        b = float(args.b)

    n = int(args.n)

    # Get mac bins
    macs = read_mac_bins(args.mac)

    num_variants = nvariants(n, omega, phi)
    rows = afs(alpha, beta, b, macs)
    write_expected_variants(args.output, num_variants, rows)

if __name__ == '__main__':
    main()
