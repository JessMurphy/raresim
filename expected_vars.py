import argparse
from nvariants import nvariants,fit_nvars
from afs import afs,fit_afs
from expected_variants import expected_variants, write_expected_variants
import os
import pandas as pd

DEFAULT_PARAMS = {
    'AFR': {"phi":0.1576, "omega":0.6247, "alpha": 1.5883, "beta": -0.3083, "b": 0.2872},
    'EAS': {"phi":0.1191, "omega":0.6369, "alpha": 1.6656, "beta": -0.2951, "b": 0.3137},
    'NFE': {"phi":0.1073, "omega":0.6539, "alpha": 1.9470, "beta": 0.1180, "b": 0.6676},
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

def check_for_stratification(file):
    with open(file) as f:
        lines = f.readlines()
        header = lines[0].split('\t')
        if len(header) > 3:
            return True
        return False

def main():
    args = get_args()
    n = int(args.n)
    macs = read_mac_bins(args.mac)

    # Validate inputs
    if args.pop:
        if args.pop.strip() not in DEFAULT_PARAMS.keys():
            raise Exception(f"{args.pop} is not a valid population. Valid populations are: {','.join(DEFAULT_PARAMS.keys())}")
    elif (args.nvar_target_data is None or args.afs_target_data is None) and (args.alpha is None or args.beta is None or args.omega is None or args.phi is None or args.b is None):
        raise Exception('Error: either a default population should be specified or alpha, beta, omega, phi, and b parameters provided or target values for nvars and afs should be provided')
    is_stratified=False
    if args.afs_target_data is not None:
        is_stratified = check_for_stratification(args.afs_target_data)

    if not is_stratified:
        # Get alpha, beta, omega, phi, b
        if args.pop:
            pop = args.pop.strip()
            alpha = DEFAULT_PARAMS[pop]['alpha']
            beta = DEFAULT_PARAMS[pop]['beta']
            omega = DEFAULT_PARAMS[pop]['omega']
            phi = DEFAULT_PARAMS[pop]['phi']
            b = DEFAULT_PARAMS[pop]['b']
        elif args.nvar_target_data is not None and args.afs_target_data is not None:
            df_afs = pd.read_csv(args.afs_target_data, delimiter='\t')
            df_nvar = pd.read_csv(args.nvar_target_data, delimiter='\t')

            alpha, beta, b = fit_afs(df_afs)
            omega, phi = fit_nvars(df_nvar)
        else:
            alpha = float(args.alpha)
            beta = float(args.beta)
            omega = float(args.omega)
            phi = float(args.phi)
            b = float(args.b)

        num_variants = nvariants(n, omega, phi) * float(args.reg_size)
        rows = afs(alpha, beta, b, macs)
        write_expected_variants(args.output, num_variants, rows)

    else:
        if args.nvar_target_data is None or args.afs_target_data is None:
            raise Exception('Error: stratification is currently only supported when target data is provided for nvars and afs')

        df_nvar_fun = pd.read_csv(args.nvar_target_data, delimiter='\t')
        df_nvar_fun.drop('syn_per_kb', axis=1, inplace=True)

        df_nvar_syn = pd.read_csv(args.nvar_target_data, delimiter='\t')
        df_nvar_syn.drop('fun_per_kb', axis=1, inplace=True)

        df_afs_fun = pd.read_csv(args.afs_target_data, delimiter='\t')
        df_afs_fun.drop('syn_prop', axis=1, inplace=True)
        df_afs_fun.rename(columns={'fun_prop' : 'Prop'}, inplace=True)

        df_afs_syn = pd.read_csv(args.afs_target_data, delimiter='\t')
        df_afs_syn.drop('fun_prop', axis=1, inplace=True)
        df_afs_syn.rename(columns={'syn_prop' : 'Prop'}, inplace=True)

        # Get values and write for Synonymous first
        alpha, beta, b = fit_afs(df_afs_syn)
        omega, phi = fit_nvars(df_nvar_syn)
        num_variants = nvariants(n, omega, phi) * float(args.reg_size)
        rows = afs(alpha, beta, b, macs)
        syn_output_file = os.path.splitext(args.output)[0] + '_syn.txt'
        write_expected_variants(syn_output_file, num_variants, rows)

        # Now do it for Functional
        alpha, beta, b = fit_afs(df_afs_fun)
        omega, phi = fit_nvars(df_nvar_fun)
        num_variants = nvariants(n, omega, phi) * float(args.reg_size)
        rows = afs(alpha, beta, b, macs)
        fun_output_file = os.path.splitext(args.output)[0] + '_fun.txt'
        write_expected_variants(fun_output_file, num_variants, rows)

if __name__ == '__main__':
    main()
