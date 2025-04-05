import argparse
import os
import pandas as pd
import numpy as np
from scipy.optimize import minimize

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
    parser.add_argument('-w',
                        dest='w',
                        default='1.0',
                        help='Weight value to multiple total number of variants by in non-stratified runs. Valid range of values is [0,2] with a default of 1')
    parser.add_argument('--w_fun',
                        dest='w_fun',
                        default='1.0',
                        help='Weight value to multiple total number of functional variants by in stratified runs. Valid range of values is [0,2] with a default of 1')
    parser.add_argument('--w_syn',
                        dest='w_syn',
                        default='1.0',
                        help='Weight value to multiple total number of synonymous variants by in stratified runs. Valid range of values is [0,2] with a default of 1')

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

def write_expected_variants(out_file, n, macs):
    with open(out_file, 'w') as output:
        output.writelines("Lower\tUpper\tExpected_var\n")
        for tup in macs:
            output.writelines(f"{tup[0]}\t{tup[1]}\t{tup[2]*n}\n")

def fit_nvars(observed_variants_per_kb):
    # Check that each column is numeric
    if not np.issubdtype(observed_variants_per_kb.iloc[:, 0].dtype, np.number) or not np.issubdtype(
            observed_variants_per_kb.iloc[:, 1].dtype, np.number):
        raise ValueError('Columns need to be numeric')

    # Check that there are not any NA values
    if observed_variants_per_kb.iloc[:, 1].isna().any():
        raise ValueError('Number of variants per Kb need to be numeric with no NA values')

    # Check that the sample sizes go from smallest to largest
    if not observed_variants_per_kb.iloc[:, 0].is_monotonic_increasing:
        raise ValueError('The sample sizes need to be ordered from smallest to largest')

    # define the least squares loss function
    def leastsquares(inner_tune):
        # calculate the expected number of variants (from the function)
        E = inner_tune[0] * (observed_variants_per_kb.iloc[:, 0] ** inner_tune[1])
        sq_err = (E - observed_variants_per_kb.iloc[:, 1]) ** 2  # calculate the squared error of expected - observed
        return np.sum(sq_err)  # return the squared error

    def hin_tune(x):  # constraints
        h = np.zeros(3)
        h[0] = x[0]  # phi greater than 0
        h[1] = x[1]  # omega greater than 0
        h[2] = 1 - x[1]  # omega less than 1
        return h

    # define the starting value for phi so the end of the function matches with omega = 0.45
    phi = observed_variants_per_kb.iloc[observed_variants_per_kb.iloc[:, 0].idxmax(), 1] / (
            observed_variants_per_kb.iloc[observed_variants_per_kb.iloc[:, 0].idxmax(), 0] ** 0.45)
    tune = np.array([phi, 0.45])  # specify the starting values

    # Use SLSQP to find phi and omega
    constraints = {'type': 'ineq', 'fun': hin_tune}
    re_LS = minimize(leastsquares, tune, constraints=constraints, options={'disp': False, 'ftol': 0.0, 'maxiter': 100})

    # If the original starting value resulted in a large loss (>1000), iterate over starting values
    if re_LS.fun > 1000:
        re_tab1 = []  # create to hold the new parameters
        for omega in np.arange(0.15, 0.66, 0.1):  # optimize with different values of omega

            # specify phi to fit the end of the function with the current value of omega
            phi = observed_variants_per_kb.iloc[observed_variants_per_kb.iloc[:, 0].idxmax(), 1] / (
                    observed_variants_per_kb.iloc[observed_variants_per_kb.iloc[:, 0].idxmax(), 0] ** omega)
            tune = np.array([phi, omega])  # updated starting values

            re_LS1 = minimize(leastsquares, tune, constraints=constraints,
                              options={'disp': False, 'ftol': 0.0, 'maxiter': 100})  # estimate parameters with SLSQP
            to_bind1 = np.concatenate((re_LS1.x, [re_LS1.fun]))  # record parameters and loss value
            re_tab1.append(to_bind1)  # bind information from each iteration together

        re_tab1 = np.array(re_tab1)
        re_fin = re_tab1[np.argmin(re_tab1[:, 2])]  # select the minimum least squared error
        phi= re_fin[0]
        omega = re_fin[1]
    else:  # if the loss was <1000, bring the parameters forward
        phi = re_LS.x[0]
        omega = re_LS.x[1]

    print(f"Calculated the following params from nvar target data. omega: {omega}, phi: {phi}")
    return omega,phi

def nvariants(n, omega, phi, reg_size, weight):
    ret = float(phi) * (int(n)**float(omega)) * reg_size * weight
    print(f"Calculated {ret} total variants (accounting for region size)")
    return ret

def fit_afs(observed_bin_props_df):
    # Check the column names of Observed_bin_props
    if list(observed_bin_props_df.columns[:3]) != ['Lower', 'Upper', 'Prop']:
        raise Exception('Observed_bin_props needs to have column names Lower, Upper, and Prop')

    # Make sure the Observed_bin_props are numeric
    observed_bin_props_df['Prop'] = pd.to_numeric(observed_bin_props_df['Prop'], errors='raise')

    # Make sure there are not any NAs in the proportions
    if observed_bin_props_df['Prop'].isnull().any():
        raise Exception('Proportions in Observed_bin_props need to be numeric with no NA values')

    if not pd.api.types.is_numeric_dtype(observed_bin_props_df['Lower']) or not pd.api.types.is_numeric_dtype(
            observed_bin_props_df['Upper']):
        raise Exception('Observed_bin_props MAC bins need to be numeric')

    # Check the order of the MAC bins
    if not observed_bin_props_df['Upper'].is_monotonic_increasing or not observed_bin_props_df[
        'Lower'].is_monotonic_increasing:
        raise Exception('The MAC bins need to be ordered from smallest to largest')

    # Set the default value for p_rv to the sum of the rare variant bins
    p_rv = observed_bin_props_df['Prop'].sum()

    # c1: individual MACs to use in the function
    upper_last = observed_bin_props_df['Upper'].iloc[-1]
    c1 = np.arange(1, int(upper_last) + 1)  # creates a sequence from 1 to the last Upper value

    # define the least squares loss function
    def prob_leastsquares(inner_tune):
        local_alpha = inner_tune[0]
        beta_val = inner_tune[1]
        # Calculate b
        # calculate the function completely without b for each individual MAC
        individual_prop_no_b = 1 / ((c1 + beta_val) ** local_alpha)
        # solve for b
        b_val = p_rv / np.sum(individual_prop_no_b)
        # calculate the function with b for each individual MAC
        individual_prop = b_val * individual_prop_no_b

        total_error = 0
        # loop over the bins
        for i, row in observed_bin_props_df.iterrows():
            # Calculate expected (from the function)
            # Adjust index by subtracting 1 because Python arrays are 0-indexed
            lower_index = int(row['Lower']) -1
            upper_index = int(row['Upper'])
            E = np.sum(individual_prop[lower_index:upper_index])
            # record the observed proportion in the target data
            O = row['Prop']
            # calculate the squared error
            c = (E - O) ** 2
            # sum the squared error over all MAC bins
            total_error += c

        # The output here is the sum of the squared error over MAC bins
        return total_error

    # start with the function 1/x (alpha = 1, beta = 0)
    tune = np.array([1.0, 0.0])

    # Minimize with the SLSQP function using the starting values (tune), the least squares loss function (calc_prob_LS), and constraints (hin_tune)
    # Constraint: x[0] > 0
    cons = {'type': 'ineq', 'fun': lambda x: x[0]}
    S = minimize(prob_leastsquares, tune, constraints=cons, options={'disp': False, 'ftol': 0.0, 'maxiter': 25})

    # back calculate b after the parameters have been solved for
    alpha_opt = S.x[0]
    beta_opt = S.x[1]
    b = p_rv / np.sum(1 / ((c1 + beta_opt) ** alpha_opt))

    # Return the parameters alpha, beta, and b, as well as the proportions as calculated by the function
    alpha = alpha_opt
    beta = beta_opt
    b = b
    print(f"Calculated the following params from AFS target data. alpha: {alpha}, beta: {beta}, b: {b}")
    return alpha, beta, b


def afs(alpha, beta, b, macs):
    lowers = []
    uppers = []
    props = []

    for mac in macs:
        lowers.append(mac[0])
        uppers.append(mac[1])
        if sorted(lowers) != lowers or sorted(uppers) != uppers:
            raise Exception("Mac bins need to be in numeric order")

        fit = [b / ((beta + i + 1) ** alpha) for i in range(uppers[-1])]

        prop = sum(fit[i - 1] for i in range(mac[0], mac[1] + 1))
        props.append(prop)

    ret = [(lowers[i], uppers[i], props[i]) for i in range(len(props))]
    return ret

def main():
    args = get_args()
    n = int(args.n)
    macs = read_mac_bins(args.mac)
    reg_size = float(args.reg_size)

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

        weight = max(0.0, min(float(args.w), 2.0)) # clamp the weight to be in the range [0.2]

        num_variants = nvariants(n, omega, phi, reg_size, weight)
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
        print("\nCalculating synonymous values")
        syn_weight = max(0.0, min(float(args.w_syn), 2.0)) # clamp the weight to be in the range [0.2]
        alpha, beta, b = fit_afs(df_afs_syn)
        omega, phi = fit_nvars(df_nvar_syn)
        num_variants = nvariants(n, omega, phi, reg_size, syn_weight)
        rows = afs(alpha, beta, b, macs)
        syn_output_file = os.path.splitext(args.output)[0] + '_syn.txt'
        write_expected_variants(syn_output_file, num_variants, rows)

        # Now do it for Functional
        print("\nCalculating functional values")
        fun_weight = max(0.0, min(float(args.w_fun), 2.0)) # clamp the weight to be in the range [0.2]
        alpha, beta, b = fit_afs(df_afs_fun)
        omega, phi = fit_nvars(df_nvar_fun)
        num_variants = nvariants(n, omega, phi, reg_size, fun_weight)
        rows = afs(alpha, beta, b, macs)
        fun_output_file = os.path.splitext(args.output)[0] + '_fun.txt'
        write_expected_variants(fun_output_file, num_variants, rows)

if __name__ == '__main__':
    main()
