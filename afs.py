import argparse
import pandas as pd
import numpy as np
from scipy.optimize import minimize, show_options

DEFAULT_PARAMS = {
    'AFR': {"alpha":1.5883, "beta":-0.3083, "b":0.2872},
    'EAS': {"alpha":1.6656, "beta":-0.2951, "b":0.3137},
    'NFE': {"alpha":1.9470, "beta":-0.1180, "b":0.6676},
    'SAS': {"alpha":1.6977, "beta":-0.2273, "b":0.3564} 
}

def get_args():
    parser = argparse.ArgumentParser()

    parser.add_argument('--pop',
                        dest='pop',
                        help='Population to use default values for')
    parser.add_argument('--mac',
                        dest='macs',
                        required=True,
                        help='csv file of MAC bins')
    parser.add_argument('--alpha',
                        dest='alpha',
                        help='Provided alpha value')
    parser.add_argument('--beta',
                        dest='beta',
                        help='Provided beta value')
    parser.add_argument('-b',
                        dest='b',
                        help='Provided b value')
    parser.add_argument('-o',
                        dest='output',
                        help='Output file to be written')

    args = parser.parse_args()

    if args.pop and args.pop not in ["EAS", "AFR", "NFE", "SAS"]:
        raise Exception(f"{args.pop} is not a valid population")

    if (not args.alpha or not args.beta or not args.b) and not args.pop:
        raise Exception('Error: either a default population should be specified or all three parameters provided')

    return args

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

    # specify the function to define the constraints
    def hin_tune(x):
        h = np.zeros(1)
        h[0] = x[0]  # the first parameter (alpha) must be > 0
        return h

    # c1: individual MACs to use in the function
    upper_last = observed_bin_props_df['Upper'].iloc[-1]
    c1 = np.arange(1, int(upper_last) + 1)  # creates a sequence from 1 to the last Upper value

    # define the least squares loss function
    def prob_leastsquares(tune):
        alpha = tune[0]
        beta_val = tune[1]
        # Calculate b
        # calculate the function completely without b for each individual MAC
        individual_prop_no_b = 1 / ((c1 + beta_val) ** alpha)
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
        fit = []
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
    macs = []
    with open(args.macs) as macs_file:
        lines = macs_file.readlines()
        header = lines[0].split(',')
        if header[0].strip() != 'Lower' or header[1].strip() != 'Upper':
            raise Exception("Mac bins file needs to have column names Lower and Upper")

        for line in lines[1:]:
            l = line.strip().split(',')
            macs.append((int(l[0]), int(l[1])))

    if args.pop:
        alpha = DEFAULT_PARAMS[args.pop]['alpha']
        beta = DEFAULT_PARAMS[args.pop]['beta']
        b = DEFAULT_PARAMS[args.pop]['b']
    else:
        alpha = int(args.alpha)
        beta = int(args.beta)
        b = int(args.b)

    rows = afs(alpha, beta, b, macs)

    with open(args.output, 'w') as f:
        f.write("Lower,Upper,Prop\n")
        for i in range(len(rows)):
            f.writelines(f"{rows[i]},{rows[i]},{rows[i]}\n")

if __name__ == '__main__':
    main()
        
