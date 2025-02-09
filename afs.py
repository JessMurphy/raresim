import argparse
import pandas as pd
from scipy.optimize import minimize

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

def fit_afs(dataframe):

    # Check the column names of dataframe
    if not all(col in dataframe.columns for col in ['Lower', 'Upper', 'Prop']):
        raise ValueError('dataframe needs to have column names Lower, Upper, and Prop')

    # Make sure the dataframe are numeric
    dataframe['Prop'] = pd.to_numeric(dataframe['Prop'])

    # Make sure there are not any NA's in the proportions
    if dataframe['Prop'].isnull().any():
        raise ValueError('Proportions in dataframe need to be numeric with no NA values')

    if not all(pd.api.types.is_numeric_dtype(dataframe[col]) for col in ['Lower', 'Upper']):
        raise ValueError('dataframe MAC bins need to be numeric')

    # Check the order of the MAC bins
    if not list(dataframe['Upper']) == sorted(dataframe['Upper']):
        raise ValueError('The MAC bins need to be ordered from smallest to largest')

    # Set the default value for p_rv to the sum of the rare variant bins
    p_rv = dataframe['Prop'].sum()

    # Define the function to calculate the least squares loss
    def calc_prob_LS(tune):
        # Calculate b
        c1 = range(1, dataframe['Upper'].max() + 1)
        indivual_prop_no_b = 1 / ((c1 + tune[1]) ** tune[0])
        b = p_rv / sum(indivual_prop_no_b)

        # Calculate the function with b for each individual MAC
        indivual_prop = b * indivual_prop_no_b

        all = 0
        for i in range(len(dataframe)):
            E = sum(indivual_prop[dataframe['Lower'].iloc[i] - 1:dataframe['Upper'].iloc[i]])
            O = dataframe['Prop'].iloc[i]
            c = (E - O) ** 2
            all += c

        return all

    # Define the constraints
    def hin_tune(x):
        return x[0]  # The first parameter (alpha) must be > 0

    # Minimize with the SLSQP function
    tune = [1, 0]  # Start with the function 1/x (alpha = 1, beta = 0)
    res = minimize(calc_prob_LS, tune, method='SLSQP', constraints={'type': 'ineq', 'fun': hin_tune})

    # Back calculate b after the parameters have been solved for
    b = p_rv / sum(1 / ((range(1, dataframe['Upper'].max() + 1) + res.x[1]) ** res.x[0]))

    # Calculate the MAC bin proportions given the parameters
    def afs_internal(alpha, beta, b, mac_bins):
        c1 = range(1, mac_bins['Upper'].max() + 1)
        indivual_prop_no_b = 1 / ((c1 + beta) ** alpha)
        indivual_prop = b * indivual_prop_no_b
        proportions = []
        for i in range(len(mac_bins)):
            proportions.append(sum(indivual_prop[mac_bins['Lower'].iloc[i] - 1:mac_bins['Upper'].iloc[i]]))
        return proportions

    re = afs_internal(res.x[0], res.x[1], b, dataframe[['Lower', 'Upper']])

    # Return the parameters alpha, beta, and b, as well as the proportions as calculated by the function
    return res.x[0], res.x[1], b



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

    return [(lowers[i], uppers[i], props[i]) for i in range(len(props))]

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
        
