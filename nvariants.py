import argparse


DEFAULT_PARAMS = {
    'AFR': {"phi":0.1576, "omega":0.6247},
    'EAS': {"phi":0.1191, "omega":0.6369},
    'NFE': {"phi":0.1073, "omega":0.6539},
    'SAS': {"phi":0.1249, "omega":0.6495} 
}

def get_args():
    parser = argparse.ArgumentParser()

    parser.add_argument('--pop',
                        dest='pop',
                        help='Population to use default values for')
    parser.add_argument('--phi',
                        dest='phi',
                        help='Provided phi value')
    parser.add_argument('--omega',
                        dest='alpha',
                        help='Provided omega value')
    parser.add_argument('-N',
                        dest='n',
                        required=True,
                        help='Provided N value')

    args = parser.parse_args()

    if args.pop and args.pop not in ["EAS", "AFR", "NFE", "SAS"]:
        raise Exception(f"{args.pop} is not a valid population")

    if (not args.phi or not args.omega) and not args.pop:
        raise Exception('Error: either a default population should be specified or all three parameters provided')

    return args

def fit_nvars_from_file(file):
    df = pd.read_csv(file, delimiter='\t')


import numpy as np
import pandas as pd
from scipy.optimize import minimize


def fit_nvars(Observed_variants_per_kb):
    # Check that each column is numeric
    if not np.issubdtype(Observed_variants_per_kb.iloc[:, 0].dtype, np.number) or not np.issubdtype(
            Observed_variants_per_kb.iloc[:, 1].dtype, np.number):
        raise ValueError('Columns need to be numeric')

    # Check that there are not any NA values
    if Observed_variants_per_kb.iloc[:, 1].isna().any():
        raise ValueError('Number of variants per Kb need to be numeric with no NA values')

    # Check that the sample sizes go from smallest to largest
    if not Observed_variants_per_kb.iloc[:, 0].is_monotonic_increasing:
        raise ValueError('The sample sizes need to be ordered from smallest to largest')

    # define the least squares loss function
    def leastsquares(tune):
        # calculate the expected number of variants (from the function)
        E = tune[0] * (Observed_variants_per_kb.iloc[:, 0] ** tune[1])
        sq_err = (E - Observed_variants_per_kb.iloc[:, 1]) ** 2  # calculate the squared error of expected - observed
        return np.sum(sq_err)  # return the squared error

    def hin_tune(x):  # constraints
        h = np.zeros(3)
        h[0] = x[0]  # phi greater than 0
        h[1] = x[1]  # omega greater than 0
        h[2] = 1 - x[1]  # omega less than 1
        return h

    # define the starting value for phi so the end of the function matches with omega = 0.45
    phi = Observed_variants_per_kb.iloc[Observed_variants_per_kb.iloc[:, 0].idxmax(), 1] / (
                Observed_variants_per_kb.iloc[Observed_variants_per_kb.iloc[:, 0].idxmax(), 0] ** 0.45)
    tune = np.array([phi, 0.45])  # specify the starting values

    # Use SLSQP to find phi and omega
    constraints = {'type': 'ineq', 'fun': hin_tune}
    re_LS = minimize(leastsquares, tune, constraints=constraints, options={'disp': False, 'ftol': 0.0, 'maxiter': 100})

    # If the original starting value resulted in a large loss (>1000), iterate over starting values
    if re_LS.fun > 1000:
        re_tab1 = []  # create to hold the new parameters
        for omega in np.arange(0.15, 0.66, 0.1):  # optimize with different values of omega

            # specify phi to fit the end of the function with the current value of omega
            phi = Observed_variants_per_kb.iloc[Observed_variants_per_kb.iloc[:, 0].idxmax(), 1] / (
                        Observed_variants_per_kb.iloc[Observed_variants_per_kb.iloc[:, 0].idxmax(), 0] ** omega)
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

def main():
    args = get_args()
    n = int(args.n)
    if args.pop:
        phi = DEFAULT_PARAMS[args.pop]["phi"]
        omega = DEFAULT_PARAMS[args.pop]["omega"]
    else:
        phi = float(args.phi)
        omega = float(args.omega)

    print(nvariants(n, omega, phi, 1))

if __name__ == '__main__':
        main()
