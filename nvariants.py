import argparse
import pandas as pd
from scipy.optimize import minimize


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

def fit(file):
    Observed_variants_per_kb = pd.read_csv(file)
    # Check that each column is numeric
    if not all(pd.api.types.is_numeric_dtype(Observed_variants_per_kb[col]) for col in Observed_variants_per_kb.columns):
        raise ValueError('Error: Columns need to be numeric')

    # Check that there are not any NA values
    if Observed_variants_per_kb.isnull().any().any():
        raise ValueError('Number of variants per Kb need to be numeric with no NA values')

    # Check that the sample sizes go from smallest to largest
    if not all(Observed_variants_per_kb.iloc[i, 0] <= Observed_variants_per_kb.iloc[i + 1, 0] for i in range(len(Observed_variants_per_kb) - 1)):
        raise ValueError('The sample sizes need to be ordered from smallest to largest')

    # Define the least squares loss function
    def leastsquares(tune):
        E = tune[0] * (Observed_variants_per_kb.iloc[:, 0] ** tune[1])
        sq_err = (E - Observed_variants_per_kb.iloc[:, 1]) ** 2
        d = sum(sq_err)
        return d

    # Define the constraints
    def hin_tune(x):
        return [x[0], x[1], 1 - x[1]]

    # Define the starting value for phi so the end of the function matches with omega = 0.45
    phi = Observed_variants_per_kb.iloc[-1, 1] / (Observed_variants_per_kb.iloc[-1, 0] ** 0.45)
    tune = [phi, 0.45]

    # Use SLSQP to find phi and omega
    res = minimize(leastsquares, tune, method='SLSQP', constraints={'type': 'ineq', 'fun': lambda x: hin_tune(x)})

    # If the original starting value resulted in a large loss (>1000), iterate over starting values
    if res.fun > 1000:
        re_tab1 = []
        for omega in [i / 10 for i in range(15, 66)]:
            phi = Observed_variants_per_kb.iloc[-1, 1] / (Observed_variants_per_kb.iloc[-1, 0] ** omega)
            tune = [phi, omega]
            res1 = minimize(leastsquares, tune, method='SLSQP', constraints={'type': 'ineq', 'fun': lambda x: hin_tune(x)})
            to_bind1 = [res1.x[0], res1.x[1], res1.fun]
            re_tab1.append(to_bind1)

        re_fin = min(re_tab1, key=lambda x: x[2])
        return {'phi': re_fin[0], 'omega': re_fin[1]}
    else:
        return {'phi': res.x[0], 'omega': res.x[1]}


def nvariants(n, omega, phi):
    return phi * (n**omega)

def main():
    args = get_args()
    n = int(args.n)
    if args.pop:
        phi = DEFAULT_PARAMS[args.pop]["phi"]
        omega = DEFAULT_PARAMS[args.pop]["omega"]
    else:
        phi = float(args.phi)
        omega = float(args.omega)

    print(nvariants(n, omega, phi))

if __name__ == '__main__':
    nvariants()
        
