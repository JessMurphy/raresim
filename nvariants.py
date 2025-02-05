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
        
