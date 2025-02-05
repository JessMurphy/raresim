import argparse

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

        for j in range(len(uppers)):
            prop = sum(fit[i - 1] for i in range(lowers[j], uppers[j] + 1))
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
        
