def loadBins(filename: str) -> list:
    bins = []
    with open(filename) as f:
        # Skip the first line since it is the header
        f.readline()

        line = f.readline()
        while line and line.strip() != "" and line.strip() != "\n":
            row = line.rstrip().split('\t')
            bins.append((int(row[0]), int(row[1]), float(row[2])))
            line = f.readline()
    return bins
