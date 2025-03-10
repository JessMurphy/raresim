def write_expected_variants(out_file, n, macs):
    with open(out_file, 'w') as output:
        output.writelines("Lower\tUpper\tExpected_var\n")
        for tup in macs:
            output.writelines(f"{tup[0]}\t{tup[1]}\t{tup[2]*n}\n")
