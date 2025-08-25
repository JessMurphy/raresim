# RAREsim2
Python interface for flexible simulation of rare-variant genetic data using real haplotypes


## Installation Steps
```
$ cd ~
$ conda activate
$ pip install --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple/ raresim==3.0.0
```

## Main Functions

### CALC
Calculate the expected number of variants per MAC bin using default population parameters, user-provided parameters, or target data.
```
usage: __main__.py calc [-h] --mac MAC -o OUTPUT -N N [--pop POP]
                        [--alpha ALPHA] [--beta BETA] [--omega OMEGA]
                        [--phi PHI] [-b B]
                        [--nvar_target_data NVAR_TARGET_DATA]
                        [--afs_target_data AFS_TARGET_DATA]
                        [--reg_size REG_SIZE] [-w W] [--w_fun W_FUN]
                        [--w_syn W_SYN]

options:
  -h, --help            show this help message and exit
  --mac MAC             MAC bin bounds (lower and upper allele counts) for the simulated sample size
  -o OUTPUT             Output file name
  -N N                  Simulated sample size
  --pop POP             Population (AFR, EAS, NFE, or SAS) to use default values for if not providing
                        alpha, beta, omega, phi, and b values or target data
  --alpha ALPHA         Shape parameter to estimate the expected AFS distribution (must be > 0)
  --beta BETA           Shape parameter to estimate the expected AFS distribution
  --omega OMEGA         Scaling parameter to estimate the expected number of variants per (Kb) for
                        sample size N (range of 0-1)
  --phi PHI             Shape parameter to estimate the expected number of variants per (Kb) for
                        sample size N (must be > 0)
  -b B                  Scale parameter to estimate the expected AFS distribution
  --nvar_target_data NVAR_TARGET_DATA
                        Target downsampling data with the number of variants per Kb to estimate the
                        expected number of variants per Kb for sample size N
  --afs_target_data AFS_TARGET_DATA
                        Target AFS data with the proportion of variants per MAC bin to estimate the
                        expected AFS distribution
  --reg_size REG_SIZE   Size of simulated genetic region in kilobases (Kb)
  -w W                  Weight to multiply the expected number of variants by in non-stratified
                        simulations (default value of 1)
  --w_fun W_FUN         Weight to multiply the expected number of functional variants by in
                        stratified simulations (default value of 1)
  --w_syn W_SYN         Weight to multiply the expected number of synonymous variants by in
                        stratified simulations (default value of 1)
```

* #### default population parameters
The expected number of functional and synonymous variants can be estimated using default paramters for the following populations: African (AFR), East Asian (EAS), Non-Finnish European (NFE), and South Asian (SAS).
```
$ python3 -m raresim calc \
    --mac data/mac_bins.csv \
    -o <output file> \
    -N 15000 \
    --pop EAS \
    --reg_size 19.029
```

* #### target data
The user can also use their own target data - this is necessary to calculate the expected number of functional and/or synonymous variants for stratified simulations. Note, the simulation paramters are output if the user wants to use them instead of target data for future simulations.
(add target data to repo)
```
$ python3 -m raresim calc
    --mac data/mac_bins.csv
    -o <output file>
    -N 15000
   --nvar_target_data chr19_block37_NFE_nvar_target_data.txt
   --afs_target_data chr19_block37_NFE_AFS_target_data.txt
   --reg_size 19.029
```

* #### user-provided parameters
If parameters are known from previous simulations, the user can provide those instead of having to provide and fit target data.
```
$ python3 -m raresim calc
    --mac data/mac_bins.csv
    -o <output file>
    -N 15000
    --alpha 1.5
    --beta -.25
    -b .25
    --omega .15
    --phi .65
    --reg_size 19.029
```

### SIM
Simulate new allele frequencies given input haplotypes, a legend file, and the expected number of variants for the simulated sample size. A list of pruned variants (.legend-pruned-variants) is also output.
```
usage: __main__.py sim [-h] -m SPARSE_MATRIX [-b EXP_BINS]
                       [--functional_bins EXP_FUN_BINS]
                       [--synonymous_bins EXP_SYN_BINS] -l INPUT_LEGEND
                       [-L OUTPUT_LEGEND] -H OUTPUT_HAP
                       [--f_only FUN_BINS_ONLY] [--s_only SYN_BINS_ONLY] [-z]
                       [-prob] [--small_sample] [--keep_protected]
                       [--stop_threshold STOP_THRESHOLD]
                       [--activation_threshold ACTIVATION_THRESHOLD]
                       [--verbose]

options:
  -h, --help            show this help message and exit
  -m SPARSE_MATRIX      Input haplotype file (can be a .haps, .sm, or .gz file)
  -b EXP_BINS           Expected number of functional and synonymous variants per MAC bin
  --functional_bins EXP_FUN_BINS
                        Expected number of variants per MAC bin for functional variants (must be used
                        with --synonymous_bins) 
  --synonymous_bins EXP_SYN_BINS
                        Expected number of variants per MAC bin for synonymous variants (must be used
                        with --functional_bins) 
  -l INPUT_LEGEND       Input legend file
  -L OUTPUT_LEGEND      Output legend file (only required when using -z)
  -H OUTPUT_HAP         Output compressed haplotype file
  --f_only FUN_BINS_ONLY
                        Expected number of variants per MAC bin for only functional variants
  --s_only SYN_BINS_ONLY
                        Expected number of variants per MAC bin for only synonymous variants
  -z                    Monomorphic and pruned variants (rows of zeros) are removed from the output
                        haplotype file
  -prob                 Variants are pruned allele by allele given a probability of removal in the
                        legend file
  --small_sample        Overrides error to allow for simulation of small sample sizes (<10,000
                        haplotypes)
  --keep_protected      Variants designated with a 1 in the protected column of the legend file will
                        not be pruned
  --stop_threshold STOP_THRESHOLD
                        Percentage threshold for stopping the pruning process (0-100). Prevents the
                        number of variants from falling below the specified percentage of the expected
                        count for any given MAC bin during pruning (default value of 20)
  --activation_threshold ACTIVATION_THRESHOLD
                        Percentage threshold for activating the pruning process (0-100). Requires that
                        the actual number of variants for a MAC bin must be more than the given
                        percentage different from the expected number to activate pruning on the bin
                        (default value of 10)
  --verbose             
```

```
$ python3 -m raresim sim \
    -m Simulated_80k_9.controls.haps.gz \
    -b data/Expected_variants_per_bin_80k.txt \
    -l data/Simulated_80k.legend \
    -L new.legend \
    -H new.hap.gz

Input allele frequency distribution:
(1, 1, 20.0) 9
(2, 2, 10.0) 5
(3, 5, 5.0) 6
(6, 10, 5.0) 7
(11, 20, 1.0) 11
(21, 1000, 0.0) 48

New allele frequency distribution:
(1, 1, 20.0) 15
(2, 2, 10.0) 11
(3, 5, 5.0) 6
(6, 10, 5.0) 3
(11, 20, 1.0) 1
(21, 1000, 0.0) 0

Writing new variant legend

Writing new haplotype file............
```

* #### stratified (functional/synonymous) pruning
To perform stratified simulations where functional and synonymous variants are pruned separately:
1. add a column to the legend file (`-l`) named "fun", where functional variants have the value "fun" and synonymous variants have the value "syn"
2. provide separate MAC bin files with the expected number of variants per bin for functional (`--functional_bins`) and synonymous (`--synonymous_bins`) variants
```
$ python3 -m raresim sim \
    -m chr19.block37.NFE.sim100.stratified.haps.gz \
    --functional_bins data/Expected_variants_functional.txt \
    --synonymous_bins data/Expected_variants_synonymous.txt \
    -l data/chr19.block37.NFE.sim100.stratified.legend \
    -L new.legend \
    -H new.hap.gz

Input allele frequency distribution:
Functional
[1,1]   610.213692400324    686
[2,2]   199.745137641156    351
[3,5]   185.434393821117    598
[6,10]  73.1664075520905    472
[11,20] 37.132127271035 432
[21,220]    34.4401706091422    768
[221,440]   1.98761248740743    10
[441, ]     30

Synonymous
[1,1]   215.389082675548    276
[2,2]   73.1166493377018    140
[3,5]   73.6972836211026    240
[6,10]  33.4315406970657    181
[11,20] 19.1432926816897    181
[21,220]    20.2848171294807    331
[221,440]   1.38678884898772    11
[441, ]     20

New allele frequency distribution:
Functional
[1,1]   610.213692400324    607
[2,2]   199.745137641156    217
[3,5]   185.434393821117    178
[6,10]  73.1664075520905    82
[11,20] 37.132127271035 40
[21,220]    34.4401706091422    41
[221,440]   1.98761248740743    1
[441, ]     30

Synonymous
[1,1]   215.389082675548    220
[2,2]   73.1166493377018    66
[3,5]   73.6972836211026    63
[6,10]  33.4315406970657    31
[11,20] 19.1432926816897    20
[21,220]    20.2848171294807    20
[221,440]   1.38678884898772    1
[441, ]     20

Writing new variant legend

Writing new haplotype file...........
```

* #### only functional/synonymous variants
To pruned only functional or only synonymous variants:
1. add a column to the legend file (`-l`) named "fun", where functional variants have the value "fun" and synonymous variants have the value "syn"
2. provide a MAC bin file with the expected number of variants per bin for only functional (`--f_only`) or only synonymous (`--s_only`) variants
```
$ python3 -m raresim sim \
    -m chr19.block37.NFE.sim100.stratified.haps.gz \
    --f_only data/Expected_variants_functional.txt \
    -l data/chr19.block37.NFE.sim100.stratified.legend \
    -L new.legend \
    -H new.hap.gz
```

* #### given probabilities
To prune variants using known or given probabilities, add a column to the legend file (`-l`) named "prob". A random number between 0 and 1 is generated for each variant, and if the number is greater than the probability, the variant is removed from the data.
```
$ python3 -m raresim sim \
    -m data/ProbExample.haps.gz \
    -H new.hap.gz \
    -l data/ProbExample.probs.legend \
    -prob
```

* #### protected status
To exclude protected variants from the pruning process, add a column to the legend file (`-l`) named "protected". Any row with a 0 in this column will be eligible for pruning while any row with a 1 will still be counted but will not be eligible for pruning.
```
$ python3 -m raresim sim \
    -m data/ProbExample.haps.gz \
    -H new.hap.gz \
    -l data/ProtectiveExample.legend \
    --keep_protected \
    -b data/fonlyBins.txt \
    --small_sample \
    -L out.test
```

### EXTRACT
Randomly extract a subset of haplotypes (.haps-sample.gz) and output the remaining haplotypes separately (.haps-remainder.gz).
```
options:
  -h, --help            show this help message and exit
  -i INPUT_FILE         Input haplotype file
  -o OUTPUT_FILE        Output haplotype file name
  -s SEED, --seed SEED  Optional seed for reproducibility
  -n NUM                Number of haplotypes to extract
```

```
$ python3 -m raresim extract \
    -i lib/raresim/test/data/Simulated_80k_9.controls.haps.gz \
    -o extracted_hap_subset.haps.gz \
    -n 20 \
    --seed 123
```
