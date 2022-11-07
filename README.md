# raresim
Python Interface for Scalable rare-variant simulations.


## To Install:  
  $ Clone git ..  <br/>
  $ cd Path to raresim/    <br/>
  $ python3 setup.py install <br/>

## Usage:

### Covert haplotype files to a sparse matrix

```
usage: convert.py [-h] -i INPUT_FILE -o OUTPUT_FILE

optional arguments:
  -h, --help      show this help message and exit
  -i INPUT_FILE   Input haplotype file path
  -o OUTPUT_FILE  Ouput sparse matrix path
```

```
$ python convert.py \
    -i lib/raresim/test/data/Simulated_80k_9.controls.haps.gz \
    -o Simulated_80k_9.controls.haps.gz.sm
```

### Extract haplotype subset

```
usage extract.py [-h] -i INPUT_FILE -o OUTPUT_FILE -n NUMBER -seed SEED
  
optional arguments:
  -h, --help        show this help message and exit
  -i INPUT_FILE     Input haplotype file path
  -o OUTPUT_FILE    Output haplotype subset file path
  -n NUM            Size of haplotype subset
  --seed SEED       Random seed for replication of random sample
```

```
$ python extract.py \
    -i lib/raresim/test/data/Simulated_80k_9.controls.haps \
    -o extracted_hap_subset.haps \
    -n 20 \
    --seed 123
```

### Simulate new allele frequencies

```
usage: sim.py [-h] -m SPARSE_MATRIX -b EXP_BINS -l INPUT_LEGEND -L
              OUTPUT_LEGEND -H OUTPUT_HAP

optional arguments:
 -h, --help        show this help message and exit
 -m SPARSE_MATRIX  Input sparse matrix path
 -b EXP_BINS       Input expected bin sizes
 -l INPUT_LEGEND   Input variant site legend
 -L OUTPUT_LEGEND  Output variant site legend
 -H OUTPUT_HAP     Output compress hap file
```

```
$ python sim.py \
    -m Simulated_80k_9.controls.haps.gz.sm \
    -b lib/raresim/test/data/Expected_variants_per_bin_80k.txt \
    -l lib/raresim/test/data/Simulated_80k.legend \
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

### Simulations that consider variant affect (functional/synonymous)

Simulations can independently be considered variants by their impact if 
1. the legend file (`-l` option) contains a column labeled fun where functional
variants have the value fun, and synonymous variants have the value syn.
2. separate expected bin size files are given for functional
(`--functional_bins`) and synonymous (`--synonymous_bins`) variants

```
$ python convert.py \
    -i lib/raresim/test/data/chr19.block37.NFE.sim100.stratified.haps.gz \
    -o chr19.block37.NFE.sim100.stratified.haps.gz.sm

$ python sim.py \
    -m chr19.block37.NFE.sim100.stratified.haps.gz.sm \
    --functional_bins lib/raresim/test/data/Expected_variants_functional.txt \
    --synonymous_bins lib/raresim/test/data/Expected_variants_synonymous.txt \
    -l lib/raresim/test/data/chr19.block37.NFE.sim100.stratified.legend \
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

### Prune only one type of variant
```
$ python convert.py \
    -i lib/raresim/test/data/chr19.block37.NFE.sim100.stratified.haps.gz \
    -o chr19.block37.NFE.sim100.stratified.haps.gz.sm

$ python sim.py
    -m chr19.block37.NFE.sim100.stratified.haps.gz.sm \
    --f_only lib/raresim/test/data/Expected_variants_functional.txt \
    -l lib/raresim/test/data/chr19.block37.NFE.sim100.stratified.legend \
    -L new.legend \
    -H new.hap.gz
```

### Prune by given probabilities

Rows can be pruned allele by allele using probabilities given in the legend file.

```
$ python sim.py
    -m testData/ProbExample.haps.sm \
    -H new.hap.gz \
    -l testData/ProbExample.probs.legend \
    -prob
```

Python code example: <br/>
```python
from rareSim import sparse

def main():
    hap = "lib/raresim/test/data/Simulated_80k_9.controls.haps.gz"
    M = sparse(hap)
    print(M.num_rows(), M.num_cols())
    for row in range( M.num_rows()):
        alts = []
        for i in range(M.row_num(row)):
            alts.append(M.get(row,i))
        print(row, alts)

    M.write('out.haps.dat')

if __name__ == '__main__': main()
  ```


## Running C code

### Build

```
cd lib/raresim/src/
make
```

### Run

```
gunzip test/data/Simulated_80k_9.controls.haps.gz
./read \
    -i ../test/data/Simulated_80k_9.controls.haps \
    -o Simulated_80k_9.controls.haps.dat \
```
