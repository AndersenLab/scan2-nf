# scan2-nf
Run scan2 (R/qtl) and permutations on cluster.

### Requirements
The scan2-nf pipeline requires installation of several R packages including [linkage mapping]("https://github.com/AndersenLab/linkagemapping"), `dplyr`, `readr`, and `qtl`.

```
# Install linkagemapping package
devtools::install_github("AndersenLab/linkagemapping")
```

### Input
For the Andersen Lab, the input file is generally going to be the output from the `easysorter` pipeline. For now, only one trait is accepted. In general, you need a tsv dataframe of the following format:

| trait | strain | phenotype |
| --- | --- | --- |
| docetaxel.mean.TOF | QX318 | -31.25 |
| docetaxel.mean.TOF | QX299 | 5.154 |
| ... | ... | ... |
| docetaxel.mean.TOF | QX306 | 2.98 |

### Usage
```
nextflow run AndersenLab/scan2-nf --in <input.tsv>
```

#### Optional parameters
| param | default | optional | explanation |
| --- | --- | --- | --- |
| `--cross` | 'N2xCB4856cross_full' | Any cross object in linkagemapping package | default cross object is generated from whole genome data for RIAILs. Other common option is `N2xCB4856cross` which is generated from ~1400 SNPs. |
| `--set` | 2 | 1, 2, 3 | Defines the set to use for RIAILs for the cross object. Default is 2. |
| `--nperm` | 1000 | Any number | Number of permutations to run for the scan2 to determine significance. Default is 1000. |
| `--out` | Directory you are running the code in | Any directory | Determines the location of the output files. | 

### Complex Usage
```
nextflow run AndersenLab/scan2-nf --in <input.tsv> --cross N2xCB4856cross --set 1 --nperm 500
```
