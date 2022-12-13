# Differential Expression Pipeline

New snakemake-based version of IMGAG pipeline, developed by Jakob

## Install packages from Renv

Before running the pipeline you need to install the R packages. Use this approach:

- Run snakemake pipeline with conda-environments. You will get an error in `rule diff_expr`. From the error log, take the location of the R conda environment
- Activate the R conda environment
-  Start R and install the missing packages
```
R
renv::install()
```

Now the dependencies should be installed, you need to repeat this step when the pipeline is run in a new conda environment.


## How to use:


### 1) Setup project

####Export metadata from NGSD into the project folder

```bash
python3 expression-analysis-smk/scripts/ngsd-project-metadata.py | cut -f1,5,6 > metadata.tsv
```

If folder name is not the project name, it can be given as parameter

```bash
python3 expression-analysis-smk/scripts/ngsd-project-metadata.py  '220100_My_Project` | cut -f1,5,6 > metadata.tsv
```

Update the group column in the `metadata.tsv` file and add additional columns for additional factors used in the analysis (age, sex, treatment ..)


#### Setup the config file

Analysis parameters are set in the `config.yml` file. Look in the `templates/config.yml` for default values and available options. Important settings are:

- Reference genome
- File names
- Filter thresholds
- Model and factor names
- Highlighted genes
- Colouring options


#### Setup the contrasts file

The contrasts file defines the group comparisons calculate in the DGE analysis. Its a tab separated text file with a row for each contrast and the following columns:

- `name`  Appears as name of the sheet in the Excel results
- `title` Title of Comparison, can be multiple words
- `description` Subtitle with description of the comparison, free text
- `coefficients` Defines which groups should be compared. Format: `groupwt-groupko` The group of the first row is implied and can be omitted. Example: if the group in first row is `wt`, it would enough to write `groupko` in the coefficients field.


### Run the pipeline

```
snakemake --profile local -s /mnt/users/ahgrosc1/pipeline/expression-analysis-smk/workflow/Snakefile --cores 5
```


