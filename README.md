# Biological QA example

## Background

This is a companion file to the tutorial 07_non_model_organism_exploratory_data_analysis from this [repos](https://github.com/Bn-Bioinformatics-Handelsbolag/RnaSeqTutorials).

In the tutorial we realised that one sample was undersampled. Unfortunately, that sample could not be recreated. Instead another replicate from an earlier time point (B3) is being used here, in an iteration of the Exploratory Data Analysis for that project.

## Set up

1. Create a new project in RStudio to clone the GitHub repository: [https://github.com/nicolasDelhomme/BioQAExample.git](https://github.com/nicolasDelhomme/BioQAExample.git)

2. Populate the UPSCb-common submodule (RStudio does not retrieve submodules by 
default). In the RStudio terminal do:

```{bash}
git submodule init
git submodule update --remote
```

3. Copy the Biological QA template from UPSCb-common/templates/R to src/R (for making the report prettier, also copy bulogo2.png,footer.html,header.html and style.css)

4. Open your copy of BiologicalQA.R

## Biological QA

1. The sample info is in `doc` and is a `tsv` file

2. The data is in data. However, I already imported the `salmon` data as a `tximport` object. Instead of reading the tx2gene file, and loading the quant.sf using tximport, simply `load(here("data/tximport.rda"))`