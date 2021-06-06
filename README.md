# Biological QA example

## Set up

1. Create a new project in RStudio to clone the GitHub repository: [https://github.com/nicolasDelhomme/BioQAExample.git](https://github.com/nicolasDelhomme/BioQAExample.git)

2. Create a link to the data and to the persisent storage

```{bash}
ln -s ../raw_data data
ln -s ../persistent analysis
```

3. Populate the UPSCb-common submodule (RStudio does not retrieve submodules by 
default). In the RStudio terminal do:

```{bash}
git submodule init
git submodule update --remote
```

4. Copy the Biological QA template from UPSCb-common/templates/R to src/R

5. Open your copy of BiologicalQA.R

## Lecture

## Homework