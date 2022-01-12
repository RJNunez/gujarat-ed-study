# All-cause excess mortality in the State of Gujarat, India, during the COVID-19 pandemic (March 2020-April2021)

This repository holds all the code and data necessary to replicate the analyses in the study _All-cause excess mortality in the State of Gujarat, India, during the COVID-19 pandemic (March 2020-April2021)_.

### Folders
- **code**: Contains all the code to replicate the analyses 
  - `code-for-publication`: This is the only script needed to replicate the analysis. It relies on the `init.R` script and the already processed data in the `data` folder.
  - `init`: Contains all the libraries, code for functions and the visualization template used for the figures.
  - `wrangling`: Contains code for the preprocessing pipeline used to wrangle the csv files.

- **data**: Contains the raw and processed data for the analyses
  - **csvs**: Municipality-specific csv files
  - **rdas**: Processed data. This object is used in the `code-for-publication` script.

- **figures**: Contains the pdf files of all figures in the paper
