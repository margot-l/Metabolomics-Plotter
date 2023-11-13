# Metabolomics-Plotter

Start by opening base_analysis.Rmd and follow the instructions therein.

## Packages

## Data Requirements

### Sample Sheet
'sample_sheet' must include the following columns:
- 'sampleName'
- 'Replicate'
- Each variable of interest (ie a drug treatment, mutant) must have one column
- Variables cannot include '-' ie 'kynu-1' should be written as 'kynu1'
- For concentrations, '0var' should be the buffer only and each concentration that follows should have the pattern of '1col1', '2col1', etc with 'col1' as the column name ie a 'KCN' column would have '0KCN' and '1KCN' while a 'Rot' column would have '0Rot' and '1Rot'. The number at the front will reflect how many '+'s will be used at the bottom of the plot. Whole numbers only.
- A 'Time' column with a number of hours the experiment was run for.
- A 'Contents' column may be helpful ie the row where 'KCN''s value is '0', 'Contents' may say 'Buffer alone'; this is purely for human aid and will not be used by the program.

### Compounds
'compounds' must include the following columns:
- 'featureName'
- 'featureClass' -> whether a feature is a 'met' or a 'ref' (an internal reference).
- 'KEGG_ID' -> the KEGG ID of each metabolite, required for pathway information.

### LC-MS Data
'objectDataFrame' should be the result of a chromeXtractor run using samples that match the 'sample_sheet' and 'compounds'. The following columns should be generated automatically but be sure to double check.

'objectDataFrame' must include the following columns:
- 'sampleName': must match the 'sampleName' from the 'sample_sheet'
- 'sampleClass': whether the sample is a 'blank' (ie a run black), a 'filter' (ie a rep blank) or a 'sample'
- 'featureName'
- 'featureClass'
- A value for each peak ie 'integratedIntensity' and/or 'maxIntensity' -> if you would like to use a 'cut_height' in dataCleanup, a 'maxIntensity' column is required.
