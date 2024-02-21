This project uses SNP data from Canada's Apple Biodiversity Collection (ABC) to examine whether the genetic diversity of apple cultivars has increased, decreased, or stayed the same over the last 200 years.

The .Rmd file can be executed as is (you'll just have to update it to set the working directory to wherever you've saved the AppleProject folder).  All necessary data files are available in the Data directory. However, some of those files are derived using secondary scripts that are available in the Code directory. If you'd like to execute those scripts, you can do so. But, for the PlinkRuns.R script you will need to download and extract some additional data files. Directions are at the top of that script.

Code: Contains .R files for exploration and some data derivation. AppleExploration.R was used for scratch code during initial data exploration. PloidismExploration.R was used to derive a data file used to filter out likely polyploid cultivars. PlinkRuns.R was used to analyze SNP data using Plink.

Data: Contains two directories. Snps contains SNP data and output of Plink runs. Phenomics contains metadata on ABC accessions.

Pngs: Contains full size .png images for plots from the most recent .Rmd run.

Please reach out to mgj23@ic.ac.uk with any questions.
