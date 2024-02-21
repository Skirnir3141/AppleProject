# To run this as is, you'll need a few things. First, of course the directory
# structure from https://github.com/Skirnir3141/AppleProject/ should be
# maintained. Second, you'll need to download data from 
# https://datadryad.org/stash/dataset/doi:10.5061/dryad.zkh1893cd and extract
# the abc_combined_maf001_sort_vineland_imputed .map and .ped files to the
# ./Data/Snps directory. These files are rather large, so I didn't want to
# upload them to github. Third, Plink v1.90 should be downloaded from the
# cog-genomics.org website (https://www.cog-genomics.org/plink2/) and installed
# in the /usr/local/bin directory on a computer running linux in a directory
# named Plink. Of course, you could install it elsewhere or on a different OS,
# but the directory path of the Plink calls will need to be updated.

# Prune SNPS data to remove SNPs in linkage disequilibrium
# Can't quite figure out how to get Plink to output files here to a given
# directory, so I'm just changing the wd back and forth as a hack.
setwd("/home/skirnir314/EECMSc/AppleProject/Data/Snps")
system("/usr/local/bin/Plink/plink --file abc_combined_maf001_sort_vineland_imputed --indep-pairwise 10 3 0.5")
setwd("/home/skirnir314/EECMSc/AppleProject")

# Run IBD analysis in Plink for the three time periods
system("/usr/local/bin/Plink/plink --file ./Data/Snps/abc_combined_maf001_sort_vineland_imputed --keep ./Data/keep.p1.tsv --extract ./Data/Snps/plink.prune.in --allow-no-sex --geno .05 --mind .1 --maf .01 --genome --out ./Data/Snps/p1_ibd")
system("/usr/local/bin/Plink/plink --file ./Data/Snps/abc_combined_maf001_sort_vineland_imputed --keep ./Data/keep.p2.tsv --extract ./Data/Snps/plink.prune.in --allow-no-sex --geno .05 --mind .1 --maf .01 --genome --out ./Data/Snps/p2_ibd")
system("/usr/local/bin/Plink/plink --file ./Data/Snps/abc_combined_maf001_sort_vineland_imputed --keep ./Data/keep.p3.tsv --extract ./Data/Snps/plink.prune.in --allow-no-sex --geno .05 --mind .1 --maf .01 --genome --out ./Data/Snps/p3_ibd")

# Get the genetic distance between accessions using Plink (this is for all
# periods combined)
system("/usr/local/bin/Plink/plink --file ./Data/Snps/abc_combined_maf001_sort_vineland_imputed --extract ./Data/Snps/plink.prune.in --allow-no-sex --geno .05 --mind .1 --maf .01 --distance-matrix --out ./Data/Snps/pca")

# Use Plink to calculate inbreeding coefficient per accession
system("/usr/local/bin/Plink/plink --file ./Data/Snps/abc_combined_maf001_sort_vineland_imputed --extract ./Data/Snps/plink.prune.in --allow-no-sex --geno .05 --mind .1 --maf .01 --het --out ./Data/Snps/ploidy-check")
