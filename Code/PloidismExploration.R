# Clear environment and set working directory
rm(list = ls())
setwd("/home/skirnir314/EECMSc/AppleProject")

# Load libraries
library(dplyr)
library(ggplot2)

# Import apple access phenotype metadata.
apple.p.metadata <- readxl::read_excel(
  "./Data/Phenomics/20201102_meta_data.xlsx")

# Combine discovery and release year.
apple.p.metadata$Year <- apply(
  apple.p.metadata[, c("Discovered/described/cultivated", "Release Year")],
  1,
  function(x) ifelse(is.na(x[1]) & is.na(x[2]), NA, min(x, na.rm = TRUE)))
apple.p.metadata$Year <- as.numeric(apple.p.metadata$Year)

# Aggregate plant id and cultivar to one field
apple.p.metadata$Cultivar.New <- apply(
  apple.p.metadata[, c("PLANTID", "Cultivar")],
  1,
  function(x) ifelse(x[1] == "NA", x[2], x[1])
)

# Filter for only malus domestica accessions with a year, drop irrelevant
# columns, and rename cultivar variable.
apple.p.metadata <- apple.p.metadata[
  apple.p.metadata$species == "domestica" &
    !is.na(apple.p.metadata$Year),
  c("apple_id", "Cultivar.New", "Year", "Country", "Use", "type")]
names(apple.p.metadata)[
  names(apple.p.metadata) == "Cultivar.New"] <- "Cultivar"

# Import pomiferous ploidism search methods, name and arrange columns
ploid.1 <- read.csv("./Data/pomiferous_pull_1.csv", header = FALSE)
ploid.2 <- read.csv("./Data/pomiferous_pull_2.csv", header = FALSE)
names(ploid.1) <- c("Cultivar", "Ploidism.1")
names(ploid.2) <- c("Cultivar", "Ploidism.2")
ploid.1 <- dplyr::arrange(ploid.1, "Cultivar")
ploid.2 <- dplyr::arrange(ploid.2, "Cultivar")

# Rank ploidism results and filter out cultivars where we have duplicates
ploid.1 <- ploid.1 %>%
  dplyr::group_by(Cultivar) %>%
  dplyr::arrange(desc(length(Cultivar))) %>%
  dplyr::mutate(rank = row_number())
ploid.2 <- ploid.2 %>%
  dplyr::group_by(Cultivar) %>%
  dplyr::arrange(desc(length(Cultivar))) %>%
  dplyr::mutate(rank = row_number())
ploid.1 <- ploid.1[ploid.1$rank == 1, c("Cultivar", "Ploidism.1")]
ploid.2 <- ploid.2[ploid.2$rank == 1, c("Cultivar", "Ploidism.2")]

# Join ploidism results with metadata
apple.p.metadata <- dplyr::left_join(apple.p.metadata, ploid.1, by = "Cultivar")
apple.p.metadata <- dplyr::left_join(apple.p.metadata, ploid.2, by = "Cultivar")

# Make ploidism string lower case for simpler string matching
apple.p.metadata$Ploidism.1 <- tolower(apple.p.metadata$Ploidism.1)
apple.p.metadata$Ploidism.2 <- tolower(apple.p.metadata$Ploidism.2)

# Create flags for references to typical apple ploidism: diploid, triploid, and
# tetraploid
apple.p.metadata$Dip.Flag.1 <- apply(
  apple.p.metadata[, c("Ploidism.1")],
  1,
  function(x) ifelse(grepl("diploid", x[1], fixed = TRUE), TRUE, FALSE)
)
apple.p.metadata$Trip.Flag.1 <- apply(
  apple.p.metadata[, c("Ploidism.1")],
  1,
  function(x) ifelse(grepl("triploid", x[1], fixed = TRUE), TRUE, FALSE)
)
apple.p.metadata$Tetra.Flag.1 <- apply(
  apple.p.metadata[, c("Ploidism.1")],
  1,
  function(x) ifelse(grepl("tetraploid", x[1], fixed = TRUE), TRUE, FALSE)
)
apple.p.metadata$Dip.Flag.2 <- apply(
  apple.p.metadata[, c("Ploidism.2")],
  1,
  function(x) ifelse(grepl("diploid", x[1], fixed = TRUE), TRUE, FALSE)
)
apple.p.metadata$Trip.Flag.2 <- apply(
  apple.p.metadata[, c("Ploidism.2")],
  1,
  function(x) ifelse(grepl("triploid", x[1], fixed = TRUE), TRUE, FALSE)
)
apple.p.metadata$Tetra.Flag.2 <- apply(
  apple.p.metadata[, c("Ploidism.2")],
  1,
  function(x) ifelse(grepl("tetraploid", x[1], fixed = TRUE), TRUE, FALSE)
)

# Confirm that multiple flags don't get return for either ploidism check. I.e.,
# this tells me that no text field contains multiple ploidism references, so
# no need to scrutinize any example.
max(apply(apple.p.metadata[, c(9:11)], 1, sum))
max(apply(apple.p.metadata[, c(12:14)], 1, sum))

# Looks like there are only 7 cases where the two search methods return
# different results.
apple.p.metadata[
  apple.p.metadata$Dip.Flag.1 != apple.p.metadata$Dip.Flag.2 |
    apple.p.metadata$Trip.Flag.1 != apple.p.metadata$Trip.Flag.2 |
    apple.p.metadata$Tetra.Flag.1 != apple.p.metadata$Tetra.Flag.2,
  c(2, 9:14)]

# Create column for ploidism
apple.p.metadata$Ploidy <- apply(
  apple.p.metadata[, c(9:14)],
  1,
  function(x) ifelse(
    x[1] == TRUE | x[4] == TRUE,
    "Diploid",
    ifelse(
      x[2] == TRUE | x[5] == TRUE,
      "Triploid",
      ifelse(
        x[3] == TRUE | x[6] == TRUE,
        "Tetraploid",
        NA)))
)

# Remove unnecessary columns
apple.p.metadata <- apple.p.metadata[, c(1:6, 15)]

# Join inbreeding coefficient to accession data
ploidy.check <- read.table("./Data/Snps/ploidy-check.het", header = TRUE)
apple.p.metadata <- dplyr::left_join(
  apple.p.metadata,
  ploidy.check[, c("IID", "F")],
  by = join_by(apple_id == IID))

# As expected, triploid accessions show up as having abnormally high
# heterozygosity as measured by the inbreeding coefficient.
boxplot(F ~ Ploidy, data = apple.p.metadata)

# We can also overlay these as histograms
ggplot2::ggplot(apple.p.metadata, aes(x = F)) +
  ggplot2::geom_histogram(
    data = subset(apple.p.metadata, Ploidy == "Diploid"),
    fill = "green") +
  ggplot2::geom_histogram(
    data = subset(apple.p.metadata, Ploidy == "Triploid"),
    fill = "blue")

# Let's take a look at some abnormally high/low F values for Diploids and
# Triploids. Cross checking some in pomiferous, they're accurate.
apple.p.metadata[
  apple.p.metadata$Ploidy %in% c("Triploid") &
    apple.p.metadata$F > 0, ]
apple.p.metadata[
  apple.p.metadata$Ploidy %in% c("Diploid") &
    apple.p.metadata$F < 0, ]

# We're missing Ploidism for most accessions, so it wouldn't make sense to
# filter by Ploidy search data from pomiferous.
apple.p.metadata %>%
  dplyr::group_by(Ploidy) %>%
  dplyr::summarise(n = n())

# Looking at all accessions, vast majority have an F-statistic greater than zero
# So, it makes sense to filter out accessions with an F-statistic less than zero
# as this likely indicates non-diploid accessions. 
ggplot2::ggplot(apple.p.metadata, aes(x = F)) +
  ggplot2::geom_histogram(fill = "green")

# write a csv of apple ids to remove from analysis due to low F-statistic
write.csv(
  dplyr::filter(apple.p.metadata, F <0 & !is.na(apple_id))[, c("apple_id")],
  "./Data/polyploid.ids.csv",
  row.names = FALSE)