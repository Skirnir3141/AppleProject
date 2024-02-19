rm(list = ls())
setwd("/home/skirnir314/EECMSc/AppleProject")
library(readxl)
library(dplyr)
library(igraph)

# Import apple phenotype data and metadata.
apple.p.metadata <- read_excel("./phenomics/20201102_meta_data.xlsx")

# Aggregate plant id and cultivar to one field
apple.p.metadata$Cultivar.New <- apply(
  apple.p.metadata[, c("PLANTID", "Cultivar")],
  1,
  function(x) ifelse(x[1] == "NA", x[2], x[1])
)

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

# Filter for only malus domestica accessions with a year and drop irrelevant
# columns
apple.p.metadata <- apple.p.metadata[
  apple.p.metadata$species == "domestica" &
  !is.na(apple.p.metadata$Year),
  c("apple_id", "Cultivar.New", "Year", "Country", "Use", "type")]
names(apple.p.metadata)[
  names(apple.p.metadata) == "Cultivar.New"] <- "Cultivar"

# Add a column for century
apple.p.metadata$Century <- apply(
  apple.p.metadata[, c("Year")],
  1,
  function(x) round(x[1] / 100) * 100
)

# Add a column for time period
apple.p.metadata$Period <- apply(
  apple.p.metadata[, c("Year")],
  1,
  function(x) ifelse(
    x[1] <= 1850,
    1,
    ifelse(x[1] > 1850 & x[1] <= 1950, 2, 3))
)

# Convert to dataframe and create TSV files for filtering in Plink
apple.p.metadata <- as.data.frame(apple.p.metadata)
keep.p1 <- data.frame(
  PID = sort(apple.p.metadata[apple.p.metadata$Period == 1, c("apple_id")]))
keep.p1$IID <- sort(
  apple.p.metadata[apple.p.metadata$Period == 1, c("apple_id")])
keep.p2 <- data.frame(
  PID = sort(apple.p.metadata[apple.p.metadata$Period == 2, c("apple_id")]))
keep.p2$IID <- sort(
  apple.p.metadata[apple.p.metadata$Period == 2, c("apple_id")])
keep.p3 <- data.frame(
  PID = sort(apple.p.metadata[apple.p.metadata$Period == 3, c("apple_id")]))
keep.p3$IID <- sort(
  apple.p.metadata[apple.p.metadata$Period == 3, c("apple_id")])

# Write to TSV
write.table(keep.p1, "keep.p1.tsv", sep = "\t", row.names = FALSE)
write.table(keep.p2, "keep.p2.tsv", sep = "\t", row.names = FALSE)
write.table(keep.p3, "keep.p3.tsv", sep = "\t", row.names = FALSE)

# Prune SNPS data to remove SNPs in linkage disequilibrium
system("/usr/local/bin/Plink/plink --file ./snps/base_snp_data --indep-pairwise 10 3 0.5")

# Run IBD analysis in Plink for the three time periods
system("/usr/local/bin/Plink/plink --file ./snps/base_snp_data --keep keep.p1.tsv --extract plink.prune.in --allow-no-sex --geno .05 --mind .1 --maf .01 --genome --out ./snps/p1_ibd")
system("/usr/local/bin/Plink/plink --file ./snps/base_snp_data --keep keep.p2.tsv --extract plink.prune.in --allow-no-sex --geno .05 --mind .1 --maf .01 --genome --out ./snps/p2_ibd")
system("/usr/local/bin/Plink/plink --file ./snps/base_snp_data --keep keep.p3.tsv --extract plink.prune.in --allow-no-sex --geno .05 --mind .1 --maf .01 --genome --out ./snps/p3_ibd")

# Import ibd results
p1.ibd <- read.table("./snps/p1_ibd.genome", header = TRUE)
p2.ibd <- read.table("./snps/p2_ibd.genome", header = TRUE)
p3.ibd <- read.table("./snps/p3_ibd.genome", header = TRUE)

# Create mappings to substitute apple id for cultivar name
map1 <- apple.p.metadata[, c("apple_id", "Cultivar")]
map2 <- apple.p.metadata[, c("apple_id", "Cultivar")]
names(map1) <- c("IID1", "Cultivar")
names(map2) <- c("IID2", "Cultivar")

# Rename cultivars
p1.ibd <- dplyr::left_join(p1.ibd, map1, by = "IID1")
names(p1.ibd)[names(p1.ibd) == "Cultivar"] <- "Cultivar1"
p1.ibd <- dplyr::left_join(p1.ibd, map2, by = "IID2")
names(p1.ibd)[names(p1.ibd) == "Cultivar"] <- "Cultivar2"
p2.ibd <- dplyr::left_join(p2.ibd, map1, by = "IID1")
names(p2.ibd)[names(p2.ibd) == "Cultivar"] <- "Cultivar1"
p2.ibd <- dplyr::left_join(p2.ibd, map2, by = "IID2")
names(p2.ibd)[names(p2.ibd) == "Cultivar"] <- "Cultivar2"
p3.ibd <- dplyr::left_join(p3.ibd, map1, by = "IID1")
names(p3.ibd)[names(p3.ibd) == "Cultivar"] <- "Cultivar1"
p3.ibd <- dplyr::left_join(p3.ibd, map2, by = "IID2")
names(p3.ibd)[names(p3.ibd) == "Cultivar"] <- "Cultivar2"

#creat node and edge objects for network graphing
p1.nodes <- data.frame(Cultivar = unique(c(p1.ibd$Cultivar1, p1.ibd$Cultivar2)))
p1.edges <- p1.ibd[p1.ibd$PI_HAT > .125, c("Cultivar1", "Cultivar2", "PI_HAT")]
names(p1.edges) <- c("from", "to", "weight")
p2.nodes <- data.frame(Cultivar = unique(c(p2.ibd$Cultivar1, p2.ibd$Cultivar2)))
p2.edges <- p2.ibd[p2.ibd$PI_HAT > .125, c("Cultivar1", "Cultivar2", "PI_HAT")]
names(p2.edges) <- c("from", "to", "weight")
p3.nodes <- data.frame(Cultivar = unique(c(p3.ibd$Cultivar1, p3.ibd$Cultivar2)))
p3.edges <- p3.ibd[p3.ibd$PI_HAT > .125, c("Cultivar1", "Cultivar2", "PI_HAT")]
names(p3.edges) <- c("from", "to", "weight")

# Create graph and plot objects and plot in a grid using cowplot.
p1.graph <- as_tbl_graph(
  p1.edges,
  directed = FALSE,
  node_key = "Cultivar",
  nodes = p1.nodes)
p2.graph <- as_tbl_graph(
  p2.edges,
  directed = FALSE,
  node_key = "Cultivar",
  nodes = p2.nodes)
p3.graph <- as_tbl_graph(
  p3.edges,
  directed = FALSE,
  node_key = "Cultivar",
  nodes = p3.nodes)
p1.plot <- ggraph(p1.graph, layout = 'fr', weights = weight) + 
  geom_edge_link(aes(width = weight)) +
  scale_edge_width(range = c(0.1, 3)) +
  geom_node_point(aes(size = centrality_pagerank(), color = "green")) +
  theme(legend.position = "none")
p2.plot <- ggraph(p2.graph, layout = 'fr', weights = weight) + 
  geom_edge_link(aes(width = weight)) +
  scale_edge_width(range = c(0.1, 3)) +
  geom_node_point(aes(size = centrality_pagerank(), color = "green")) +
  theme(legend.position = "none")
p3.plot <- ggraph(p3.graph, layout = 'fr', weights = weight) + 
  geom_edge_link(aes(width = weight)) +
  scale_edge_width(range = c(0.1, 3)) +
  geom_node_point(aes(size = centrality_pagerank(), color = "green")) +
  theme(legend.position = "none")
cowplot::plot_grid(
  p1.plot,
  p2.plot,
  p3.plot,
  nrow = 1,
  ncol = 3,
  labels = c("<= 1850", "1851 - 1950", "1951-Present"))

p1.graph %>%
  activate(nodes) %>%
  as_tibble() %>%
  summarise(n = n())

nrow(p2.nodes) - nrow(as_tibble(activate(p2.graph, nodes)))

test <- ggraph(p2.graph, layout = 'fr', weights = weight) + 
  geom_edge_link(aes(width = weight)) +
  scale_edge_width(range = c(0.1, 3)) +
  geom_node_point(
    aes(size = centrality_pagerank(), color = "green"),
    show.legend = FALSE) +
  labs(edge_width = "Pi Hat") +
  theme()
test2 <- get_legend(test)

cowplot::plot_grid(
  p1.plot,
  p2.plot,
  p3.plot,
  test2,
  nrow = 2,
  ncol = 2,)

print(
  apple.p.metadata %>%
    group_by(Year) %>%
    summarise(n = n()) %>%
    arrange(Year),
  n = 200)


rm(list = ls())

# Run PCA analysis in Plink for the three time periods
system("/usr/local/bin/Plink/plink --file ./snps/abc_combined_maf001_sort_vineland_imputed --extract plink.prune.in --allow-no-sex --geno .05 --mind .1 --maf .01 --distance-matrix --out ./snps/pca-test")

# Load genetic distance data
dist <- read.table("./snps/pca-test.mdist", header = F)
PID <- data.frame(PID = read.table("./snps/pca-test.mdist.id", header = F)[, 1])
IID <- data.frame(IID = read.table("./snps/pca-test.mdist.id", header = F)[, 2])

# Perform PCA
mds <- cmdscale(dist, eig = T, 5)

# Extract eigen vectors and bind with IDs and period
eigenvec <- cbind(PID, IID, mds$points)
temp <- apple.p.metadata[, c("apple_id", "Period")]
names(temp) <- c("IID", "Period")
eigenvec <- dplyr::left_join(eigenvec, temp, by = "IID")

# Calculate proportion of variation captured by each eigenvector
eigen.perc <- round(((mds$eig) / sum(mds$eig)) * 100, 2)

eigenvec$P

# Graph
ggplot(data = eigenvec[!is.na(eigenvec$Period), ]) +
  geom_point(
    mapping = aes(x = `1`, y = `2`, color = as.factor(Period)),
    show.legend = TRUE) +
  scale_color_manual("Period", values = c("1" = "yellow", "2" = "red", "3" = "blue")) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  labs(
    title = "PCA of apple cultivar SNPs",
    x = paste0("Principle Component 1 (", eigen.perc[1], " %)"),
    y = paste0("Principle Component 2 (", eigen.perc[2], " %)")) +
  theme_minimal()


ggraph::ggraph(p2.graph, layout = 'fr', weights = weight) + 
  ggraph::geom_edge_link(aes(width = weight)) +
  ggraph::scale_edge_width(range = c(0.1, 3)) +
  ggraph::geom_node_point(
    aes(size = tidygraph::centrality_pagerank()),
    shape = 21,
    fill = "green",
    color = "black") +
  geom_node_text(
    aes(filter = tidygraph::centrality_pagerank() > .025, label = name),
    colour = 'white',
    size=8,
    family = "serif") +
  ggplot2::theme(legend.position = "none")


ggplot2::ggplot(data = eigenvec[!is.na(eigenvec$Period), ]) +
  ggplot2::geom_point(
    mapping = aes(x = `1`, y = `2`, color = as.factor(Period)),
    show.legend = TRUE) +
  ggplot2::scale_color_manual(
    "Period",
    values = c("1" = "yellow", "2" = "red", "3" = "blue")) +
  ggplot2::geom_hline(yintercept = 0, linetype = "dotted") +
  ggplot2::geom_vline(xintercept = 0, linetype = "dotted") +
  ggplot2::labs(
    title = "PCA of apple cultivar SNPs",
    x = paste0("Principle Component 1 (", eigen.perc[1], " %)"),
    y = paste0("Principle Component 2 (", eigen.perc[2], " %)")) +
  ggplot2::theme(legend.position = "top")


p1.ibd <- read.table("./Data/Snps/p1_ibd.genome", header = TRUE)
p2.ibd <- read.table("./Data/Snps/p2_ibd.genome", header = TRUE)
p3.ibd <- read.table("./Data/Snps/p3_ibd.genome", header = TRUE)


test <- cbind(
  size    = vcount(p1.graph),
  nedges  = ecount(p1.graph),
  density = edge_density(p1.graph),
  recip   = reciprocity(p1.graph),
  centr   = centr_betw(p1.graph)$centralization,
  pathLen = mean_distance(p1.graph)
)

test <- rbind(
  test,
  c(
    vcount(p3.graph),
    ecount(p3.graph),
    edge_density(p3.graph),
    reciprocity(p3.graph),
    centr_betw(p3.graph)$centralization,
    mean_distance(p3.graph)))

nodes
edges
nodes <- data.frame(id = c(1, 2, 3))
edges <- data.frame(to = c(1), from = c(2))
net <- graph_from_data_frame(d=edges, vertices=nodes, directed=FALSE)
plot(net)



p3.test <- bind_nodes(
  p3.graph,
  data.frame(
    name = p3.nodes[
      !(p3.nodes$Cultivar %in% unlist(as_tibble(activate(p3.graph, nodes)))),
      ]))

ggraph::ggraph(p1.test, layout = 'fr', weights = weight) + 
  ggraph::geom_edge_link(aes(width = weight)) +
  ggraph::scale_edge_width(range = c(0.1, 3)) +
  ggraph::geom_node_point(
    aes(size = tidygraph::centrality_pagerank()),
    shape = 21,
    fill = "green",
    color = "black") +
  ggplot2::theme(legend.position = "none")



cowplot::plot_grid(
  bar1,
  p1.plot,
  p2.plot,
  p3.plot,
  #cowplot::get_legend(dummy.plot),
  nrow = 2,
  ncol = 2,
  labels = c("", "P1", "P2", "P3"),
  label_size = 14,
  label_fontfamily = "serif",
  hjust = -.7,
  vjust = 2.0,
  scale = 1.02)

# Average degree of network
mean(degree(p3.graph))

# Clustering coefficient
transitivity(p1.graph)

# Size of the largest component
igraph::components(p2.graph, mode = "weak")

igraph::degree(p1.graph)
igraph::degree_distribution(p1.graph)

hist(igraph::degree(p1.graph))

mean(igraph::closeness(p1.graph))
mean(igraph::betweenness(p3.graph))


ggraph::ggraph(p3.graph, layout = "fr", weights = weight) + 
  ggplot2::ggtitle("P3: ") +
  ggraph::geom_edge_link(aes(width = weight)) +
  ggraph::scale_edge_width(range = c(0.1, 3)) +
  ggraph::geom_node_point(
    aes(size = tidygraph::centrality_pagerank()),
    shape = 21,
    fill = "green",
    color = "black") +
  ggplot2::theme(
    plot.title = element_text(hjust = .5),
    legend.position = "none")
