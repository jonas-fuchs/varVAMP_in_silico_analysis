# libraries
list.of.packages <- c("ggplot2", "ggtree", "treeio", "data.table", "dplyr", "optparse", "viridis")

invisible(lapply(list.of.packages, function(sm) suppressMessages(require(sm, character.only = TRUE,quietly=TRUE,warn.conflicts = FALSE))))

# option parser
option_list = list(
  make_option(c("-t", "--tree"), type="character", default= "test",
              help="treefile location", metavar="character"),
  make_option(c("-c", "--cluster"), type="character", default= NULL,
              help="cluster assignments", metavar="numeric"),
  make_option(c("-o", "--output"), type="character", default= NULL,
              help="output of pdf", metavar="character"),
  make_option(c("-e", "--exclude"), type="numeric", default= 6,
              help="exclude clusters unter n sequences (default = 6)", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$tree)|is.null(opt$cluster)|is.null(opt$output)){
  print_help(opt_parser)
  stop("options not defined", call.=FALSE)
}

# load tree
tree <- read.iqtree(opt$tree)
tree <- as.phylo(tree)

# load cluster assignment
metadata <- fread(opt$cluster)

if (isTRUE(!colnames(metadata) == c("seq", "cluster"))) {
  colnames(metadata) <- c("seq", "cluster")
}

cluster_counts <- metadata %>% group_by(cluster) %>% summarize(count=n())

for (i in 1:nrow(cluster_counts)) {
  if (cluster_counts$count[i] < opt$exclude) {
    metadata$cluster[metadata$cluster == cluster_counts$cluster[i]] <- "none"
  }
}

print(paste0("Showing only sequences with > ", opt$exclude-1, " sequences"))

# vizualization
tree_vis <- ggtree(tree, layout="daylight", show.legend = T) %<+% metadata +
  geom_tippoint(aes(fill = cluster),
                shape = 21,
                size = 5,
                color = "black") +
  scale_fill_viridis(option="H", discrete = T, begin = 0.1, end = 1)

pdf(opt$output, height = 20, width = 20)

plot(tree_vis)

dev.off()