# The clustering function we use throughout the script.
cluster_fun <- function(d) hclust(d,method="average")

# Spearman correlation distance function.
dist_fun <- function(x) as.dist(1-cor(t(x), method="spearman"))

# Function to generate a blue to red to white color palette from.
colMap <- div_gradient_pal(muted("red"), mid="#ededed", high=muted("blue"), space="Lab")
redBlueMap <- function(n) colMap(seq(0,1, length.out=n))

# Function to generate palette from
colMap <- div_gradient_pal(muted("white"), high=muted("black"), space="Lab")
whiteBlackMap <- function(n) colMap(seq(0,1, length.out=n))

# Function to transform train and test data based on clusters obtained by cutting a 
# hclust object at cut-off thresh.
clusterFeatures <- function(hcd, orig_data, thresh) {
  # This function cuts the hierarchical tree at asked for threshold. It returns a list with as the
  # first element data merged according to clusters and second argument a list with the 
  # chromosome locations that are within each cluster.
  
  # Do a cut of the hierarchical tree. Use this to generate a list with the number
  # of chromosomes per cluster.
  hcd.cut <- data.frame(cluster = cutree(hcd, h=thresh))
  hcd.cut$chrom_loc <- rownames(hcd.cut) 
  scores_per_cluster <- hcd.cut %>% merge(metadata, by = "chrom_loc") %>% as.data.frame() %>% group_by(cluster) %>% mutate(numInCluster = n())
  not_in_clusters <- scores_per_cluster$chrom_loc[scores_per_cluster$numInCluster == 1]
  scores_per_cluster <- scores_per_cluster %>% filter(numInCluster > 1)
  
  # Append the original data to the tree and derive some summary statistics per cluster.
  transposed.data <- t(orig_data[-1]) %>% melt(varnames = c("chrom_loc", "sample"), value.name = "measurement") %>% mutate(measurement = as.integer(measurement))
  cluster.data <- merge(transposed.data, scores_per_cluster, by = "chrom_loc") %>% mutate(cluster = paste0("clus",cluster)) %>% group_by(cluster, sample) %>% summarize(sum_measurement = sum(measurement)) %>% spread(cluster, sum_measurement)
  cluster.data$sample <- NULL
  merged.data <- data.frame(lapply(orig_data[c("Subgroup",not_in_clusters)],factor), cluster.data)
  
  # Make a dictionary of clusters.
  clusterlist <- scores_per_cluster %>% group_by(cluster) %>% summarize(elements = paste(chrom_loc,collapse=" ")) %>%
    mutate(cluster = paste0("clus",cluster))
  
  result <- list(merged.data, clusterlist)
  return(result)
}


# Function to calculate node purity and number of clusters.
calculateNodeImpurityAndNumClusters <- function(hcd,thresh) {
  # Input of this function is an object of class "hclust" and a threshold for the Spearman correlation distance.
  # The function returns the purity for retrieved clusters, i.e. the fraction of clusters that are confined to 
  # one chromosome and the total number of clusters for the threshold.
  
  # Construct a dataframe with clusters and chromosome locations per segment.
  hcd.cut <- data.frame(cluster = cutree(hcd, h=thresh))
  hcd.cut$chrom_loc <- rownames(hcd.cut)
  
  # Merge this dataframe to the metadata. Then summarize per cluster how many chromosomes are in there.
  purity.df <- merge(hcd.cut, metadata) %>%
    arrange(cluster, Chromosome) %>%
    group_by(cluster, Chromosome) %>%
    filter(row_number() == 1) %>%   # At most one entry per cluster chromosome combination
    group_by(cluster) %>%
    summarize(numberDistinctChromosomes = n())
  
  numClusters <- length(unique(purity.df$cluster))
  purity <- sum(purity.df$numberDistinctChromosomes == 1)/length(purity.df$numberDistinctChromosomes)
  
  return(list(purity, numClusters))
}

# Function to do the cross-validation procedure for gbm.
cv.hcgbm.informal <- function(train, cvfolds, rho) {
  
  # First cluster all training data. Note that this is somewhat problematic as we split the training data
  # by cross-validation, so this allows leakage of information from train to validation. However, we want
  # to use the caret package to tune parameters of gbm. This package does not allow for custom pre-processing
  # so we cannot do hierarchical clustering on the (k-1) folds that are the true training data. Note that with
  # 10-fold cross-validation, this may not be so problematic practically as we could expect the correlations to 
  # be relatively stable for 90% of data.
  
  # Do hierarchical clustering for all training data. Cluster.result is a list with as first element
  # transformed data, second element untransformed data.
  m.train <- train[-1] %>% as.matrix()
  hcd <- cluster_fun(dist_fun(t(m.train)))
  cluster.result <- clusterFeatures(hcd, train, rho)
  
  # Now that we have clustered the data, do training with caret.
  clustered.data <- cluster.result[[1]]
  tc <- trainControl(index = cvfolds, method = 'cv', number = 10, allowParallel = T)
  gbm.caret <- caret::train(Subgroup ~ ., data=clustered.data, method="gbm", 
                            distribution="multinomial", bag.fraction = .75, 
                            trControl=tc, verbose=FALSE, tuneGrid=caretGrid, 
                            metric="Accuracy")
  
  # Return the caret results.
  gbm.caret$results %>% mutate(rho = rho) %>% return()
}