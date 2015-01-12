start_parallel <- function(source_func) {
  require(parallel)
  require(foreach)
  require(iterators)
  require(doParallel)
  if (!is.na(pmatch("Windows", Sys.getenv("OS")))) {
    cluster <- makeCluster(detectCores() - 1, type = "SOCK")
    registerDoParallel(cluster)
    clusterExport(cluster, source_func) 
  } else {
    registerDoParallel(cores=detectCores() - 1)
  }
  cluster
}