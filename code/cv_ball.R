library(dplyr)
library(cluster)
options(java.parameters = "-Xmx30g")
library(bartMachine)
set_bart_machine_num_cores(5)
library(reticulate)
python3_path = system('which python3', intern = TRUE)
#use_python(python3_path, required = TRUE)
use_python('/Library/Frameworks/Python.framework/Versions/3.7/bin/python3', required = TRUE)
pkg = import('SciGlass_pkg')


###############################################################################
#:------------------- Parameters setting in BartMachine ---------------------:#
###############################################################################
args_list = list(
  cluster = 1,
  seed = 22,
  beta = 0.5,
  alpha = 0.9,
  k = 2,
  q = 0.9,
  num_trees = 30,
  num_burn_in = 3000,
  num_iterations_after_burn_in = 10000
)

## get parameters from command line
args = commandArgs(trailingOnly = TRUE)
print(args)

# functions for finding named arguments
args_to_list = function(args){
  ind = grep('=', args)  
  args_list = strsplit(args[ind], '=')
  names(args_list) = sapply(args_list, function(x) x[1])
  
  args_list = lapply(args_list, function(x) as.numeric(x[2]))
  args_list
}

# get named arguments
args_list_in = args_to_list(args)

# update non default arguments
ignored = c()
for ( arg in names(args_list_in) ) {
  # Check for unknown argument
  if ( is.null(args_list[[arg]]) ) {
    ignored = c(ignored, arg)
  } else{
    # update if known
    args_list[[arg]] = args_list_in[[arg]]
  }
}

# Print warning message about unknown arguments
if ( length(ignored) > 0 ) {
  cat('Ignoring unkown arguments:',paste(ignored,collapse=', '), '\n')
}



###############################################################################
#:------------------------------- Import datasets ---------------------------:#
###############################################################################
data_file = 'datasets/4_clear_data_overall.csv'
info_file = 'datasets/4_Ca-Mg-Al-Si-O_info.csv'
property_name = c('formation_energy', 'band_gap', 'atom_volume', 'density', 'G', 'K')

ball_obj = pkg$NewDescriptors(data_file, info_file, transform = TRUE)
#comp_obj = pkg$NewDescriptors(data_file, info_file, transform = FALSE)

glass = read.csv(data_file, header = TRUE)
Temperature = glass$Tliq.C
print(nrow(glass))
################################################################################
#: --------------------------- Clustering Splitting ------------------------- :#
################################################################################
ncluster = 10
data_hc = glass[1:4]
dissimilar = daisy(data_hc)
hc_complete_raw = hclust(dissimilar, method = "complete")
cluster_idx = cutree(hc_complete_raw, ncluster)
print(table(cluster_idx))

sil = silhouette (cluster_idx, dissimilar)

#########################################################################
intercept_seq = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)
width_seq = c(0.05, 0.07, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45)
para.grid = expand.grid('width' = width_seq ,'intercept' = intercept_seq)
##########################################################################
model_list = list()
test_pred_overall_df = c()
pred_list = list()
for(i in 1:nrow(para.grid)){
  # Apply one combination of bandwidth & intercept to 'width' and 'intercept'
  width = para.grid$width[i]
  intercept = para.grid$intercept[i]
  
  # Construct new descriptor with the bandwidth & intercept we choose
  ball_data = ball_obj$all_df(width = width, intercept = intercept, property_name = property_name)
  
  # Check if the construction is successful
  print(width)
  print(intercept)
  
  # Obtain the indices for test set
  test_idx = (1:nrow(glass))[(cluster_idx==args_list$cluster)]
  print(length(test_idx))
  
  # Training set
  x.train = ball_data[-test_idx,]
  y.train = Temperature[-test_idx]
  print(nrow(x.train))
  
  # Test set
  x.test = ball_data[test_idx,]
  y.test = Temperature[test_idx]
  print(nrow(x.test))
  
  # Fit a BART model
  bart_fit = bartMachine(x.train, y.train, seed = args_list$seed,
                         beta = args_list$beta,
                         alpha = args_list$alpha,
                         k = args_list$k,
                         q = args_list$q,
                         num_trees = args_list$num_trees,
                         num_burn_in = args_list$num_burn_in,
                         num_iterations_after_burn_in = args_list$num_iterations_after_burn_in,
                         verbose = TRUE, serialize = TRUE)
  model_list[[i]] = bart_fit
  # Predict temperatures in test set
  pred = calc_prediction_intervals(bart_fit, x.test)
  pred_list[[i]] = pred$all_prediction_samples
  pred_df_test = data.frame('pred' = predict(bart_fit, x.test), 'lb' = pred$interval[,1], 
                            'ub' = pred$interval[,2], 'ytrue' = y.test, 'width' = width, 
                            'intercept' = intercept, 'cluster' = args_list$cluster)
  
  test_pred_overall_df = rbind(test_pred_overall_df, pred_df_test)
  
  print(paste0('---------- cluster ', args_list$cluster, '; Bandwidth: ', width, 
               '; Intercept: ', intercept, ' ----------'))
  
  # Save results
  save(test_pred_overall_df, para.grid, model_list, file = paste0('OUTPUT/cv_ball/cluster', args_list$cluster,'.RData'))
}




