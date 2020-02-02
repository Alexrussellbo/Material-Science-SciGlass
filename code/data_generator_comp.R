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


###############################################################################
#:------------------------------- Import datasets ---------------------------:#
###############################################################################
data_file = 'datasets/4_clear_data_overall.csv'
info_file = 'datasets/4_Ca-Mg-Al-Si-O_info.csv'
property_name = c('formation_energy', 'band_gap', 'atom_volume', 'density', 'G', 'K')

comp_obj = pkg$NewDescriptors(data_file, info_file, transform = FALSE)

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
width_seq = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)
para.grid = expand.grid('width' = width_seq ,'intercept' = intercept_seq)

data_list = list()
for(i in 1:nrow(para.grid)){
  # Apply one combination of bandwidth & intercept to 'width' and 'intercept'
  width = para.grid$width[i]
  intercept = para.grid$intercept[i]
  
  # Construct new descriptor with the bandwidth & intercept we choose
  comp_data = comp_obj$all_df(width = width, intercept = intercept, property_name = property_name)
  data_list[[i]] = comp_data
  
  # Check if the construction is successful
  print(width)
  print(intercept)
  print(any(is.na(comp_data)))
  
}
save(data_list, file = 'new_data_comp.RData')
