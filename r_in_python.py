import rpy2.robjects as robjects
from rpy2.robjects.packages import importr, data
utils = importr('utils')

utils.install_packages('stats')

r_code = 'write_clustering = function(outdir, tag, label_df, cell_col, cluster_col) {\
  to_write = label_df[c(cell_col, cluster_col)]\
  colnames(to_write)[colnames(to_write) == cell_col] = "cell"\
  colnames(to_write)[colnames(to_write) == cluster_col] = "cluster"\
  write.csv(to_write, paste(outdir, "/", "clustering_", tag, ".csv", sep=""), row.names = FALSE)\
}'



# translate the string r_code to python code





r_vector = robjects.r(r_code)

print(r_vector)