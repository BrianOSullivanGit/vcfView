
#
# Extensible hook for "neutralitytestr"
# Check for neutral evolution within this range..
# (TODO!!!! put in checks to ensure range make sence with plody assumptions etc..)
#

if (!requireNamespace("neutralitytestr", quietly = TRUE))
  BiocManager::install("neutralitytestr")

library(neutralitytestr)

exten_str_Nevol="Neutrality"

exten_fn_Nevol <- function(result, lower_frange, upper_frange)
{ 
  # Subset out all the SNVs within required allele freq. range.
  snvIdx = which(result$AF > lower_frange &
                 result$AF <= upper_frange)
  
  n <- neutralitytest(result$AF[snvIdx], fmin = lower_frange, fmax = upper_frange)
  
  return(lsq_plot(n))
}

