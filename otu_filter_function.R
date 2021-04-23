
otu_filter_function <- function(data, readsum, freqfilter, readcount){
  # readsum: threshold of the sum of reads that has to be present in an otu to be kept. 
  data_t <- t(data)
  data_F <- rowSums(data_t)
  data_new <- cbind(data_t, data_F)
  data_new_da <- as.data.frame(data_new)
  colnames(data_new_da)[ncol(data_new_da)] <- "Sum"
  data_ord <- data_new_da[order(-data_new_da$Sum), ]
  if (readsum > 0){
  index <- which(data_ord$Sum < readsum)
  data_wo <- data_ord[-index, -ncol(data_ord)]}
  # set to -3 for single test samples: needs fixing (ncol to delete last column will also delete colnames, which are needed for comparison)
  else {data_wo <- data_ord[,-3]}
  
  # freqfilter: threshold of how many samples have to have read counts present in the otu to be kept. 
  if(freqfilter > 0){
    data_com <- data_wo[apply(data_wo == 0,1, sum) <= freqfilter,]
  } else {
    data_com <- data_wo
  }
  
  # readcount: sets readcounts below a certain threshold to zero, to ignore low read count data 
  data_com_t <- t(data_com)
  if (readcount > 0){
    data_fin <- ifelse(data_com_t < readcount, 0, data)
  } else {
    data_fin <- data_com_t
  }
  
  return(data_fin)
}


