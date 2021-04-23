  # function to prepare the data set:
  # including normalisation (two different methods are available (specify: "relative abundance" or "NMDS")) and factorisation of the wanted classes 
  
  model_input <- function(data, class_name, class_data, method){
    
    # sorting of the data according the sample site numbers
    dataS <- data[order(as.numeric(substr(rownames(data),3,5))),]
    
    #Normalisation step (here: relative abundance)
    if (method == "relative abundance"){
    data_pre <- sweep(dataS, 1, rowSums(dataS), "/")
    data_pre <- as.data.frame(data_pre^0.5)
    # Factorisation of the wanted class 
    data_wc <- cbind(data_pre, class_data)
    data_wc[class_name] <- lapply(data_wc[class_name], as.factor)
    }
    if (method == "NMDS"){
    # NMDS 
    data_pre <- metaMDS(dataS, distance = "bray", k = 4, try = 200, trymax = 2000)
    # Factorisation of the wanted class
    data_wc <- cbind(data_pre$points, class_data)
    data_wc[class_name] <- lapply(data_wc[class_name], as.factor)
    }
    
    if (method == "nothing") {
      # Factorisation of the wanted class 
      data_wc <- cbind(dataS, class_data)
      data_wc[class_name] <- lapply(data_wc[class_name], as.factor)
    }
    
    
    
    return(data_wc)
  }


