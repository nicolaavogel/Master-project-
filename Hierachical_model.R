library(here)
library(dplyr)
library(caret)
library(stringr)
library(vegan)
library(cvms)
library(ggrepel)
library(data.tree)
source("model_input_function.R")
source("otu_filter_function.R")



nested.classification <- function(trainingData, testData, classData, testData_colnames){
  # Adjusting classes for the first round of the model: 
  # classification on the big three classes (Agriculture, Forest and Grassland) 
  fi <- classData
  fi$Forest <- ifelse(fi$Forest == 1, "Forest", NA)
  fi$Agriculture <- ifelse(fi$Agriculture == 1, "Agriculture", NA)
  # seperate adjusting of the Heathland class: Dwarfshrub will be added to the "Grassland" class making it a more general class
  fi$Heathland <- ifelse(fi$Heathland == 1, "Grassland", NA)
  fi$Dwarfshrubs <- ifelse(fi$Dwarfshrubs == 1, "Grassland", NA)
  
  fi$testClass <- apply(fi[,c(2,3,8,9)], 1, function(x) x[!is.na(x)][1])
  fi$testClass[is.na(fi$testClass)] <- "No Habitat"
  # Adding two "No Habitat" sample to the new Heathland class, based on the MDS plots. 
  fi$testClass[c(61, 63)] <- "Grassland"
  fi$testClass[67] <- "Forest"
  
  
  # Adjusting testData to the needed format: input can be single data or a loop through a dataset
  testData_frame <- as.data.frame(testData, col.names = testData_colnames)
  testData_t <- t(testData_frame)
  testData_f <- otu_filter_function(testData_t, 0, 0, 0)
  # Relative abundance normalisation for the test Data 
  testData_rel <- sweep(testData_f, 1, rowSums(testData_f), "/")
  testData_rel_f<- as.data.frame(testData_rel^0.5)
  
  # Indexing and finding matching column name (sequence tags) for classification of the test Data. 
  index <- which(colnames(trainingData) %in% colnames(testData_rel_f))
  index2 <- which(colnames(testData_rel_f) %in% colnames(trainingData))
  
  tar_new <- trainingData[,index]
  tar_test <- testData_rel_f[,index2]
  
  
  # model_input function; here: set to a relative abundance normalisation 
  tar_pre <- model_input(tar_new, "testClass", fi, "relative abundance")
  
  # First KNN model training, the Classes are Agriculture, Forest, Grassland and No Habitat
  tar_model <- train(testClass ~. , data = tar_pre[,c(1:length(colnames(tar_new)), length(colnames(tar_pre)))],
                     method = "knn", metric = "Kappa", 
                     trControl = trainControl("LOOCV"),
                     tuneGrid = expand.grid(k = c(3:15)))
  print(tar_model)
  
  # Adding a plot for the kNN algorithm to show the bes t selection of k 
  knn1 <- plot(tar_model)
  
  # creating a confusion matrix for the best k parameter
  best.k <- tar_model$bestTune[[1]]
  best.knn <- tar_model$pred[tar_model$pred$k == best.k, ]
  cm1 <- confusionMatrix(best.knn$pred, best.knn$obs)
  cm_df <- as.data.frame(cm1$table)
  CMP1 <- plot_confusion_matrix(cm_df,
                        place_x_axis_above = F, targets_col = "Reference", predictions_col = "Prediction",
                        counts_col = "Freq", palette = "Greens", add_normalized = F,
                        add_row_percentages = T, add_col_percentages = T, counts_on_top = T)
  print(CMP1)
  
  # First prediction output
  pred <- predict(tar_model, newdata = tar_test)
  testPred <- predict(tar_model, newdata =  tar_test, type = "prob")
  print(testPred) # type give the proportions of the vote as Probability
  
  # Adding the unknown sample to the training Data set
  tar_new_plus <- plyr::rbind.fill(tar_pre, tar_test[2,])
  rownames(tar_new_plus)[1:130] <- rownames(tar_pre)
  rownames(tar_new_plus)[131] <- "Unknown sample"
  tar_new_plus[131,ncol(tar_new_plus)] <- pred[1]
  
  # getting the distances for a visualisation including the unknown sample 
  unknown.tar <- metaMDS(tar_new_plus[,1:length(colnames(tar_new))], distance = "bray", k = 6, try = 200, trymax = 2000)
  # Stress lowers down to 0.0832 (good fit)
  
  # NMDS data mostly for visualisation 
  nmds.data.test <- data.frame(Sample = rownames(tar_new_plus), 
                               color = tar_new_plus$testClass,
                               X = unknown.tar$points[,1],
                               Y = unknown.tar$points[,2],
                               X2 = unknown.tar$points[,3],
                               Y2 = unknown.tar$points[,4],
                               X3 = unknown.tar$points[,5],
                               Y3 = unknown.tar$points[,6])
  
  
  # Calculating the centroids of the 4 classes 
  class.dist <- cbind.data.frame(unknown.tar$points, tar_new_plus$testClass)
  colnames(class.dist)[7] <- "testClass"
  
  centroids <- aggregate(. ~ testClass, data = class.dist, mean)
  
  # NMDS plot with the centroids of the class and the unknown sample 
  NMDS_all <- ggplot()+
    geom_point(data = nmds.data.test, aes(x = X3, y = Y3, group = color, shape = color, color = color, size = color)) +
    scale_shape_manual(name = "Habitat classes", values = c("Agriculture" = 16, "Forest" = 16, "Grassland" = 16, "No Habitat" = 16)) +
    scale_colour_manual(name = "Habitat classes", values = c("Agriculture" = "#D16103",
                                                             "Forest" = "#52854C",
                                                             "Grassland" = "#4E84C4", 
                                                             "No Habitat" = "#C4961A")) +
    scale_size_manual(name = "Habitat classes", values = c("Agriculture" = 4, "Forest" = 4, "Grassland" = 4, "No Habitat" = 4)) +
    geom_point(data = centroids, aes(x = MDS5, y = MDS6, color = centroids$testClass), size = 5, shape = 15) +
    theme_bw() +
    ylab("NMDS 6") +
    xlab("NMDS 5") +
    theme(axis.title = element_text(size = 18)) +
    theme(axis.text = element_text(size = 18)) +
    theme(legend.title = element_text(size = 23)) +
    theme(legend.text = element_text(size = 23)) +
    geom_text_repel(data = nmds.data.test[131,], aes(x = X3, y = Y3, label = Sample), size = 5, 
                    box.padding = unit(0.35, "lines"), point.padding = unit(0.3, "lines"))
  
  print(NMDS_all)
  
  # computing the euclidean distance between every in class and out of class data point
  # to each centroid of the class for a likelihood estimation 
  # defining the euclidean distance function manually (could also be done with dist(rbind(x1, x2)))
  euc.dist <- function(x1, x2) {sqrt(sum((x1 - x2) ^2))}
  euc.dist2 <- function(xa, xb, x1a, x1b, x2a, x2b, x3a, x3b, x4a, x4b, x5a, x5b)
  {sqrt(sum((xa - xb) ^2,(x1a - x1b)^2,(x2a - x2b)^2,(x3a - x3b)^2,(x4a - x4b)^2,(x5a - x5b)^2))}
  
  dist1 <- NULL 
  for (i in 1:nrow(unknown.tar$points)){
    dist1[i] <- euc.dist2(unknown.tar$points[i,1],centroids[1,2],
                          unknown.tar$points[i,2],centroids[1,3],
                          unknown.tar$points[i,3],centroids[1,4],
                          unknown.tar$points[i,4],centroids[1,5],
                          unknown.tar$points[i,5],centroids[1,6],
                          unknown.tar$points[i,6],centroids[1,7])
  }
  
  dist2 <- NULL
  for (i in 1:nrow(unknown.tar$points)){
    dist2[i] <- euc.dist2(unknown.tar$points[i,1],centroids[2,2],
                          unknown.tar$points[i,2],centroids[2,3],
                          unknown.tar$points[i,3],centroids[2,4],
                          unknown.tar$points[i,4],centroids[2,5],
                          unknown.tar$points[i,5],centroids[2,6],
                          unknown.tar$points[i,6],centroids[2,7])
  }
  
  dist3 <- NULL
  for (i in 1:nrow(unknown.tar$points)){
    dist3[i] <- euc.dist2(unknown.tar$points[i,1],centroids[3,2],
                          unknown.tar$points[i,2],centroids[3,3],
                          unknown.tar$points[i,3],centroids[3,4],
                          unknown.tar$points[i,4],centroids[3,5],
                          unknown.tar$points[i,5],centroids[3,6],
                          unknown.tar$points[i,6],centroids[3,7])
  }
  
  dist4 <- NULL
  for (i in 1:nrow(unknown.tar$points)){
    dist4[i] <- euc.dist2(unknown.tar$points[i,1],centroids[4,2],
                          unknown.tar$points[i,2],centroids[4,3],
                          unknown.tar$points[i,3],centroids[4,4],
                          unknown.tar$points[i,4],centroids[4,5],
                          unknown.tar$points[i,5],centroids[4,6],
                          unknown.tar$points[i,6],centroids[4,7])
  }
  
  distMatrix_all_w <- cbind.data.frame(dist1, dist2, dist3, dist4, nmds.data.test$color)
  rownames(distMatrix_all_w) <- nmds.data.test$Sample
  colnames(distMatrix_all_w) <- c("Agriculture", "Forest", "Grassland", "No Habitat", "Class")
  
  distMatrix_all <- distMatrix_all_w[-131,]
  
  clSubset <- subset(distMatrix_all[,which(colnames(distMatrix_all) == pred[1])], distMatrix_all$Class == pred[1])
  mean_F <- mean(clSubset)
  sd_F <- sd(clSubset)
  
  ncSubsetwo <- subset(distMatrix_all[,c(which(colnames(distMatrix_all) == pred[1]))], distMatrix_all$Class != pred[1])
  mean_nF <- mean(ncSubsetwo)
  sd_nF <- sd(ncSubsetwo)
  
  unknown_dist <- distMatrix_all_w[nrow(distMatrix_all_w),which(colnames(distMatrix_all_w) == pred[1])]
  
  # Probability of the unknown_dist begin from the distribution (from the wanted centroid class)
  A_B <- pnorm(unknown_dist, mean_F, sd_F) 
  
  # Probability for event 2 the unknown samples does not come from the predicted class
  A_B2 <-  pnorm(unknown_dist, mean_nF, sd_nF)
  
  prior <- 0.25 #The chosen prior here is 0.25 based on the four centroids. This suggest that each unknown samples's distance can be equally belong to one of the four classes.
  
  # Likelihood function for event A 
  Bayestest <- A_B * prior 
  
  #Likelihood function for event B 
  Bayestest2 <- A_B2 * prior 
  
  LR <- Bayestest / Bayestest2 
  print(LR)
  
  
  p1 <- ggplot()+ 
    geom_histogram(aes(x = clSubset,y = ..density.., fill = "In class distance"), binwidth = 0.1, alpha = .8, position = "identity") + 
    geom_histogram(aes(x = ncSubsetwo, y = ..density.., fill = "Out of class distance"), binwidth = 0.1,alpha = .8, position = "identity") +
    labs(fill = "Distributions") +
    ylab("Density")+
    xlab("Histogram of in-class distances and out-of class distances") + 
    theme_bw() + 
    stat_function(aes(x = clSubset), fun = dnorm, args = list(mean = mean(clSubset), sd = sd(clSubset)), color = "red") +
    stat_function(aes(x = ncSubsetwo),fun = dnorm, args = list(mean = mean(ncSubsetwo), sd = sd(ncSubsetwo)), color = "blue") +
    geom_vline(xintercept = unknown_dist, linetype = "dotted", color = "black") +
    geom_text(aes(x = unknown_dist + 0.03, label = "distance of the unknown test sample", y = 1), color = "black", 
              angle = 90) +
    theme(axis.title = element_text(size = 18)) +
    theme(axis.text = element_text(size = 18)) +
    theme(legend.text = element_text(size = 23)) +
    theme(legend.title = element_text(size = 23))

  print(p1)
  
  # Adding the best matching site based on the multidimensionl euclidean distance between the points (comparison of picture)
  
  nmds.values.c <- subset(nmds.data.test[-131,], color == pred[1])[,3:8]
  
  
  matching <- NULL
  for (i in 1:nrow(nmds.values.c)){
    
    matching[i] <- euc.dist2(nmds.values.c[i, 1], unknown.tar$points[131,1], nmds.values.c[i, 2],unknown.tar$points[131,2],
                             nmds.values.c[i, 3],unknown.tar$points[131,3], nmds.values.c[i, 4], unknown.tar$points[131,4],
                             nmds.values.c[i, 5], unknown.tar$points[131,5], nmds.values.c[i, 6], unknown.tar$points[131,6])
    
    
  }
  
  min.dist.id <- which.min(matching)
  best.match <- rownames(nmds.values.c)[min.dist.id]
  print(best.match)
  
  list1 <- list(tar_model, pred, testPred,knn1,CMP1, NMDS_all, LR, p1, best.match)
  
  if(pred[1] == "Forest"){
    # classData will be adjusted once more for the classification of Forest samples: Coniferous or Deciduous
    fifi <- classData
    fifi$Coniferous <- ifelse(fifi$Coniferous == 1, "Coniferous", NA)
    fifi$Beech <- ifelse(fifi$Beech == 1, "Deciduous", NA)
    fifi$Oak <- ifelse(fifi$Oak == 1, "Deciduous", NA)
    fifi$Willow <- ifelse(fifi$Willow == 1, "Deciduous", NA)
    fifi$Alnus <- ifelse(fifi$Alnus == 1, "Deciduous", NA)
    
    fifi$testClass2 <- apply(fifi[,c(4:7,12)], 1, function(x) x[!is.na(x)][1])
    fifi$testClass2[is.na(fifi$testClass2)] <- "No Specification"
    fifi$Forest[67] <- 1
    
    new_data2 <- model_input(tar_new, "testClass2", fifi, "relative abundance")
    new_input_data <- subset(new_data2, Forest == 1)
    
    
    
    
    nested.model1 <- train(testClass2 ~ . , new_input_data[,c(1:dim(tar_new)[2], ncol(new_input_data))], 
                           method = "knn", metric = "Kappa", trControl = trainControl("LOOCV"), tuneGrid = expand.grid(k = c(3:10)))
    print(nested.model1)
    
    knn2 <- plot(nested.model1)
    print(knn2)
    
    pred2 <- predict(nested.model1, newdata = tar_test)
    testPred2 <- predict(nested.model1, newdata =  tar_test, type = "prob")
    print(testPred2)
    
    best.k <- nested.model1$bestTune[[1]]
    best.knn <- nested.model1$pred[nested.model1$pred$k == best.k, ]
    cm2 <- confusionMatrix(best.knn$pred, best.knn$obs)
    cm2_df <- as.data.frame(cm2$table)
    CMP2 <- plot_confusion_matrix(cm2_df,
                                  place_x_axis_above = F, targets_col = "Reference", predictions_col = "Prediction",
                                  counts_col = "Freq", palette = "Greens", add_normalized = F,
                                  add_row_percentages = T, add_col_percentages = T, counts_on_top = T)
    print(CMP2)
    
    tar_new_plus <- plyr::rbind.fill(new_input_data, tar_test[2,])
    rownames(tar_new_plus)[1:53] <- rownames(new_input_data)
    rownames(tar_new_plus)[54] <- "Unknown sample"
    tar_new_plus[54,ncol(tar_new_plus)] <- pred2[1]
    
    # getting the distances for a visualisation including the unknown sample 
    unknown.tar <- metaMDS(tar_new_plus[,1:length(colnames(tar_new))], distance = "bray", k = 6, try = 200, trymax = 2000)
    # Stress lowers down to 0.05674773 (good fit)
    
    # NMDS data mostly for visualisation 
    nmds.data.test <- data.frame(Sample = rownames(tar_new_plus), 
                                 color = tar_new_plus$testClass2,
                                 X = unknown.tar$points[,1],
                                 Y = unknown.tar$points[,2],
                                 X2 = unknown.tar$points[,3],
                                 Y2 = unknown.tar$points[,4],
                                 X3 = unknown.tar$points[,5],
                                 Y3 = unknown.tar$points[,6])
    print(nmds.data.test)
    
    
    # Calculating the centroids of the 4 classes 
    class.dist <- cbind.data.frame(unknown.tar$points, tar_new_plus$testClass2)
    colnames(class.dist)[7] <- "testClass"
    
    centroids <- aggregate(. ~ testClass, data = class.dist, mean)
    
    # NMDS plot with the centroids of the class and the unknown sample 
    NMDS_Forest <- ggplot()+
      geom_point(data = nmds.data.test, aes(x = X, y = Y, group = color, shape = color, color = color, size = color)) +
      scale_shape_manual(name = "Habitat classes", values = c("Coniferous" = 16, "Deciduous" = 16, "No Specification" = 16)) +
      scale_colour_manual(name = "Habitat classes", values = c("Coniferous" = "#52854C",
                                                               "Deciduous" = "#D16103",
                                                               "No Specification" = "#4E84C4")) +
      scale_size_manual(name = "Habitat classes", values = c("Coniferous" = 4, "Deciduous" = 4, "No Specification" = 4)) +
      geom_point(data = centroids, aes(x = MDS1, y = MDS2, color = centroids$testClass), size = 5, shape = 15) +
      theme_bw() +
      ylab("NMDS 2") +
      xlab("NMDS 1") +
      theme(axis.title = element_text(size = 18)) +
      theme(axis.text = element_text(size = 18)) +
      theme(legend.title = element_text(size = 23)) +
      theme(legend.text = element_text(size = 23)) +
      geom_text_repel(data = nmds.data.test[54,], aes(x = X, y = Y, label = Sample), size = 5, 
                      box.padding = unit(0.35, "lines"), point.padding = unit(0.3, "lines"))
    print(NMDS_Forest)
    
    # computing the euclidean distance between every in class and out of class data point
    # to each centroid of the class for a likelihood estimation
    
    dist1 <- NULL 
    for (i in 1:nrow(unknown.tar$points)){
      dist1[i] <- euc.dist2(unknown.tar$points[i,1],centroids[1,2],
                            unknown.tar$points[i,2],centroids[1,3],
                            unknown.tar$points[i,3],centroids[1,4],
                            unknown.tar$points[i,4],centroids[1,5],
                            unknown.tar$points[i,5],centroids[1,6],
                            unknown.tar$points[i,6],centroids[1,7])
    }
    
    dist2 <- NULL
    for (i in 1:nrow(unknown.tar$points)){
      dist2[i] <- euc.dist2(unknown.tar$points[i,1],centroids[2,2],
                            unknown.tar$points[i,2],centroids[2,3],
                            unknown.tar$points[i,3],centroids[2,4],
                            unknown.tar$points[i,4],centroids[2,5],
                            unknown.tar$points[i,5],centroids[2,6],
                            unknown.tar$points[i,6],centroids[2,7])
    }
    
    dist3 <- NULL
    for (i in 1:nrow(unknown.tar$points)){
      dist3[i] <- euc.dist2(unknown.tar$points[i,1],centroids[3,2],
                            unknown.tar$points[i,2],centroids[3,3],
                            unknown.tar$points[i,3],centroids[3,4],
                            unknown.tar$points[i,4],centroids[3,5],
                            unknown.tar$points[i,5],centroids[3,6],
                            unknown.tar$points[i,6],centroids[3,7])
    }
    
    
    
    distMatrix_Forest_w <- cbind.data.frame(dist1, dist2, dist3, nmds.data.test$color)
    rownames(distMatrix_Forest_w) <- nmds.data.test$Sample
    colnames(distMatrix_Forest_w) <- c("Coniferous", "Deciduous", "No Specification", "Class")
    
    distMatrix_Forest <- distMatrix_Forest_w[-nrow(distMatrix_Forest_w),]
    
    clSubset <- subset(distMatrix_Forest[,which(colnames(distMatrix_Forest) == pred2[1])], distMatrix_Forest$Class == pred2[1])
    mean_F <- mean(clSubset)
    sd_F <- sd(clSubset)
    
    ncSubsetwo <- subset(distMatrix_Forest[,c(which(colnames(distMatrix_Forest) == pred2[1]))], distMatrix_Forest$Class != pred2[1])
    mean_nF <- mean(ncSubsetwo)
    sd_nF <- sd(ncSubsetwo)
    
    unknown_dist <- distMatrix_Forest_w[nrow(distMatrix_Forest_w),which(colnames(distMatrix_Forest_w) == pred2[1])]
    
    A_B <- pnorm(unknown_dist, mean_F, sd_F) # Probability of the unknown_dist begin from the distribution (from the wanted centroid class)
    
    # Probability for event 2 the unknown samples does not come from the predicted class
    A_B2 <-  pnorm(unknown_dist, mean_nF, sd_nF)
    
    prior <- 0.25 #The chosen prior here is 0.25 based on the four centroids. This suggest that each unknown samples's distance can be equally belong to one of the four classes.
    
    Bayestest <- A_B * prior 
    Bayestest2 <- A_B2 * prior 
    
    
    LR2 <- Bayestest / Bayestest2
    print(LR2)
    
    p2 <- ggplot()+ 
      geom_histogram(aes(x = clSubset,y = ..density.., fill = "In class distance"), binwidth = 0.1, alpha = .8, position = "identity") + 
      geom_histogram(aes(x = ncSubsetwo, y = ..density.., fill = "Out of class distance"), binwidth = 0.1,alpha = .8, position = "identity") +
      labs(fill = "Distributions") +
      ylab("Density")+
      xlab("Histogram of in-class distances and out-of class distances") + 
      theme_bw() + 
      stat_function(aes(x = clSubset), fun = dnorm, args = list(mean = mean(clSubset), sd = sd(clSubset)), color = "red") +
      stat_function(aes(x = ncSubsetwo),fun = dnorm, args = list(mean = mean(ncSubsetwo), sd = sd(ncSubsetwo)), color = "blue") +
      geom_vline(xintercept = unknown_dist, linetype = "dotted", color = "black") +
      geom_text(aes(x = unknown_dist + 0.03, label = "distance of the unknown test sample", y = 1), color = "black", 
                angle = 90) +
      theme(axis.title = element_text(size = 18)) +
      theme(axis.text = element_text(size = 18)) +
      theme(legend.text = element_text(size = 23)) +
      theme(legend.title = element_text(size = 23))
    print(p2)
    
    # Adding the best matching site based on the euclidean distance between the points (comparison of picture)
    
    nmds.values.c <- subset(nmds.data.test[-54,], color == pred2[1])[,3:8]
    
    matching2 <- NULL
    for (i in 1:nrow(nmds.values.c)){
      
      matching2[i] <- euc.dist2(nmds.values.c[i, 1], unknown.tar$points[54,1], nmds.values.c[i, 2],unknown.tar$points[54,2],
                                nmds.values.c[i, 3],unknown.tar$points[54,3], nmds.values.c[i, 4], unknown.tar$points[54,4],
                                nmds.values.c[i, 5], unknown.tar$points[54,5], nmds.values.c[i, 6], unknown.tar$points[54,6])
      
      
    }
    
    min.dist.id <- which.min(matching2)
    best.match2 <- rownames(nmds.values.c)[min.dist.id]
    print(best.match2)
    
    list2 <- list(list1, nested.model1, knn2, CMP2, pred2, testPred2, NMDS_Forest, LR2, p2, best.match2)
  }
  else { return(list1)}
  
  if (pred2[1] == "Deciduous"){
    # ClassDara from the previous classification will be used and altered for all Deciduous samples: Oak, Willow, Alnus, Beech
    fifi$Beech <- ifelse(fifi$Beech == "Deciduous", "Beech", NA)
    fifi$Oak <- ifelse(fifi$Oak == "Deciduous", "Oak", NA)
    fifi$Willow <- ifelse(fifi$Willow == "Deciduous", "Willow", NA)
    fifi$Alnus <- ifelse(fifi$Alnus == "Deciduous", "Alnus", NA)
    fifi$testClass3 <- apply(fifi[,c(5:7,12)], 1, function(x) x[!is.na(x)][1])
    #fifi$testClass3[is.na(fifi$testClass3)] <- "No Specification"
    
    new_data3 <-model_input(tar_new, "testClass3", fifi, "relative abundance")
    input_data <- subset(new_data3, testClass2 == "Deciduous" & Forest == 1)
    input_data$testClass3
    
    nested.model2 <- train(testClass3 ~ ., input_data[,c(1:dim(tar_new)[2], ncol(input_data))],
                           method = "knn", metric = "Kappa", trControl = trainControl("LOOCV"), tuneGrid = expand.grid(k = c(3:7)))
    print(nested.model2)
    
    knn3 <- plot(nested.model2)
    print(knn3)
    
    pred3 <- predict(nested.model2, newdata = tar_test)
    print(pred3)
    testPred3 <- predict(nested.model2, newdata =  tar_test, type = "prob")
    print(testPred3)
    
    best.k <- nested.model2$bestTune[[1]]
    best.knn <- nested.model2$pred[nested.model2$pred$k == best.k, ]
    cm3 <- confusionMatrix(best.knn$pred, best.knn$obs)
    cm3_df <- as.data.frame(cm3$table)
    CMP3 <- plot_confusion_matrix(cm3_df,
                                  place_x_axis_above = F, targets_col = "Reference", predictions_col = "Prediction",
                                  counts_col = "Freq", palette = "Greens", add_normalized = F,
                                  add_row_percentages = T, add_col_percentages = T, counts_on_top = T)
    print(CMP3)
   
    tar_new_plus <- plyr::rbind.fill(input_data, tar_test[2,])
    rownames(tar_new_plus)[1:41] <- rownames(input_data)
    rownames(tar_new_plus)[42] <- "Unknown sample"
    tar_new_plus[42,ncol(tar_new_plus)] <- pred3[1]
    
    # getting the distances for a visualisation including the unknown sample 
    unknown.tar <- metaMDS(tar_new_plus[,1:length(colnames(tar_new))], distance = "bray", k = 6, try = 200, trymax = 2000)
      # stress = 0.0547 (good fit)
    
    # NMDS data mostly for visualisation 
    nmds.data.test <- data.frame(Sample = rownames(tar_new_plus), 
                                 color = tar_new_plus$testClass3,
                                 X = unknown.tar$points[,1],
                                 Y = unknown.tar$points[,2],
                                 X2 = unknown.tar$points[,3],
                                 Y2 = unknown.tar$points[,4],
                                 X3 = unknown.tar$points[,5],
                                 Y3 = unknown.tar$points[,6])
    
    
    # Calculating the centroids of the 4 classes 
    class.dist <- cbind.data.frame(unknown.tar$points, tar_new_plus$testClass3)
    colnames(class.dist)[7] <- "testClass"
    
    centroids <- aggregate(. ~ testClass, data = class.dist, mean)
    
    # NMDS plot with the centroids of the class and the unknown sample 
    NMDS_Dec <- ggplot()+
      geom_point(data = nmds.data.test, aes(x = X, y = Y, group = color, shape = color, color = color, size = color)) +
      scale_shape_manual(name = "Tree classes", values = c("Alnus" = 16, "Beech" = 16, "Oak" = 16, "Willow" = 16)) +
      scale_colour_manual(name = "Tree classes", values = c("Alnus" = "#D16103",
                                                               "Beech" = "#52854C",
                                                               "Oak" = "#4E84C4", 
                                                               "Willow" = "#C4961A")) +
      scale_size_manual(name = "Tree classes", values = c("Alnus" = 4, "Beech" = 4, "Oak" = 4, "Willow" = 4)) +
      geom_point(data = centroids, aes(x = MDS1, y = MDS2, color = centroids$testClass), size = 5, shape = 15) +
      theme_bw() +
      ylab("NMDS 2") +
      xlab("NMDS 1") +
      theme(axis.title = element_text(size = 18)) +
      theme(axis.text = element_text(size = 18)) +
      theme(legend.title = element_text(size = 23)) +
      theme(legend.text = element_text(size = 23)) +
      geom_text_repel(data = nmds.data.test[42,], aes(x = X, y = Y, label = Sample), size = 5, 
                      box.padding = unit(0.35, "lines"), point.padding = unit(0.3, "lines"))
    print(NMDS_Dec)
    
    # computing the euclidean distance between every in class and out of class data point
    # to each centroid of the class for a likelihood estimation 
    
    dist1 <- NULL 
    for (i in 1:nrow(unknown.tar$points)){
      dist1[i] <- euc.dist2(unknown.tar$points[i,1],centroids[1,2],
                            unknown.tar$points[i,2],centroids[1,3],
                            unknown.tar$points[i,3],centroids[1,4],
                            unknown.tar$points[i,4],centroids[1,5],
                            unknown.tar$points[i,5],centroids[1,6],
                            unknown.tar$points[i,6],centroids[1,7])
    }
    
    dist2 <- NULL
    for (i in 1:nrow(unknown.tar$points)){
      dist2[i] <- euc.dist2(unknown.tar$points[i,1],centroids[2,2],
                            unknown.tar$points[i,2],centroids[2,3],
                            unknown.tar$points[i,3],centroids[2,4],
                            unknown.tar$points[i,4],centroids[2,5],
                            unknown.tar$points[i,5],centroids[2,6],
                            unknown.tar$points[i,6],centroids[2,7])
    }
    
    dist3 <- NULL
    for (i in 1:nrow(unknown.tar$points)){
      dist3[i] <- euc.dist2(unknown.tar$points[i,1],centroids[3,2],
                            unknown.tar$points[i,2],centroids[3,3],
                            unknown.tar$points[i,3],centroids[3,4],
                            unknown.tar$points[i,4],centroids[3,5],
                            unknown.tar$points[i,5],centroids[3,6],
                            unknown.tar$points[i,6],centroids[3,7])
    }
    
    dist4 <- NULL
    for (i in 1:nrow(unknown.tar$points)){
      dist4[i] <- euc.dist2(unknown.tar$points[i,1],centroids[4,2],
                            unknown.tar$points[i,2],centroids[4,3],
                            unknown.tar$points[i,3],centroids[4,4],
                            unknown.tar$points[i,4],centroids[4,5],
                            unknown.tar$points[i,5],centroids[4,6],
                            unknown.tar$points[i,6],centroids[4,7])
    }
    
    distMatrix_Dec_w <- cbind.data.frame(dist1, dist2, dist3, dist4, nmds.data.test$color)
    rownames(distMatrix_Dec_w) <- nmds.data.test$Sample
    colnames(distMatrix_Dec_w) <- c("Alnus", "Beech", "Oak", "Willow", "Class")
    
    distMatrix_Dec <- distMatrix_Dec_w[-nrow(distMatrix_Dec_w),]
    
    clSubset <- subset(distMatrix_Dec[,which(colnames(distMatrix_Dec) == pred3[1])], distMatrix_Dec$Class == pred3[1])
    mean_F <- mean(clSubset)
    sd_F <- sd(clSubset)
    
    ncSubsetwo <- subset(distMatrix_Dec[,c(which(colnames(distMatrix_Dec) == pred3[1]))], distMatrix_Dec$Class != pred3[1])
    mean_nF <- mean(ncSubsetwo)
    sd_nF <- sd(ncSubsetwo)
    
    unknown_dist <- distMatrix_Dec_w[nrow(distMatrix_Dec_w),which(colnames(distMatrix_Dec_w) == pred3[1])]
    
    A_B <- pnorm(unknown_dist, mean_F, sd_F) # Probability of the unknown_dist begin from the distribution (from the wanted centroid class)

    # Probability for event 2 the unknown samples does not come from the predicted class
    A_B2 <-  pnorm(unknown_dist, mean_nF, sd_nF)
    
    
    prior <- 0.25 #The chosen prior here is 0.25 based on the four centroids. This suggest that each unknown samples's distance can be equally belong to one of the four classes.
    
    
    Bayestest <- A_B * prior / (A_B2 * prior)
    
    LR3 <- Bayestest 
    print(LR3)
    
    p3 <- ggplot()+ 
      geom_histogram(aes(x = clSubset,y = ..density.., fill = "In class distance"), binwidth = 0.1, alpha = .8, position = "identity") + 
      geom_histogram(aes(x = ncSubsetwo, y = ..density.., fill = "Out of class distance"), binwidth = 0.1,alpha = .8, position = "identity") +
      labs(fill = "Distributions") +
      ylab("Density")+
      xlab("Histogram of in-class distances and out-of class distances") + 
      theme_bw() + 
      stat_function(aes(x = clSubset), fun = dnorm, args = list(mean = mean(clSubset), sd = sd(clSubset)), color = "red") +
      stat_function(aes(x = ncSubsetwo),fun = dnorm, args = list(mean = mean(ncSubsetwo), sd = sd(ncSubsetwo)), color = "blue") +
      geom_vline(xintercept = unknown_dist, linetype = "dotted", color = "black") +
      geom_text(aes(x = unknown_dist + 0.03, label = "distance of the unknown test sample", y = 1), color = "black", 
                angle = 90) +
      theme(axis.title = element_text(size = 18)) +
      theme(axis.text = element_text(size = 18)) +
      theme(legend.text = element_text(size = 23)) +
      theme(legend.title = element_text(size = 23))
    print(p3)
    
    # Adding the best matching site based on the euclidean distance between the points (comparison of picture)
    
    nmds.values.c <- subset(nmds.data.test[-42,], color == pred3[1])[,3:8]
    
    matching3 <- NULL
    for (i in 1:nrow(nmds.values.c)){
      
      matching3[i] <- euc.dist2(nmds.values.c[i, 1], unknown.tar$points[42,1], nmds.values.c[i, 2],unknown.tar$points[42,2],
                                nmds.values.c[i, 3],unknown.tar$points[42,3], nmds.values.c[i, 4], unknown.tar$points[42,4],
                                nmds.values.c[i, 5], unknown.tar$points[42,5], nmds.values.c[i, 6], unknown.tar$points[42,6])
      
      
    }
    
    min.dist.id <- which.min(matching3)
    best.match3 <- rownames(nmds.values.c)[min.dist.id]
    print(best.match3)
    
    list3 <- list(list1, list2, nested.model2, knn3, CMP3, pred3, testPred3, NMDS_Dec, LR3, p3, best.match3)
    
    return(list3)
  }
  else { return(list2)}
  
}

tar <- readRDS(here::here("bw", "tar/bw_tar_seqtab.nochim_Both_RDS"))
fun_otu_tab <- readRDS(here::here("data","samples_focussed_otu_tab.RDS"))
Data <- read.table(here::here("data","SoilTrackerData.txt"), header = TRUE)
tar_urb <- readRDS(here::here("st _urban", "tar/st_urb_tar_seqtab.nochim_Both_RDS"))
urb_colnames <- colnames(tar_urb)

  
tar_f <- otu_filter_function(tar, 50,0,50)


firstrun <- nested.classification(tar_f, tar_urb[23,], Data[,32:43], urb_colnames)

blind_tar <- readRDS(here::here("st_mix", "tar/tar_mix_seqtab.nochim_Both_RDS"))
blind_colnames <- colnames(blind_tar)

# Forest - Beech dominated
firstrun_blind <- nested.classification(tar_f, blind_tar[52,], Data[,32:43], blind_colnames)
firstrun_blind2 <- nested.classification(tar_f, blind_tar[53,], Data[,32:43], blind_colnames)
firstrun_blind3 <- nested.classification(tar_f, blind_tar[54,], Data[,32:43], blind_colnames)

# Grassland (Dwarfshurbs)
blind_mean <- rowMeans(cbind(blind_tar[3,], blind_tar[2,], blind_tar[1,]))
blind_mean_81 <- rowMeans(cbind(blind_tar[25,], blind_tar[26,], blind_tar[27,]))
secondrun_blind <- nested.classification(tar_f, blind_mean, Data[,32:43], blind_colnames)
secondrun_blind2 <- nested.classification(tar_f, blind_mean_81, Data[,32:43], blind_colnames)
  secondrun_blind3 <- nested.classification(tar_f, blind_tar[1,], Data[,32:43], blind_colnames)

# Forest - Alnus dominated 
blind_mean2 <- rowMeans(cbind(blind_tar[55,], blind_tar[56,], blind_tar[57,]))
thirdrun_blind3 <- nested.classification(tar_f, blind_tar[56,], Data[,32:43], blind_colnames)
thirdrun_blind2 <- nested.classification(tar_f, blind_tar[55,], Data[,32:43], blind_colnames)
thirdrun_blind <- nested.classification(tar_f, blind_mean2, Data[,32:43], blind_colnames)

# No Habitat 
fourthrun_blind <- nested.classification(tar_f, blind_tar[82,], Data[,32:43], blind_colnames)
fourthrun_blind2 <- nested.classification(tar_f, blind_tar[83,], Data[,32:43], blind_colnames)
fourthrun_blind3 <- nested.classification(tar_f, blind_tar[84,], Data[,32:43], blind_colnames)

# Agriculture
blind_mean3 <- rowMeans(cbind(blind_tar[109,], blind_tar[110,], blind_tar[111,]))

fifthrun_blind <- nested.classification(tar_f, blind_mean3, Data[,32:43], blind_colnames)
fifthrun_blind2 <- nested.classification(tar_f, blind_tar[110,], Data[,32:43], blind_colnames)
fifthrun_blind3 <- nested.classification(tar_f, blind_tar[111,], Data[,32:43], blind_colnames)

  # 
pla <- readRDS(here::here("bw", "pla/bw_pla_seqtab.nochim_Both_RDS"))
pla_info <- read.table("bw/pla/Plants_ITS2_sample_sites.csv", header = T, sep = ";" )
pla <- prep(pla_info,pla)

pla_f <- otu_filter_function(pla, 50, 0, 20)
blind_pla <- readRDS(here::here("st_mix", "pla/pla_mix_seqtab.nochim_Both_RDS"))
blind_colnames <- colnames(blind_pla)

pla_run <- nested.classification(pla_f, blind_pla[27,], Data[,32:43], blind_colnames)

fun <- readRDS(here::here("bw", "fun/bw_fun_seqtab.nochim_Both_RDS"))
fun_info <- read.table("bw/fun300/Sample_site_replicate_soilFungi.txt", header = T, sep = "\t", stringsAsFactors = F)
fun <- prep(fun_info,fun)

fun_f <- otu_filter_function(fun, 50, 0, 50)
blind_fun <- readRDS(here::here("st_mix", "fun/fun_mix_seqtab.nochim_Both_RDS"))
blind_colnames <- colnames(blind_fun)

fun_run <- nested.classification(fun_f, blind_fun[27,], Data[,32:43], blind_colnames)
