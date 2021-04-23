library(here)
library(dplyr)
library(caret)
library(stringr)
library(vegan)
library(cvms)
library(ggrepel)
source("model_input_function.R")
source("otu_filter_function.R")



# Adjusting classes for the first round of the model: 
# classification on the big three classes (Agriculture, Forest and Grassland) 
fi <- Data[,32:43]
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
testData_colnames <- blind_colnames
testData <- blind_mean
testData_frame <- as.data.frame(testData, col.names = testData_colnames)
testData_t <- t(testData_frame)
testData_f <- otu_filter_function(testData_t, 0, 0, 0)
# Relative abundance normalisation for the test Data 
testData_rel <- sweep(testData_f, 1, rowSums(testData_f), "/")
testData_rel_f<- as.data.frame(testData_rel^0.5)
  
# Indexing and finding matching column name (sequence tags) for classification of the test Data. 
trainingData <- otu_filter_function(tar, 50, 0, 50)
trainingData <- model_input(trainingData, "testClass", fi, "relative abundance")
index <- which(colnames(trainingData) %in% colnames(testData_rel_f))
index2 <- which(colnames(testData_rel_f) %in% colnames(trainingData))
  
tar_new <- trainingData[,index]
tar_test <- testData_rel_f[,index2]
 
# Knn classification based on Bray-Curtis dissimilarity 
brayknn <- function(trainingData, testData, classData, k){
  
  # Bray-Curtis dissimilarity formula 
  bcd <- function(f,p){ 
    1-2*sum(pmin(f, p))/ (sum(f) + sum(p))}
  
  # loop to get every dissimilarity 
  dist <- NULL
  for (i in 1:nrow(tar_new)){
      dist[i] <- bcd(tar_new[i,], tar_test[1,])}
 
 tar_new_w <- cbind(tar_new, fi$testClass)
    
 knn <- order(dist)[1:k]
 c <- rownames(tar_new_w)[knn]
 
 res <- NULL
 for (i in c){
 res[i] <- tar_new_w[which(rownames(tar_new_w) == i), ncol(tar_new_w)]} 
 
}
 

bc_knn_pred <- function(testData, trainingData, k){
  pred <- c()
  
    bc_dis <- c()
    bc_char <- c()
    good <- 0
    bad <- 0 
    
    for (j in c(1:nrow(trainingData))){
      
      bc_dis <- c(bc_dis, bcd(trainingData[j,-131], testData))
      
      bc_char <- c(bc_char, as.character(trainingData[j,][[nrow(trainingData)]]))
    }
    
    bc <- data.frame(bc_char, bc_dis)
    bc <- bc[order(bc$bc_dis),]
    bc <- bc[1:k,]
    
    for (k in c(1:nrow(bc))){
      if(as.character(bc[k, "bc_char"]) == "g"){
        good = good + 1 
      }
      else 
        bad = bad + 1
    }
    
    if (good > bad){
      pred <- c(pred, "g")
    }
    else if(good < bad){
      pred <- c(pred, "b")
    }
  
  return(pred)
}

tar_new_ww <- tar_new_w
colnames(tar_new_ww) <- NULL

bc_knn_pred(tar_test[1,], tar_new_ww, 10)

 
 
  
  
 
 
 for (i in 1:130){
   best.qda <- qda(as.matrix(Data[-i,paste(unlist(strsplit(pred.vector, split=" ")))]), Data$Forest[-i])
   Pred <- predict(best.qda, newdata=Data[i,paste(unlist(strsplit(pred.vector, split=" ")))])
   PredictionBi$Foclass[i] <- Pred$class
   PredictionBi$Fopost[i] <- Pred$posterior[2]
 }
 
 print(res)
 res2 <-sort(table(res), decreasing = T)
 
  
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
  unknown.tar <- vegdist(tar_new_plus[,1:length(colnames(tar_new))], distance = "bray")
  unknown.tar <- as.matrix(unknown.tar)
  
  #CAP analysis approach 
  testset <- cbind(testClass = tar_new_plus$testClass, tar_new_plus[,1:length(colnames(tar_new))])
  pred <-  tar_new_plus[,1:length(colnames(tar_new))]
  
  test2 <- BiodiversityR::CAPdiscrim(pred ~ testClass, data = testset, dist = "bray",
                                     axes = 6)
  nmds.data.test <- cbind.data.frame(test2$PCoA, tar_new_plus$testClass)
  colnames(nmds.data.test) <- c("X", "Y", "X1", "Y1", "X2", "Y2", "color")
  
  NMDS_all <- ggplot()+
    geom_point(data = nmds.data.test, aes(x = X, y = Y, group = color, shape = color, color = color, size = color)) +
    scale_shape_manual(name = "Habitat classes", values = c("Agriculture" = 16, "Forest" = 16, "Grassland" = 16, "No Habitat" = 16)) +
    scale_colour_manual(name = "Habitat classes", values = c("Agriculture" = "#D16103",
                                                             "Forest" = "#52854C",
                                                             "Grassland" = "#4E84C4", 
                                                             "No Habitat" = "#C4961A")) +
    scale_size_manual(name = "Habitat classes", values = c("Agriculture" = 4, "Forest" = 4, "Grassland" = 4, "No Habitat" = 4)) +
    #geom_point(data = centroids, aes(x = MDS1, y = MDS2, color = centroids$testClass), size = 5, shape = 15) +
    theme_bw() +
    ylab("PCoA 2") +
    xlab("PCoA 1") +
    theme(axis.title = element_text(size = 18)) +
    theme(axis.text = element_text(size = 18)) +
    theme(legend.title = element_text(size = 23)) +
    theme(legend.text = element_text(size = 23)) 
    #geom_text_repel(data = nmds.data.test[131,], aes(x = X, y = Y, label = rownames(nmds.data.test), size = 5, 
                    #box.padding = unit(0.35, "lines"), point.padding = unit(0.3, "lines")))
  
  print(NMDS_all)

  
  # Calculating the centroids of the 4 classes
  class.dist <- cbind.data.frame(unknown.tar, tar_new_plus$testClass)
  colnames(class.dist)[132] <- "testClass"
                               
  centroids <- aggregate(. ~ testClass, data = class.dist, mean)
  centroids_raw <- aggregate(. ~ testClass, data = tar_new_plus[,c(1:length(colnames(tar_new)),ncol(tar_new_plus))], mean)
  dim(centroids_raw)
  mat_all <- rbind.data.frame(centroids_raw[,-1], tar_new_plus[,1:length(colnames(tar_new))])
  dis_all <- vegdist(mat_all, distance = "bray")
  dis_all <- as.matrix(dis_all)

  disMatrix <- as.data.frame(dis_all[5:135,1:4])
  disMatrix <- cbind(disMatrix, tar_new_plus$testClass)
  colnames(disMatrix) <- c("Agriculture", "Forest", "Grassland", "No Habitat", "Class")
  
  tar_model2 <- train(Class ~ ., data = disMatrix, method = "knn", metric = "Kappa", 
                      trControl = trainControl("LOOCV"),
                      tuneGrid = expand.grid(k = c(3:15)))
  
  
  clSubset <- subset(disMatrix[,which(colnames(disMatrix) == pred[1])], disMatrix$Class == pred[1])
  mean_F <- mean(clSubset)
  sd_F <- sd(clSubset)
  
  ncSubsetwo <- subset(disMatrix[,c(which(colnames(disMatrix) == pred[1]))], disMatrix$Class != pred[1])
  mean_nF <- mean(ncSubsetwo)
  sd_nF <- sd(ncSubsetwo)
  
  unknown_dist <- disMatrix[nrow(disMatrix),which(colnames(disMatrix) == pred[1])]
  
  # Probability of the unknown_dist begin from the distribution (from the wanted centroid class)
  A_B <- dnorm(unknown_dist, mean_F, sd_F) 
  
  # Probability for event 2 the unknown samples does not come from the predicted class
  A_B2 <-  dnorm(unknown_dist, mean_nF, sd_nF)
  
  prior <- 0.25 #The chosen prior here is 0.25 based on the four centroids. This suggest that each unknown samples's distance can be equally belong to one of the four classes.
  
  # Likelihood function for event A 
  Bayestest <- A_B * prior 
  
  #Likelihood function for event B 
  Bayestest2 <- A_B2 * prior 
  
  LR <- Bayestest / Bayestest2 
  print(LR)
  
  
    p1 <- ggplot()+ 
    geom_histogram(aes(x = clSubset,y = ..density.., fill = "In class distance"), binwidth = 0.1, alpha = .8, position = "identity") + 
    stat_function(mapping = NULL, fun = dnorm, args = list(mean = mean(clSubset), sd = sd(clSubset)), color = "red") +
    geom_histogram(aes(x = ncSubsetwo, y = ..density.., fill = "Out of class distance"), binwidth = 0.1,alpha = .8, position = "identity") +
    stat_function(mapping = NULL ,fun = dnorm, args = list(mean = mean(ncSubsetwo), sd = sd(ncSubsetwo)), color = "blue") +
    labs(fill = "Distributions") +
    ylab("Density")+
    xlab("Histogram of in-class distances and out-of class distances") + 
    theme_bw() +
    geom_vline(xintercept = unknown_dist, linetype = "dotted", color = "black") +
    geom_text(aes(x = unknown_dist + 0.03, label = "distance of the unknown test sample", y = 1), color = "black", 
              angle = 90) +
    theme(axis.title = element_text(size = 18)) +
    theme(axis.text = element_text(size = 18)) +
    theme(legend.text = element_text(size = 23)) +
    theme(legend.title = element_text(size = 23))
  
  print(p1)
  
  # Adding the best matching site based on the multidimensionl euclidean distance between the points (comparison of picture)
  
  dist1 <- vegdist(tar_new_plus[,1:length(colnames(tar_new))], method = "bray")
  dist1 <- as.matrix(dist1)
  dist1[,131]
  best.ind <- which.min(dist1[1:130,131])
  best.match <- rownames(tar_new_plus)[best.ind]
  
  print(best.match)
  
  list1 <- list(tar_model, pred, testPred,knn1,CMP1, NMDS_all, LR, p1, best.match)
  
  
  return(list1)}
  


tar <- readRDS(here::here("bw", "tar/bw_tar_seqtab.nochim_Both_RDS"))
fun_otu_tab <- readRDS(here::here("data","samples_focussed_otu_tab.RDS"))
Data <- read.table(here::here("data","SoilTrackerData.txt"), header = TRUE)
tar_urb <- readRDS(here::here("st _urban", "tar/st_urb_tar_seqtab.nochim_Both_RDS"))
urb_colnames <- colnames(tar_urb)


tar_f <- otu_filter_function(tar, 50,0,50)


firstrun <- nested.classification2(tar_f, tar_urb[23,], Data[,32:43], urb_colnames)

blind_tar <- readRDS(here::here("st_mix", "tar/tar_mix_seqtab.nochim_Both_RDS"))
blind_colnames <- colnames(blind_tar)

# Forest - Beech dominated
firstrun_blind <- nested.classification(tar_f, blind_tar[52,], Data[,32:43], blind_colnames)
firstrun_blind2 <- nested.classification(tar_f, blind_tar[53,], Data[,32:43], blind_colnames)
firstrun_blind3 <- nested.classification(tar_f, blind_tar[54,], Data[,32:43], blind_colnames)

# Grassland (Dwarfshurbs)
blind_mean <- rowMeans(cbind(blind_tar[3,], blind_tar[2,], blind_tar[1,]))
blind_mean_81 <- rowMeans(cbind(blind_tar[25,], blind_tar[26,], blind_tar[27,]))
secondrun_blind <- nested.classification2(tar_f, blind_mean, Data[,32:43], blind_colnames)
secondrun_blind2 <- nested.classification2(tar_f, blind_mean_81, Data[,32:43], blind_colnames)


# Forest - Alnus dominated 
blind_mean2 <- rowMeans(cbind(blind_tar[55,], blind_tar[56,], blind_tar[57,]))
thirdrun_blind <- nested.classification2(tar_f, blind_mean2, Data[,32:43], blind_colnames)

# Agriculture
blind_mean3 <- rowMeans(cbind(blind_tar[109,], blind_tar[110,], blind_tar[111,]))

fifthrun_blind <- nested.classification2(tar_f, blind_mean3, Data[,32:43], blind_colnames)


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
