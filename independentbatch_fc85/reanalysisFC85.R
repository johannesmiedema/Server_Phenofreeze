#Analysis of batch fc85 animals - using ML models with phenoFreeze

library(PhenoFreeze)
library(openxlsx)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(emmeans)
library(lme4)
library(mclust)

source("functions.R")
################################################################################
############# CLASSIFICATION ############# 

#MR1 classification
mr1 <- openxlsx::read.xlsx("FC85_F_master table_freezing.xlsx", 
                           sheet = "R1 0.1_1 sec", colNames = T)

#Get freezing data and id´s as rownames
mr1_freezing <- mr1
rownames(mr1_freezing) <- mr1$X40
mr1_freezing <- mr1_freezing[2:nrow(mr1_freezing),16:27]


#Classify freezers and convert factors to character
results_MR1 <- PhenoFreeze::classify_freezer(data_MR1 = mr1_freezing, sex = "female", MR = 1, shifter = F)

mr2 <- openxlsx::read.xlsx("FC85_F_master table_freezing.xlsx", 
                           sheet = "R2 0.1_1 sec", colNames = T)
mr2_freezing <- mr2
rownames(mr2_freezing) <- mr2$X40
mr2_freezing <- mr2_freezing[3:nrow(mr2_freezing),16:27]

#Get only animals which are both in MR1 and MR2 present for shifter model
mr2_filtered <- mr2_freezing[rownames(mr2_freezing) %in% rownames(mr1_freezing),]
mr1_filtered <- mr1_freezing[rownames(mr1_freezing) %in% rownames(mr2_freezing),]

#Check if filtered animals are in correct order for consensus modelling
identical(rownames(mr1$filtered), rownames(mr2$filtered))

#MR2 classification 3 class model 
results_MR2 <- PhenoFreeze::classify_freezer(data_MR1 = mr1_filtered, data_MR2 = mr2_filtered,
                                            sex = "female", MR = 2, shifter = T)

#MR2 classification 2 class consensus model
results_MR2_noshifter <- PhenoFreeze::classify_freezer(data_MR2 = mr2_filtered,
                                                       sex = "female", MR = 2, shifter = F)

results_MR2_consensus <- c()
#Classify shifter using consensus model
for (i in 1:length(results_MR2)){
  if (as.character(results_MR1[i]) == as.character(results_MR2_noshifter[i])){
    results_MR2_consensus[i] <- as.character(results_MR1[i])
  } else {
    results_MR2_consensus[i] <- "shifter"
  }
}

#combine all results
MR1 <- data.frame(id = rownames(mr1_freezing), results_MR1)
MR2 <- data.frame(id = rownames(mr2_freezing), results_MR2)
MR2_consensus <- data.frame(id = rownames(mr2_freezing), results_MR2_consensus)
MR2_noshifter <- data.frame(id = rownames(mr2_freezing), results_MR2_noshifter)

################################################################################
############# Plots of average freezing curves ############# 
#Get freezing values and IDs

#MR1
mr1_freezing_all <- mr1[2:nrow(mr1),4:39]
df_MR1 <- mutate_all(mr1_freezing_all, function(x) as.numeric(as.character(x)))

mr1_animal <- mr1$X40[2:length(mr1$X40)]

#MR2
mr2_freezing_all <- mr2[3:nrow(mr2),4:39]
df_MR2 <- mutate_all(mr2_freezing_all, function(x) as.numeric(as.character(x)))

mr2_animal <- mr2$X40[3:length(mr2$X40)]

#Start PDF to collect all plots
pdf("average_freezing_fc85.pdf",
    height= 4, width = 5)

#MR1 classification results
avg.MR1 <- freezing_curves(df = df_MR1, class = results_MR1,
                           animal = mr1_animal, retrieval = "MR1", stat = TRUE)

avg.MR1 = avg.MR1 + ggtitle("MR1 Classification")
print(avg.MR1)

#MR2 classification results without any shifter
avg.MR2.noshifter <- freezing_curves(df = df_MR2, class = results_MR2_noshifter,
                           animal = mr2_animal, retrieval = "MR2", stat = TRUE)

avg.MR2.noshifter = avg.MR2.noshifter + ggtitle("MR2 Classification")
print(avg.MR2.noshifter)

#MR2 classification results shifter consensus model
avg.MR2.consensus.shifter <- freezing_curves(df = df_MR2, class = results_MR2_consensus,
                                     animal = mr2_animal, retrieval = "MR2", stat = TRUE)

avg.MR2.consensus.shifter = avg.MR2.consensus.shifter + ggtitle("MR2 Classification consensus shifter")
print(avg.MR2.consensus.shifter)

#MR2 classification results shifter 3class model
avg.MR2.3class.shifter <- freezing_curves(df = df_MR2, class = results_MR2,
                                            animal = mr2_animal, retrieval = "MR2", stat = TRUE)

avg.MR2.3class.shifter <- avg.MR2.3class.shifter + ggtitle("MR2 Classification 3 class model")
print(avg.MR2.3class.shifter)

dev.off()

################################################################################
######### Bivariate Plots ############

pdf("bivariate_plots_fc85.pdf", height = 4, width = 3)

#First perform regression analysis manually
data_MR1 <- mr1_freezing

# 12 time bins during tone was played
#Initialize variables
bins <- 1:12
beta_MR1 <- 0
int_MR1 <- 0
#Calculate average freezing of each animal
freeze_MR1 <- rowMeans(data_MR1[,6:12])

#Fit regression model for each animal
for (i in 1:nrow(data_MR1)){
  #transpose x and y
  y <- t(data_MR1[i,])
  #pseudocount if y is zero to allow log operation
  if (0 %in% y){y<-y+1}
  #fit loglinear model
  mod.loglinear <- stats::lm(log(y)~bins)
  #obtain beta coefficients
  beta_MR1[i] <- mod.loglinear$coefficients[2]
  #obtain intercepts
  int_MR1[i] <- mod.loglinear$coefficients[1]
}

#Initialize data frame for regression parameters
params_MR1 <- data.frame(freeze=freeze_MR1, beta=beta_MR1, int=int_MR1)

bi_mr1 <- data.frame(params_MR1, animal = rownames(data_MR1), class = MR1$results_MR1)

bivariate.mr1 <-ggplot(bi_mr1, aes(x=freeze, y = -beta, col = class)) +
  geom_point() + theme_test() + theme(axis.text=element_text(size=12, color="black")) +
  labs(x = "Freezing MR1 bins 18-24 (%)", y ="Decay rate MR1", col="") +
  scale_color_manual(values=c("magenta" ,"green2")) + 
  ggtitle("MR1 Classification") + theme(legend.position = "none") +
  theme(plot.title=element_text(hjust=0.5)) +
  stat_ellipse()


#MR2 regression parameters
data_MR2 <- mr2_freezing
# 12 time bins during tone was played
#Initialize variables
bins <- 1:12
beta_MR2 <- 0
int_MR2 <- 0
#Calculate average freezing of each animal
freeze_MR2 <- rowMeans(data_MR2[,6:12])

#Fit regression model for each animal
for (i in 1:nrow(data_MR2)){
  #transpose x and y
  y <- t(data_MR2[i,])
  #pseudocount if y is zero to allow log operation
  if (0 %in% y){y<-y+1}
  #fit loglinear model
  mod.loglinear <- stats::lm(log(y)~bins)
  #obtain beta coefficients
  beta_MR2[i] <- mod.loglinear$coefficients[2]
  #obtain intercepts
  int_MR2[i] <- mod.loglinear$coefficients[1]
}

#Initialize data frame for regression parameters
params_MR2 <- data.frame(freeze=freeze_MR2, beta=beta_MR2, int=int_MR2)

bi_mr2_noshifter <- data.frame(params_MR2, animal = MR2_noshifter$id, class = results_MR2_noshifter)

bivariate.mr2 <-ggplot(bi_mr2_noshifter, aes(x=freeze, y = -beta, col = class)) +
  geom_point() + theme_test() + theme(axis.text=element_text(size=12, color="black")) +
  labs(x = "Freezing MR2 bins 18-24 (%)", y ="Decay rate MR2", col="") +
  scale_color_manual(values=c("magenta","green2")) + 
  ggtitle("MR2 Classification") + theme(legend.position = "none") +
  theme(plot.title=element_text(hjust=0.5)) +
  stat_ellipse()

bi_mr2_consensus <- data.frame(params_MR2, animal = MR2_consensus$id, class = results_MR2_consensus)

bivariate.mr2.consensus <-ggplot(bi_mr2_consensus, aes(x=freeze, y = -beta, col = class)) +
  geom_point() + theme_test() + theme(axis.text=element_text(size=12, color="black")) +
  labs(x = "Freezing MR2 bins 18-24 (%)", y ="Decay rate MR2", col="") +
  scale_color_manual(values=c("green2","gray50", "magenta")) + 
  ggtitle("MR2 Classification consensus") + theme(legend.position = "none") +
  theme(plot.title=element_text(hjust=0.5)) +
  stat_ellipse()

bi_mr2_3class <- data.frame(params_MR2, animal = MR2$id, class = results_MR2)

bivariate.mr2.3class <-ggplot(bi_mr2_3class, aes(x=freeze, y = -beta, col = class)) +
  geom_point() + theme_test() + theme(axis.text=element_text(size=12, color="black")) +
  labs(x = "Freezing MR2 bins 18-24 (%)", y ="Decay rate MR2", col="") +
  scale_color_manual(values=c("green2","gray50", "magenta")) + 
  ggtitle("MR2 Classification 3class") + theme(legend.position = "none") +
  theme(plot.title=element_text(hjust=0.5)) +
  stat_ellipse()

print(bivariate.mr1)
print(bivariate.mr2)
print(bivariate.mr2.consensus)
print(bivariate.mr2.3class)

dev.off()

#MR1 filtered regression parameters
data <- mr1_filtered
# 12 time bins during tone was played
#Initialize variables
bins <- 1:12
beta <- 0
int <- 0
#Calculate average freezing of each animal
freeze <- rowMeans(data)

#Fit regression model for each animal
for (i in 1:nrow(data)){
  #transpose x and y
  y <- t(data[i,])
  #pseudocount if y is zero to allow log operation
  if (0 %in% y){y<-y+1}
  #fit loglinear model
  mod.loglinear <- stats::lm(log(y)~bins)
  #obtain beta coefficients
  beta[i] <- mod.loglinear$coefficients[2]
  #obtain intercepts
  int[i] <- mod.loglinear$coefficients[1]
}
params.mr1.filtered <- data.frame(freeze, beta, int)
################################################################################
#Comparison of classified bivariate plots against conventionally clustered animals

assign.phenotype <- function(freeze_vector, cluster_vector){
  data <- data.frame(freeze_vector, cluster_vector)
  colnames(data) <- c("freeze", "cluster")
  mean1 <- mean(data[data$cluster==1,]$freeze)
  mean2 <- mean(data[data$cluster==2,]$freeze)
  if (mean1 < mean2){
    data$cluster <-c("phasic", "sustained")[ match( data$cluster, c(1,2))]
    
  } else {
    
    data$cluster <-c("sustained", "phasic")[ match( data$cluster, c(1,2))]
  }
  return(data$cluster)
}

clustered <- Mclust(params_MR1, G = 2)
clustered <- assign.phenotype(params_MR1$freeze, clustered$classification)

clustered.filtered <- Mclust(params.mr1.filtered, G = 2)
clustered.filtered <- assign.phenotype(params.mr1.filtered$freeze, clustered.filtered$classification)

clustered.MR2 <- Mclust(params_MR2, G = 2)
clustered.MR2 <- assign.phenotype(params_MR2$freeze, clustered.MR2$classification)

clustered.MR2.shifter <- c()
for (i in 1:length(clustered.filtered)){
  if (clustered.filtered[i] == clustered.MR2[i]){
    clustered.MR2.shifter[i] <- clustered.filtered[i]
  } else {
    clustered.MR2.shifter[i] <- "shifter"
  }
}

#Plot bivariate comparison MR1
plot <- data.frame(params_MR1, clustered = clustered, classified = MR1$results_MR1)

classified.MR1 <- ggplot(plot, aes(x=freeze, y = -beta, col = classified)) +
  geom_point() + theme_test() + theme(axis.text=element_text(size=12, color="black")) +
  labs(x = "Freezing MR1 bins 18-24 (%)", y ="Decay rate MR1", col="") +
  scale_color_manual(values=c("magenta", "green2")) + 
  ggtitle("MR1 Classified") + theme(legend.position = "none") +
  theme(plot.title=element_text(hjust=0.5)) +
  stat_ellipse()

clustered.MR1 <- ggplot(plot, aes(x=freeze, y = -beta, col = clustered)) +
  geom_point() + theme_test() + theme(axis.text=element_text(size=12, color="black")) +
  labs(x = "Freezing MR1 bins 18-24 (%)", y ="Decay rate MR1", col="") +
  scale_color_manual(values=c("green2", "magenta")) + 
  ggtitle("MR1 Clustered") + theme(legend.position = "none") +
  theme(plot.title=element_text(hjust=0.5)) +
  stat_ellipse()

comparison.MR1 <- ggpubr::ggarrange(clustered.MR1, classified.MR1)

#Plot bivariate comparison MR2 without shifter
plot.MR2 <- data.frame(params_MR2, clustered = clustered.MR2, classified = MR2_noshifter$results_MR2_noshifter)

bi.classified.MR2 <- ggplot(plot.MR2, aes(x=freeze, y = -beta, col = classified)) +
  geom_point() + theme_test() + theme(axis.text=element_text(size=12, color="black")) +
  labs(x = "Freezing MR2 bins 18-24 (%)", y ="Decay rate MR2", col="") +
  scale_color_manual(values=c("magenta", "green2")) + 
  ggtitle("MR2 Classified") + theme(legend.position = "none") +
  theme(plot.title=element_text(hjust=0.5)) +
  stat_ellipse()

bi.clustered.MR2 <- ggplot(plot.MR2, aes(x=freeze, y = -beta, col = clustered)) +
  geom_point() + theme_test() + theme(axis.text=element_text(size=12, color="black")) +
  labs(x = "Freezing MR2 bins 18-24 (%)", y ="Decay rate MR2", col="") +
  scale_color_manual(values=c("green2", "magenta")) + 
  ggtitle("MR2 Clustered") + theme(legend.position = "none") +
  theme(plot.title=element_text(hjust=0.5)) +
  stat_ellipse()

comparison.MR2 <- ggpubr::ggarrange(bi.clustered.MR2, bi.classified.MR2)

#Plot bivariate comparison MR2 with shifter consensus
plot.MR2.consensus <- data.frame(params_MR2, clustered = clustered.MR2.shifter, classified = MR2_consensus$results_MR2_consensus)

bi.classified.MR2.consensus <- ggplot(plot.MR2.consensus, aes(x=freeze, y = -beta, col = classified)) +
  geom_point() + theme_test() + theme(axis.text=element_text(size=12, color="black")) +
  labs(x = "Freezing MR2 bins 18-24 (%)", y ="Decay rate MR2", col="") +
  scale_color_manual(values=c("green2", "gray38", "magenta")) + 
  ggtitle("MR2 Classified consensus") + theme(legend.position = "none") +
  theme(plot.title=element_text(hjust=0.5)) +
  stat_ellipse()

bi.clustered.MR2.consensus <- ggplot(plot.MR2.consensus, aes(x=freeze, y = -beta, col = clustered)) +
  geom_point() + theme_test() + theme(axis.text=element_text(size=12, color="black")) +
  labs(x = "Freezing MR2 bins 18-24 (%)", y ="Decay rate MR2", col="") +
  scale_color_manual(values=c("green2", "gray38", "magenta")) + 
  ggtitle("MR2 Clustered consensus") + theme(legend.position = "none") +
  theme(plot.title=element_text(hjust=0.5)) +
  stat_ellipse()

comparison.MR2.consensus <- ggpubr::ggarrange(bi.clustered.MR2.consensus, bi.classified.MR2.consensus)

#Plot bivariate comparison MR2 with shifter 3class model
plot.MR2.3class <- data.frame(params_MR2, clustered = clustered.MR2.shifter, classified = MR2$results_MR2)

bi.classified.MR2.3class <- ggplot(plot.MR2.3class, aes(x=freeze, y = -beta, col = classified)) +
  geom_point() + theme_test() + theme(axis.text=element_text(size=12, color="black")) +
  labs(x = "Freezing MR2 bins 18-24 (%)", y ="Decay rate MR2", col="") +
  scale_color_manual(values=c("green2", "gray38", "magenta")) + 
  ggtitle("MR2 Classified 3class")  + theme(legend.position = "none") +
  theme(plot.title=element_text(hjust=0.5)) +
  stat_ellipse()

bi.clustered.MR2.consensus <- ggplot(plot.MR2.3class, aes(x=freeze, y = -beta, col = clustered)) +
  geom_point() + theme_test() + theme(axis.text=element_text(size=12, color="black")) +
  labs(x = "Freezing MR2 bins 18-24 (%)", y ="Decay rate MR2", col="") +
  scale_color_manual(values=c("green2", "gray38", "magenta")) + 
  ggtitle("MR2 Clustered consensus")  + theme(legend.position = "none") +
  theme(plot.title=element_text(hjust=0.5)) +
  stat_ellipse()

comparison.MR2.3class <- ggpubr::ggarrange(bi.clustered.MR2.consensus, bi.classified.MR2.3class)

pdf("FC85_Clustered_vs_Classified.pdf", width = 7, height = 4)
print(comparison.MR1)
print(comparison.MR2)
print(comparison.MR2.3class)
print(comparison.MR2.consensus)
dev.off()

################################################################################
#Calculate Concordance measure of Classified ~ Clustered for MR1

#Prepare confusion matrix
data <- data.frame(Clustered = clustered, Classified = MR1$results_MR1)
data$Clustered <- factor(data$Clustered, levels = c("phasic", "sustained"))
data$Classified <- factor(data$Classified, levels = c("phasic", "sustained"))
cf_matrix <- table(data)

#Calculate Cohen´s Kappa
DescTools::CohenKappa(cf_matrix, conf.level = 0.95)

#Visualize results

cm <- caret::confusionMatrix(data$Classified, data$Clustered, dnn = c("Prediction", "Reference"))

plt <- as.data.frame(cm$table)
plt$Prediction <- factor(plt$Prediction, levels=rev(levels(plt$Prediction)))

plot <- ggplot(plt, aes(Prediction,Reference, fill= Freq)) +
        geom_tile(color = "black",
                  lwd = 0.4,
                  linetype = 1) +
        geom_text(aes(label=Freq)) +
        scale_fill_gradient(low="white", high="blue") +
        labs(x = "Clustered",y = "Classified") + theme_minimal() +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
        ggtitle("Cohen`s Kappa = 0.8872")

pdf("CF_matrix_MR1.pdf", width = 3.5, height = 2.5)
print(plot)
dev.off()

################################################################################
#Calculate Concordance measure of Classified ~ Clustered for MR2

#Prepare confusion matrix
data <- data.frame(Clustered = clustered.MR2, Classified = MR2_noshifter$results_MR2_noshifter)
data$Clustered <- factor(data$Clustered, levels = c("phasic", "sustained"))
data$Classified <- factor(data$Classified, levels = c("phasic", "sustained"))
cf_matrix <- table(data)

#Calculate Cohen´s Kappa
DescTools::CohenKappa(cf_matrix, conf.level = 0.95)

#Visualize results

cm <- caret::confusionMatrix(data$Classified, data$Clustered, dnn = c("Prediction", "Reference"))

plt <- as.data.frame(cm$table)
plt$Prediction <- factor(plt$Prediction, levels=rev(levels(plt$Prediction)))

plot <- ggplot(plt, aes(Prediction,Reference, fill= Freq)) +
  geom_tile(color = "black",
            lwd = 0.4,
            linetype = 1) +
  geom_text(aes(label=Freq)) +
  scale_fill_gradient(low="white", high="blue") +
  labs(x = "Clustered",y = "Classified") + theme_minimal() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ggtitle("Cohen`s Kappa = 0.8852")

pdf("CF_matrix_MR2.pdf", width = 3.5, height = 2.5)
print(plot)
dev.off()

################################################################################
#Subsampling and then plotting Clustered VS Classified bootstrap accuracy compared to initial clustering
s.e. <- function(x)sd(x)/sqrt(length(x))
#-----------------------------------
#----------------MR1----------------
set.seed(42)
groundtruth <- mclust::Mclust(params_MR1, G=2)
groundtruth <- assign.phenotype(params_MR1$freeze, groundtruth$classification)

params_MR1$groundtruth <- groundtruth
mr1_freezing$groundtruth <- groundtruth

mean.classified <- c()
se.classified <- c()
mean.clustered <- c()
se.clustered <- c()

animal_numbers <- c(38:3)
iterations <- 1000

for (animal in 1:length(animal_numbers)){
  Accuracy.clust <- c()
  Accuracy.classified <- c()
  for (i in 1:iterations){
    print(as.numeric(animal_numbers[animal]))
    print(i)
    #draw random animals from dataframe
    samplenum <- sample(nrow(params_MR1), as.numeric(animal_numbers[animal]))
    
    #Subset each df
    current_df <- params_MR1[samplenum, ]
    classify_df <- mr1_freezing[samplenum,]
    
    #Cluster the animals using mclust
    clust.current <- invisible(mclust::Mclust(current_df[,1:3], G=2))
    cluster_res <- assign.phenotype(current_df$freeze, clust.current$classification)
    current_df$clustered <- cluster_res
    
    #Classify animals with phenoFreeze
    classify_res <- PhenoFreeze::classify_freezer(data_MR1 = classify_df[,1:12], sex = "female", MR = 1, shifter = F)
    classify_df$classified <- classify_res
    
    #Declare phenotypes as factors
    current_df$groundtruth <- factor(current_df$groundtruth, levels = c("phasic", "sustained"))
    current_df$clustered <- factor(current_df$clustered, levels = c("phasic", "sustained"))
    
    classify_df$groundtruth <- factor(current_df$groundtruth, levels = c("phasic", "sustained"))
    classify_df$classified <- factor(classify_df$classified, levels = c("phasic", "sustained"))
    
    #Calculate Accuracy of clustered vs original clusterdataset
    cm <- table(current_df$clustered, current_df$groundtruth)
    cmMatrix <- caret::confusionMatrix(cm, mode="everything")
    Accuracy.clust[i] <- cmMatrix$overall["Accuracy"]
    print(cmMatrix$overall["Accuracy"])
    
    cm.classified <- table(classify_df$classified, classify_df$groundtruth)
    cmMatrix.classified <- caret::confusionMatrix(cm.classified, mode="everything")
    Accuracy.classified[i] <- cmMatrix.classified$overall["Accuracy"]
  }
  #Calculate means and standard errors
  mean.clustered[animal] <- mean(Accuracy.clust)
  se.clustered[animal] <- s.e.(Accuracy.clust)
  mean.classified[animal] <- mean(Accuracy.classified)
  se.classified[animal] <- s.e.(Accuracy.classified)
}

clustered.data <- data.frame(Average.accuracy=mean.clustered, se=se.clustered, Subsample.size=animal_numbers)
clustered.data$prediction <- "Clustering approach"

classified.data <- data.frame(Average.accuracy=mean.classified, se=se.classified, Subsample.size=animal_numbers)
classified.data$prediction <- "ML-Model"

data <- rbind(clustered.data, classified.data)
data$prediction <- factor(data$prediction, levels = c("Clustering approach", "ML-Model"))

lineplot <- ggplot(data, aes(x=Subsample.size, y=Average.accuracy, group=prediction, color=prediction)) + 
  geom_errorbar(aes(ymin=Average.accuracy-se, ymax=Average.accuracy+se), width=.2) +
  geom_line() + geom_point()+ theme_test(base_size=14) + theme(axis.text = element_text(colour="black", size=12)) + 
  ylab("Average Bootstrap Accuracy") + xlab("Subsample Size") + scale_color_manual(values = c("#005f73", "#fb8500")) +
  scale_x_continuous(trans='reverse', breaks = c(38, 35, 30, 25, 20, 15, 10, 5)) + theme(legend.title = element_blank())

lineplot

pdf("../Figures/Downsampling/Cluster.vs.ML.fc85.Downsampling.pdf", width = 6, height = 4)
print(lineplot)
dev.off()
