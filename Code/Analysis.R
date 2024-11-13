################################################################################
#PhenoFreeze Manuscript analysis code
#Clustering and analysis of clustering results

#Load required libraries 
library(openxlsx)
library(mclust)
library(emmeans)
library(ggplot2)
library(PhenoFreeze)

#Source function script
source("Functions.R")

#Set seed for reproducibility
set.seed(42)

#lmer and emmeans error caused by exceeding observation limit of 3000
#Increase them manually
emm_options(pbkrtest.limit = 14112)

################################################################################
#Data preprocessing
females_MR1 <- openxlsx::read.xlsx("../Data/all_females_MR1.xlsx")
females_MR2 <- openxlsx::read.xlsx("../Data/all_females_MR2.xlsx")

males_MR1 <- openxlsx::read.xlsx("../Data/all_males_MR1.xlsx")

#Two different Protocols for different batches
#Normal protocol with CS 13-24 bins
males_MR2 <- openxlsx::read.xlsx("../Data/males_MR2.xlsx")

#Altered protocol with CS 5-16 bins
males_MR2_diff <- openxlsx::read.xlsx("../Data/males_differentCS.MR2.xlsx")


#---------------------------------------
#Merge Batch and IDs to construct unique IDs

#Females MR1
females_MR1.id <- c()
for (i in 1:nrow(females_MR1)){females_MR1.id[i] <- paste(females_MR1$batch[i], 
                                                          females_MR1$id[i],
                                                          sep = "-")}
rownames(females_MR1) <- females_MR1.id

#Females MR2
females_MR2.id <- c()
for (i in 1:nrow(females_MR2)){females_MR2.id[i] <- paste(females_MR2$batch[i], 
                                                          females_MR2$id[i],
                                                          sep = "-")}
rownames(females_MR2) <- females_MR2.id

#Males MR1
males_MR1.id <- c()
for (i in 1:nrow(males_MR1)){males_MR1.id[i] <- paste(males_MR1$batch[i], 
                                                          males_MR1$id[i],
                                                          sep = "-")}
rownames(males_MR1) <- males_MR1.id

#Males MR2 normal CS protocol
males_MR2.id <- c()
for (i in 1:nrow(males_MR2)){males_MR2.id[i] <- paste(males_MR2$batch[i], 
                                                      males_MR2$id[i],
                                                      sep = "-")}
rownames(males_MR2) <- males_MR2.id

#Males MR2 - different CS protocol
males_MR2_diff.id <- c()
for (i in 1:nrow(males_MR2_diff)){males_MR2_diff.id[i] <- paste(males_MR2_diff$batch[i], 
                                                      males_MR2_diff$id[i],
                                                      sep = "-")}
rownames(males_MR2_diff) <- males_MR2_diff.id

################################################################################
#Perform loglinear regression for bins 13:24 (first two columns are batch and ID

#Females MR1
females_MR1_param <- perform_loglinear(females_MR1[,15:26])

#Add batch and id information
females_MR1_param$id <- females_MR1$id
females_MR1_param$batch <- females_MR1$batch

#Females MR2
females_MR2_param <- perform_loglinear(females_MR2[,15:26])

#Add batch and id information
females_MR2_param$id <- females_MR2$id
females_MR2_param$batch <- females_MR2$batch

#Males MR1
males_MR1_param <- perform_loglinear(males_MR1[,15:26])

#Add batch and id information
males_MR1_param$id <- males_MR1$id
males_MR1_param$batch <- males_MR1$batch

#Males MR2 
#Merge CS bins from both protocols
CS.normal <- males_MR2[,c(1,2,15:26)]
CS.diff <- males_MR2_diff[,c(2,1,7:18)]

#Rename CS bins to 1-12 for both protocols
colnames(CS.normal)[3:14] <- c(1:12)
colnames(CS.diff)[3:14] <- c(1:12)
males_MR2_CS <- rbind(CS.normal, CS.diff)

#Perform loglinear regression for males MR2 for both protocols merged together
males_MR2_CS_param <- perform_loglinear(males_MR2_CS[,3:14])
males_MR2_CS_param$id <- males_MR2_CS$id
males_MR2_CS_param$batch <- males_MR2_CS$batch

################################################################################
#Perform two component Gaussian Mixture Model clustering
#Using average freezing of bins 6-12 during CS and beta and int parameters from
#The loglinear regression model

#Females MR1
mclust.f.MR1 <- mclust::Mclust(females_MR1_param[,1:3], G=2)
females_MR1_param$class <- assign.phenotype(females_MR1_param$freeze,
                                            mclust.f.MR1$classification)

#Females MR2
mclust.f.MR2 <- mclust::Mclust(females_MR2_param[,1:3], G=2)
females_MR2_param$class <- assign.phenotype(females_MR2_param$freeze,
                                            mclust.f.MR2$classification)

#Males MR1
mclust.m.MR1 <- mclust::Mclust(males_MR1_param[,1:3], G=2)
males_MR1_param$class <- assign.phenotype(males_MR1_param$freeze,
                                          mclust.m.MR1$classification)

#Males MR2
mclust.m.MR2 <- mclust::Mclust(males_MR2_CS_param[,1:3], G=2)
males_MR2_CS_param$class <- assign.phenotype(males_MR2_CS_param$freeze,
                                          mclust.m.MR2$classification)

#----------------------------------------------------
#Identify Shifter between MR1 and MR2
#First, create unique combination of batch and id as there are duplicate IDs for 
#some batches

#Female MR1
f.mr1.id <- c()
for (i in 1:nrow(females_MR1_param)){f.mr1.id[i] <- paste(females_MR1_param$batch[i], 
                                                          females_MR1_param$id[i],
                                                          sep = "-")}
females_MR1_param$B_ID <- f.mr1.id

#Reorganize dataframe 
females_MR1_param <- females_MR1_param[,c(8,1,2,3,7)]

#Female MR2
f.mr2.id <- c()
for (i in 1:nrow(females_MR2_param)){f.mr2.id[i] <- paste(females_MR2_param$batch[i], 
                                                          females_MR2_param$id[i],
                                                          sep = "-")}
females_MR2_param$B_ID <- f.mr2.id

#Reorganize dataframe 
females_MR2_param <- females_MR2_param[,c(8,1,2,3,7)]

#Male MR1
m.mr1.id <- c()
for (i in 1:nrow(males_MR1_param)){m.mr1.id[i] <- paste(males_MR1_param$batch[i], 
                                                        males_MR1_param$id[i],
                                                        sep = "-")}
males_MR1_param$B_ID <- m.mr1.id

#Reorganize dataframe 
males_MR1_param <- males_MR1_param[,c(8,1,2,3,7)]

#Male MR2
m.mr2.id <- c()
for (i in 1:nrow(males_MR2_CS_param)){m.mr2.id[i] <- paste(males_MR2_CS_param$batch[i], 
                                                          males_MR2_CS_param$id[i],
                                                          sep = "-")}
males_MR2_CS_param$B_ID <- m.mr2.id

#Reorganize dataframe 
males_MR2_CS_param <- males_MR2_CS_param[,c(8,1,2,3,7)]

#Merge MR1 and MR2 param dataframes
females.shifter <- merge(females_MR1_param, females_MR2_param, by = "B_ID")
males.shifter <- merge(males_MR1_param, males_MR2_CS_param, by = "B_ID")

#Determine shifter for females
shifter.f <- c()
for (i in 1:nrow(females.shifter)){
  if(females.shifter$class.x[i] == females.shifter$class.y[i]){
    shifter.f[i] <- females.shifter$class.y[i]
  } else {
    shifter.f[i] <- "shifter"
  }
}

females.shifter$class <- shifter.f

#Determine shifter for males
shifter.m <- c()
for (i in 1:nrow(males.shifter)){
  if(males.shifter$class.x[i] == males.shifter$class.y[i]){
    shifter.m[i] <- males.shifter$class.y[i]
  } else {
    shifter.m[i] <- "shifter"
  }
}

males.shifter$class <- shifter.m


#----------------------------------------------------
#Identification of shifter subtypes S-P and P-S

#Determine shifter for females
shifter.f.sub <- c()
for (i in 1:nrow(females.shifter)){
  if(females.shifter$class.x[i] == "phasic" && females.shifter$class.y[i] == "sustained"){
    shifter.f.sub[i] <- "P-S"
  }
  if(females.shifter$class.x[i] == "sustained" && females.shifter$class.y[i] == "phasic"){
    shifter.f.sub[i] <- "S-P"
  }
  if(females.shifter$class.x[i] == females.shifter$class.y[i]){
    shifter.f.sub[i] <- females.shifter$class.x[i]
  }
}

females.shifter$subclass <- shifter.f.sub

#Determine shifter for males
shifter.m.sub <- c()
for (i in 1:nrow(males.shifter)){
  if(males.shifter$class.x[i] == "phasic" && males.shifter$class.y[i] == "sustained"){
    shifter.m.sub[i] <- "P-S"
  }
  if(males.shifter$class.x[i] == "sustained" && males.shifter$class.y[i] == "phasic"){
    shifter.m.sub[i] <- "S-P"
  }
  if(males.shifter$class.x[i] == males.shifter$class.y[i]){
    shifter.m.sub[i] <- males.shifter$class.x[i]
  }
}


males.shifter$subclass <- shifter.m.sub

#----------------------------------------------------
#Save labelled datasets for model training 

write.xlsx(females_MR1_param, "../Data/labelled.f.MR1.xlsx")
write.xlsx(females_MR2_param, "../Data/labelled.f.MR2.xlsx")
write.xlsx(males_MR1_param, "../Data/labelled.m.MR1.xlsx")
write.xlsx(males_MR2_CS_param, "../Data/labelled.m.MR2.xlsx")
write.xlsx(females.shifter, "../Data/labelled.f.shifter.xlsx")
write.xlsx(males.shifter, "../Data/labelled.m.shifter.xlsx")

################################################################################
#Plot average freezing curves as well as individual freezing curves

#----------------------------------------------------
#Actual Freezing curves

#Females MR1
actual.freezing.f.MR1 <- freezing_curves(df = females_MR1[,3:38],
                                         class = females_MR1_param$class, 
                                         animal = females_MR1$id,
                                         retrieval = "MR1",
                                         stat = T, individual = TRUE)


pdf("../Figures/Freezing/Figure.1.freezing.f.MR1.pdf", height = 4, width = 5)
print(actual.freezing.f.MR1)
dev.off()

#Females MR2
actual.freezing.f.MR2 <- freezing_curves(df = females_MR2[,3:38],
                                         class = females_MR2_param$class, 
                                         animal = females_MR2$id,
                                         retrieval = "MR2", stat = T,
                                         individual = TRUE)

pdf("../Figures/Freezing/Figure.2.freezing.f.MR2.pdf", height = 4, width = 5)
print(actual.freezing.f.MR2)
dev.off()

#Males MR1
actual.freezing.m.MR1 <- freezing_curves(df = males_MR1[,3:38],
                                         class = males_MR1_param$class, 
                                         animal = males_MR1$id,
                                         retrieval = "MR1", stat = T, 
                                         individual = TRUE, shape = 15)


pdf("../Figures/Freezing/Figure.3.freezing.m.MR1.pdf", height = 4, width = 5)
print(actual.freezing.m.MR1)
dev.off()

#Males MR2 - merge different protocols of CS of different mails 

#Reformat batches so the columns match
males_MR2_diff_formatted <- data.frame(id = males_MR2_diff$id, batch = males_MR2_diff$batch,
                             "1" = NA, "2" = NA, "3" = NA, "4" = NA, "5" = NA,
                             "6" = NA, "7" = NA, "8" = NA, males_MR2_diff[,3:24],
                            "31" = NA, "32" = NA, "33" = NA,
                             "34" = NA, "35" = NA, "36" = NA)

colnames(males_MR2_diff_formatted) <- c("id", "batch", 1:36)

males_MR2_freezing <- rbind(males_MR2, males_MR2_diff_formatted)


actual.freezing.m.MR2 <- freezing_curves(df = males_MR2_freezing[,3:38],
                                         class = males_MR2_CS_param$class, 
                                         animal = males_MR2_freezing$id,
                                         retrieval = "MR2", stat = T,
                                         individual = TRUE, shape = 15)

pdf("../Figures/Freezing/Figure.4.freezing.m.MR2.pdf", height = 4, width = 5)
print(actual.freezing.m.MR2)
dev.off()

#----------------------------------------------------
#Fitted Freezing curves

#Females MR1
fitted.freezing.f.MR1 <- fitted_freezing_curves(df = females_MR1[,15:26],
                                                class = females_MR1_param$class,
                                                retrieval = "MR1", stat = T,
                                                individual = TRUE)

pdf("../Figures/Freezing/Figure.5.fitted.f.MR1.pdf", height = 4, width = 5)
print(fitted.freezing.f.MR1)
dev.off()

#Females MR2
fitted.freezing.f.MR2 <- fitted_freezing_curves(df = females_MR2[,15:26],
                                                class = females_MR2_param$class,
                                                retrieval = "MR2", stat = T,
                                                individual = TRUE)

pdf("../Figures/Freezing/Figure.6.fitted.f.MR2.pdf", height = 4, width = 5)
print(fitted.freezing.f.MR2)
dev.off()

#Males MR1
fitted.freezing.m.MR1 <- fitted_freezing_curves(df = males_MR1[,15:26],
                                                class = males_MR1_param$class,
                                                retrieval = "MR1", stat = T,
                                                individual = TRUE)

pdf("../Figures/Freezing/Figure.7.fitted.m.MR1.pdf", height = 4, width = 5)
print(fitted.freezing.m.MR1)
dev.off()

#Males MR2 - use both protocols for plotting
fitted.freezing.m.MR2 <- fitted_freezing_curves(df = males_MR2_CS[,3:14],
                                                class = males_MR2_CS_param$class,
                                                retrieval = "MR2", stat = T,
                                                individual = TRUE)

pdf("../Figures/Freezing/Figure.8.fitted.m.MR2.pdf", height = 4, width = 5)
print(fitted.freezing.m.MR2)
dev.off()

################################################################################
#Bivariate Scatterplots of sustained and phasic freezers using decay rate and 
#Average freezing bins 18-24

scatter.f.mr1 <-ggplot(females_MR1_param, aes(x=freeze, y = -beta, col = class)) +
  geom_point() + theme_test() + theme(axis.text=element_text(size=12, color="black")) +
  labs(x = "Freezing MR1 bins 18-24 (%)", y ="Decay rate MR1", col="") +
  scale_color_manual(values=c("green2", "magenta")) + 
  ggtitle("MR1 Classification") +
  stat_ellipse() 

pdf("../Figures/Freezing/Figure.9.bivariate.f.MR1.pdf", height = 4, width = 4.5)
print(scatter.f.mr1)
dev.off()

scatter.f.mr2 <- ggplot(females_MR2_param, aes(x=freeze, y = -beta, col = class)) +
  geom_point() + theme_test() + theme(axis.text=element_text(size=12, color="black")) +
  labs(x = "Freezing MR2 bins 18-24 (%)", y ="Decay rate MR2", col="") +
  scale_color_manual(values=c("green2", "magenta")) + 
  ggtitle("MR2 Classification") +
  stat_ellipse()

pdf("../Figures/Freezing/Figure.10.bivariate.f.MR2.pdf", height = 4, width = 4.5)
print(scatter.f.mr2)
dev.off()

scatter.m.mr1 <-ggplot(males_MR1_param, aes(x=freeze, y = -beta, col = class)) +
  geom_point(shape=15) + theme_test() + theme(axis.text=element_text(size=12, color="black")) +
  labs(x = "Freezing MR1 bins 18-24 (%)", y ="Decay rate MR1", col="") +
  scale_color_manual(values=c("green2", "magenta")) + 
  ggtitle("MR1 Classification") +
  stat_ellipse() 

pdf("../Figures/Freezing/Figure.11.bivariate.m.MR1.pdf", height = 4, width = 4.5)
print(scatter.m.mr1)
dev.off()

scatter.m.mr2 <- ggplot(males_MR2_CS_param, aes(x=freeze, y = -beta, col = class)) +
  geom_point(shape=15) + theme_test() + theme(axis.text=element_text(size=12, color="black")) +
  labs(x = "Freezing MR2 bins 18-24 (%)", y ="Decay rate MR2", col="") +
  scale_color_manual(values=c("green2", "magenta")) + 
  ggtitle("MR2 Classification") +
  stat_ellipse() 

pdf("../Figures/Freezing/Figure.12.bivariate.m.MR2.pdf", height = 4, width = 4.5)
print(scatter.m.mr2)
dev.off()

#----------------------------------------------------
#Shifter subtype plots 

#TBD

################################################################################
#Average freezing values female vs male

#MR1
females_MR1$sex = "female"
males_MR1$sex = "male"
freezing.all.mr1 <- rbind(females_MR1, males_MR1)


freezing.m.vs.f.mr1 <- freezing_curves_general(df = freezing.all.mr1[,3:38],
                                               class = freezing.all.mr1$sex,
                                               animal = freezing.all.mr1$id,
                                               retrieval = "MR1",  
                                               individual = F, stat = T,
                                               shape = c(16,15),
                                               color = c("#ea7eb2", "#2397e6"))
#Rescale theme sizes
freezing.m.vs.f.mr1 <- freezing.m.vs.f.mr1 + 
  theme(axis.text = element_text(colour="black", size=12)) + theme_test(base_size=14) 

pdf("../Figures/Freezing/Figure.13.freezing.sex.MR1.pdf", height = 4, width = 5)
print(freezing.m.vs.f.mr1)
dev.off()

#MR2
females_MR2$sex = "female"
males_MR2$sex = "male"
freezing.all.mr2 <- rbind(females_MR2, males_MR2)

freezing.m.vs.f.mr2 <- freezing_curves_general(df = freezing.all.mr2[,3:38],
                                               class = freezing.all.mr2$sex,
                                               animal = freezing.all.mr2$id,
                                               retrieval = "MR2",  
                                               individual = F, stat = T,
                                               shape = c(16,15),
                                               color = c("#ea7eb2", "#2397e6"))

#Rescale theme sizes
freezing.m.vs.f.mr2 <- freezing.m.vs.f.mr2 + 
  theme(axis.text = element_text(colour="black", size=12)) + theme_test(base_size=14) 

pdf("../Figures/Freezing/Figure.14.freezing.sex.MR2.pdf", height = 4, width = 5)
print(freezing.m.vs.f.mr2)
dev.off()

################################################################################
#Parameter Distribution females vs males:

#Plot Density female vs male MR1
females_MR1_param.hist <- females_MR1_param
males_MR1_param.hist <- males_MR1_param

females_MR1_param.hist$class <- "sustained"
males_MR1_param.hist$class <- "phasic"

data <- rbind(females_MR1_param.hist, males_MR1_param.hist)

all.density.mr1 <- plot_density(variable = data$freeze, 
                            class = data$class, xlab = "Average freeze MR1 % Bins 18-24", title = "", breaks = c(10, 15, 20, 25), n=300)
all.density.mr1 <- all.density.mr1 +
  scale_color_manual(values=c("steelblue2", "violet"), labels = c("males", "female")) 



pdf("../Figures/Freezing/Figure.15.sex.dist.MR1.pdf", height = 4, width = 5)
print(all.density.mr1)
dev.off()


#Plot Density female vs male MR2
females_MR2_param.hist <- females_MR2_param
males_MR2_param.hist <- males_MR2_CS_param

females_MR2_param.hist$class <- "sustained"
males_MR2_param.hist$class <- "phasic"

data <- rbind(females_MR2_param.hist, males_MR2_param.hist)

all.density.mr2 <- plot_density(variable = data$freeze, 
                            class = data$class, xlab = "Average freeze MR1 % Bins 18-24", title = "", breaks = c(10, 15, 20, 25),
                            n = 392)
all.density.mr2 <- all.density.mr2 +
  scale_color_manual(values=c("steelblue2", "violet"), labels = c("males", "female")) 

pdf("../Figures/Freezing/Figure.16.sex.dist.MR2.pdf", height = 4, width = 5)
print(all.density.mr2)
dev.off()

################################################################################
#Independent Batch Plots 

#Code is currently in separate directory/script- tidy up later :) 

################################################################################
#Down sampling Analysis - Clustering VS pretrained ML-Models

female.mr1 <- openxlsx::read.xlsx("../Data/labelled.f.MR1.xlsx")
female.mr2 <- openxlsx::read.xlsx("../Data/labelled.f.MR2.xlsx")
male.mr1 <- openxlsx::read.xlsx("../Data/labelled.m.MR1.xlsx")
male.mr2 <- openxlsx::read.xlsx("../Data/labelled.m.MR2.xlsx")

#Obtain CS freezing values of all conditions
f.mr1.cs <- females_MR1[,15:26]
f.mr2.cs <- females_MR2[,15:26]
m.mr1.cs <- males_MR1[,15:26]

#Males MR2 was merged before due to different CS protocols
m.mr2.cs <- males_MR2_CS[,3:14]

#Determine Accuracy of downsampled cluster results compared to original cluster dataset

animal_numbers <- c(100, 95, 90, 85, 80, 75, 70, 65, 60, 55, 50, 45, 40, 35, 30, 25, 20, 15, 10, 5)
iterations <- 1000

#--------------------------------------
#Female MR1
results.f.mr1.clust <- data.frame(matrix(ncol=length(animal_numbers),nrow=iterations, dimnames=list(NULL, animal_numbers)))
results.f.mr1.ML <- data.frame(matrix(ncol=length(animal_numbers),nrow=iterations, dimnames=list(NULL, animal_numbers)))

data.f.mr1.clust <- data.frame(matrix(ncol=4,nrow=length(animal_numbers), 
                                dimnames=list(NULL, c("Dataset", "Subsample size", "Average accuracy", "se"))))

data.f.mr1.ML <- data.frame(matrix(ncol=4,nrow=length(animal_numbers), 
                                dimnames=list(NULL, c("Dataset", "Subsample size", "Average accuracy", "se"))))

for (animal in 1:length(animal_numbers)){
  Accuracy.clust <- c()
  Accuracy.ML <- c()
  for (i in 1:iterations){
    print(as.numeric(animal_numbers[animal]))
    #draw 40 random animals from dataframe
    current_df <- female.mr1[sample(nrow(female.mr1), as.numeric(animal_numbers[animal])), ]
    
    #Cluster the animals using mclust
    clust.current <- mclust::Mclust(current_df[,2:4], G=2)
    cluster_res <- assign.phenotype(current_df$freeze, clust.current$classification)
    current_df$clustered <- cluster_res
    
    #Classify animals using PhenoFreeze - get animal IDs
    current_freeze <- f.mr1.cs[current_df$B_ID,]
    print(all.equal.character(rownames(current_freeze), current_df$B_ID))
    
    current_df$ML <- PhenoFreeze::classify_freezer(data_MR1 = current_freeze, sex = "female", 
                                                   MR = 1, shifter = FALSE)
    
    #Declare phenotypes as factors
    current_df$class <- factor(current_df$class, levels = c("phasic", "sustained"))
    current_df$clustered <- factor(current_df$clustered, levels = c("phasic", "sustained"))
    current_df$ML <- factor(current_df$ML, levels = c("phasic", "sustained"))
    
    #Calculate Accuracy of clustered vs original clusterdataset
    cm.clust <- table(current_df$clustered, current_df$class)
    cmMatrix.clust <- caret::confusionMatrix(cm.clust, mode="everything")
    Accuracy.clust[i] <- cmMatrix.clust$overall["Accuracy"]
    
    cm.ML <- table(current_df$ML, current_df$class)
    cmMatrix.ML <- caret::confusionMatrix(cm.ML, mode="everything")
    Accuracy.ML[i] <- cmMatrix.ML$overall["Accuracy"]
    
    print(Accuracy.ML[i])
    print(cm.ML)
  }
  results.f.mr1.clust[,animal] <- Accuracy.clust
  results.f.mr1.ML[,animal] <- Accuracy.ML
}

colnames(results.f.mr1.clust) <- animal_numbers
colnames(results.f.mr1.ML) <- animal_numbers

#Get Average Accuracies and SD
mean.f.mr1.clust = base::colMeans(results.f.mr1.clust)
se.f.mr1.clust = base::apply(results.f.mr1.clust, 2, function(x)sd(x)/sqrt(length(x)))

mean.f.mr1.ML = base::colMeans(results.f.mr1.ML)
se.f.mr1.ML = base::apply(results.f.mr1.ML, 2, function(x)sd(x)/sqrt(length(x)))

#Append results to DF
data.f.mr1.clust$Dataset <- "F MR1"
data.f.mr1.clust$Subsample.size <- animal_numbers
data.f.mr1.clust$Average.accuracy <- mean.f.mr1.clust
data.f.mr1.clust$se <- se.f.mr1.clust

data.f.mr1.ML$Dataset <- "F MR1"
data.f.mr1.ML$Subsample.size <- animal_numbers
data.f.mr1.ML$Average.accuracy <- mean.f.mr1.ML
data.f.mr1.ML$se <- se.f.mr1.ML

#--------------------------------------
#Female MR2
results.f.mr2.clust <- data.frame(matrix(ncol=length(animal_numbers),nrow=iterations, dimnames=list(NULL, animal_numbers)))
results.f.mr2.ML <- data.frame(matrix(ncol=length(animal_numbers),nrow=iterations, dimnames=list(NULL, animal_numbers)))

data.f.mr2.clust <- data.frame(matrix(ncol=4,nrow=length(animal_numbers), 
                                      dimnames=list(NULL, c("Dataset", "Subsample size", "Average accuracy", "se"))))

data.f.mr2.ML <- data.frame(matrix(ncol=4,nrow=length(animal_numbers), 
                                   dimnames=list(NULL, c("Dataset", "Subsample size", "Average accuracy", "se"))))

for (animal in 1:length(animal_numbers)){
  Accuracy.clust <- c()
  Accuracy.ML <- c()
  for (i in 1:iterations){
    print(as.numeric(animal_numbers[animal]))
    #draw 40 random animals from dataframe
    current_df <- female.mr2[sample(nrow(female.mr2), as.numeric(animal_numbers[animal])), ]
    
    #Cluster the animals using mclust
    clust.current <- invisible(mclust::Mclust(current_df[,2:4], G=2))
    cluster_res <- assign.phenotype(current_df$freeze, clust.current$classification)
    current_df$clustered <- cluster_res
    
    #Classify animals using PhenoFreeze - get animal IDs
    current_freeze <- f.mr2.cs[current_df$B_ID,]
    
    current_df$ML <- PhenoFreeze::classify_freezer(data_MR2 = current_freeze, sex = "female", 
                                                   MR = 2, shifter = FALSE)
    print(current_df$ML)
    #Declare phenotypes as factors
    current_df$class <- factor(current_df$class, levels = c("phasic", "sustained"))
    current_df$clustered <- factor(current_df$clustered, levels = c("phasic", "sustained"))
    current_df$ML <- factor(current_df$ML, levels = c("phasic", "sustained"))
    
    #Calculate Accuracy of clustered vs original clusterdataset
    cm.clust <- table(current_df$clustered, current_df$class)
    cmMatrix.clust <- caret::confusionMatrix(cm.clust, mode="everything")
    Accuracy.clust[i] <- cmMatrix.clust$overall["Accuracy"]
    
    cm.ML <- table(current_df$ML, current_df$class)
    cmMatrix.ML <- caret::confusionMatrix(cm.ML, mode="everything")
    Accuracy.ML[i] <- cmMatrix.ML$overall["Accuracy"]
    print(Accuracy.ML[i])
    
  }
  results.f.mr2.clust[,animal] <- Accuracy.clust
  results.f.mr2.ML[,animal] <- Accuracy.ML
}

colnames(results.f.mr2.clust) <- animal_numbers
colnames(results.f.mr2.ML) <- animal_numbers

#Get Average Accuracies and SD
mean.f.mr2.clust = base::colMeans(results.f.mr2.clust)
se.f.mr2.clust = base::apply(results.f.mr2.clust, 2, function(x)sd(x)/sqrt(length(x)))

mean.f.mr2.ML = base::colMeans(results.f.mr2.ML)
se.f.mr2.ML = base::apply(results.f.mr2.ML, 2, function(x)sd(x)/sqrt(length(x)))

#Append results to DF
data.f.mr2.clust$Dataset <- "F MR2"
data.f.mr2.clust$Subsample.size <- animal_numbers
data.f.mr2.clust$Average.accuracy <- mean.f.mr2.clust
data.f.mr2.clust$se <- se.f.mr2.clust

data.f.mr2.ML$Dataset <- "F MR2"
data.f.mr2.ML$Subsample.size <- animal_numbers
data.f.mr2.ML$Average.accuracy <- mean.f.mr2.ML
data.f.mr2.ML$se <- se.f.mr2.ML


#--------------------------------------
#Male MR1
results.m.mr1.clust <- data.frame(matrix(ncol=length(animal_numbers),nrow=iterations, dimnames=list(NULL, animal_numbers)))
results.m.mr1.ML <- data.frame(matrix(ncol=length(animal_numbers),nrow=iterations, dimnames=list(NULL, animal_numbers)))

data.m.mr1.clust <- data.frame(matrix(ncol=4,nrow=length(animal_numbers), 
                                      dimnames=list(NULL, c("Dataset", "Subsample size", "Average accuracy", "se"))))

data.m.mr1.ML <- data.frame(matrix(ncol=4,nrow=length(animal_numbers), 
                                   dimnames=list(NULL, c("Dataset", "Subsample size", "Average accuracy", "se"))))


for (animal in 1:length(animal_numbers)){
  Accuracy.clust <- c()
  Accuracy.ML <- c()
  for (i in 1:iterations){
    print(as.numeric(animal_numbers[animal]))
    current_df <- male.mr1[sample(nrow(male.mr1), as.numeric(animal_numbers[animal])), ]
    
    #Cluster the animals using mclust
    clust.current <- invisible(mclust::Mclust(current_df[,2:4], G=2))
    cluster_res <- assign.phenotype(current_df$freeze, clust.current$classification)
    current_df$clustered <- cluster_res
    
    #Classify animals using PhenoFreeze - get animal IDs
    current_freeze <- m.mr1.cs[current_df$B_ID,]
    
    current_df$ML <- PhenoFreeze::classify_freezer(data_MR1 = current_freeze, sex = "male", 
                                                   MR = 1, shifter = FALSE)
    
    #Declare phenotypes as factors
    current_df$class <- factor(current_df$class, levels = c("phasic", "sustained"))
    current_df$clustered <- factor(current_df$clustered, levels = c("phasic", "sustained"))
    current_df$ML <- factor(current_df$ML, levels = c("phasic", "sustained"))
    
    #Calculate Accuracy of clustered vs original clusterdataset
    cm.clust <- table(current_df$clustered, current_df$class)
    cmMatrix.clust <- caret::confusionMatrix(cm.clust, mode="everything")
    Accuracy.clust[i] <- cmMatrix.clust$overall["Accuracy"]
    
    cm.ML <- table(current_df$ML, current_df$class)
    cmMatrix.ML <- caret::confusionMatrix(cm.ML, mode="everything")
    Accuracy.ML[i] <- cmMatrix.ML$overall["Accuracy"]
    print(Accuracy.ML[i])
    
  }
  results.m.mr1.clust[,animal] <- Accuracy.clust
  results.m.mr1.ML[,animal] <- Accuracy.ML
}

colnames(results.m.mr1.clust) <- animal_numbers
colnames(results.m.mr1.ML) <- animal_numbers

#Get Average Accuracies and SD
mean.m.mr1.clust = base::colMeans(results.m.mr1.clust)
se.m.mr1.clust = base::apply(results.m.mr1.clust, 2, function(x)sd(x)/sqrt(length(x)))

mean.m.mr1.ML = base::colMeans(results.m.mr1.ML)
se.m.mr1.ML = base::apply(results.m.mr1.ML, 2, function(x)sd(x)/sqrt(length(x)))

#Append results to DF
data.m.mr1.clust$Dataset <- "M MR1"
data.m.mr1.clust$Subsample.size <- animal_numbers
data.m.mr1.clust$Average.accuracy <- mean.m.mr1.clust
data.m.mr1.clust$se <- se.m.mr1.clust

data.m.mr1.ML$Dataset <- "M MR1"
data.m.mr1.ML$Subsample.size <- animal_numbers
data.m.mr1.ML$Average.accuracy <- mean.m.mr1.ML
data.m.mr1.ML$se <- se.m.mr1.ML


#--------------------------------------
#Male MR2
results.m.mr2.clust <- data.frame(matrix(ncol=length(animal_numbers),nrow=iterations, dimnames=list(NULL, animal_numbers)))
results.m.mr2.ML <- data.frame(matrix(ncol=length(animal_numbers),nrow=iterations, dimnames=list(NULL, animal_numbers)))

data.m.mr2.clust <- data.frame(matrix(ncol=4,nrow=length(animal_numbers), 
                                      dimnames=list(NULL, c("Dataset", "Subsample size", "Average accuracy", "se"))))

data.m.mr2.ML <- data.frame(matrix(ncol=4,nrow=length(animal_numbers), 
                                   dimnames=list(NULL, c("Dataset", "Subsample size", "Average accuracy", "se"))))

for (animal in 1:length(animal_numbers)){
  Accuracy.clust <- c()
  Accuracy.ML <- c()
  for (i in 1:iterations){
    print(as.numeric(animal_numbers[animal]))
    current_df <- male.mr2[sample(nrow(male.mr2), as.numeric(animal_numbers[animal])), ]
    
    #Cluster the animals using mclust
    clust.current <- invisible(mclust::Mclust(current_df[,2:4], G=2))
    cluster_res <- assign.phenotype(current_df$freeze, clust.current$classification)
    current_df$clustered <- cluster_res
    
    #Classify animals using PhenoFreeze - get animal IDs
    current_freeze <- m.mr2.cs[current_df$B_ID,]
    print(all.equal.character(rownames(current_freeze), current_df$B_ID))
    
    current_df$ML <- PhenoFreeze::classify_freezer(data_MR2 = current_freeze, sex = "male", 
                                                   MR = 2, shifter = FALSE)
    
    #Declare phenotypes as factors
    current_df$class <- factor(current_df$class, levels = c("phasic", "sustained"))
    current_df$clustered <- factor(current_df$clustered, levels = c("phasic", "sustained"))
    current_df$ML <- factor(current_df$ML, levels = c("phasic", "sustained"))
    
    #Calculate Accuracy of clustered vs original clusterdataset
    cm.clust <- table(current_df$clustered, current_df$class)
    cmMatrix.clust <- caret::confusionMatrix(cm.clust, mode="everything")
    Accuracy.clust[i] <- cmMatrix.clust$overall["Accuracy"]
    
    cm.ML <- table(current_df$ML, current_df$class)
    cmMatrix.ML <- caret::confusionMatrix(cm.ML, mode="everything")
    Accuracy.ML[i] <- cmMatrix.ML$overall["Accuracy"]
    print(Accuracy.ML[i])
    
  }
  results.m.mr2.clust[,animal] <- Accuracy.clust
  results.m.mr2.ML[,animal] <- Accuracy.ML
}

colnames(results.m.mr2.clust) <- animal_numbers
colnames(results.m.mr2.ML) <- animal_numbers

#Get Average Accuracies and SD
mean.m.mr2.clust = base::colMeans(results.m.mr2.clust)
se.m.mr2.clust = base::apply(results.m.mr2.clust, 2, function(x)sd(x)/sqrt(length(x)))

mean.m.mr2.ML = base::colMeans(results.m.mr2.ML)
se.m.mr2.ML = base::apply(results.m.mr2.ML, 2, function(x)sd(x)/sqrt(length(x)))

#Append results to DF
data.m.mr2.clust$Dataset <- "M MR2"
data.m.mr2.clust$Subsample.size <- animal_numbers
data.m.mr2.clust$Average.accuracy <- mean.m.mr2.clust
data.m.mr2.clust$se <- se.m.mr2.clust

data.m.mr2.ML$Dataset <- "M MR2"
data.m.mr2.ML$Subsample.size <- animal_numbers
data.m.mr2.ML$Average.accuracy <- mean.m.mr2.ML
data.m.mr2.ML$se <- se.m.mr2.ML


#-------------------------------------------------------------------------------
#Plotting

#Plot cluster downsampling
lineplot.data.clust <- rbind(data.f.mr1.clust, data.f.mr2.clust, data.m.mr1.clust, data.m.mr2.clust)

lineplot.clust <- ggplot(lineplot.data.clust, aes(x=Subsample.size, y=Average.accuracy, group=Dataset, color=Dataset)) + 
  geom_errorbar(aes(ymin=Average.accuracy-se, ymax=Average.accuracy+se), width=.2) +
  geom_line() + geom_point(aes(shape=Dataset))+ theme_test(base_size=14) + theme(axis.text = element_text(colour="black", size=12)) + 
  ylab("Average Accuracy") + xlab("Subsample Size") + scale_color_manual(values = c("red1", "red4", "royalblue1", "royalblue4")) +
  scale_shape_manual(values = c(16, 16, 15, 15)) +
  scale_x_continuous(trans='reverse', breaks = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100)) + theme(legend.title = element_blank())

pdf("../Figures/Downsampling/Accuracies.Cluster.Downsampling.pdf", width = 5, height = 4)
print(lineplot.clust)
dev.off()

#Plot PhenoFreeze downsampling
lineplot.data.ML <- rbind(data.f.mr1.ML, data.f.mr2.ML, data.m.mr1.ML, data.m.mr2.ML)

lineplot.ML <- ggplot(lineplot.data.ML, aes(x=Subsample.size, y=Average.accuracy, group=Dataset, color=Dataset)) + 
  geom_errorbar(aes(ymin=Average.accuracy-se, ymax=Average.accuracy+se), width=.2) +
  geom_line() + geom_point(aes(shape=Dataset))+ theme_test(base_size=14) + theme(axis.text = element_text(colour="black", size=12)) + 
  ylab("Average Accuracy") + xlab("Subsample Size") + scale_color_manual(values = c("#ea7eb2ff", "#aa0044ff", "#2395e5ff", "#27408bff")) +
  scale_shape_manual(values = c(16, 16, 15, 15)) +
  scale_x_continuous(trans='reverse', breaks = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100)) + theme(legend.title = element_blank())

pdf("../Figures/Downsampling/Accuracies.ML.Downsampling.pdf", width = 5, height = 4)
print(lineplot.ML)
dev.off()


################################################################################
#-------------------------------------------------------------------------------

#Cluster Accuracy vs phenoFreeze Accuracy, clustered/predicted vs original training cluster dataset

#Check if IDs are in the same order compared to the raw freezing value table, 
#both are used in the next step
f.mr1.raw <- read.xlsx("../Data/all_females_MR1.xlsx")
all(f.mr1.raw$id == female.mr1$id)

f.mr1.raw <- f.mr1.raw[, c(1,2,15:26)]
results.f.mr1 <- data.frame(matrix(ncol=length(animal_numbers),nrow=iterations, dimnames=list(NULL, animal_numbers)))

results.f.mr1.classify <- data.frame(matrix(ncol=length(animal_numbers),nrow=iterations, dimnames=list(NULL, animal_numbers)))


data.f.mr1 <- data.frame(matrix(ncol=4,nrow=length(animal_numbers), 
                                dimnames=list(NULL, c("Dataset", "Subsample size", "Average accuracy", "se"))))

s.e. <- function(x)sd(x)/sqrt(length(x))

mean.classified <- c()
se.classified <- c()
mean.clustered <- c()
se.clustered <- c()

animal_numbers <- c(100, 95, 90, 85, 80, 75, 70, 65, 60, 55, 50, 45, 40, 35, 30, 25, 20, 15, 10, 5)
iterations <- 1000

for (animal in 1:length(animal_numbers)){
  Accuracy.clust <- c()
  Accuracy.classified <- c()
  for (i in 1:iterations){
    print(as.numeric(animal_numbers[animal]))
    #draw random animals from dataframe
    samplenum <- sample(nrow(female.mr1), as.numeric(animal_numbers[animal]))
    
    #Subset each df
    current_df <- female.mr1[samplenum, ]
    classify_df <- f.mr1.raw[samplenum,]
    
    #Cluster the animals using mclust
    clust.current <- invisible(mclust::Mclust(current_df[,2:4], G=2))
    cluster_res <- assign.phenotype(current_df$freeze, clust.current$classification)
    current_df$clustered <- cluster_res
    
    #Classify animals with phenoFreeze
    classify_res <- PhenoFreeze::classify_freezer(data_MR1 = classify_df[,3:14], sex = "female", MR = 1, shifter = F)
    classify_df$classified <- classify_res
    
    #Declare phenotypes as factors
    current_df$class <- factor(current_df$class, levels = c("phasic", "sustained"))
    current_df$clustered <- factor(current_df$clustered, levels = c("phasic", "sustained"))
    
    classify_df$class <- factor(current_df$class, levels = c("phasic", "sustained"))
    classify_df$classified <- factor(classify_df$classified, levels = c("phasic", "sustained"))
    
    #Calculate Accuracy of clustered vs original clusterdataset
    cm <- table(current_df$clustered, current_df$class)
    cmMatrix <- caret::confusionMatrix(cm, mode="everything")
    Accuracy.clust[i] <- cmMatrix$overall["Accuracy"]
    
    cm.classified <- table(classify_df$classified, classify_df$class)
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
  ylab("Average Accuracy") + xlab("Subsample Size") + scale_color_manual(values = c("#005f73", "#fb8500")) +
  scale_x_continuous(trans='reverse', breaks = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100)) + theme(legend.title = element_blank())

lineplot

pdf("../Figures/Downsampling/Cluster.vs.ML.Downsampling.pdf", width = 6, height = 4)
print(lineplot)
dev.off()

################################################################################
#Cluster Stability Analysis with Bootstrapping

#-------------------------------------------------------------------------------
#Females MR1
f.mr1.bootstrap <- bootstrapping(data = females_MR1_param, data.freezing = females_MR1,
                                 iterations = 200, retrieval = "MR1", shape = 16)

#Save Plots
pdf("../Figures/Bootstrapping/f.mr1.jaccard.pdf", height = 4, width = 4.5)
print(f.mr1.bootstrap$Jaccard)
dev.off()

pdf("../Figures/Bootstrapping/f.mr1.sustained.pdf", height = 4, width = 5)
print(f.mr1.bootstrap$Sustained.freezing)
dev.off()

pdf("../Figures/Bootstrapping/f.mr1.phasic.pdf", height = 4, width = 5)
print(f.mr1.bootstrap$Phasic.freezing)
dev.off()

#-------------------------------------------------------------------------------
#Females MR2
f.mr2.bootstrap <- bootstrapping(data = females_MR2_param, data.freezing = females_MR2,
                                 iterations = 200, retrieval = "MR2", shape = 16)

#Save Plots
pdf("../Figures/Bootstrapping/f.mr2.jaccard.pdf", height = 4, width = 4.5)
print(f.mr2.bootstrap$Jaccard)
dev.off()

pdf("../Figures/Bootstrapping/f.mr2.sustained.pdf", height = 4, width = 5)
print(f.mr2.bootstrap$Sustained.freezing)
dev.off()

pdf("../Figures/Bootstrapping/f.mr2.phasic.pdf", height = 4, width = 5)
print(f.mr2.bootstrap$Phasic.freezing)
dev.off()


#-------------------------------------------------------------------------------
#Males MR1
m.mr1.bootstrap <- bootstrapping(data = males_MR1_param, data.freezing = males_MR1,
                                 iterations = 200, retrieval = "MR1", shape = 15)

#Save Plots
pdf("../Figures/Bootstrapping/m.mr1.jaccard.pdf", height = 4, width = 4.5)
print(m.mr1.bootstrap$Jaccard)
dev.off()

pdf("../Figures/Bootstrapping/m.mr1.sustained.pdf", height = 4, width = 5)
print(m.mr1.bootstrap$Sustained.freezing)
dev.off()

pdf("../Figures/Bootstrapping/m.mr1.phasic.pdf", height = 4, width = 5)
print(m.mr1.bootstrap$Phasic.freezing)
dev.off()

#-------------------------------------------------------------------------------
#Males MR2
m.mr2.bootstrap <- bootstrapping(data = males_MR2_CS_param, data.freezing = males_MR2_freezing,
                                 iterations = 200, retrieval = "MR2", shape = 15)

#Save Plots
pdf("../Figures/Bootstrapping/m.mr2.jaccard.pdf", height = 4, width = 4.5)
print(m.mr2.bootstrap$Jaccard)
dev.off()

pdf("../Figures/Bootstrapping/m.mr2.sustained.pdf", height = 4, width = 5)
print(m.mr2.bootstrap$Sustained.freezing)
dev.off()

pdf("../Figures/Bootstrapping/m.mr2.phasic.pdf", height = 4, width = 5)
print(m.mr2.bootstrap$Phasic.freezing)
dev.off()

################################################################################
#Training and Envaluation of Machine-learning Models

#Load labelled data sets for binary training
female.MR1 <- openxlsx::read.xlsx("../Data/labelled.f.MR1.xlsx")
female.MR2 <- openxlsx::read.xlsx("../Data/labelled.f.MR2.xlsx")
male.MR1 <- openxlsx::read.xlsx("../Data/labelled.m.MR1.xlsx")
male.MR2 <- openxlsx::read.xlsx("../Data/labelled.m.MR2.xlsx")
female.shifter <- openxlsx::read.xlsx("../Data/labelled.f.shifter.xlsx")
male.shifter <- openxlsx::read.xlsx("../Data/labelled.m.shifter.xlsx")

#-------------------------------------------------------------------------------

#Train binary models
#Exclude batch and id column during training

females_MR1_models <- train_binary(data=female.MR1[,c(2:5)],
                                   steps = 1000, title="", shape = 16)

#Save all model metrics
females_MR1_models_stats <- females_MR1_models$Stats
openxlsx::write.xlsx(females_MR1_models_stats, "../ML_metrics/metrics.xlsx", sheetName="females MR1", rowNames = TRUE)

#### REOPEN EXISTING EXCEL FILE as wb
wb <- loadWorkbook("../ML_metrics/metrics.xlsx")


pdf("../Figures/ML/Figure.17.females.MR1.acc.pdf", width = 4, height = 4)
print(females_MR1_models$Accuracy)
dev.off()

females_MR2_models <- train_binary(data=female.MR2[,c(2:5)],
                                   steps = 1000, title="", shape = 16)


#Save all model metrics
females_MR2_models_stats <- females_MR2_models$Stats
addWorksheet(wb,"females MR2")
writeData(wb,"females MR2",females_MR2_models_stats, rowNames = TRUE)


pdf("../Figures/ML/Figure.18.females.MR2.acc.pdf", width = 4, height = 4)
print(females_MR2_models$Accuracy)
dev.off()

males_MR1_models <- train_binary(data=male.MR1[,c(2:5)],
                                 steps = 1000, title="", shape = 15)

#Save all model metrics
males_MR1_models_stats <- males_MR1_models$Stats
addWorksheet(wb,"males MR1")
writeData(wb,"males MR1",males_MR1_models_stats, rowNames = TRUE)

pdf("../Figures/ML/Figure.19.males.MR1.acc.pdf", width = 4, height = 4)
print(males_MR1_models$Accuracy)
dev.off()


males_MR2_models <- train_binary(data=male.MR2[,c(2:5)],
                                 steps = 1000, title="", shape = 15)

#Save all model metrics
males_MR2_models_stats <- males_MR2_models$Stats
addWorksheet(wb,"males MR2")
writeData(wb,"males MR2",males_MR2_models_stats, rowNames = TRUE)

pdf("../Figures/ML/Figure.20.males.MR2.acc.pdf", width = 4, height = 4)
print(males_MR2_models$Accuracy)
dev.off()

#-------------------------------------------------------------------------------
#Train shifter models

#Declare class as factor
female.shifter$class <- as.factor(female.shifter$class)
male.shifter$class <- as.factor(male.shifter$class)

females_shifter.MR2 <- train_shifter(data=female.shifter[,c(6:8, 10)],
                                     steps = 1000, title="", shape = 16)

#Save all model metrics
females_MR2_shifter_stats <- females_shifter.MR2$Stats
addWorksheet(wb,"females shifter MR2")
writeData(wb,"females shifter MR2",females_MR2_shifter_stats, rowNames = TRUE)


pdf("../Figures/ML/Figure.21.females.shifter.MR2.acc.pdf", width = 4, height = 4)
print(females_shifter.MR2$Accuracy)
dev.off()

females_shifter.MR1.MR2 <- train_shifter(data=female.shifter[,c(2:4, 6:8, 10)],
                                         steps = 1000, title="", shape = 16)

#Save all model metrics
females_MR1_MR2_shifter_stats <- females_shifter.MR1.MR2$Stats
addWorksheet(wb,"females shifter MR1 MR2")
writeData(wb,"females shifter MR1 MR2",females_MR1_MR2_shifter_stats, rowNames = TRUE)

pdf("../Figures/ML/Figure.22.females.shifter.MR2.MR1.acc.pdf", width = 4, height = 4)
print(females_shifter.MR1.MR2$Accuracy)
dev.off()

males_shifter.MR2 <- train_shifter(data=male.shifter[,c(6:8, 10)],
                                   steps = 1000, title="", shape = 15)

#Save all model metrics
males_MR2_shifter_stats <- males_shifter.MR2$Stats
addWorksheet(wb,"males shifter MR2")
writeData(wb,"males shifter MR2", males_MR2_shifter_stats, rowNames = TRUE)

pdf("../Figures/ML/Figure.23.males.shifter.MR2.acc.pdf", width = 4, height = 4)
print(males_shifter.MR2$Accuracy)
dev.off()

males_shifter.MR1.MR2 <- train_shifter(data=male.shifter[,c(2:4, 6:8, 10)],
                                       steps = 1000, title="", shape = 15)

males_MR1_MR2_shifter_stats <- males_shifter.MR1.MR2$Stats
addWorksheet(wb,"males shifter MR1 MR2")
writeData(wb,"males shifter MR1 MR2",males_MR1_MR2_shifter_stats, rowNames = TRUE)


pdf("../Figures/ML/Figure.24.males.shifter.MR1.MR2.acc.pdf", width = 4, height = 4)
print(males_shifter.MR1.MR2$Accuracy)
dev.off()

#Save all ML-metrics as xlsx file 
saveWorkbook(wb,"../ML_metrics/metrics.xlsx",overwrite = TRUE)
