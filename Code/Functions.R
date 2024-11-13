################################################################################
#PhenoFreeze analysis code
#Script with various functions used for the analysis 

library(ggplot2)
library(mclust)
library(lme4)
library(emmeans)
library(e1071)
library(MASS)
library(randomForest)
library(caret) 
library(pROC)
library(reshape2)
library(ggplot2)
library(ggsci)
library(nnet)

#-------------------------------------------------------------------------------
#This function performs loglinear regression
#param df_CS: Dataframe with 12 time bins where CS was played

perform_loglinear <- function(df_CS){
  #Initialize bin variable
  bins <- 1:12
  
  #Initialize R squared, decay rate (beta) and intercept (int) variables
  r2 <- 0
  beta <- 0
  int <- 0 
  
  #Perform loglinear regression for each row/animal
  for (i in 1:nrow(df_CS)){
    #Transpose matrix
    y <- t(df_CS[i,])
    #Impute value 1 in value is 0 to allow logarithmic transformation
    if (0 %in% y){y<-y+1}
    #Fit loglinear model
    mod <- lm(log(y)~bins)
    #Obtain R2 value
    r2[i] <- summary(mod)$r.squared
    #Otain beta coefficient
    beta[i] <- mod$coefficients[2]
    #Obtain intercept
    int[i] <- mod$coefficients[1]
  }
  
  #Output result dataframe with int, beta, r2 and average freezing of CS bins of
  #7-12 as additional columns
  freeze <- rowMeans(df_CS[,6:12])
  result <- data.frame(freeze, beta, int, r2)
  
  return(result)
}

#-------------------------------------------------------------------------------
#This function assigns phenotype to clusters based on average freezing
#Sustained gets assigned to the cluster with higher average freezing

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

#-------------------------------------------------------------------------------
#Function to plot freezing curves of sustained and phasic freezers (as well as
#shifters if available)

freezing_curves = function(df, class, animal, retrieval,  individual = FALSE, 
                           stat = FALSE,
                           shape = 16){
  
  #Time bins (30 sec)    
  bins_all = 1:36
  
  #Calculate mean freezing values per time bin and phenotype
  mean_ph = base::colMeans(df[class=="phasic",], na.rm = TRUE)
  df_ph = df[class=="phasic",]
  se_ph = c()
  for (i in 1:ncol(df_ph)){
    col = as.vector(df_ph[,i])
    col = col[!is.na(col)]
    se_ph[i] = sd(col)/sqrt(length(col))
  }
  
  mean_sus = base::colMeans(df[class=="sustained",], na.rm = TRUE)
  df_sus = df[class=="sustained",]
  se_sus = c()
  for (i in 1:ncol(df_sus)){
    col = as.vector(df_sus[,i])
    col = col[!is.na(col)]
    se_sus[i] = sd(col)/sqrt(length(col))
  }
  
  if ("shifter" %in% class){
    mean_shift = base::colMeans(df[class=="shifter",])
    se_shift = base::apply(df[class=="shifter",], 2, sd(x)/sqrt(length(x)))
  }
  
  #-----
  data.plot = data.frame(mean_ph, se_ph, bins_all)
  colnames(data.plot) = c("freeze", "se", "bins")
  
  
  p1 = ggplot2::ggplot(data.plot, aes(bins, freeze))  + theme_test(base_size=14) +
    labs(x = "Bin", y = paste0("Average freezing per bin ", retrieval,  " [%]")) + 
    geom_line(col =   "#006400ff", size=1) + geom_point(size=2.5, col=  "#006400ff", shape=shape) + 
    geom_errorbar(aes(ymin=freeze-se, ymax=freeze+se), width=.2,
                  position=position_dodge(0.05), col =   "#006400ff") +
    scale_x_continuous(breaks = c(1,13,18,24,36), labels = c(1,13,18,24,36)) + 
    scale_y_continuous(expand = c(0,1), limits = c(-4,100)) +
    theme(axis.text = element_text(colour="black", size=12)) 
  
  
  
  #Add individual curves if individual == TRUE
  
  if (individual == TRUE){
    
    df_sus = df[class == "sustained",]
    df_ph = df[class == "phasic",]
    
    
    #Add individual curves sustained
    for (i in 1:nrow(df_sus)){
      data = data.frame(t(df_sus[i,]), bins_all)
      colnames(data) = c("freeze", "bins")
      
      p1 = p1 + geom_line(col="magenta", data=data, alpha = 0.1) + 
        geom_point(size=1.5, col="magenta", data=data, alpha = 0.1, shape=shape)
    }
    
    #Add individual curves phasic
    for (i in 1:nrow(df_ph)){
      data = data.frame(t(df_ph[i,]), bins_all)
      colnames(data) = c("freeze", "bins")
      
      p1 = p1 + geom_line(col="green2", data=data, alpha = 0.1) + 
        geom_point(size=1.5, col="green2", data=data, alpha = 0.1, shape=shape)
    }
    
    
    #Add shifters if available
    if("shifter" %in% class){
      df_shift = df[class == "shifter",]
      
      #Add individual curves shifters
      for (i in 1:nrow(df_shift)){
        data = data.frame(t(df_shift[i,]), bins_all)
        colnames(data) <- c("freeze", "bins")
        
        p1 = p1 + geom_line(col="gray50", data=data, alpha = 0.1) + 
          geom_point(size=1.5, col="gray50", data=data, alpha = 0.1, shape=shape)
      }
    }
    
    #Replot average curve for phasic animals
    p1 = p1 + geom_line(col= "#006400ff", data=data.plot, size=1) + 
      geom_point(size=2.5, col= "#006400ff", data=data.plot, shape=shape) + 
      geom_errorbar(data= data.plot, mapping=aes(ymin=freeze-se, ymax=freeze+se), width=.2,
                    position=position_dodge(0.05), col =  "#006400ff") 
  }
  
  #Add sustained animals
  data.plot.sus = data.frame(mean_sus, se_sus, bins_all)
  colnames(data.plot.sus) = c("freeze", "se", "bins")
  p1 = p1 + geom_line(col="#c300b4ff", data=data.plot.sus, size=1) + 
    geom_point(size=2.5, col="#c300b4ff", data=data.plot.sus, shape=shape) + 
    geom_errorbar(data= data.plot.sus, mapping=aes(ymin=freeze-se, ymax=freeze+se), width=.2,
                  position=position_dodge(0.05), col = "#c300b4ff") 
  
  #Add shifter animals if available 
  if("shifter" %in% class){
    data.plot.shift = data.frame(mean_shift, se_shift, bins_all)
    colnames(data.plot.shift) = c("freeze", "se", "bins")
    p1 = p1 + geom_line(col="gray50", data=data.plot.shift, size=1) + 
      geom_point(size=2.5, col="gray50", data=data.plot.shift, shape=shape) + 
      geom_errorbar(data= data.plot.shift, mapping=aes(ymin=freeze-se, ymax=freeze+se), 
                    width=.2,
                    position=position_dodge(0.05), col = "gray50")
  }
  
  
  #Add a line for region designating the tone
  p1 = p1 + ggplot2::geom_rect(aes(xmin=13, xmax = 18, ymin = -4, ymax = -2), 
                               color= "green2", fill = "green2")
  
  p1 = p1 + ggplot2::geom_rect(aes(xmin=18, xmax = 24, ymin = -4, ymax = -2), 
                               color= "magenta", fill = "magenta")
  p1 = p1 + annotate(geom = "text", x = c(17.9, 18.1), y = -3, label= c("C", "S"),
                     fontface = 2, size =3, hjust = c(1,0))
  
  #############################################
  #Perform linear mixed effects analysis comparing sustained versus phasic animals
  
  if (stat == TRUE){
    df$animal = animal
    df$phenotype = class
    
    #Keep only sustained and phasic phenotypes
    df1 = df[df$phenotype != "shifter",]
    
    #Melt to long format
    df_stat = reshape2::melt(df1, id.vars = c("animal", "phenotype"))
    colnames(df_stat)[3:4] = c("bin", "freezing")
    
    #Fit a linear model
    lmm = lme4::lmer(freezing~phenotype*bin + (1|animal), data=df_stat)
    #Compare marginal means
    em = emmeans::emmeans(lmm, pairwise~phenotype|bin)
    
    #Extract p-values
    p.vals = summary(em$contrasts)[7]
    
    #Plot significant p-values
    if(all(p.vals$p.value>0.05)==F){
      x = which(p.vals$p.value<0.05)
      
      p.vals$symbol = ifelse(p.vals$p.value<0.05, "*", "-")
      p.vals$symbol = ifelse(p.vals$p.value<0.01, "**", p.vals$symbol)
      p.vals$symbol = ifelse(p.vals$p.value<0.001, "***", p.vals$symbol)
      
      label = p.vals$symbol[x]
      
      #Get y coordinate for the statistics
      
      y = mean_sus[x] + se_sus[x] + 1
      
      p1 = p1 + annotate(geom="text", x = x, y = y, label = label,
                         size=3.5, fontface=2)
      
    }
  }
  
  
  return(p1)
}

#-------------------------------------------------------------------------------

#Function to plot the fitted freezing curves obtained with a log-linear fitting strategy

fitted_freezing_curves <- function(df, class, 
                                   retrieval,  individual = FALSE, stat = FALSE){
  
  
  time <- seq(1,12,0.012)
  bins <- 1:12
  
  #Average freezing curve sustained animals
  data_sus <- colMeans(df[class =="sustained",])
  
  data.plot <- data.frame(data_sus, bins)
  colnames(data.plot) <- c("freeze", "bins")
  
  if (0 %in% data.plot$freeze) {data.plot$freeze <- data.plot$freeze + 1}
  mod.sus <- lm(log(freeze)~bins, data=data.plot)
  y_new <- exp(coefficients(mod.sus)[1] + coefficients(mod.sus)[2]*time)
  p <- ggplot(data.plot, aes(bins, freeze))  + theme_test(base_size=14) +
    labs(x = "Bin", y = paste0("Fitted freezing curve ", retrieval, " [%]")) + 
    geom_line(data=data.frame(time, y_new), 
              mapping=aes(x = time, y = y_new), col = "#c300b4ff", size=1) +
    scale_x_continuous(breaks = 1:12, labels = 13:24) +
    theme(axis.text = element_text(colour="black", size=12))
  
  #Add fitted curves for individual animals if individual = T
  
  if (individual == T){
    cols <- data.frame(pheno = c("sustained", "phasic", "shifter"),
                       col = c("magenta", "green2", "gray50"))
    for (i in 1:nrow(df)){
      col <- cols$col[cols$pheno == class[i]]
      data <- data.frame(t(df[i,]), bins)
      colnames(data) <- c("freeze", "bins")
      if (0 %in% data$freeze) {data$freeze <- data$freeze + 1}
      mod <- lm(log(freeze)~bins, data=data)
      y_new <- exp(coefficients(mod)[1] + coefficients(mod)[2]*time) 
      p <- p + geom_line(data=data.frame(time, y_new), mapping=aes(x = time, y = y_new),
                         col = col, alpha = 0.2)
    }
    
    #Replot the sustained animals average cruve so that it is overlayed
    y_new <- exp(coefficients(mod.sus)[1] + coefficients(mod.sus)[2]*time)
    p <- p + geom_line(data=data.frame(time, y_new), mapping=aes(x = time, y = y_new),
                       col = "#c300b4ff", size=1)
  }
  
  
  #Add phasic animals
  data_ph <- colMeans(df[class =="phasic", ])
  data.plot.ph <- data.frame(data_ph, bins)
  colnames(data.plot.ph) <- c("freeze", "bins")
  if (0 %in% data.plot.ph$freeze) {data.plot.ph$freeze <- data.plot.ph$freeze + 1}
  mod.ph <- lm(log(freeze)~bins, data=data.plot.ph)
  y_new <- exp(coefficients(mod.ph)[1] + coefficients(mod.ph)[2]*time)
  
  p <- p + geom_line(data=data.frame(time, y_new), 
                     mapping=aes(x = time, y = y_new), col =  "#006400ff", size = 1)
  
  if ("shifter" %in% class){
    data_shift <- colMeans(df[class =="shifter", ])
    data.plot.shift <- data.frame(data_shift, bins)
    colnames(data.plot.shift) <- c("freeze", "bins")
    if (0 %in% data.plot.shift$freeze) {data.plot.shift$freeze <- data.plot.shift$freeze + 1}
    mod.shift <- lm(log(freeze)~bins, data=data.plot.shift)
    y_new <- exp(coefficients(mod.shift)[1] + coefficients(mod.shift)[2]*time)
    
    p <- p + geom_line(data=data.frame(time, y_new), 
                       mapping=aes(x = time, y = y_new), col = "gray50", size=1)
  }
  
  
  ymax = layer_scales(p)$y$range$range[2]
  
  
  #Add goodness of fit statistics for the average freezing curves
  if (stat == TRUE){
    r2.sus <- round(summary(mod.sus)[9]$adj.r.squared,2)
    p.sus <- summary(mod.sus)$coefficients[2,4]
    
    if (p.sus <0.001){
      p.sus <- "p<0.001"
    }else{p.sus <- paste0("p=", round(p.sus),3)}
    
    #Phasic
    r2.ph <- round(summary(mod.ph)[9]$adj.r.squared,2)
    p.ph <- summary(mod.ph)$coefficients[2,4]
    
    if (p.ph <0.001){
      p.ph <- "p<0.001"
    }else{p.ph <- paste0("p=", round(p.ph),3)}
    
    
    lab1 <- bquote(paste("R"^2*"=", .(r2.sus), ", ", .(p.sus)))
    
    p <- p + annotate(geom="text", x = 3.5, y = 1.1*ymax, label = lab1, color = "magenta", 
                      size = 4, fontface=3)
    
    lab2 <-  bquote(paste("R"^2*"=", .(r2.ph), ", ", .(p.ph)))
    
    p <- p + annotate(geom="text", x = 9.5, y = 1.1*ymax, label = lab2, color = "green3", 
                      size = 4, fontface=3)
    
    
    if ("shifter" %in% class){
      
      r2.shift <- round(summary(mod.shift)[9]$adj.r.squared,2)
      p.shift <- summary(mod.shift)$coefficients[2,4]
      
      if (p.shift <0.001){
        p.shift <- "p<0.001"
      }else{p.shift <- paste0("p=", round(p.shift),3)}
      
      lab3 <-  bquote(paste("R"^2*"=", .(r2.shift), ", ", .(p.shift)))
      
      p <- p + annotate(geom="text", x = 6.5, y = 1.2*ymax, label = lab3, color = "gray50", 
                        size = 4, fontface=3)
      
    }
    
    p = p + ggplot2::scale_y_continuous(expand = c(0.005,1), limits = c(-5,1.3*ymax))
  }else(p = p + ggplot2::scale_y_continuous(expand = c(0.005,1), limits = c(-5, ymax)))
  
  p = p + ggplot2::geom_rect(aes(xmin=1, xmax = 6, ymin = -5, ymax = -1), 
                             color= "green2", fill = "green2")
  
  p = p + ggplot2::geom_rect(aes(xmin=6, xmax = 12, ymin = -5, ymax = -1), 
                             color= "magenta", fill = "magenta")
  p = p + annotate(geom = "text", x = 6, y = -2.5, label= c("C", "S"),
                   fontface = 2, size =3, hjust = c(1,0))
  return(p)
}

#-------------------------------------------------------------------------------

#A more general version of the function to plot freezing curves that should world
#For any 2-group assignment
freezing_curves_general = function(df, class, animal, retrieval,  
                                   individual = FALSE, stat = FALSE,
                                   shape = c(16,16), color = c("magenta", "green2")){
  
  bins_all = 1:36
  
  #Calculate average freezing per time bin and group
  mean1 = base::colMeans(df[class==unique(class)[1],], na.rm = TRUE)
  se1 = base::apply(df[class==unique(class)[1], ], 2, function(x)sd(na.omit(x))/sqrt(length(na.omit(x))))
  
  mean2 = base::colMeans(df[class==unique(class)[2],], na.rm = TRUE)
  se2 = base::apply(df[class==unique(class)[2], ], 2, function(x)sd(na.omit(x))/sqrt(length(na.omit(x))))
  
  
  #-----
  data.plot = data.frame(mean1, se1, bins_all)
  colnames(data.plot) = c("freeze", "se", "bins")
  
  
  p1 = ggplot2::ggplot(data.plot, aes(bins, freeze))  + theme_test(base_size=19) +
    labs(x = "Bin", y = paste0("Average freezing per bin ", retrieval,  " [%]")) + 
    geom_line(col = color[1], linewidth=1) + geom_point(size=2.5, col=color[1], shape=shape[1]) + 
    geom_errorbar(aes(ymin=freeze-se, ymax=freeze+se), width=.2,
                  position=position_dodge(0.05), col = color[1]) +
    scale_x_continuous(breaks = c(1,13,18,24,36), labels = c(1,13,18,24,36)) + 
    scale_y_continuous(expand = c(0,1), limits = c(-4,100)) +
    theme(axis.text = element_text(colour="black")) 
  
  
  
  #Add individual curves if individual == TRUE
  
  if (individual == TRUE){
    
    df1 = df[class == unique(class)[1],]
    df2 = df[class == unique(class)[2],]
    
    
    #Add individual curves sustained
    for (i in 1:nrow(df1)){
      data = data.frame(t(df1[i,]), bins_all)
      colnames(data) = c("freeze", "bins")
      
      p1 = p1 + geom_line(col=color[1], data=data, alpha = 0.1) + 
        geom_point(size=1.5, col=color[1], data=data, alpha = 0.1, shape=shape[1])
    }
    
    #Add individual curves phasic
    for (i in 1:nrow(df2)){
      data = data.frame(t(df2[i,]), bins_all)
      colnames(data) = c("freeze", "bins")
      
      p1 = p1 + geom_line(col=color[2], data=data, alpha = 0.1) + 
        geom_point(size=1.5, col=color[2], data=data, alpha = 0.1, shape=shape[2])
    }
    
    
    
    #Replot average curve for group 1
    p1 = p1 + geom_line(col=color[1], data=data.plot, linewidth=1) + 
      geom_point(size=2.5, col=color[1], data=data.plot, shape=shape[1]) + 
      geom_errorbar(data= data.plot, mapping=aes(ymin=freeze-se, ymax=freeze+se), width=.2,
                    position=position_dodge(0.05), col = color[1]) 
  }
  
  #Add group 2 animals
  data.plot2 = data.frame(mean2, se2, bins_all)
  colnames(data.plot2) = c("freeze", "se", "bins")
  p1 = p1 + geom_line(col=color[2], data=data.plot2, linewidth=1) + 
    geom_point(size=2.5, col=color[2], data=data.plot2, shape=shape[2]) + 
    geom_errorbar(data= data.plot2, mapping=aes(ymin=freeze-se, ymax=freeze+se), width=.2,
                  position=position_dodge(0.05), col = color[2]) 
  
  
  #Add a line for region designating the tone
  p1 = p1 + ggplot2::geom_rect(aes(xmin=13, xmax = 18, ymin = -4, ymax = -2), 
                               color= "green2", fill = "green2")
  
  p1 = p1 + ggplot2::geom_rect(aes(xmin=18, xmax = 24, ymin = -4, ymax = -2), 
                               color= "magenta", fill = "magenta")
  p1 = p1 + annotate(geom = "text", x = c(17.9, 18.1), y = -3, label= c("C", "S"),
                     fontface = 2, size =3, hjust = c(1,0))
  
  #############################################
  #Perform linear mixed effects analysis comparing sustained versus phasic animals
  
  if (stat == TRUE){
    df$animal = animal
    df$group = class
    
    
    #Melt to long format
    df_stat = reshape2::melt(df, id.vars = c("animal", "group"))
    colnames(df_stat)[3:4] = c("bin", "freezing")
    
    #Fit a linear model
    lmm = lme4::lmer(freezing~group*bin + (1|animal), data=df_stat)
    #Compare marginal means
    em = emmeans::emmeans(lmm, pairwise~group|bin)
    
    #Extract p-values
    p.vals = summary(em$contrasts)[7]
    
    #Plot significant p-values
    if(all(p.vals$p.value>0.05)==F){
      x = which(p.vals$p.value<0.05)
      
      p.vals$symbol = ifelse(p.vals$p.value<0.05, "*", "-")
      p.vals$symbol = ifelse(p.vals$p.value<0.01, "**", p.vals$symbol)
      p.vals$symbol = ifelse(p.vals$p.value<0.001, "***", p.vals$symbol)
      
      label = p.vals$symbol[x]
      
      #Get y coordinate for the statistics
      
      y = mean1[x] + se1[x] + 1
      
      p1 = p1 + annotate(geom="text", x = x, y = y, label = label,
                         size=3.5, fontface=2)
      
    }
  }
  
  
  return(p1)
}


########################################################################################
# FUNCTION TO PLOT GAUSSIAN DISTRIBUTIONS OF PARAMETERS FOR SUSTAINED AND PHASIC FREEZERS
#---------------------------------------------------------------------------------------
#The function creates histograms and density plots of freezing parameters for 
#Sustained and phasic freezers

plot_density <- function(variable, class, xlab, title, n, breaks){
  data <- data.frame(variable, class)
  colnames(data) <- c("var", "class")
  data$class <- factor(data$class)
  bw <- round(max(variable),1)/30
  p <- ggplot(data, aes(x=var, y = ..density..))  + 
    geom_histogram(fill = "gray87", col = "black", binwidth = bw, show.legend=F) + 
    theme_test(base_size = 14) + labs(x=xlab, y = "Density", col = "") + 
    geom_density(aes(col = class), size=1.5) +
    theme(axis.text = element_text(colour="black", size=12)) +
    scale_color_manual(values=c("green3", "magenta"), labels = c("phasic", "sustained")) + 
    ggtitle(title) + 
    theme(plot.title = element_text(hjust = 0.5)) + 
    scale_y_continuous(expand=expansion(mult = c(0,0.05)), sec.axis = sec_axis(trans=~.*bw*n, name = "Counts", breaks = breaks)) +
    scale_x_continuous(expand=(c(0,0)))
  
  return(p)
  
}


#-------------------------------------------------------------------------------

#Function to train several ML Models and validate them via Monte-Carlo Cross validation
#BINARY CLASSIFICATION

#Input params:
#data: Dataframe to train with - prediction column should be named "class"
#ALL other columns will be used to train, exclude batch or id column before
train_binary <- function(data, steps=1000, title, shape){
  message("Training models ...")
  
  #Init progress bar for console
  pb = txtProgressBar(min = 0, max = steps, initial = 0) 
  
  #Declare target variable as factor
  data$class <- factor(data$class, levels = c("sustained", "phasic"))
  
  #Initialize empty data.frames to collect all performance metrics
  svm <- data.frame()
  rsvm <- data.frame()
  log <- data.frame()
  lda <- data.frame()
  rf <- data.frame()
  
  #Monte Carlo Cross Validation: Use Random Subsampling - 500 iterations as default
  iterations <- steps
  
  suppressMessages(for (i in 1:iterations){
    
    #Update Progress Bar
    setTxtProgressBar(pb,i)
    
    #Randomly split data into training and test set
    #use 70% of data set as training set and 30% as test set
    split <- base::sample(c(TRUE, FALSE), nrow(data), replace=TRUE, prob=c(0.7,0.3))
    training  <- data[split, ]
    testing   <- data[!split, ]
    
    metrics <- c()
    model_names <- c("SVM _linear", "SVM_radial", "logistic", "lda", "rf")
    
    #-------------------------
    #Fit SVM Model using linear kernel and no scaling
    svm_model <- e1071::svm(class~., data = training, kernel = "linear")
    
    #Obtain predictions
    svm_pred <- stats::predict(svm_model, newdata = testing)
    svm_cm <- table(testing$class, svm_pred)
    
    #Calculate performance metrics 
    svm_metrics <- caret::confusionMatrix(svm_cm, mode = "everything")
    svm_auroc <- pROC::auc(unclass(testing$class), unclass(svm_pred))
    
    #Get metrics
    svm_all_metrics <- c()
    svm_all_metrics[1] <- svm_metrics$overall["Accuracy"]
    svm_all_metrics[2] <- svm_metrics$byClass["Sensitivity"]
    svm_all_metrics[3] <- svm_metrics$byClass["Specificity"]
    svm_all_metrics[4] <- svm_metrics$byClass["F1"]
    svm_all_metrics[5] <- svm_auroc
    
    #-------------------------
    #Fit SVM Model using radial kernel and no scaling
    rsvm_model <- e1071::svm(class~., data = training, kernel = "radial")
    
    #Obtain predictions
    rsvm_pred <- stats::predict(rsvm_model, newdata = testing)
    rsvm_cm <- table(testing$class, rsvm_pred)
    
    #Calcuate performance metrics 
    rsvm_metrics <- caret::confusionMatrix(rsvm_cm, mode = "everything")
    rsvm_auroc <- pROC::auc(unclass(testing$class), unclass(rsvm_pred))
    
    #Get metrics
    rsvm_all_metrics <- c()
    rsvm_all_metrics[1] <- rsvm_metrics$overall["Accuracy"]
    rsvm_all_metrics[2] <- rsvm_metrics$byClass["Sensitivity"]
    rsvm_all_metrics[3] <- rsvm_metrics$byClass["Specificity"]
    rsvm_all_metrics[4] <- rsvm_metrics$byClass["F1"]
    rsvm_all_metrics[5] <- rsvm_auroc
    
    #-------------------------
    #Fit GLM Model for logistic Regression
    log_model <- suppressWarnings(stats::glm(class~., data =training, family = "binomial"))
    
    #Obtain predictions
    log_pred <- stats::predict(log_model, newdata = testing)
    
    #Convert values to phenotypes -> is there a more generalized approach?
    log_pred <- ifelse(log_pred>=0.5, "phasic", "sustained")
    log_pred <- factor(log_pred, levels = c("sustained", "phasic"))
    #Calcuate performance metrics 
    log_metrics <- caret::confusionMatrix(data=log_pred, reference = testing$class,mode = "everything")
    log_auroc <- pROC::auc(unclass(testing$class), unclass(factor(log_pred)))
    
    #Get metrics
    log_all_metrics <- c()
    log_all_metrics[1] <- log_metrics$overall["Accuracy"]
    log_all_metrics[2] <- log_metrics$byClass["Sensitivity"]
    log_all_metrics[3] <- log_metrics$byClass["Specificity"]
    log_all_metrics[4] <- log_metrics$byClass["F1"]
    log_all_metrics[5] <- log_auroc
    
    #-------------------------
    #Fit linear discriminant analysis models
    lda_model <- MASS::lda(class~., data =training)
    
    #Obtain predictions
    lda_pred <- stats::predict(lda_model, newdata = testing)
    #Get Confusion matrix of the lda model
    lda_cm <- table(testing$class, lda_pred$class)
    
    #Calcuate performance metrics 
    lda_metrics <- caret::confusionMatrix(lda_cm, mode = "everything")
    lda_auroc <- pROC::auc(unclass(testing$class), unclass(factor(lda_pred$class)))
    
    
    #Get metrics
    lda_all_metrics <- c()
    lda_all_metrics[1] <- lda_metrics$overall["Accuracy"]
    lda_all_metrics[2] <- lda_metrics$byClass["Sensitivity"]
    lda_all_metrics[3] <- lda_metrics$byClass["Specificity"]
    lda_all_metrics[4] <- lda_metrics$byClass["F1"]
    lda_all_metrics[5] <- lda_auroc
    
    #-------------------------
    #Fit random Forests
    rf_model <- randomForest::randomForest(class~., data =training)
    
    #Obtain predictions
    rf_pred <- stats::predict(rf_model, newdata = testing)
    
    #Get Confusion matrix of the random forest model
    rf_cm <- table(testing$class, rf_pred)
    
    #Calcuate performance metrics 
    rf_metrics <- caret::confusionMatrix(rf_cm, mode = "everything")
    rf_auroc <- pROC::auc(unclass(testing$class), unclass(rf_pred))
    
    #Get metrics
    rf_all_metrics <- c()
    rf_all_metrics[1] <- rf_metrics$overall["Accuracy"]
    rf_all_metrics[2] <- rf_metrics$byClass["Sensitivity"]
    rf_all_metrics[3] <- rf_metrics$byClass["Specificity"]
    rf_all_metrics[4] <- rf_metrics$byClass["F1"]
    rf_all_metrics[5] <- rf_auroc
    
    #Add new metrics to dfs
    svm <- rbind(svm, svm_all_metrics)
    rsvm <- rbind(rsvm, rsvm_all_metrics)
    log <- rbind(log, log_all_metrics)
    lda <- rbind(lda, lda_all_metrics)
    rf <- rbind(rf, rf_all_metrics)
  })
  close(pb)
  
  #Assign Columnnames of Metric dataframes
  colnames(svm) <- c("Accuracy", "Sensitivity", "Specificity", "F1", "Auroc")
  colnames(rsvm) <- c("Accuracy", "Sensitivity", "Specificity", "F1", "Auroc")
  colnames(log) <- c("Accuracy", "Sensitivity", "Specificity", "F1", "Auroc")
  colnames(lda) <- c("Accuracy", "Sensitivity", "Specificity", "F1", "Auroc")
  colnames(rf) <- c("Accuracy", "Sensitivity", "Specificity", "F1", "Auroc")
  
  #----------------------------
  #Gather mean and standard deviation of all metrics
  acc.stats <- data.frame(SVM=paste(round(mean(svm$Accuracy), digits = 6), "±", round(sd(svm$Accuracy), digits = 6)), 
                          RSVM=paste(round(mean(rsvm$Accuracy), digits = 6), "±", round(sd(rsvm$Accuracy), digits = 6)),
                          LOG=paste(round(mean(log$Accuracy), digits = 6), "±", round(sd(log$Accuracy), digits = 6)),
                          LDA=paste(round(mean(lda$Accuracy), digits = 6), "±", round(sd(lda$Accuracy), digits = 6)),
                          RandomForest=paste(round(mean(rf$Accuracy), digits = 6), "±", round(sd(rf$Accuracy), digits = 6)))
  
  
  
  sens.stats <- data.frame(SVM=paste(round(mean(svm$Sensitivity), digits = 6), "±", round(sd(svm$Sensitivity), digits = 6)), 
                          RSVM=paste(round(mean(rsvm$Sensitivity), digits = 6), "±", round(sd(rsvm$Sensitivity), digits = 6)),
                          LOG=paste(round(mean(log$Sensitivity), digits = 6), "±", round(sd(log$Sensitivity), digits = 6)),
                          LDA=paste(round(mean(lda$Sensitivity), digits = 6), "±", round(sd(lda$Sensitivity), digits = 6)),
                          RandomForest=paste(round(mean(rf$Sensitivity), digits = 6), "±", round(sd(rf$Sensitivity), digits = 6)))
  
  spec.stats <- data.frame(SVM=paste(round(mean(svm$Specificity), digits = 6), "±", round(sd(svm$Specificity), digits = 6)), 
                           RSVM=paste(round(mean(rsvm$Specificity), digits = 6), "±", round(sd(rsvm$Specificity), digits = 6)),
                           LOG=paste(round(mean(log$Specificity), digits = 6), "±", round(sd(log$Specificity), digits = 6)),
                           LDA=paste(round(mean(lda$Specificity), digits = 6), "±", round(sd(lda$Specificity), digits = 6)),
                           RandomForest=paste(round(mean(rf$Specificity), digits = 6), "±", round(sd(rf$Specificity), digits = 6)))

  f1.stats <- data.frame(SVM=paste(round(mean(svm$F1), digits = 6), "±", round(sd(svm$F1), digits = 6)), 
                           RSVM=paste(round(mean(rsvm$F1), digits = 6), "±", round(sd(rsvm$F1), digits = 6)),
                           LOG=paste(round(mean(log$F1), digits = 6), "±", round(sd(log$F1), digits = 6)),
                           LDA=paste(round(mean(lda$F1), digits = 6), "±", round(sd(lda$F1), digits = 6)),
                           RandomForest=paste(round(mean(rf$F1), digits = 6), "±", round(sd(rf$F1), digits = 6)))  
  
  auroc.stats <- data.frame(SVM=paste(round(mean(svm$Auroc), digits = 6), "±", round(sd(svm$Auroc), digits = 6)), 
                         RSVM=paste(round(mean(rsvm$Auroc), digits = 6), "±", round(sd(rsvm$Auroc), digits = 6)),
                         LOG=paste(round(mean(log$Auroc), digits = 6), "±", round(sd(log$Auroc), digits = 6)),
                         LDA=paste(round(mean(lda$Auroc), digits = 6), "±", round(sd(lda$Auroc), digits = 6)),
                         RandomForest=paste(round(mean(rf$Auroc), digits = 6), "±", round(sd(rf$Auroc), digits = 6)))  
  
  metrics <- rbind(acc.stats, sens.stats, spec.stats, f1.stats, auroc.stats)
  rownames(metrics) <- c("Accuracy", "Sensitivity", "Specificity", "F1", "Auroc")
  
  #----------------------------
  #Plot Accuracies
  accuracies <- data.frame(SVM=svm$Accuracy, RSVM=rsvm$Accuracy, LogisticRegression=log$Accuracy,
                           LDA=lda$Accuracy, randomForest=rf$Accuracy)
  
  accuracies.long <- suppressMessages(reshape2::melt(accuracies))
  
  plot.acc <- ggplot(data=accuracies.long, aes(x = variable, y = value, colour = variable)) +
    geom_boxplot(outlier.shape = NA)  + 
    theme(axis.text = element_text(colour="black", size=12)) + theme_test(base_size=14) +
    labs(y= "Accuracy", x="") + 
    scale_x_discrete(labels=c("linear SVM", "radial SVM", "logistic regression", "LDA", "random Forest")) + 
    geom_jitter(size=0.5, width = 0.2, show.legend = F, shape = shape) + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") +
    scale_color_aaas() + scale_fill_aaas() + 
    theme(axis.text.x = element_text(angle=35, hjust=1))
  
  #----------------------------
  #Plot Sensitivity
  
  sensitivities <- data.frame(SVM=svm$Sensitivity, RSVM=rsvm$Sensitivity, LogisticRegression=log$Sensitivity,
                              LDA=lda$Sensitivity, randomForest=rf$Sensitivity)
  
  sensitivities.long <- suppressMessages(reshape2::melt(sensitivities))
  
  plot.sens <- ggplot(data=sensitivities.long, aes(x = variable, y = value, colour = variable)) +
    geom_boxplot(outlier.shape = NA) + theme_test(base_size = 14) +
    labs(y= "Sensitivity", x="") + 
    scale_x_discrete(labels=c("linear SVM", "radial SVM", "logistic regression", "LDA", "random Forest")) + 
    ggtitle(label=paste("Sensitivity", title, sep = " ")) + 
    geom_jitter(size=1, width = 0.2, show.legend = F) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") +
    scale_color_aaas() + scale_fill_aaas() + 
    theme(axis.text.x = element_text(angle=35, hjust=1))
  
  #----------------------------
  #Plot Specificity 
  
  specificities <- data.frame(SVM=svm$Specificity, RSVM=rsvm$Specificity, LogisticRegression=log$Specificity,
                              LDA=lda$Specificity, randomForest=rf$Specificity)
  
  specificities.long <- suppressMessages(reshape2::melt(specificities))
  
  plot.spec <- ggplot(data=specificities.long, aes(x = variable, y = value, colour = variable)) +
    geom_boxplot(outlier.shape = NA)  + theme_test(base_size = 14) +
    labs(y= "Specificity", x="") + scale_x_discrete(labels=c("linear SVM", "radial SVM", "logistic regression", "LDA", "random Forest")) + 
    ggtitle(label=paste("Specificity", title, sep = " ")) +
    geom_jitter(size=1, width = 0.2, show.legend = F) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") +
    scale_color_aaas() + scale_fill_aaas() + 
    theme(axis.text.x = element_text(angle=35, hjust=1))  
  
  #----------------------------
  #Plot F1 Scores
  
  f1 <- data.frame(SVM=svm$F1, RSVM=rsvm$F1, LogisticRegression=log$F1,
                   LDA=lda$F1, randomForest=rf$F1)
  
  f1.long <- suppressMessages(reshape2::melt(f1))
  
  plot.f1 <- ggplot(data=f1.long, aes(x = variable, y = value, colour = variable)) +
    geom_boxplot(outlier.shape = NA)  + theme_test(base_size = 14) +
    labs(y= "F1", x="") + scale_x_discrete(labels=c("linear SVM", "radial SVM", "logistic regression", "LDA", "random Forest")) + 
    ggtitle(label=paste("F1", title, sep = " ")) + 
    geom_jitter(size=1, width = 0.2, show.legend = F) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") +
    scale_color_aaas() + scale_fill_aaas() + 
    theme(axis.text.x = element_text(angle=35, hjust=1)) 
  
  #----------------------------
  #Plot Auroc Scores
  aurocs <- data.frame(SVM=svm$Auroc, RSVM=rsvm$Auroc, LogisticRegression=log$Auroc,
                       LDA=lda$Auroc, randomForest=rf$Auroc)
  
  aurocs.long <- suppressMessages(reshape2::melt(aurocs))
  
  plot.aurocs <- ggplot(data=aurocs.long, aes(x = variable, y = value, colour = variable)) +
    geom_boxplot(outlier.shape = NA)  + theme_test(base_size = 14) +
    labs(y= "Auroc", x="") + 
    scale_x_discrete(labels=c("linear SVM", "radial SVM", "logistic regression", "LDA", "random Forest")) + 
    ggtitle(label=paste("Auroc", title, sep = " "))+ geom_jitter(size=1, width = 0.2, show.legend = F) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") +
    scale_color_aaas() + scale_fill_aaas() + 
    theme(axis.text.x = element_text(angle=35, hjust=1))
  
  ##############################################################################
  #Train final models on all data
  SVM <- suppressMessages(e1071::svm(class~., data = data, kernel = "linear"))
  RSVM <- suppressMessages(e1071::svm(class~., data = data, kernel = "radial"))
  Logistic <- suppressWarnings(stats::glm(class~., data = data, family = "binomial"))
  LDA <- suppressMessages(MASS::lda(class~., data = data))
  RF <- suppressMessages(randomForest::randomForest(class~., data = data))
  
  message("Finished Training successfully")
  #Return plots and Models
  return(list(Accuracy=plot.acc, Sensitivity=plot.sens, Specificity=plot.spec,
              F1=plot.f1, Auroc=plot.aurocs, SVM=SVM, RSVM=RSVM, Logistic=Logistic,
              LDA=LDA, RF=RF, Stats=metrics))
}

################################################################################
#Function to train several ML Models and validate them via Monte-Carlo Crossvalidation
#NON BINARY MULTINOMIAL CLASSIFICATION

#Input params:
#data: Dataframe to train with - prediction column should be named "class"
#Class column should be a factor
#ALL other columns will be used to train, exclude batch or id column before

train_shifter <- function(data, steps=500, title, shape){
  
  shifter <- data
  
  #Initialize empty data.frames to collect all performance metrics
  svm.shifter <- c()
  rsvm.shifter <- c()
  log.shifter <- c()
  lda.shifter <- c()
  rf.shifter <- c()
  
  invisible(for (i in 1:steps){
    #Randomly split data into training and test set
    #use 70% of data set as training set and 30% as test set
    split <- base::sample(c(TRUE, FALSE), nrow(shifter), replace=TRUE, prob=c(0.7,0.3))
    training  <- shifter[split, ]
    testing   <- shifter[!split, ]
    
    metrics <- c()
    model_names <- c("SVM _linear", "SVM_radial", "logistic", "lda", "rf")
    
    #Fit SVM Model using linear kernel and no scaling
    svm_model <- e1071::svm(class~., data = training, kernel = "linear")
    
    #Obtain predictions
    svm_pred <- stats::predict(svm_model, newdata = testing)
    svm_cm <- table(testing$class, svm_pred)
    
    #Calculate accuracy
    svm_metrics <- caret::confusionMatrix(svm_cm, mode = "everything")
    
    #Get metrics
    svm_accuracy <- svm_metrics$overall["Accuracy"]
    
    #-------------------------
    #Fit SVM Model using radial kernel and no scaling
    rsvm_model <- e1071::svm(class~., data = training, kernel = "radial")
    
    #Obtain predictions
    rsvm_pred <- stats::predict(rsvm_model, newdata = testing)
    rsvm_cm <- table(testing$class, rsvm_pred)
    
    #Calcuate performance metrics 
    rsvm_metrics <- caret::confusionMatrix(rsvm_cm, mode = "everything")
    rsvm_auroc <- pROC::auc(unclass(testing$class), unclass(rsvm_pred))
    
    #Get metrics
    rsvm_accuracy <- rsvm_metrics$overall["Accuracy"]
    
    #-------------------------
    #Fit multinomial logistic regression Model using neural net package
    log_model <- suppressMessages(nnet::multinom(class ~., data = training))
    
    #Obtain predictions
    log_pred <- stats::predict(log_model, newdata = testing)
    
    #Get Confusion matrix of the logisticc regression model
    log_cm <- table(testing$class, log_pred)
    
    #Calcuate performance metrics 
    log_metrics <- caret::confusionMatrix(data=log_pred, reference = testing$class,mode = "everything")
    
    log_accuracy <- log_metrics$overall["Accuracy"]
    
    #-------------------------
    #Fit linear discriminant analysis models
    lda_model <- MASS::lda(class~., data =training)
    
    #Obtain predictions
    lda_pred <- stats::predict(lda_model, newdata = testing)
    #Get Confusion matrix of the lda model
    lda_cm <- table(testing$class, lda_pred$class)
    
    #Calcuate performance metrics 
    lda_metrics <- caret::confusionMatrix(lda_cm, mode = "everything")
    
    lda_accuracy <- lda_metrics$overall["Accuracy"]
    
    #-------------------------
    #Fit random Forests
    rf_model <- randomForest::randomForest(class~., data =training)
    
    #Obtain predictions
    rf_pred <- stats::predict(rf_model, newdata = testing)
    
    #Get Confusion matrix of the random forest model
    rf_cm <- table(testing$class, rf_pred)
    
    #Calcuate performance metrics 
    rf_metrics <- caret::confusionMatrix(rf_cm, mode = "everything")
    
    #Get metrics
    rf_accuracy <- rf_metrics$overall["Accuracy"]
    
    #add accuracies to total vectors
    svm.shifter[i] <- svm_accuracy
    rsvm.shifter[i] <- rsvm_accuracy
    log.shifter[i] <- log_accuracy
    lda.shifter[i] <- lda_accuracy
    rf.shifter[i] <- rf_accuracy
  })
  
  #Plot Accuracies
  accuracies <- data.frame(SVM=svm.shifter, RSVM=rsvm.shifter, LogisticRegression=log.shifter,
                           LDA=lda.shifter, randomForest=rf.shifter)
  
  acc.stats <- data.frame(SVM=paste(round(mean(svm.shifter), digits = 6), "±", round(sd(svm.shifter), digits = 6)), 
                          RSVM=paste(round(mean(rsvm.shifter), digits = 6), "±", round(sd(rsvm.shifter), digits = 6)),
                          LOG=paste(round(mean(log.shifter), digits = 6), "±", round(sd(log.shifter), digits = 6)),
                          LDA=paste(round(mean(lda.shifter), digits = 6), "±", round(sd(lda.shifter), digits = 6)),
                          RandomForest=paste(round(mean(rf.shifter), digits = 6), "±", round(sd(rf.shifter), digits = 6)))
  
  rownames(acc.stats) <- c("Accuracy")
  
  accuracies.long <- reshape2::melt(accuracies)
  
  plot.acc <- ggplot(data=accuracies.long, aes(x = variable, y = value, colour = variable)) +
    geom_boxplot(outlier.shape = NA)  + 
    theme(axis.text = element_text(colour="black", size=12)) + theme_test(base_size=14) +
    labs(y= "Accuracy", x="") + scale_x_discrete(labels=c("linear SVM", "radial SVM", "logistic regression", "LDA", "random Forest")) +
    geom_jitter(size=1, width = 0.2, show.legend = F, shape = shape) + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") +
    scale_color_aaas() + scale_fill_aaas() + 
    theme(axis.text.x = element_text(angle=35, hjust=1))  
  
  ##############################################################################
  #Train final models on all data
  SVM <- suppressMessages(e1071::svm(class~., data = shifter, kernel = "linear"))
  RSVM <- suppressMessages(e1071::svm(class~., data = shifter, kernel = "radial"))
  Logistic <- suppressWarnings(nnet::multinom(class ~., data = shifter))
  LDA <- suppressMessages(MASS::lda(class~., data = shifter))
  RF <- suppressMessages(randomForest::randomForest(class~., data = shifter))
  
  message("Finished Training successfully")
  #Return plots and Models
  return(list(Accuracy=plot.acc, SVM=SVM, RSVM=RSVM, Logistic=Logistic,
              LDA=LDA, RF=RF, Stats=acc.stats))
}

#-------------------------------------------------------------------------------
#Function for bootstrapping analysis

bootstrapping <- function(data, data.freezing, iterations=200, retrieval, shape){
  
  #Merge Batch and ID for freezing df
  B_ID <- c()
  for (i in 1:nrow(data.freezing)){
    B_ID[i] <- paste(data.freezing$batch[i], data.freezing$id[i], sep = "-")
  }
  rownames(data.freezing) <- B_ID
  
  #Get original cluster assignments
  original_p <- data[data$class=="phasic",1]
  original_s <- data[data$class=="sustained",1]
  
  jaccards_p <- c()
  jaccards_s <- c()
  
  #Get animal B-IDs
  animals <- data$B_ID
  
  #Initialize iteration dataframe
  bootstrap_results <- data.frame(matrix(nrow = length(animals), ncol = iterations)) 
  bootstrap_results[is.na(bootstrap_results)] <- "no"
  
  #Initialize dataframes for freezing plots
  cluster_p <-  data.frame(matrix(nrow = iterations, ncol = 36)) 
  cluster_s <-  data.frame(matrix(nrow = iterations, ncol = 36)) 
  
  rownames(bootstrap_results) <- animals
  colnames(bootstrap_results) <- c(1:iterations)
  
  #Collect Average Freezing of Original cluster group
  freezing_p_original <- data.freezing[original_p,]
  freezing_s_original <- data.freezing[original_s,]
  
  average_p_original <- colMeans(freezing_p_original[,3:38], na.rm = TRUE)
  average_s_original <- colMeans(freezing_s_original[,3:38], na.rm = TRUE)
  
  for (step in 1:iterations){
    
    #sample the dataset with replacement
    rand_sample <- sample(data$B_ID , size = nrow(data), replace = TRUE)
    
    #Initialize new dataframe for the subsampled animals
    rand_data <- data.frame(B_ID = c(), freeze = c(), beta = c(), int = c(), class = c())
    
    #Add all sampled animals to the data frame
    for (i in 1:length(rand_sample)){
      current_row <- data[data$B_ID==rand_sample[i],]
      rand_data <- rbind(rand_data, current_row)
    }
    
    #perform GMM clustering on subsampled set
    mclust <- mclust::Mclust(data = rand_data[,2:4], G=2)
    rand_data$class_bootstrap <- assign.phenotype(rand_data$freeze,
                                                  mclust$classification)
    
    
    #Subset to phasic and sustained, collect B_ID
    bootstrap_p <- rand_data[rand_data$class_bootstrap=="phasic",1]
    bootstrap_s <- rand_data[rand_data$class_bootstrap=="sustained",1]
    
    #Collect Average Freezing of bootstrap cluster group
    freezing_p <- data.freezing[bootstrap_p,]
    freezing_s <- data.freezing[bootstrap_s,]
    
    average_p <- colMeans(freezing_p[,3:38], na.rm = TRUE)
    average_s <- colMeans(freezing_s[,3:38], na.rm = TRUE)
    
    cluster_p[step,] <-  average_p
    cluster_s[step,] <-  average_s
    
    
    for (i in bootstrap_p){
      if (i %in% original_p){
        bootstrap_results[i,step] <- "retained"
      } else {
        bootstrap_results[i,step] <- "swapped"
      }
    }
    
    for (i in bootstrap_s){
      if (i %in% original_s){
        bootstrap_results[i,step] <- "retained"
      } else {
        bootstrap_results[i,step] <- "swapped"
      }
    }
  }
  
  #Initialize iteration dataframe
  jaccard_s <- c()
  jaccard_p <- c()
  
  for (i in 1:length(original_p)){
    #Get row of current animals as a character vector
    individual <- as.character(bootstrap_results[original_p[i],])
    #Count all instances
    counted <- table(individual)
    no <- as.numeric(counted["no"])
    retained <- as.numeric(counted["retained"])
    swapped <- as.numeric(counted["swapped"])
    #Calculate Jaccard Index
    jaccard_p[i] <- retained/(iterations-no)
  }
  
  for (i in 1:length(original_s)){
    #Get row of current animals as a character vector
    individual <- as.character(bootstrap_results[original_s[i],])
    #Count all instances
    counted <- table(individual)
    no <- as.numeric(counted["no"])
    retained <- as.numeric(counted["retained"])
    swapped <- as.numeric(counted["swapped"])
    #Calculate Jaccard Index
    jaccard_s[i] <- retained/(iterations-no)
  }
  
  #Prepare data for plotting
  res_p <- data.frame(Jaccard=jaccard_p)
  res_s <- data.frame(Jaccard=jaccard_s)
  
  mean_p <- colMeans(res_p)
  mean_s <- colMeans(res_s)
  
  res_p$Cluster <- "phasic"
  res_s$Cluster <- "sustained"
  
  res_all <- rbind(res_p, res_s)
  
  #Plot distributions of Jaccard indices per original cluster
  plot.jac <- ggplot(data=res_all, aes(x = Cluster, y = Jaccard, colour = Cluster)) +
    theme_test(base_size = 14) + geom_boxplot() +
    labs(y= "Jaccard index", x="") + scale_y_continuous(limits = c(0, 1.01), breaks = seq(0, 1, by = 0.2)) +
    geom_jitter(size=1.5, width = 0.2, show.legend = F, shape = shape) + 
    scale_color_manual(values = c("green2", "magenta")) +
    theme(axis.text = element_text(colour="black", size=12)) 
    
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") 
  
  #Plot average bootstrap sample freezing per cluster
  bins_all = 1:36
  
  bootstrap_mean_p <- colMeans(cluster_p)
  bootstrap_mean_s <- colMeans(cluster_s)
  
  #########PHASIC PLOT
  data.phasic = data.frame(bootstrap_mean_p, bins_all)
  colnames(data.phasic) = c("freeze", "bins")
  
  plot.freezing.phasic = ggplot2::ggplot()
  for (i in 1:nrow(cluster_p)){
    data = data.frame(t(cluster_p[i,]), bins_all)
    colnames(data) = c("freeze", "bins")
    
    plot.freezing.phasic = plot.freezing.phasic + geom_line(aes(x=bins, y=freeze), col="gray70", alpha = 0.5,data=data) 
  }
  
  #Add average bootstrap trajectory and add Plot asthetics
  plot.freezing.phasic = plot.freezing.phasic + geom_line(aes(x=bins, y=freeze), col="black", data=data.phasic, linewidth=1, linetype = "longdash") +
    theme_test(base_size=19) +
    labs(x = "Bin", y = paste0("Average freezing per bin ", retrieval,  " [%]")) +
    scale_x_continuous(breaks = c(1,13,18,24,36), labels = c(1,13,18,24,36)) + 
    scale_y_continuous(expand = c(0,1), limits = c(-4,100)) +
    theme(axis.text = element_text(colour="black")) 
  
  data.phasic.original = data.frame(average_p_original, bins_all)
  colnames(data.phasic.original) = c("freeze", "bins")
  
  plot.freezing.phasic = plot.freezing.phasic + geom_line(aes(x=bins, y=freeze), col="green2", data=data.phasic.original, linewidth=1)
  
  #########SUSTAINED PLOT
  data.sus = data.frame(bootstrap_mean_s, bins_all)
  colnames(data.sus) = c("freeze", "bins")
  
  plot.freezing.sus = ggplot2::ggplot()
  for (i in 1:nrow(cluster_s)){
    data = data.frame(t(cluster_s[i,]), bins_all)
    colnames(data) = c("freeze", "bins")
    
    plot.freezing.sus = plot.freezing.sus + geom_line(aes(x=bins, y=freeze), col="gray70", alpha = 0.5,data=data) 
  }
  
  #Add average bootstrap trajectory and add Plot asthetics
  plot.freezing.sus = plot.freezing.sus + geom_line(aes(x=bins, y=freeze), col="black", data=data.sus, linewidth=1, linetype = "longdash") +
    theme_test(base_size=19) +
    labs(x = "Bin", y = paste0("Average freezing per bin ", retrieval,  " [%]")) +
    scale_x_continuous(breaks = c(1,13,18,24,36), labels = c(1,13,18,24,36)) + 
    scale_y_continuous(expand = c(0,1), limits = c(-4,100)) +
    theme(axis.text = element_text(colour="black")) 
  
  data.sus.original = data.frame(average_s_original, bins_all)
  colnames(data.sus.original) = c("freeze", "bins")
  
  plot.freezing.sus = plot.freezing.sus + geom_line(aes(x=bins, y=freeze), col="magenta", data=data.sus.original, linewidth=1)
  
  
  #RESIZE
  plot.freezing.phasic = plot.freezing.phasic + 
    theme(axis.text = element_text(colour="black", size=12)) + theme_test(base_size=14)
  
  plot.freezing.sus= plot.freezing.sus + 
    theme(axis.text = element_text(colour="black", size=12)) + theme_test(base_size=14)
  
  return(list(Jaccard=plot.jac, Phasic.freezing=plot.freezing.phasic, Sustained.freezing=plot.freezing.sus))
}
