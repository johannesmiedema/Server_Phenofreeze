freezing_curves = function(df, class, animal, retrieval,  individual = FALSE, stat = FALSE,
                           shape = 16){
  
  bins_all = 1:36
  #Average freezing R1 
  mean_ph = base::colMeans(df[class=="phasic",])
  se_ph = base::apply(df[class=="phasic", ], 2, function(x)sd(x)/sqrt(length(x)))
  
  mean_sus = base::colMeans(df[class=="sustained",])
  se_sus = base::apply(df[class=="sustained", ], 2, function(x)sd(x)/sqrt(length(x)))
  
  if ("shifter" %in% class){
    mean_shift = base::colMeans(df[class=="shifter",])
    se_shift = base::apply(df[class=="shifter",], 2, function(x)sd(x)/sqrt(length(x)))
  }
  
  #-----
  data.plot = data.frame(mean_ph, se_ph, bins_all)
  colnames(data.plot) = c("freeze", "se", "bins")
  
  
  p1 = ggplot2::ggplot(data.plot, aes(bins, freeze))  + theme_test(base_size=14) +
    labs(x = "Bin", y = paste0("Average freezing per bin ", retrieval,  " [%]")) + 
    geom_line(col = "green2", size=1) + geom_point(size=2.5, col="green2", shape=shape) + 
    geom_errorbar(aes(ymin=freeze-se, ymax=freeze+se), width=.2,
                  position=position_dodge(0.05), col = "green2") +
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
    p1 = p1 + geom_line(col="green2", data=data.plot, size=1) + 
      geom_point(size=2.5, col="green2", data=data.plot, shape=shape) + 
      geom_errorbar(data= data.plot, mapping=aes(ymin=freeze-se, ymax=freeze+se), width=.2,
                    position=position_dodge(0.05), col = "green2") 
  }
  
  #Add sustained animals
  data.plot.sus = data.frame(mean_sus, se_sus, bins_all)
  colnames(data.plot.sus) = c("freeze", "se", "bins")
  p1 = p1 + geom_line(col="magenta", data=data.plot.sus, size=1) + 
    geom_point(size=2.5, col="magenta", data=data.plot.sus, shape=shape) + 
    geom_errorbar(data= data.plot.sus, mapping=aes(ymin=freeze-se, ymax=freeze+se), width=.2,
                  position=position_dodge(0.05), col = "magenta") 
  
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