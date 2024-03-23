#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 0-FUNCTIONS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Functions for testing zero sum dynamics
#
# Author: Haoran Wu (haoran.wu@wolfson.ox.ac.uk)
# Institute: Environmental Change Institute, University of Oxford
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#remove invalid matrices
compadre_filter <- function(compadre){
  #preserve annual matrix
  is_annual_matrix <- !(compadre$metadata$MatrixSeasonal=="Yes")
  is_annual_matrix[is.na(is_annual_matrix)] <- FALSE
  #remove matrices with NA values 
  has_na <- do.call(c, lapply(compadre$mat, function(mat) any(is.na(mat$matA))))
  #update data
  compadre$metadata <- compadre$metadata[is_annual_matrix & !has_na,]
  compadre$matrixClass <- compadre$matrixClass[which(is_annual_matrix & !has_na)]
  compadre$mat <- compadre$mat[which(is_annual_matrix & !has_na)]
  compadre
}

#only preserve specified organism type
compadre_filter_organismType <- function(compadre, organs){
  filter <- compadre$metadata$OrganismType %in% organs
  
  compadre$metadata <- compadre$metadata[filter,]
  compadre$matrixClass <- compadre$matrixClass[filter]
  compadre$mat <- compadre$mat[filter]
  compadre
}

#calculate dominant eigenvalues
#  NOTE: for cases with equal absolute values of negative and positive eigenvalues,
#  I take the positive ones for simplification.
compadre_calc_domin_eigen <- function(compadre){
  dat <- compadre
  as.data.frame(do.call(rbind, lapply(1:nrow(dat$metadata), function(i){
    ii<- i
    #estimate eigenvalue
    eigens_mat <- eigen(dat$mat[[i]]$matA)
    eigens_values <- eigens_mat$values
    
    #special cases... all values are complex numbers
    #  but no cases.
    if(all(Im(eigens_values)!=0)) stop("all eigen values are complex numbers.") 
    
    #special cases... all real eigenvalues are negative
    #  but no cases.
    real_values <- as.numeric(eigens_values[Im(eigens_values)==0])
    if(all(real_values<0)) stop("all real values are negative.") #special cases...
    
    #special cases... there are complex eigenvalues having the largest real part
    #  but no cases...
    if(max(Re(eigens_values))>max(abs(real_values))) stop("complex eigenvalues have larger real part than real eigenvalues.")
    
    #special cases... one positive and one negative eigenvalues with equal absolute values
    #  use the maximum absolute value
    #  this would demonstrate oscillate dynamics with a period T = 2yr but I just
    #    simplify it and take a dummy dominant eigenvalue as an average for two years.
    if(any(real_values<0)){
      if(abs(abs(min(real_values[real_values<0])) - max(real_values[real_values>=0])) < 0.0000001){
        #...just go next
      }
      
      #special cases... a negative eigenvalue having the highest absolute value
      #  but no cases.
      if(abs(min(real_values[real_values<0])) - max(real_values[real_values>=0]) > 0.000001) stop("a negative eigenvalue having the highest absolute value")
    }
    
    #calculate dominate eigenvalue
    domin <- max(real_values)
    if(length(eigens_values)==1){
      damp <- Inf
    } else if(max(real_values) == 0 & sort(abs(Re(eigens_values)))[2] == 0){
      damp <- 0
    } else{
      damp <- max(real_values)/sort(abs(Re(eigens_values)))[2]
    }
    
    c(domin = domin, ratio = damp)
  })))
}

#Figure 1AB. asymptotic behavior
eigen_plot <- function(eigens){
  theme1 <- theme_bw()+
    theme(axis.text.x=element_text(size=16,angle=0,colour="black"),
          axis.text.y=element_text(size=16,angle=0,colour="black"),
          axis.title=element_text(size=18),
          axis.line=element_line(linetype=1,color="black",size=0.1),
          axis.ticks = element_line(colour="black"),
          panel.grid.major = element_blank(), #change the major and minor grid lines, 
          panel.grid.minor = element_blank(), #if want to change, check this parameters, I think it's easier to dao that
          #strip.background = element_rect(colour = "black",size = 0.8),
          #panel.background = element_rect(colour="black", fill="white"),
          #panel.border = element_blank(),
          panel.border = element_rect(colour = "black",fill=NA,size = 1.2),
          plot.title=element_text(size=14,angle=0,colour="black", face = "italic"),
          plot.tag=element_text(size=14,angle=0,colour="black", face = "bold"),
          plot.caption=element_text(size=14,angle=0,colour="black",face = "italic"),
          axis.title.y=element_text(vjust=1.9),
          axis.title.x=element_text(vjust=0.5),
          legend.text=element_text(colour="black",size=14),
          legend.background= element_rect(fill = "transparent", color = NA),
          legend.position = "none",
          legend.title = element_text(colour="black", size=14,angle=0))
  
  eigens$color <- ifelse(log10(eigens$domin)>=0,"blue","red")
  ggplot() + geom_point(data = eigens, aes(x=log10(domin), y=log10(ratio), color=color), alpha = 0.3) + theme1 +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.9) +
    scale_color_manual(values = c(red="#ff0000",blue="#0000ff")) + 
    xlab("Dominant eigenvalue") + ylab("Damping Ratio")
}

#Figure 1AB-Signif.Test. Geometric mean deviates from 1
geometric_mean_test <- function(vect, delta = 0.01){
  vect[vect==0] <- delta
  cat("Geometric Mean =", 10^mean(log10(vect)), "\n")
  cat("Test of geometric mean = 1:\n")
  print(wilcox.test(log10(vect)))
  NULL
}

#Figure 1CD. Species Abundance Distribution
SAC_plot <- function(lambda_vect){
  lambda_vect <- lambda_vect[lambda_vect!=0]
  X <- 1:length(lambda_vect)
  
  theme1 <- theme_bw()+
    theme(axis.text.x=element_text(size=16,angle=0,colour="black"),
          axis.text.y=element_text(size=16,angle=0,colour="black"),
          axis.title=element_text(size=18),
          axis.line=element_line(linetype=1,color="black",size=0.1),
          axis.ticks = element_line(colour="black"),
          panel.grid.major = element_blank(), #change the major and minor grid lines, 
          panel.grid.minor = element_blank(), #if want to change, check this parameters, I think it's easier to dao that
          #strip.background = element_rect(colour = "black",size = 0.8),
          #panel.background = element_rect(colour="black", fill="white"),
          #panel.border = element_blank(),
          panel.border = element_rect(colour = "black",fill=NA,size = 1.2),
          plot.title=element_text(size=14,angle=0,colour="black", face = "italic"),
          plot.tag=element_text(size=14,angle=0,colour="black", face = "bold"),
          plot.caption=element_text(size=14,angle=0,colour="black",face = "italic"),
          axis.title.y=element_text(vjust=1.9),
          axis.title.x=element_text(vjust=0.5),
          legend.text=element_text(colour="black",size=14),
          legend.background= element_rect(fill = "transparent", color = NA),
          legend.title = element_blank())
  
  ggplot() + 
    geom_line(aes(x=X, y=log10(dmzsm(X, J = sum(lambda_vect), theta = 1)), color="θ = 1"), size=1) +
    geom_line(aes(x=X, y=log10(dmzsm(X, J = sum(lambda_vect), theta = 10)), color="θ = 10"),size=1) +
    geom_line(aes(x=X, y=log10(dmzsm(X, J = sum(lambda_vect), theta = 20)), color="θ = 20"),size=1) +
    geom_line(aes(x=X, y=log10(dmzsm(X, J = sum(lambda_vect), theta = 50)), color="θ = 50"), size=1) +
    geom_line(aes(x=X, y=sort(log10(lambda_vect/sum(lambda_vect)), decreasing = TRUE), color="observed"),linetype="dashed", size=1) +
    scale_color_manual(values = c("θ = 1" = "#9D9DF6", "θ = 10" = "#7475F6",
                                  "θ = 20" = "#4747F6", "θ = 50" = "#1E1EF6",
                                  "observed" = "black"),
                       breaks = c("observed","θ = 1","θ = 10","θ = 20","θ = 50")) +
    theme1 +xlab("Species rank") + ylab(expression("log"[10]*"(Abundance)"))
  
}

#Figure 2AB. eigenvalues by growth from
eigen_plot_Type <- function(eigen_vect, organ_type){
  for(organ in unique(organ_type)){
    if(sum(organ_type==organ)<20){
      eigen_vect <- eigen_vect[organ_type!=organ]
      organ_type <- organ_type[organ_type!=organ]
    }
  }
  
  data <- data.frame(
    group = names(tapply(eigen_vect, organ_type, median)),
    median = as.numeric(tapply(eigen_vect, organ_type, median)),
    max = as.numeric(tapply(eigen_vect, organ_type, quantile, 0.75)),
    min = as.numeric(tapply(eigen_vect, organ_type, quantile, 0.25))
  )
  
  theme1 <- theme_bw()+
    theme(axis.text.x=element_text(size=16,angle=30,colour="black",hjust = 1),
          axis.text.y=element_text(size=16,angle=0,colour="black"),
          axis.title=element_text(size=18),
          axis.line=element_line(linetype=1,color="black",size=0.1),
          axis.ticks = element_line(colour="black"),
          panel.grid.major = element_blank(), #change the major and minor grid lines, 
          panel.grid.minor = element_blank(), #if want to change, check this parameters, I think it's easier to dao that
          #strip.background = element_rect(colour = "black",size = 0.8),
          #panel.background = element_rect(colour="black", fill="white"),
          #panel.border = element_blank(),
          panel.border = element_rect(colour = "black",fill=NA,size = 1.2),
          plot.title=element_text(size=14,angle=0,colour="black", face = "italic"),
          plot.tag=element_text(size=14,angle=0,colour="black", face = "bold"),
          plot.caption=element_text(size=14,angle=0,colour="black",face = "italic"),
          axis.title.y=element_text(vjust=1.9),
          axis.title.x=element_text(vjust=0.5),
          legend.text=element_text(colour="black",size=14),
          legend.background= element_rect(fill = "transparent", color = NA),
          legend.title = element_blank())
  
  ggplot(data, aes(x = group, y = median)) +
    geom_point(size=3) +
    geom_hline(yintercept = 1, linetype = "dashed", size = 1, color = "blue") +
    geom_errorbar(aes(ymin = min, ymax = max), width = 0.2, size=1) + theme1 +
    ylab("Dominant eigenvalue") +xlab(NULL)
}

#Figure 2AB-Signif.Test. Geometric mean deviates from 1
geometric_mean_test_Type <- function(eigen_vect, organ_type, delta = 0.01){
  for(organ in unique(organ_type)){
    if(sum(organ_type==organ)<20){
      eigen_vect <- eigen_vect[organ_type!=organ]
      organ_type <- organ_type[organ_type!=organ]
    }
  }
  tapply(eigen_vect, organ_type, function(vect){
    vect[vect==0] <- delta
    wilcox.test(log10(vect))$p.value
  })
}

#Figure 3AB. Species Abundance Distribution by Class
SAC_plotType <- function(eigen_vect, organ_type){
  for(organ in unique(organ_type)){
    if(sum(organ_type==organ)<20){
      eigen_vect <- eigen_vect[organ_type!=organ]
      organ_type <- organ_type[organ_type!=organ]
    }
  }
  
  data <- tapply(eigen_vect, organ_type, function(vect){
    ret <- sort(vect, decreasing = TRUE)
    ret[ret!=0]
  })
  data <- do.call(rbind, lapply(1:length(data), function(ii){
    data.frame(
      x = 1:length(data[[ii]]),
      value = data[[ii]],
      group = names(data)[ii]
    )
  }))
  
  theme1 <- theme_bw()+
    theme(axis.text.x=element_text(size=16,angle=0,colour="black"),
          axis.text.y=element_text(size=16,angle=0,colour="black"),
          axis.title=element_text(size=18),
          axis.line=element_line(linetype=1,color="black",size=0.1),
          axis.ticks = element_line(colour="black"),
          panel.grid.major = element_blank(), #change the major and minor grid lines, 
          panel.grid.minor = element_blank(), #if want to change, check this parameters, I think it's easier to dao that
          #strip.background = element_rect(colour = "black",size = 0.8),
          #panel.background = element_rect(colour="black", fill="white"),
          #panel.border = element_blank(),
          panel.border = element_rect(colour = "black",fill=NA,size = 1.2),
          plot.title=element_text(size=14,angle=0,colour="black", face = "italic"),
          plot.tag=element_text(size=14,angle=0,colour="black", face = "bold"),
          plot.caption=element_text(size=14,angle=0,colour="black",face = "italic"),
          axis.title.y=element_text(vjust=1.9),
          axis.title.x=element_text(vjust=0.5),
          legend.text=element_text(colour="black",size=14),
          legend.background= element_rect(fill = "transparent", color = NA),
          legend.title = element_blank())
  
  ggplot(data, aes(x=x, y=log10(value/sum(value)), group=group, color=group))+geom_path(size=1)+theme1+
    scale_color_manual(values = colorRampPalette(c("lightblue", "darkblue"))(length(unique(data$group)))) +
    xlab("Species rank") + ylab(expression("log"[10]*"(Abundance)"))
}

#numerical simulation
compadre_calc_simulation <- function(dat){
  for(i in 1:length(dat$mat)){
    mat <- dat$mat[[i]]$matA
    init <- runif(nrow(mat), 0, 1)
    init <- init/sum(init)
    current <- t(init)
    total_individual <- 1
    for(year in 1:10){
      current <- current %*% mat  
      total_individual <- c(total_individual, sum(current))
    }
    if(i==1){
      result <- t(total_individual)
    } else{
      result <- rbind(result, t(total_individual))
    }
  }
  
  do.call(rbind, lapply(2:ncol(result), function(ii){
    data.frame(
      x=1:nrow(result),
      y=sort(result[,ii], decreasing = TRUE)/sum(result[,ii]),
      year=ii
    )  
  })) 
}

#FIGURE 4. plot of numerical simulation
simulation_plot <- function(simu){
  theme1 <- theme_bw()+
    theme(axis.text.x=element_text(size=16,angle=0,colour="black"),
          axis.text.y=element_text(size=16,angle=0,colour="black"),
          axis.title=element_text(size=18),
          axis.line=element_line(linetype=1,color="black",size=0.1),
          axis.ticks = element_line(colour="black"),
          panel.grid.major = element_blank(), #change the major and minor grid lines, 
          panel.grid.minor = element_blank(), #if want to change, check this parameters, I think it's easier to dao that
          #strip.background = element_rect(colour = "black",size = 0.8),
          #panel.background = element_rect(colour="black", fill="white"),
          #panel.border = element_blank(),
          panel.border = element_rect(colour = "black",fill=NA,size = 1.2),
          plot.title=element_text(size=14,angle=0,colour="black", face = "italic"),
          plot.tag=element_text(size=14,angle=0,colour="black", face = "bold"),
          plot.caption=element_text(size=14,angle=0,colour="black",face = "italic"),
          axis.title.y=element_text(vjust=1.9),
          axis.title.x=element_text(vjust=0.5),
          legend.text=element_text(colour="black",size=14),
          legend.background= element_rect(fill = "transparent", color = NA),
          legend.title = element_blank())
  
  ggplot(simu, aes(x=x, y=log10(y), alpha=year, group=year))+geom_path(size=1, color="#121292")+
    theme1+xlab("Species Rank")+ylab(expression("log"[10]*"(Abundace)"))
}

#Downsampling and estimating eigenvalue
eigen_downsampling <- function(eigen_vect, delta = 0.01){
  eigen_vect[eigen_vect==0] <- delta
  dsam <- do.call(rbind, lapply(round(seq(50, length(eigen_vect), (length(eigen_vect)-50)/100)), function(size){
    cat("size = ", size, "\n")
    sample_100 <- do.call(rbind,lapply(1:100, function(i){
      sample_result <- sample(eigen_vect,size)
      data.frame(
        lambda=10^mean(log10(sample_result)),
        p=wilcox.test(log10(sample_result))$p.value
      )
    }))
    data.frame(size=size,max=max(sample_100$lambda),
               mean=mean(sample_100$lambda),min=min(sample_100$lambda),
               p=sum(sample_100$p<=0.05))
  }))
  cat("size = ")
  print(dsam$size)
  cat("\n")
  cat("p values = ")
  print(dsam$p)
  cat("\n")
  theme1 <- theme_bw()+
    theme(axis.text.x=element_text(size=16,angle=0,colour="black"),
          axis.text.y=element_text(size=16,angle=0,colour="black"),
          axis.title=element_text(size=18),
          axis.line=element_line(linetype=1,color="black",size=0.1),
          axis.ticks = element_line(colour="black"),
          panel.grid.major = element_blank(), #change the major and minor grid lines, 
          panel.grid.minor = element_blank(), #if want to change, check this parameters, I think it's easier to dao that
          #strip.background = element_rect(colour = "black",size = 0.8),
          #panel.background = element_rect(colour="black", fill="white"),
          #panel.border = element_blank(),
          panel.border = element_rect(colour = "black",fill=NA,size = 1.2),
          plot.title=element_text(size=14,angle=0,colour="black", face = "italic"),
          plot.tag=element_text(size=14,angle=0,colour="black", face = "bold"),
          plot.caption=element_text(size=14,angle=0,colour="black",face = "italic"),
          axis.title.y=element_text(vjust=1.9),
          axis.title.x=element_text(vjust=0.5),
          legend.text=element_text(colour="black",size=14),
          legend.background= element_rect(fill = "transparent", color = NA),
          legend.title = element_blank())
  ggplot(data = dsam, aes(x=size,y=mean)) + 
    geom_line(size=1)+
    geom_ribbon(aes(ymax=max, ymin=min), alpha=0.4) + 
    xlab("Sample size") + ylab("Dominant eigenvalue") +theme1
}
