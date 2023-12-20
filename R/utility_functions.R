# Generate Cluster Statistics
# @author jbard
#

PlotACED_Original <- function(df,xaxis="resolution"){
  if (xaxis == "cluster"){
    optimal_cluster = order(df$ACE_SCORE,decreasing = T)[1]
    optimal_cluster = df$Clusters[optimal_cluster]
    plot(df$ACE_Random~df$Clusters,col="red",
         ylim=c(min(c(df$ACE_Random,df$ACE_SCORE),0),max(df$ACE_Random) + (max(df$ACE_Random) * .5)),ylab = "ACE",xlab="Resolution")
    arrows(x0=df$Clusters, y0=df$BGM-df$BGSD,
           x1=df$Clusters, y1=df$BGM+df$BGSD,
           code=3, angle=90, length=0.05,col="red",lty=2)
    points(df$ACE_SCORE~df$Clusters,col="darkgreen")
    points(df$ACE~df$Clusters,col="blue")
    abline(v=optimal_cluster,h=max(df$ACE_SCORE),lty=3,col="orange")
    legend("top",inset=c(-0.2,0),legend=c("Background ACE", "Prediction ACE","ACED Score"),fill = c("red","blue","darkgreen"),xpd=TRUE)
  } else {
    optimal_cluster = order(df$ACE_SCORE,decreasing = T)[1]
    optimal_cluster = df$Resolution[optimal_cluster]
    plot(df$ACE_Random~df$Resolution,col="red",
         ylim=c(min(c(df$ACE_Random,df$ACE_SCORE),0),max(df$ACE_Random) + (max(df$ACE_Random) * .5)),ylab = "ACE",xlab="Resolution")
    arrows(x0=df$Resolution, y0=df$BGM-df$BGSD,
           x1=df$Resolution, y1=df$BGM+df$BGSD,
           code=3, angle=90, length=0.05,col="red",lty=2)
    points(df$ACE_SCORE~df$Resolution,col="darkgreen")
    points(df$ACE~df$Resolution,col="blue")
    abline(v=optimal_cluster,h=max(df$ACE_SCORE),lty=3,col="orange")
    legend("top",inset=c(-0.2,0),legend=c("Background ACE", "Prediction ACE","ACE Score"),fill = c("red","blue","darkgreen"),xpd=TRUE)
  }

}

### https://rdrr.io/github/statlab/permuter/src/R/utils.R
permute_within_rows <- function(x) {
  for (row in seq_len(nrow(x))) {
    x[row, ] <- sample(x[row, ])
  }
  return(x)
}

get_optimal_resolution <- function(aced_object){
  return(aced_object[head(n=1,order(aced_object$ACE_SCORE,decreasing = T)),]$Resolution)
}

get_top_resolutions <- function(aced_object,n=3){
  return(head(aced_object[order(aced_object$ACE_SCORE,decreasing = T),]$Resolution,n=n))
}

get_least_optimal_resolution <- function(aced_object){
  return(aced_object[head(n=1,order(aced_object$ACE_SCORE,decreasing = F)),]$Resolution)
}


PlotClusterBreakpoints <- function(aced_object){
  ggplot(aced_object,aes(x=Resolution,y=Clusters,colour=log(ACE_SCORE,2)))+ geom_point() +
    viridis::scale_colour_viridis(direction = -1) + theme_minimal() +ylab("Clusters") + labs(colour="Log2 Ace")
}


validate_bulk_gedit <- function(bulk_tsv=NULL,refObj=NULL,start=0.1,stop=1,step=0.25){
  values_ACE = c()
  values_res = c()
  values_BACE =c()
  for (res in c(seq(from=start,to=stop,by=step))){
    message("Running GEDIT3 against BULK Truth using the following settings:")
    system(paste0(py_config()$python," ",package_info(include_base = F)$library[1],"/ACED/GEDIT3.py -mix $PWD/",bulk_tsv," -ref $PWD/RefObj.",res,".csv -outFile $PWD/GEDIT_Deconv"),TRUE)
    #message(paste0(py_config()$python," ",package_info("ACED")$path,"/GEDIT3.py -mix $PWD/",bulk_tsv," -ref $PWD/RefObj.",res,".csv -outFile $PWD/GEDIT_Deconv_Bulk.",res))
    #system(paste0(py_config()$python," ",package_info("ACED")$path,"/GEDIT3.py -mix $PWD/",bulk_tsv," -ref $PWD/RefObj.",res,".csv -outFile $PWD/GEDIT_Deconv_Bulk.",res),TRUE)
    predictions = read.table(file=paste0("GEDIT_Deconv_Bulk.",res,"_CTPredictions.tsv"),header = TRUE, row.names = 1, sep = "\t")
    actual = read.table(file=paste0("ACED_ActProp",res,".csv"),header=T,row.names = 1,sep = ",")

    actual <- actual[rownames(actual) %in% rownames(predictions),]
    predictions <- predictions[rownames(actual),]
    print("Pred:")
    print(predictions)
    print("Act:")
    print(actual)

    tab <- table(refObj$orig.ident)
    tab <- as.matrix(tab)
    tab <- tab[rownames(predictions),]
    cells_per_sample <- as.vector(tab)
    print(cells_per_sample)

    actual <- as.matrix(actual)
    predictions <- as.matrix(predictions)

    print("Act:")
    print(actual)
    print("Pred:")
    print(predictions)
    print("Calculating ACE")

    a.mat <- actual ### convert actual porpotion out of table object type
    b.mat <- predictions ## gather estimated cells per cluster up
    a <- as.vector(a.mat * cells_per_sample) # multiply the actual proportion table against the total cells to get cells-per-cluster
    b <- as.vector(b.mat * cells_per_sample) # multiply the estimated proportion table against the total cells to get cells-per-cluster

    print("Bootstrapping the background")
    ACE_Boot <- c()
    for (boot in seq(1,500)){
    samples <- nrow(actual)
    # next get the number of clusters in the object
    clusts <- ncol(actual)
    m <- matrix(rnorm(samples * clusts,mean(actual),sd = 1), nrow=samples)
    #m <- matrix(rnorm(samples * clusts,mean(actual_proportion),sd = sd(actual_proportion)), nrow=samples)
    prob <- exp(m)/rowSums(exp(m))
    back <- as.vector(prob * cells_per_sample) ## background error rate
    ACE_Boot <- rbind(ACE_Boot,sum(abs(a-back)))
    }

    values_ACE = c(values_ACE,sum(abs(a-b)))
    values_BACE = c(values_BACE,mean(ACE_Boot))
    values_res = c(values_res,res)

  }
  return(data.frame("ACE"=values_ACE,"Resolution"=values_res,"ACE_Random"=values_BACE,"ACE_SCORE"=(values_BACE-values_ACE)))
}

### Updating ACE plotting to ggplot style
PlotACED <- function(ace.results,title=""){
  ggplot(ace.results,aes(x=Resolution,y=ACE_SCORE)) + geom_point() + geom_smooth() +
    geom_point(x=get_optimal_resolution(ace.results),
               y=ace.results[ace.results$Resolution == get_optimal_resolution(ace.results),]$ACE_SCORE,
               shape=3,color="red") + theme_minimal(base_size = 12) + ylab("ACE")+
    labs(subtitle = paste0("ACED res :",get_optimal_resolution(ace.results))) + ggtitle(title)
}

