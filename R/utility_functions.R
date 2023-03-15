# Generate Cluster Statistics
# @author jbard
#

PlotACED <- function(df,xaxis="resolution"){
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
    legend("top",inset=c(-0.2,0),legend=c("Background ACE", "Prediction ACE","DRRSD Score"),fill = c("red","blue","darkgreen"),xpd=TRUE)
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

get_optimal_resolution <- function(drrsd_object){
  return(drrsd_object[head(n=1,order(drrsd_object$ACE_SCORE,decreasing = T)),]$Resolution)
}

get_top_resolutions <- function(drrsd_object,n=3){
  return(head(drrsd_object[order(drrsd_object$ACE_SCORE,decreasing = T),]$Resolution,n=n))
}

get_least_optimal_resolution <- function(drrsd_object){
  return(drrsd_object[head(n=1,order(drrsd_object$ACE_SCORE,decreasing = F)),]$Resolution)
}


PlotClusterBreakpoints <- function(drrsd_object){
  ggplot(drrsd_object,aes(x=Resolution,y=Clusters,colour=log(ACE_SCORE,2)))+ geom_point() +
    viridis::scale_colour_viridis(direction = -1) + theme_minimal() +ylab("Clusters") + labs(colour="Log2 Ace")
}


validate_bulk_gedit <- function(bulk_tsv=NULL,start=0.1,stop=1,step=0.25){
  values_mae = c()
  values_res = c()
  for (res in c(seq(from=start,to=stop,by=step))){
    message("Running GEDIT3 against BULK Truth using the following settings:")
    message(paste0(py_config()$python," ",package_info("DRRSD")$path,"/GEDIT3.py -mix $PWD/",bulk_tsv," -ref $PWD/RefObj.",res,".csv -outFile $PWD/GEDIT_Deconv_Bulk.",res))
    system(paste0(py_config()$python," ",package_info("DRRSD")$path,"/GEDIT3.py -mix $PWD/",bulk_tsv," -ref $PWD/RefObj.",res,".csv -outFile $PWD/GEDIT_Deconv.",res),TRUE)
    predictions = read.table(file=paste0("GEDIT_Deconv",res,"_CTPredictions.tsv"),header = TRUE, row.names = 1, sep = "\t")
    actual = read.table(file=paste0("ACED_ActProp",res,".csv"),header = TRUE, row.names = 1, sep = "\t")
    calculate_mean_absolute_error(actual,predictions)
    values_mae = c(values_mae,calculate_mean_absolute_error(actual,predictions))
    values_res = c(values_res,res)
  }
  return(data.frame("MAE"=values_mae,"RES"=values_res))
}
