# Generate Cluster Statistics
# @author jbard
#

PlotDRRSD <- function(df,xaxis="cluster"){
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

get_least_optimal_resolution <- function(drrsd_object){
  return(drrsd_object[head(n=1,order(drrsd_object$ACE_SCORE,decreasing = F)),]$Resolution)
}


PlotClusterBreakpoints <- function(drrsd_object){
  ggplot(drrsd_object,aes(x=Resolution,y=ACE_SCORE,color=as.factor(Clusters))) + geom_point() +
    geom_point(x=get_optimal_resolution(drrsd_object),
               y=drrsd_object[drrsd_object$Resolution == get_optimal_resolution(drrsd_object),]$ACE_SCORE,
               shape=3,color="black")+theme_minimal()
}
