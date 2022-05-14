# GSEA --------------------------------------------------------------------
GSEA.EnrichmentScore <- function(gene.list, gene.set, weighted.score.type = 1, correl.vector = NULL) {

  tag.indicator <- sign(match(gene.list, gene.set, nomatch = 0))  # notice that the sign is 0 (no tag) or 1 (tag)
  no.tag.indicator <- 1 - tag.indicator
  N <- length(gene.list)
  Nh <- length(gene.set)
  Nm <- N - Nh
  if (weighted.score.type == 0) {
    correl.vector <- rep(1, N)
  }
  alpha <- weighted.score.type
  correl.vector <- abs(correl.vector^alpha)
  sum.correl.tag <- sum(correl.vector[tag.indicator == 1])
  norm.tag <- 1/sum.correl.tag
  norm.no.tag <- 1/Nm
  RES <- cumsum(tag.indicator * correl.vector * norm.tag - no.tag.indicator * norm.no.tag)
  max.ES <- max(RES)
  min.ES <- min(RES)
  if (max.ES > -min.ES) {
    # ES <- max.ES
    ES <- signif(max.ES, digits = 5)
    arg.ES <- which.max(RES)
  } else {
    # ES <- min.ES
    ES <- signif(min.ES, digits = 5)
    arg.ES <- which.min(RES)
  }
  return(list(ES = ES, arg.ES = arg.ES, RES = RES, indicator = tag.indicator))
}

# GSEA for one variable across all hallway

#' GSEA_var
#'
#' @param emm_result emmreml result data frame
#' @param varname variable name for calculation, var_varname, beta_varname
#' @param hallway the list containing info about 50 hallways
#' @param numcore number of cores to use, default 3
#' @param sim number rounds of simulation, default 10000
#'
#' @return a dataframe containing ES and ES_p
#' @export
#'
#' @examples GSEA_var(emm,"rank_lps")
GSEA_var<-function(emm_result,varname,hallway,numcore=3,sim=10000){
  emma_noper<-emm_result
  varlist<-c(grep("var_",colnames(emma_noper)))
  print(varlist)
  keep<-which(apply(emma_noper[,varlist],1,function(x){length(which(x>0))})==length(varlist))
  hall_enrich<-matrix(NA,nrow=50,ncol=2)
  hall_enrich<-as.data.frame(hall_enrich)
  for(f in 1:50){
    #Retain genes that converged in our model (indicated by genes with reasonable SE(betas))

    # rank effect in lps
    betacv<-grep(paste0("beta_",varname),colnames(emma_noper))
    varcv<-grep(paste0("var_",varname),colnames(emma_noper))
    o<-emma_noper[keep,betacv]/sqrt(emma_noper[keep,varcv])
    names(o)<-rownames(emma_noper)[keep]
    o<-o[order(o)]
    hall_enrich[f,1]<-GSEA.EnrichmentScore(gene.list = names(o),
                                           gene.set = unlist(hallmark[[f]]),
                                           correl.vector = o,weighted.score.type = 1)$ES

    temp<-1:sim
    clus <- makeCluster(numcore)
    registerDoParallel(cores=numcore)
    clusterExport(clus,varlist =c("hallmark","keep",'o',"GSEA.EnrichmentScore","f","temp"),envir=environment())
    temp2<-parSapply(cl = clus,X = 1:sim,FUN = function(i){
      #Randomly shuffle standardized betas across genes
      rand<-sample(1:length(keep),size= length(keep),replace=FALSE)
      o_perm<-o[rand]
      names(o_perm)<-names(o)
      o_perm<-o_perm[order(o_perm)]
      temp[i]<-GSEA.EnrichmentScore(gene.list = names(o_perm), gene.set = unlist(hallmark[[f]]),correl.vector = o_perm,weighted.score.type = 1)$ES
      return(temp[i])
    })
    stopCluster(clus)

    #The pvalue is the % of abs(observed ES)< abs(permuted ES)
    hall_enrich[f,2]<-length(which(abs(hall_enrich[f,1]) < abs(temp2)))/10000
    colnames(hall_enrich)<-c(paste0("ES_",varname),paste0("ES_p_",varname))
    print(paste("Finished for hallway",f,"of",varname))
  }
  row.names(hall_enrich)<-names(hallway)
  return(hall_enrich)
}
