#' gsea_pip
#'
#' @param emm_result The emm result dataframe
#' @param varname A vector of varname, var_varname, beta_varname, default to all
#' @param hallway The list containing all hallway info
#' @param numcore Number of core, default to 3
#' @param sim number of simulation round
#' @param type add a column stating the info, default to NULL
#'
#' @return a data frame containing all ES_p and ES
#' @export
#' @import dplyr
#' @import parallel
#' @import doParallel

gsea_pip<-function(emm_result,varname="all",hallway,numcore=3,sim=10000,type=NULL){
  if(varname=="all"){
    cn<-colnames(emm_result%>%select(contains("var_")))
    varname<-gsub("var_","",cn)
  }
  gsea_result<-GSEA_var(emm_result,varname[1],hallway,numcore,sim)
  for(i in 2:length(varname)){
    tmp<-GSEA_var(emm_result,varname[i],hallway,numcore,sim)
    gsea_result<-gsea_result%>%bind_cols(tmp)
  }
  if(!is.null(type)){
    gsea_result$type<-type
  }
  return(gsea_result)
}
