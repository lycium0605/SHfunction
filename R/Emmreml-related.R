# Run Emmereml ------------------------------------------------------------

#' create_design
#'
#' @param data dataframe containing meta info
#' @param nested_name var to be nested within LPS vs. Null
#' @param lps_name name of the column containing "LPS"/"NULL", default to "treatment"
#'
#' @return a design matrix
#' @export
#'
#' @examples create_design(dataframe,"rank","treatment")
create_design<-function(data,nested_name,lps_name="treatment"){
  ## With association structure without permutation
  int<-rep(1,nrow(data))
  lps<-subset(data,select=lps_name) #remember to switch it to residuals
  nested_lps<-as.data.frame(scale(subset(data,select=nested_name)))
  nested_lps[lps=="NULL",]<-0
  colnames(nested_lps)<-paste(colnames(nested_lps),"lps",sep = "_")
  nested_null<-as.data.frame(scale(subset(data,select=nested_name)))
  nested_null[lps=="LPS",]<-0
  colnames(nested_null)<-paste(colnames(nested_null),"null",sep = "_")
  lps_num<-rep(0.5,nrow(data))
  lps_num[lps=="NULL"]<--0.5
  design<-as.data.frame(int)
  design<-cbind(design,lps_num)
  colnames(design)<-c("intercept","LPS")
  for(i in 1:ncol(nested_lps)){
    tmp<-colnames(design)
    design<-cbind(design,nested_null[,i])
    design<-cbind(design,nested_lps[,i])
    colnames(design)<-c(tmp,colnames(nested_null)[i],colnames(nested_lps)[i])
  }
  return(design)
}


#' run_emm
#'
#' @param z completed cases
#' @param resids GE residual matrix
#' @param rel relatedness matrix
#' @param design design matrix
#' @param numcore number of core to use, default to 3
#'
#' @return a data frame containing emmreml results
#' @import parallel
#' @import EMMREML
#' @import doParallel
#' @export

run_emm<-function(z,resids,rel,design,numcore=3){
  # Run EMMREML model
  iclus <- makeCluster(numcore)
  registerDoParallel(cores=numcore)
  design<-as.matrix(design)
  resids_female<-as.matrix(resids)
  rel_female<-as.matrix(rel)
  clusterExport(iclus,varlist=c("resids_female","design",'z',"rel_female"),envir=environment())
  EMMA_RNA_nested=t(parApply(iclus,resids_female[,z],1,function(y){
    requireNamespace("EMMREML",quietly = T)
    emma=EMMREML::emmreml(y=y,X=design[z,],Z=diag(length(z)),K=rel_female[z,z],varbetahat=T,varuhat=T,PEVuhat=T,test=T)
    p=emma$pvalbeta
    varb=emma$varbetahat
    b=emma$betahat
    return(c(b,varb,p[,"none"]))
  }))
  EMMA_RNA_nested<-as.data.frame(EMMA_RNA_nested)
  colnames(EMMA_RNA_nested)<-c(paste("beta",colnames(design),sep = "_"),
                               paste("var",colnames(design),sep = "_"),
                               paste("p",colnames(design),sep = "_"))
  stopCluster(iclus)
  return(EMMA_RNA_nested)
}

#' batch_emm
#' @description run emmreml for multiple variable in one time
#' @param info meta info data frame containing lps and other variable
#' @param variable a vector of variables to be included in emmreml
#' @param resids_female GE residuals data frame, will automatically be transformed to matrix
#' @param rel_female relatedness data frame, same as above
#' @param lps_name name of "LPS"/"NULL" column, default to "treatment"
#' @param numcore default to 3
#' @param sname default to "sname.x"
#' @param plot plot p histogram or not, default to F
#' @param savedir if you want to save the picture, specify a dir
#' @param x_axis_size font size of the x axis title, default to 20
#' @return a combined result datarame
#' @import ggplot2
#' @export
#'
#' @examples batch_emm(meta_info,c("age","rank"),resids,rel)
batch_emm<-function(info,variable,resids,rel,
                    lps_name="treatment",numcore=3,
                    sname="sname.x",plot=F,savedir="NA",x_axis_size=20){
  message(paste("Running model with",paste(variable,collapse = " and ")))
  design<-create_design(info,variable,lps_name)
  info_sname<-subset(info,select = sname)
  z<-which(complete.cases(design) & info_sname!="SAD"& info_sname!="ACI") #the column has to be sname.x
  if(length(z)<=10){
    message(paste("Only",length(z),"samples have all information, skip."))
  }else{
    message(paste(length(z),"out of",nrow(info),"samples have all information, running emmreml"))
    result_full<-run_emm(z,resids,rel,design,numcore = numcore)
    pname<-colnames(result_full)[c(grep("p_",colnames(result_full)))]
    if(plot==T){
      plot_p(result_full,pname,pname,row = length(variable)+1,col = 2,savedir=savedir,x_axis_size=x_axis_size)
    }
    return(result_full)
  }
}

# Calculating fdr given permutation rate ---------------------------------------------------------------------

#' perm.fdr
#' @description Function for calculating permutation-based false discovery rates based on a generalization of the false discovery rate method of Storey and Tibshirani (2003)
#' @param input_df no-permutation, anything
#' @param perm_df permutation p, only contain p of interesting covariate
#' @param Pvals_col_name in input_df, only one
#' @param name just for final output, could be same with Pvals_col_name
#'
#' @return
#' @export

perm.fdr=function(input_df,perm_df,Pvals_col_name,name){ #name?

  pvals_index=which(colnames(input_df)==Pvals_col_name)
  ro<-input_df[order(input_df[,pvals_index]),] #order input df by the p value
  p_obs <- data.frame(pvalue=ro[,pvals_index]) #extract the p value
  p_vector<-matrix(as.matrix(perm_df),ncol=1) #Flatten the p value
  p_vector=data.frame(p_vector[order(p_vector)])

  F<-p_obs[,1]
  F_o<-p_obs[,1]
  pi_hat<-p_obs[,1]

  j=1
  observed<-length(p_obs[,1])
  randoms<-length(p_vector[,1])

  for(i in 1:observed)
  {
    repeat
    {
      if((p_vector[j,1]<p_obs[i,1])&j<randoms){j<-j+1}else{break}
    }
    F[i]=i/observed
    F_o[i]=(j-1)/randoms
    if(F_o[i]<1){pi_hat[i]=(1-F[i])/(1-F_o[i])}else{pi_hat[i]=1}
  }
  tabla <-data.frame(pi_hat,pval=p_obs[,1])

  tabla[1,]=c(1,0)
  last_percentile_average=mean(tabla$pi_hat[as.integer(min((length(tabla[,1])*0.99),(nrow(tabla)-1)):length(tabla[,1]))])
  tabla[nrow(tabla),]=c(last_percentile_average,1)
  constraint_matrix=as.matrix(data.frame(c(0,2),c(0,1),c(1,0)))
  f_hat<-suppressWarnings(cobs(tabla$pval,tabla$pi_hat,constraint="convex",pointwise=constraint_matrix,maxiter=1000,print.warn=FALSE,print.mesg=FALSE))

  f_hat_serie=f_hat$fitted
  pi_o=f_hat_serie[length(f_hat_serie)]
  pi_o=min(pi_o,1)
  pi_o=max(pi_o,0)

  Fdr_ST_perm=pi_o*F_o/F

  for(i in 1:length(p_obs[,1]))
  {
    Fdr_ST_perm[i]=pi_o*F_o[i]/F[i]
    if(i>1)
    {
      for(j in 1:(i-1))
      {
        if(Fdr_ST_perm[i-j]>Fdr_ST_perm[i]){Fdr_ST_perm[i-j]=Fdr_ST_perm[i]}else{break}
      }
    }
    if(Fdr_ST_perm[i]>1)  Fdr_ST_perm[i]=1
  }

  fdrs_df <-data.frame(ro,q_ST_perm=Fdr_ST_perm)
  rownames(fdrs_df)=rownames(ro)
  colnames(fdrs_df)[ncol(fdrs_df)]=paste0("fdr_",name)

  return(fdrs_df)
}


# Marker analysis ---------------------------------------------------------

single_marker<-function(gene_name,emm_merge){
  IL6<-emm_merge[grep(gene_name,emm_merge$gene),]%>%
    select(contains("lps"))%>%
    select(-contains("M"))
  # Get var list
  varlist<-unique(gsub("beta_","",
                       gsub("var_","",
                            gsub("p_","",colnames(IL6)))))
  # Process data
  IL6_effect<-data.frame(matrix(ncol=4,nrow=length(varlist)))
  colnames(IL6_effect)<-c("Beta","Var","P","Variable")
  for(i in 1:length(varlist)){ #LPS everywhere
    tmp<-IL6%>%select(contains(varlist[i],ignore.case = F))
    tmp$Variable<-varlist[i]
    IL6_effect[i,]<-tmp[1,]
  }
   IL6_effect$SD<-sqrt(IL6_effect$Var)
   IL6_effect<-IL6_effect%>%arrange(Beta)
  return(IL6_effect)
  #return(IL6)
}
