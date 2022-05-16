# Emmreml Pipeline ----------------------------------------------

#' emm_pip
#' @description Pipeline for running emmreml
#' @param ge_meta should contain sname,dart_date,treatment,group,sex,rank,age
#' @param sci sci input from Ramboseli::sci
#' @param agi agi input from Ramboseli::agi
#' @param dsi dsi_pop input from Ramboseli::dsi_pop
#' @param dsi_sum dsi_pop_sumamry input from Ramboseli::dsi_pop_summary
#' @param type "wt_avg"(default),"dir_avg", or "grp_only". calculating weighted average by days_present across group, direct average, or keep only the group observed on the darting date
#' @param res full GE residual dataframe
#' @param rel full K matrix
#' @param sex The sex to look at, "F" (default) or "M"
#' @param fix_var The variables to be included in every emmreml model, defaut to c("age","rank")
#' @param change_var The variables to be added with replacement to the basic model. Default to all SR.
#' @param lps_name The colname indicating LPS/NULL in ge_meta, default to treatment
#' @param numcore The number of core to be used for EMMREML, default to 3
#' @param sname The colname indicaiting sname, default to sname, could be sname.x
#' @param outputdir The folder to store output, end with /
#' @param emm_prefix The prefix for emmresult RDS and p histogram png
#' @param x_axis_size The font size of x axis of p histogram, default to 20
#' @return a data frame with merged emmreml result
#' @export
#' @return Nothing, results will be stored in the specified output directory
#'
#' @importFrom stats weighted.mean
#' @import EMMREML
#' @import dplyr

emm_pip<-function(
  ge_meta,sci,agi,dsi,dsi_sum,type="wt_avg",
  res,rel,sex="F",
  fix_var=c("age","rank"),
  change_var=c("eigen_wt",paste0(c("SCI","AGI","DSI","SumBond"),"_F"),paste0(c("SCI","AGI","DSI","SumBond"),"_M")),
  lps_name = "treatment",numcore = 3,sname="sname",
  outputdir,emm_prefix,x_axis_size=20

){
  # Merge ge and sr meta info
  meta<-merge_ge_sr(ge_meta,sci,agi,dsi,dsi_sum,"wt_avg")

  # Extract info for a single sex
  if(sex!="F" & sex!="M"){
    stop("The sex parameter has to be F or M")
  }else{
    message(paste("Extracting sex_specific data for",sex))
    pos_F<-which(ge_meta$sex==sex)
    res_sub<-res[,pos_F]
    K_sub<-rel[pos_F,pos_F]
    merge_meta_sub<-meta[pos_F,]
  }

  # Run emmreml model
  # 1. Run basic emmreml model with fix_var
  emm_result<-batch_emm(info=merge_meta_sub,
                        variable=fix_var,
                        resids=res_sub,rel=K_sub,
                        lps_name = lps_name,numcore = numcore,
                        sname=sname,plot = F)
  # 2. Running models adding change_var
  for(i in 1:length(change_var)){
    tmp<-batch_emm(info=merge_meta_sub,
                   variable=c(fix_var,change_var[i]),
                   resids=res_sub,rel=K_sub,
                   lps_name = lps_name,numcore = numcore,
                   sname=sname,plot = F)
    if(!is.na(tmp)){
      tmp_sr<-tmp%>%select(contains(change_var[i]))
      emm_result<-emm_result%>%bind_cols(tmp_sr)
    }
  }
  saveRDS(emm_result,paste0(outputdir,emm_prefix,"_",Sys.Date(),"_emm_result.RDS"))

  # 3. Plot p histogram for all variable presented in the giant data frame
  message("Plotting p histogram")
  pname<-colnames(emm_result)[c(grep("p_",colnames(emm_result)))]
  plot_p(emm_result,pname,pname,
         row = ceiling(length(pname)/4),col = 4,
         savedir=paste0(outputdir,emm_prefix,"_",Sys.Date(),"_p_hist.png"),
         x_axis_size=x_axis_size)
}
