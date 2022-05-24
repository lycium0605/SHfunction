# Basic visualization -----------------------------------------------------

#' define_region
#'
#' @param row no default
#' @param col no default
#'
#' @return
#' @export

define_region <- function(row, col){
  viewport(layout.pos.row = row, layout.pos.col = col)}

#' plot_gg
#' @description plot multi-panel ggplot
#' @param gg a vector containing all the ggplot object
#' @param row number of row
#' @param col number of column
#'
#' @return
#' @export

plot_gg<-function(gg,row,col){
  if(length(gg)>row*col){
    print("More ggplot than arranged")
  }
  else{
    plot_index<-NULL
    for (i in 1:row){
      for (j in 1:col){
        plot_index<-c(plot_index,c(i,j))
      }
    }
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow = row, ncol = col)))
    for(i in 1:length(gg)){
      print(gg[i], vp = define_region(row = plot_index[i*2-1], col = plot_index[i*2]))
    }
  }
}

#' save_gg_group
#'
#' @param gg a list of gg objects
#' @param row nrow
#' @param col ncol
#' @param savedir outputdir
#' @param unitwidth default 1440
#' @param unitheight default 600
#'
#' @return nothing, save plot
#' @export

save_gg_group<-function(gg,row,col,savedir,
                        unitwidth=1440,unitheight=600){
  png(filename = savedir,
      width = unitwidth*col,height = unitheight*row,unit = "px")
  plot_gg(gg,row,col)
  dev.off()
}
# Emmreml related -----------------------------------------------------------
# plot p distribution
# dat: data frame
# p_name: a vector containing the name for the p value
# x_name: a vector containing the label for x axis
# bw: binwidth, default as 0.005
# plim: p value color upper limit: default 0.05
# row: row number
# col: column number

#' plot_p
#' @description plot and save p histogram for Emmreml
#' @param dat data frame, emmreml result
#' @param p_name a vector containing the name for the p value
#' @param x_name a vector containing the label for x axis
#' @param bw binwidth, default as 0.005
#' @param plim p value color upper limit: default 0.05
#' @param row row number
#' @param col column number
#' @param x_axis_size letter size, default to 50
#' @param savedir output png dir, default to "NA"
#' @param unitwidth width of a single p histogram, default to 1440
#' @param unitheight height of a single p histogram, default to 600
#'
#' @return
#' @import ggplot2
#' @import grid
#' @export

plot_p<-function(dat,p_name,x_name,bw=0.005,plim=0.05,
                 row,col,x_axis_size=50,
                 savedir="NA",unitwidth=1440,unitheight=600){
  if(length(p_name)!=length(x_name)){
    print("Mismatch of p name and xlab name. Quit.")
  }
  else if(length(p_name)>row*col){
    print("More figure than the plot can hold, check dimension. Quit.")
  }
  else{
    if(savedir!="NA"){
      png(filename = savedir,
          width = unitwidth*col,height = unitheight*row,unit = "px")
    }

    plot_index<-NULL
    for (i in 1:row){
      for (j in 1:col){
        plot_index<-c(plot_index,c(i,j))
      }
    }
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow = row, ncol = col))) # Change to 4 for assoc
    for (i in 1:length(p_name)){
      pvals_index=which(colnames(dat)==p_name[i])
      tmp<-ggplot(data = dat)+
        geom_histogram(mapping = aes(x=dat[,pvals_index],fill = ..x..),binwidth = bw)+
        scale_fill_gradientn(colors=c("dark blue"),limits=c(0,plim))+guides(fill="none")+
        labs(x=x_name[i]) +
        theme(axis.title.x=element_text(size = x_axis_size))

      print(tmp, vp = define_region(row = plot_index[i*2-1], col = plot_index[i*2]))
    }

    if(savedir!="NA"){
      dev.off()
    }
  }
}

#' Plot same effect from two emm result
#'
#' @param emm1 emm1
#' @param emm2 emm2
#' @param var_name var_var_name
#' @param suffix var_var_name_suffix
#' @param pcut default 0.05
#' @param ptsize default 0.5
#' @param color_lab default specificity
#' @import dplyr
#' @return
#' @export

stdbeta_plot_twoemm<-function(emm1,emm2,var_name,suffix=c("X","Y"),
                               pcut=0.05,ptsize=0.5,color_lab="specificity"){
  if(length(var_name)==1){
    var_name=c(var_name,var_name)
  }
  emm1<-emm1%>%select(contains(var_name[1]))%>%
    rename_all(function(x){paste0(x,"_",suffix[1])})
  emm2<-emm2%>%select(contains(var_name[2]))%>%
    rename_all(function(x){paste0(x,"_",suffix[2])})
  merge<-bind_cols(emm1,emm2)%>%
    filter_at(vars(starts_with("var_")),all_vars(.>0))
  g<-stdbeta_plot_simple(x=paste0(var_name[1],"_",suffix[1]),
                         y=paste0(var_name[2],"_",suffix[2]),
                         merge,pcut,suffix,ptsize)
  g2<-g+labs(x=paste("Stdbeta of",var_name[1], "effect in", suffix[1]),
             y=paste("Stdbeta of",var_name[2], "effect in", suffix[2]),
             color=color_lab)
  return(g2)
}

#' stdbeta_plot_simple
#' @description std_beta calculaiton and comparison based on emmreml result.
#' @param x parameter 1
#' @param y parameter 2
#' @param data emmreml result, colname in form beta_x, var_x
#' @param p p cut-off to mark red
#' @param color_lab default c("X","Y")
#' @param ptsize default 0.5
#'
#' @return a ggplot object
#' @import ggplot2
#' @importFrom stats cor
#' @export

stdbeta_plot_simple<-function(x,y,data,p=0.1,color_lab=c("X","Y"),ptsize=0.5){
  #for variable X
  rb<-data[,which(colnames(data)==paste("beta",x,sep = "_"))]/
    sqrt(data[,which(colnames(data)==paste("var",x,sep = "_"))])

  #for variable Y
  bb<-data[,which(colnames(data)==paste("beta",y,sep = "_"))]/
    sqrt(data[,which(colnames(data)==paste("var",y,sep = "_"))])

  #print correlation
  print(paste("In",nrow(data),"genes, correlation between", x,"and", y,":",cor(rb,bb,use = "complete.obs")))

  #filter for p
  sigposX<-which(data[,which(colnames(data)==paste("p",x,sep = "_"))]<p)
  sigposY<-which(data[,which(colnames(data)==paste("p",y,sep = "_"))]<p)
  sigpos<-intersect(sigposX,sigposY)
  print(paste("After filtering for p <",p,
              ", in",length(sigpos),"genes, correlation:",
              cor(rb[sigpos],bb[sigpos],use = "complete.obs")))
  print(paste("in non-significant genes, correlation:",
              cor(rb[-sigpos],bb[-sigpos],use = "complete.obs")))

  #ggplot
  type<-rep("Non-sig",length(rb))
  type[sigposX]<-paste0(color_lab[1],"-specific")
  type[sigposY]<-paste0(color_lab[2],"-specific")
  type[sigpos]<-"Both"
  gg<-ggplot()+
    geom_point(aes(x=scale(rb), y=scale(bb),color = type),size = ptsize)+
    #geom_point(aes(x=scale(rb), y=scale(bb)),color = "grey",size = 0.5)+
    #geom_point(aes(x=scale(rb)[sigpos], y=scale(bb)[sigpos]),size = 0.7,color = type)+
    labs(x=x,y=y)
  return(gg)
}

#' stdbeta_plot
#' @description std_beta calculaiton and comparison based on emmreml result.
#' @param x parameter 1
#' @param y parameter 2
#' @param condition usually lps or null
#' @param data emmreml result, colname in form beta_x_condition, var_x_condition
#' @param p p cut-off to mark red
#'
#' @return a ggplot object
#' @import ggplot2
#' @importFrom  stats cor
#' @export


stdbeta_plot<-function(x,y,condition="lps",data,p=0.1){
  #for variable X
  rb<-data[,which(colnames(data)==paste("beta",x,condition,sep = "_"))]/
    sqrt(data[,which(colnames(data)==paste("var",x,condition,sep = "_"))])

  #for variable Y
  bb<-data[,which(colnames(data)==paste("beta",y,condition,sep = "_"))]/
    sqrt(data[,which(colnames(data)==paste("var",y,condition,sep = "_"))])

  #print correlation
  print(paste("In",nrow(data),"genes, correlation between", x,"and", y,"under condition",condition,":",cor(rb,bb,use = "complete.obs")))

  #filter for p
  sigpos<-which(data[,which(colnames(data)==paste("p",x,condition,sep = "_"))]<p &
                  data[,which(colnames(data)==paste("p",y,condition,sep = "_"))]<p  )
  print(paste("After filtering for p <",p,
              ", in",length(sigpos),"genes, correlation:",
              cor(rb[sigpos],bb[sigpos])))

  #ggplot
  gg<-ggplot()+geom_point(aes(x=scale(rb), y=scale(bb)),color = "grey",size = 0.5)+
    geom_point(aes(x=scale(rb)[sigpos], y=scale(bb)[sigpos]),size = 0.7,color = "red")+
    labs(x=paste(x,condition,sep = "_" ),y=paste(y,condition,sep = "_"))
  return(gg)
}

#' log10qq
#'
#' @param data emm result frame
#' @param var_list a vector of varname, p_varname
#' @param colorlab label for color legend, default to variable
#' @param qlength resolution for quantile, default to 1000
#' @param ptsize point size, default to 1
#' @param labsize lab text size
#' @param alpha default to 0.7
#'
#' @return p a ggplot object
#' @export
#' @import dplyr

log10qq<-function(data,var_list,colorlab="variable",qlength=1000,ptsize=1,labsize=5,alpha=0.7){

  p_theo<-runif(1:nrow(data))
  quantiles = seq(0, 1, length.out=qlength)
  quantiles<-quantiles[2:qlength]

  econ<-list()
  for(i in 1:length(var_list)){
    varname<-var_list[i]
    p_index<-which(colnames(data)==paste0("p_",varname))
    tmp<-cbind(quantile(-log10(p_theo),quantiles),quantile(-log10(data[,p_index]), quantiles))
    colnames(tmp)<-c("theo","emp")
    econ[[i]]<-tmp
    names(econ)[i]<-varname
  }

  df<- cbind(cat=rep(names(econ),sapply(econ,nrow)),do.call(rbind,econ))
  df<-as.data.frame(df)%>%
    mutate(cat=as.factor(cat))%>%
    mutate(theo=as.numeric(theo))%>%
    mutate(emp=as.numeric(emp))
  g<-ggplot(df, aes(theo,emp, group=cat,color=cat))+
    geom_abline(slope = 1,intercept = 0,lty=2,size=ptsize)+
    geom_point(size=ptsize,alpha=alpha)+
    labs(x="-log10(p) dnormal",
         y="-log10(p) empirical",
         color=colorlab)+
    theme(axis.title=element_text(size = labsize),
          legend.text = element_text(size = labsize),
          legend.title = element_text(size = labsize))
  return(g)
}

#' Title
#'
#' @param emmlist a list of emm result
#' @param suffix_list a list of suffix lps_null_suffix
#' @param var_name a pattern to grep "lps_null", an element or a list length = emmlist
#' @param colorlab label for color legend, default to variable
#' @param qlength resolution for quantile, default to 1000
#' @param ptsize point size, default to 1
#' @param labsize lab text size
#' @param alpha default to 0.7
#'
#' @return p a ggplot object
#' @export
#' @import dplyr
log10qq_acrossemm<-function(emmlist,suffix_list,var_name,
                  colorlab="variable",
                  qlength=1000,ptsize=1,labsize=5,alpha=0.7){
  suffix<-suffix_list
  if(length(var_name)==1 & length(emmlist)>1){
    var_name<-rep(var_name,length(emmlist))
  }
  if(length(emmlist)==1 & length(var_name)>1){
    emmlist<-rep(emmlist,length(var_name))
  }
  emm<-emmlist[[1]]%>%select(contains(var_name[1]))%>%
    rename_all(function(x){paste0(x,"_",suffix[1])})
  for(i in 2:length(emmlist)){
    tmp<-emmlist[[i]]%>%select(contains(var_name[i]))%>%
      rename_all(function(x){paste0(x,"_",suffix[i])})
    emm<-bind_cols(emm,tmp)
  }
    emm<-emm%>%filter_at(vars(starts_with("var_")),all_vars(.>0))
    varlist<-colnames(emm%>%select(contains("p_"))) #paste0("p_",var_name))))
    varlist<-gsub("p_","",varlist)
    p<-log10qq(emm,varlist,colorlab,qlength,ptsize,labsize,alpha=alpha)
    return(p)


}
# GSEA related ----------------------------------------------------------
#' plot_gsea_comp
#' @description compare the ES score for two variables
#' @param data GSEA result
#' @param var1 in form ES_p_var1, ES_var1
#' @param var2 in form ES_p_var2, ES_var2
#' @param pcut cut-off to label in color
#' @param lab_size default 2.5, decrease if too crowded
#' @param ptsize point size
#'
#' @return a gg object
#' @import ggplot2
#' @importFrom  stats cor
#' @import dplyr
#' @export

plot_gsea_comp<-function(data,var1,var2,pcut=0.05,lab_size=2.5,ptsize=0.8){
  gsea_result<-data
  g<-ggplot()+
    # not significant
    geom_point(data=gsea_result%>%filter((!!as.name(paste0("ES_p_",var1)))>=pcut&(!!as.name(paste0("ES_p_",var2)))>=pcut),
               aes_string(x=paste0("ES_",var1),y=paste0("ES_",var2)),color="grey",size=ptsize) +
    # significant for first variable - green
    geom_point(data=gsea_result%>%filter((!!as.name(paste0("ES_p_",var1)))<pcut&(!!as.name(paste0("ES_p_",var2)))>=pcut),
               aes_string(x=paste0("ES_",var1),y=paste0("ES_",var2)),color="green",size=ptsize) +
    # significant for second variable - blue
    geom_point(data=gsea_result%>%filter((!!as.name(paste0("ES_p_",var1)))>=pcut&(!!as.name(paste0("ES_p_",var2)))<pcut),
               aes_string(x=paste0("ES_",var1),y=paste0("ES_",var2)),color="blue",size=ptsize) +
    geom_point(data=gsea_result%>%filter((!!as.name(paste0("ES_p_",var1)))<pcut&(!!as.name(paste0("ES_p_",var2)))<pcut),
               aes_string(x=paste0("ES_",var1),y=paste0("ES_",var2)),color="red",size=ptsize) +
    geom_text_repel(data=gsea_result%>%filter((!!as.name(paste0("ES_p_",var1)))<pcut|(!!as.name(paste0("ES_p_",var2)))<pcut),
                    aes_string(x=paste0("ES_",var1),y=paste0("ES_",var2),label="hw"),size=lab_size,
                    max.overlaps = 20)
  return(g)
}


#' plot_gsea_comp_twogsea
#'
#' @param gsealist a list, one or two gsea data frame
#' @param varvec a vector, one or two varname
#' @param suffix a vector, suffix for var1 and var2
#' @param pcut default to 0.05
#' @param lab_size label size, default to 2.5
#' @param ptsize point size, default to 0.8
#' @param max_overlap default to 20
#'
#' @return
#' @export
#' @import ggrepel
#' @import ggplot2
#' @import dplyr

plot_gsea_comp_twogsea<-function(gsealist,varvec,suffix=c("X","Y"),pcut=0.05,lab_size=2.5,ptsize=0.8,max_overlap=20){
  # Check input
  if(length(gsealist)==1){
    g1<-g2<-gsealist[[1]]
  }else if(length(gsealist)==2){
    g1<-gsealist[[1]]
    g2<-gsealist[[2]]
  }else(stop("GSEA list can contain only one or two dataframe"))

  if(length(varvec)==1){
    var1<-var2<-varvec
  }else if(length(varvec)==2){
    var1<-varvec[1]
    var2<-varvec[2]
  }else(stop("Var vector can contain only one or two item"))

  # Merge data frame
  gsea_result<-bind_cols(g1%>%select(contains(paste0("ES_",var1)),contains(paste0("ES_p_",var1)))%>%rename_all(function(x){paste0(x,"_",suffix[1])}),
                         g2%>%select(contains(paste0("ES_",var2)),contains(paste0("ES_p_",var2)))%>%rename_all(function(x){paste0(x,"_",suffix[2])}))%>%
    tibble::rownames_to_column()%>%
    rename(hallway=rowname)
  gsea_result$hallway<-gsub(pattern="HALLMARK_",replacement = "",gsea_result$hallway)
  gsea_result$hallway<-gsub(pattern="_",replacement = " ",gsea_result$hallway)

  # weird Jordan formula for Bonferonni correction
  sig1<-which(gsea_result[,3]<(pcut/50))
  print(length(sig1))
  sig2<-which(gsea_result[,5]<(pcut/50))
  print(length(sig2))
  sig_both<-intersect(sig1,sig2)
  sig_either<-union(sig1,sig2)
  print(length(sig_both))
  gsea_result[,c(3,5)]<-(-1)*log10((gsea_result[,c(3,5)]+0.00001)*50)/3
  gsea_result$specific<-"Non-sig"
  gsea_result$specific[sig1]<-paste0(var1,"_",suffix[1],"-specific")
  gsea_result$specific[sig2]<-paste0(var2,"_",suffix[2],"-specific")
  gsea_result$specific[sig_both]<-"Both"
  gsea_result<-gsea_result%>%filter(specific!="Non-sig")

  g<-ggplot(gsea_result,aes_string(x=paste("ES",var1,suffix[1],sep = "_"),y=paste("ES",var2,suffix[2],sep = "_"),color="specific"))+
    geom_hline(yintercept = 0,lty=2,alpha=0.5)+
    geom_vline(xintercept = 0,lty=2,alpha=0.5)+
    geom_point(size=ptsize) +
    geom_text_repel(aes_string(x=paste("ES",var1,suffix[1],sep = "_"),
                               y=paste("ES",var2,suffix[2],sep = "_"),
                               label="hallway"),size=lab_size,max.overlaps=max_overlap)

  return(g)
#   g<-ggplot()+
#     # not significant
#     geom_point(data=gsea_result%>%filter((!!as.name(paste0("ES_p_",var1)))>=pcut&(!!as.name(paste0("ES_p_",var2)))>=pcut),
#                aes_string(x=paste0("ES_",var1),y=paste0("ES_",var2)),color="grey",size=ptsize) +
#     # significant for first variable - green
#     geom_point(data=gsea_result%>%filter((!!as.name(paste0("ES_p_",var1)))<pcut&(!!as.name(paste0("ES_p_",var2)))>=pcut),
#                aes_string(x=paste0("ES_",var1),y=paste0("ES_",var2)),color="green",size=ptsize) +
#     # significant for second variable - blue
#     geom_point(data=gsea_result%>%filter((!!as.name(paste0("ES_p_",var1)))>=pcut&(!!as.name(paste0("ES_p_",var2)))<pcut),
#                aes_string(x=paste0("ES_",var1),y=paste0("ES_",var2)),color="blue",size=ptsize) +
#     geom_point(data=gsea_result%>%filter((!!as.name(paste0("ES_p_",var1)))<pcut&(!!as.name(paste0("ES_p_",var2)))<pcut),
#                aes_string(x=paste0("ES_",var1),y=paste0("ES_",var2)),color="red",size=ptsize) +
#     geom_text_repel(data=gsea_result%>%filter((!!as.name(paste0("ES_p_",var1)))<pcut|(!!as.name(paste0("ES_p_",var2)))<pcut),
#                     aes_string(x=paste0("ES_",var1),y=paste0("ES_",var2),label="hw"),size=lab_size,
#                     max.overlaps = 20)

}

#' Title
#'
#' @param gsea_list a list of gsea results, 1 element or many
#' @param var_vec a vector of var, 1 element of many
#' @param suffix_vec the suffix, has to be multiple, will directly be x-axis-title
#' @param pathway the index of pathway to keep if only a subset
#' @param labsize the size of axis title, default to 1
#' @return a gg object
#' @export
#' @import ggplot2
#' @import dplyr
GSEA_dot_plot<-function(gsea_list,var_vec,suffix_vec,pathway=c(1:50),labsize=1){
  #Make sure the gsea_list and var_vec are of same length
  if(length(gsea_list)==1 & length(var_vec)>1){
    gsea_list<-rep(gsea_list,length(var_vec))
  }else if(length(gsea_list)>1 & length(var_vec)==1){
    var_vec<-rep(var_vec,length(gsea_list))
  }

  #generate first data frame
  es<-gsea_list[[1]]%>%select(contains(paste0("ES_",var_vec[1])),contains(paste0("ES_p_",var_vec[1])))
  es<-es[pathway,]#select specific pathways, if not then default all pathways
  #es$var<-paste(var_vec[1],suffix_vec[1],sep = "_")
  es$var<-suffix_vec[1]
  es$hallway<-gsub("HALLMARK_","",rownames(es))
  colnames(es)<-c("ES","ES_p","var","hallway")
  #add the rest
  if(length(gsea_list)>1){
    for (i in 2:length(gsea_list)){
      tmp<-gsea_list[[i]]%>%select(contains(paste0("ES_",var_vec[i])),contains(paste0("ES_p_",var_vec[i])))
      tmp<-tmp[pathway,]#select specific pathways, if not then default all pathways
      #tmp$var<-paste(var_vec[i],suffix_vec[i],sep = "_")
      tmp$var<-suffix_vec[i]
      tmp$hallway<-gsub("HALLMARK_","",rownames(tmp))
      colnames(tmp)<-c("ES","ES_p","var","hallway")
      es<-bind_rows(es,tmp)
    }
  }
  colnames(es)<-c("ES","ES_p","var","hallway")
  #es$hallway<-gsub("HALLMARK_","",rownames(es))

  #draw dotplot
  g<-ggplot(es)+geom_point(aes(x=var,y=hallway,color=ES,size=-log10(ES_p+0.00001)),alpha=0.7) +
    theme(axis.title=element_text(size = labsize))
  return(g)

}


# Gene-num related --------------------------------------------------------


#' Title
#'
#' @param single_hallway a list element, a single hallway
#' @param hallway_name a string
#' @param emm_list a emmresult data frame list, list(1) or length=var_list
#' @param var_list a vector of variable name beta_var
#' @param suffix_list a vector of suffix
#'
#'
#' @return a gg object
#' @export
#' @import ggplot2


gene_num_comparison<-function(single_hallway,hallway_name,emm_list,var_list,suffix_list){

  if(length(emm_list)==1 & length(var_list)>1){
    emm_list<-rep(emm_list,length(var_list))
  }else if(length(emm_list)>1 & length(var_list)==1){
    var_list<-rep(var_list,length(emm_list))
  }


  #For the first variable
  keep<-rownames(emm_list[[1]])%in%single_hallway
  emm_pathway<-emm_list[[1]][keep,]
  beta<-select(emm_pathway,contains(paste0("beta_",var_list[1])))[,1]
  gene_num<-table(sign(beta))
  re<-c(paste(var_list[1],suffix_list[1],sep = "_"),unname(gene_num[1])*(-1),unname(gene_num[2]))

  #For the rest variable if any
  if(length(var_list)>1){
    for ( i in 2:length(var_list)){
      keep<-rownames(emm_list[[i]])%in%single_hallway
      emm_pathway<-emm_list[[i]][keep,]
      beta<-select(emm_pathway,contains(paste0("beta_",var_list[i])))[,1]
      gene_num<-table(sign(beta))
      re<-rbind(re,c(paste(var_list[i],suffix_list[i],sep = "_"),unname(gene_num[1])*(-1),unname(gene_num[2])))
    }
  }

  re<-as.data.frame(re)
  colnames(re)<-c("var","neg","pos")
  re$neg<-as.numeric(re$neg)
  re$pos<-as.numeric(re$pos)
  #Draw figure
  g<-ggplot(re)+
    theme_bw()+
    xlab("")+ggtitle(hallway_name)+
    geom_point(aes(x=var,y=neg,fill=var),size=5,pch=21,alpha=0.8)+
    geom_point(aes(x=var,y=pos,fill=var),size=5,pch=21,alpha=0.8)+
    geom_segment(aes(x=var,xend=var,y=neg,yend=pos),lty=3)+
    geom_hline(yintercept = 0)+ylab("# of genes")+ylim(c(-170,160))+theme(text=element_text(size=10))
  return(g)

}



