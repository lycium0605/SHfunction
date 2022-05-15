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
#' @param x_axis_size letter size, default to 12
#' @param savedir output png dir, default to "NA"
#' @param unitwidth width of a single p histogram, default to 1440
#' @param unitheight height of a single p histogram, default to 600
#'
#' @return
#' @import ggplot2
#' @import grid
#' @export

plot_p<-function(dat,p_name,x_name,bw=0.005,plim=0.05,
                 row,col,x_axis_size=12,
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

#' stdbeta_plot_simple
#' @description std_beta calculaiton and comparison based on emmreml result.
#' @param x parameter 1
#' @param y parameter 2
#' @param data emmreml result, colname in form beta_x, var_x
#' @param p p cut-off to mark red
#'
#' @return a ggplot object
#' @import ggplot2
#' @importFrom stats cor
#' @export

stdbeta_plot_simple<-function(x,y,data,p=0.1){
  #for variable X
  rb<-data[,which(colnames(data)==paste("beta",x,sep = "_"))]/
    sqrt(data[,which(colnames(data)==paste("var",x,sep = "_"))])

  #for variable Y
  bb<-data[,which(colnames(data)==paste("beta",y,sep = "_"))]/
    sqrt(data[,which(colnames(data)==paste("var",y,sep = "_"))])

  #print correlation
  print(paste("In",nrow(data),"genes, correlation between", x,"and", y,":",cor(rb,bb,use = "complete.obs")))

  #filter for p
  sigpos<-which(data[,which(colnames(data)==paste("p",x,sep = "_"))]<p &
                  data[,which(colnames(data)==paste("p",y,sep = "_"))]<p  )
  print(paste("After filtering for p <",p,
              ", in",length(sigpos),"genes, correlation:",
              cor(rb[sigpos],bb[sigpos],use = "complete.obs")))
  print(paste("in non-significant genes, correlation:",
              cor(rb[-sigpos],bb[-sigpos],use = "complete.obs")))

  #ggplot
  gg<-ggplot()+geom_point(aes(x=scale(rb), y=scale(bb)),color = "grey",size = 0.5)+
    geom_point(aes(x=scale(rb)[sigpos], y=scale(bb)[sigpos]),size = 0.7,color = "red")+
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


