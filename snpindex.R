#一个根据snp个数计算snp-index的软件包
rm(list = ls())
setwd("/Users/jack/r")
#计算snp-index的选项calsnpinde，输入calsnpindex(vcffile)，vcffile是您的vcf文件名
calsnpindex<-function(vcffile){
  library("vcfR")
  library(stringr)
  vcf<-read.vcfR(vcffile)
  gt<-extract.gt(vcf)
  ad<-extract.gt(vcf,"AD")
  ad_dataframe<-as.data.frame(ad)
  pop1<-limma::strsplit2(ad_dataframe$celldeath,split = ",")
  pop2<-limma::strsplit2(ad_dataframe$nocelldeath,split = ",")
  pop1<-pop1[,c(1,2)]
  pop2<-pop2[,c(1,2)]
  all_pop<-cbind(pop1,pop2)
  all_pop_nona<-all_pop[rowSums(is.na(all_pop))==0,]
  snp_index_pop1<-as.numeric(all_pop_nona[,2])/(as.numeric(all_pop_nona[,1])+as.numeric(all_pop_nona[,2]))
  snp_index_pop2<-as.numeric(all_pop_nona[,4])/(as.numeric(all_pop_nona[,3])+as.numeric(all_pop_nona[,4]))
  daltsnp_index=snp_index_pop2 - snp_index_pop1
  contig_site<-limma::strsplit2(row.names(ad),split = "_")
  snp_index<-data.frame(contig_site,snp_index_pop1,snp_index_pop2,daltsnp_index)
  row.names(snp_index)<-row.names(ad)
  colnames(snp_index)=c("contig","pos","snp_index_pop1","snp_index_pop2","daltsnp_index")
  contigid<-as.vector(as.numeric(str_sub(snp_index$contig,6,9)))
  snp_index<-cbind(contigid,snp_index)
  snp_index<-na.omit(snp_index)
  as.numeric(snp_index$pos)
  snp_index$pos<-as.numeric(as.vector(snp_index$pos))
  snp_index<-snp_index[,-2]
  return(snp_index)
}

#计算滑窗calculcvaluentwindow,可以选择步长和窗口,参数1是计算好的snp-index，2是窗口，3是步长
calculcvaluentwindow<-function(snp_index,win,st){
  snp_number=as.numeric(dim(snp_index)[1])
  wid_snp_index <- data.frame(win_ctg="",win_pos="",pop1_wid_snpindex="",pop2_wid_snpindex="")[-1,]
  i=1
  while ((i+(win-1))<=snp_number){
    win_pos<-as.numeric(as.vector(snp_index$pos)[i])
    win_ctg<-as.vector(snp_index$contigid)[i]
    win_snpindex_ve<-c(win_ctg,win_pos,sum(snp_index[i:(i+(win-1)),3])/win,sum(snp_index[i:(i+(win-1)),4]/win))
    wid_snp_index<-rbind(wid_snp_index,win_snpindex_ve)
    i=i+st
  }
  win_daltsnp_index = wid_snp_index[,4] -  wid_snp_index[,3]
  wid_snp_index<-cbind(wid_snp_index,win_daltsnp_index)
  colnames(wid_snp_index)<-c("contigid","win_pos","pop1_wid_snpindex","pop2_wid_snpindex","wid_daltsnpindex")
  return(wid_snp_index)
}



#绘图,第一个参数是snp_index，第二个是win_snpindex
draw_snpindex_picture<-function(snp_index,win_snpindex){
  library(ggplot2)  # 加载绘图包ggplot2
  p1 <- ggplot()+ 
    geom_point(data = snp_index, aes(x =pos, y =daltsnp_index,color = factor(contigid)),size = 1)+   
    geom_line(data = win_snpindex, aes(x = win_pos, y = wid_daltsnpindex), size = 1)+ 
    ylim(0,1)+ 
    facet_wrap( ~ contigid,ncol = 10,strip.position = "bottom",scales = "free_x")+ 
    geom_hline(yintercept=0.55, colour="#CD2626", linetype="dashed", size = 0.8)+
    theme_classic()+
    theme(axis.text.x = element_blank(), strip.background = element_blank(), strip.placement = "outside",axis.line.x=element_line(linetype=1,color="black",size=1),axis.line.y=element_line(linetype=1,color="black",size=1),axis.title.x = element_text(size = 15),axis.title.y = element_text(size = 15), legend.position = 'none')+
    labs(x="Contig ID",y="dalt_SNP_index")
  return(p1)
}

#去除我们不想呈现的contig,remove_contig,第一个参数是去除contig的向量,第二个是文件
remove_contig<-function(contig,snp_index){
  for (i in contig) {
    clean_snpindex<-subset(snp_index,snp_index[,1] !=i)
    snp_index=clean_snpindex
  }
  return(clean_snpindex)
}
#选择我们想呈现的contig,remove_contig,第一个参数是去除contig的向量,第二个是文件
choose_contig<-function(contig,snp_index){
  a <- data.frame(win_ctg="",win_pos="",pop1_wid_snpindex="",pop2_wid_snpindex="",b="")[-1,]
  for (i in contig) {
    clean_snpindex1<-subset(snp_index,snp_index[,1]==i)
    a<-rbind(a,clean_snpindex1)
  }
  colnames(a)<-colnames(snp_index)
  return(a)
}
#选择一个contig
choose_snpindex<-function(snp_index,ctg){
  choosed_snpindex<-subset(snp_index,snp_index[,1]==ctg)
  return(choosed_snpindex)
}
#绘制我们想呈现的一个contig,第一个参数是选择的contig,第二个是文件
draw_onecontig<-function(snp_index,cho_cotig_win,contig){
  cho_cotig<-choose_snpindex(snp_index,contig)
  cho_cotig_win<-choose_snpindex(cho_cotig_win,contig)
  picture_title<-paste("Contig",contig,sep = "_")
  p1 <- ggplot()+
    geom_point(data = cho_cotig,aes(x = pos/1000000, y = daltsnp_index),size = 2,color = "#FF3030")+
    ylim(0,1)+
    geom_line(data =cho_cotig_win , aes(x = win_pos/1000000, y =wid_daltsnpindex), size = 1.5)+
    geom_hline(yintercept=0.55, colour="#CD2626", linetype="dashed", size = 0.8)+
    theme_classic()+
    theme(axis.line.x=element_line(linetype=1,color="black",size=1),axis.line.y=element_line(linetype=1,color="black",size=1),axis.title.x = element_text(size = 14),axis.title.y = element_text(size = 14), legend.position = 'none')+
    labs(x="Position(Mb)",y="SNP_index",title = picture_title)
  return(p1)
}
#step1计算snpindex
snp_index<-calsnpindex("allsample_snp_flt.vcf")
#step2计算滑窗snpindex
win_snpindex<-calculcvaluentwindow(snp_index,8,2)
#step3通过滑窗snpindex和snp-index绘制图片，通过图片选择需要去除的contig
p1<-draw_snpindex_picture(snp_index,win_snpindex)
#step4去除contig,绘制图片
contig<-c(38,61,96,73,104,113,129,198,237,277,283,357,395,352,351,420,430,529,537,522,495,842,472,683,649,1049,1026,1370,508,528,317,396,567,117,147,140,47,49,54,53,50,55,56,972,230,43,1298,585,157,44)
clean_snpindex<-remove_contig(contig,snp_index)
clean_winsnpindex<-remove_contig(contig,win_snpindex)  
p2<-draw_snpindex_picture(clean_snpindex,clean_winsnpindex)
#step5选取多个contig,绘制图片
contig<-c(15,32,59,109)
choosed_contig<-choose_contig(contig,snp_index)
choosed_contig2<-choose_contig(contig,win_snpindex)
p3<-draw_snpindex_picture(choosed_contig,choosed_contig2)
#step6绘制单个contig的图片
p15<-draw_onecontig(snp_index,win_snpindex,15)
p109<-draw_onecontig(snp_index,win_snpindex,109)
p32<-draw_onecontig(snp_index,win_snpindex,32)
p34<-draw_onecontig(snp_index,win_snpindex,34)
p39<-draw_onecontig(snp_index,win_snpindex,39)
p59<-draw_onecontig(snp_index,win_snpindex,59)
p266<-draw_onecontig(snp_index,win_snpindex,266)




  