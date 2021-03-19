rm(list = ls(all = TRUE))
stringsAsFactors = FALSE

DIR <- getwd()
   DIR <- "C:/Users/Takahiro Konuma/Documents/Trans-Phar/Output/Schizo"
   setwd("C:/Users/Takahiro Konuma/Documents/Trans-Phar/Output/Schizo")


args1 = commandArgs(trailingOnly=TRUE)[1]
print(args1)

###make dir
PATH_eachpair_results <- paste(DIR,"/Spearmanresults/spearman_eachpair_results",sep="")
dir.create(PATH_eachpair_results, showWarnings = F, recursive = T)
PATH_eachpair_coplots <- paste(DIR,"/Spearmanresults/spearman_eachpair_coplots",sep="")
dir.create(PATH_eachpair_coplots, showWarnings = F, recursive = T)
PATH_total_results <- paste(DIR,"/Spearmanresults/spearman_totalresults",sep="")
dir.create(PATH_total_results, showWarnings = F, recursive = T)

###TWASdata shaping and loading
DIR_TWASfiles <- paste(getwd(),"/TWASresults/ALLTISSUE",sep="")
setwd(DIR_TWASfiles)
files <- list.files()
filter <- grep("\\chr_all.focus.tsv$", files) #「chr_all.focus.tsv」listup
TWASfiles <- files[filter]
TWASfilename <- gsub(".tsv", "", TWASfiles)
#TWASfilename <- gsub(".chr_all.focus_upperlower", "", TWASfilename)
TWASfile_number <- length(TWASfiles)

#shaping (delete rows including NULL)
for (h in 1:TWASfile_number){
  TWASdf <- read.delim(TWASfiles[h],sep="\t", header=T)
  filter<- grep("NULL",TWASdf$ens_gene_id)
  TWASdf2 <- TWASdf[-filter,]
  #sort by Z-score
  order <- order(TWASdf2$twas_z)
  TWASdf3 <- TWASdf2[order(TWASdf2$twas_z, decreasing=T),]
  #save
  write.table(TWASdf3,paste(TWASfilename[h],"_shaped.tsv",sep=""),quote=FALSE, sep="\t", col.names=T) 
}



###Cmapdata loading
setwd(DIR)
DIR_Cmapfiles <- "../../Cmap_QCeddata"
files <- list.files(path = DIR_Cmapfiles, full.names = F, pattern=".txt")

###Categorynamefile_loadning
Categorynamefile <- read.table("../../listfile/TissueCelltypeCategory.csv",header=T, sep=",")

###Spearman's correlation analysis (per tissue/cell-type category)

for (i in 1:13){
  print(paste("Spearman's correlation analysis for tissue/cell-type category No.",i))
  #celltypes_files_loading
  setwd(DIR)
  DIR_Cmapfiles <- "../../Cmap_QCeddata"
  celltypefiles <- list.files(path = DIR_Cmapfiles, full.names = F, pattern=as.vector(Categorynamefile$Celltype[i]),ignore.case = TRUE)
  celltypefilesname <- gsub("_QCed.txt", "", celltypefiles)

  #GTExtissue_files_loading
  GTExtissuefiles_tmp<- list.files(path = DIR_TWASfiles, full.names = F, pattern=as.vector(Categorynamefile$Tissue[i]),ignore.case = TRUE)
  filter <- grep("shaped", GTExtissuefiles_tmp)  
  GTExtissuefiles <- GTExtissuefiles_tmp[filter]
  GTExtissuefilesname <- gsub(".chr_all.focus_shaped.tsv", "", GTExtissuefiles)

  for (j in 1:length(celltypefiles)) {
      setwd(DIR)
      DIR_Cmapfiles <- "../../Cmap_QCeddata"
      setwd(DIR_Cmapfiles)
      Celldf <- read.delim(celltypefiles[j],sep="\t", header=T,fileEncoding="cp932")
       sortlist <- order(Celldf$Gene)
       Celldf <- Celldf[sortlist,]

      for (k in 1:length(GTExtissuefiles)) {
        setwd(DIR_TWASfiles)
        TWASdf <- read.delim(GTExtissuefiles[k],sep="\t", header=T)
         TWASdf <- data.frame(Gene =TWASdf$mol_name, TWASZ = TWASdf$twas_z,pip=TWASdf$pip)

         #shaping
         TWASdf_2 <- na.omit(TWASdf)
         TWASdf_2 <- TWASdf_2[TWASdf_2$Gene %in% Celldf$Gene, ]
         TWASdf_2 <- TWASdf_2[order(abs(TWASdf_2$TWASZ),decreasing=T),]

         #ABS(TWAS_Z) TOP10%_get
         TWASdf_3 <- TWASdf_2[1:round((nrow(TWASdf_2)/10)),]
         TWASdf_3 <- aggregate(TWASZ~Gene,data=TWASdf_3,FUN=mean)

         #filtering genes which are ranked both in +(Z score) TOP and -(Z score) TOP
         filter <- abs(TWASdf_3$TWASZ) >= TWASdf_2$TWASZ[round((nrow(TWASdf_2)/10))]
         TWASdf_4<- TWASdf_3[filter,]

         #grep same genes in TWASdf_4 to get Celldf_4
          sortlist <- order(TWASdf_4$Gene)
          TWASdf_4 <- TWASdf_4[sortlist,]
         Celldf_4 <- Celldf[Celldf$Gene %in% TWASdf_4$Gene, ]
          sortlist <- order(Celldf_4$Gene)
          Celldf_4 <- Celldf_4[sortlist,]

        #Spearman correlation analysis (for detecting negative correlation) of TWASdf and Celldf
         createEmptyDf = function( nrow, ncol, colnames = c() ){
          data.frame( matrix( vector(), nrow, ncol, dimnames = list( c(), colnames ) ) )
         }

        binddat <- createEmptyDf( ncol(Celldf)-1, 3, colnames = c("Rho","Pval","Ngenes") )
        for (l in 2:ncol(Celldf)) {
            cortest <- cor.test(TWASdf_4$TWASZ, Celldf_4[,l], method="s",alternative = c("less"))
            Rho <- cortest$estimate[[1]]
            Pval <- cortest$p.value
            Ngenes <- nrow(TWASdf_4)
            binddat[l-1,1] <- Rho
            binddat[l-1,2] <- Pval
            binddat[l-1,3] <- Ngenes

            rownames(binddat)[l-1] <- paste(GTExtissuefilesname[k],colnames(Celldf)[l],sep="_")
            #Co-plot if Pval < 0.0001 
             if (Pval < 0.0001){
             setwd(PATH_eachpair_coplots)
             pdf(paste(GTExtissuefilesname[k],"_",colnames(Celldf)[l],"_Rho_",Rho,"_Pval_",Pval,".pdf",sep=""))
             plot(TWASdf_4$TWASZ, Celldf_4[,l], xlim=c(min(TWASdf_4$TWASZ)-1, max(TWASdf_4$TWASZ)+1),ylim=c(min(Celldf_4[,l])-1, max(Celldf_4[,l])+1),pch=19,
                  main=paste(GTExtissuefilesname[k],"_",colnames(Celldf)[l],sep=""),xlab="TWAS-Z", ylab="CMAP-Z")
             dev.off()
              #Also save dataframe (Z-scores in TWASdf and Celldf)
              df <- cbind(TWASdf_4, Celldf_4[,l])
              colnames(df) <- c("Gene", "TWAS_Z", "CMAP_Z")
              write.table(df, paste(GTExtissuefilesname[k],"_",colnames(Celldf)[l],"_Rho_",Rho,"_Pval_",Pval,".csv",sep=""),sep=",", quote=F,row.names=F,fileEncoding="cp932")
              }
         }

        setwd(PATH_eachpair_results)
        colnames(binddat) <-c("Rho","Pval","Ngenes")
        write.table(binddat,paste(celltypefilesname[j],GTExtissuefilesname[k],"spearmanresults.txt",sep="_"),quote=FALSE, sep="\t", col.names=T,fileEncoding="cp932")   
      }
    }
}



#bind ALL pair data
setwd(PATH_eachpair_results)

resultfiles <- list.files(path = PATH_eachpair_results, full.names = F, pattern="spearmanresults.txt")
resultfilenumber <- length(resultfiles)
result_binddat <- data.frame(Rho=character(),
                             Pval=character(),
                             Ngenes=character(),
                             stringsAsFactors=FALSE)

for (m in 1:resultfilenumber){
      resultfile<- read.delim(resultfiles[m], header=T, sep="\t",fileEncoding="cp932")
      colnames(resultfile)= c("Rho","Pval","Ngenes")
      result_binddat <- rbind(result_binddat, resultfile)
}

#calculate FDR-q
result_binddat$FDR_q <- p.adjust(result_binddat$Pval, method = "BH")

#save all-pair data
setwd(PATH_total_results)
write.table(result_binddat,file="Allpairs_spearmanresults.txt",quote=FALSE, sep="\t", col.names=T,fileEncoding="cp932")


###QQ-plot for p-values of ALL pair data
dotforplot_tmp=c(0:1000,seq(1001,10000,20),seq(10001,20000,500),seq(20001,100000,1000),seq(100001,nrow(result_binddat),10000))
dotforplot=nrow(result_binddat)-dotforplot_tmp

pdf("QQplot_allpairs.pdf",width=7, height=7)
par(pty="s")
par(mar = c(5,4,4,2)+0.1)
par(oma = c(0,0,0,0))

pvals_tmp <- as.data.frame(result_binddat$Pval)
pvals <- na.omit(pvals_tmp)
colnames(pvals) <- "Pval"
distribution_tmp=seq(0,1,length=nrow(pvals)+2)
distribution=distribution_tmp[-(length(distribution_tmp))]
distribution=distribution[-1]

plot(sort(-log10(distribution))[dotforplot],sort(-log10(pvals$Pval))[dotforplot],
     xlab=expression(Expected -log[10](P-value)), ylab=expression(Obserbed -log[10](P-value)),
     xlim=(c(0,7)),ylim=(c(0,7)),
     cex.main=0.9,
     col="red",pch=19,cex=0.7,
     xaxt="n",yaxt="n")
axis(side=1, at=c(0,2,4,6))
axis(side=2, at=c(0,2,4,6))

#y=x line
par(new=T)
abline(0,1,lty=3,lwd=3)

###QQ-plot for p-values of approved drug pair data
 setwd("../../../../listfile")
approveddruglist <- read.delim("ICD_approveddruglist.csv",sep=",",header=T,na.strings = c(""))
traitname <- args1
filter <- grep(traitname,colnames(approveddruglist),ignore.case = T,value = F)
if(length(filter)>0){
  drugs <- approveddruglist[,match(traitname,colnames(approveddruglist))]
  drugs <- as.data.frame(drugs)
  drugs <- na.omit(drugs)
  drugsfilter<-grep(paste(drugs$drugs,collapse="|"),rownames(result_binddat),ignore.case = T,value=FALSE)
  
  data_p_approvedsignal <- result_binddat[drugsfilter,]
  #save all approved signal data
  rownames(data_p_approvedsignal) <- rownames(result_binddat)[drugsfilter]
  setwd(PATH_total_results)
  write.table(data_p_approvedsignal,"Approveddrugpairs_spearmanresults.txt",quote = F,sep=",")
  #QQplot for approved drug
   
  pvals_p_tmp <- as.data.frame(data_p_approvedsignal$Pval)
  pvals_p <- na.omit(pvals_p_tmp)
  distribution_tmp=seq(0,1,length=nrow(pvals_p)+2)
  distribution=distribution_tmp[-(length(distribution_tmp))]
  distribution=distribution[-1]
  par(new=T)
  plot(sort(-log10(distribution)),sort(-log10(pvals_p$`data_p_approvedsignal$Pval`)),
       xlab=expression(Expected -log[10](P-value)), ylab=expression(Obserbed -log[10](P-value)),
       xlim=(c(0,7)),ylim=(c(0,7)),
       cex.main=0.9,
       col="steelblue",pch=19,cex=0.7,
       xaxt="n",yaxt="n")
  axis(side=1, at=c(0,2,4,6))
  axis(side=2, at=c(0,2,4,6))

  #
  cols=c("red","Steelblue")
  labels=c("All drug-cell pair","Approved drug-cell pair")
  pchs=c(19,19)
  legend("bottomright", legend=labels, col=cols,pch=pchs,cex=0.75)
  
  dev.off()
} else {
  print("no approved drugs from TTD and ChEMBL are registered")
  #
  cols=c("red")
  labels=c("All drug-cell pair")
  pchs=c(19)
  legend("bottomright", legend=labels, col=cols,pch=pchs,cex=0.75)
  
  dev.off()
}



###QQ-plot with tissue category

#install.packages("dichromat")
library(dichromat)

col_vector_tmp <- c(colorschemes$Categorical.12, 
                dichromat(colorschemes$Categorical.12)[6],dichromat(colorschemes$Categorical.12)[7],
                dichromat(colorschemes$Categorical.12)[11],dichromat(colorschemes$Categorical.12)[12])
col_vector <- c(col_vector_tmp[2], "darkgoldenrod4", "seagreen4","azure4", "firebrick", 
                "chartreuse2", "cyan2", col_vector_tmp[8], "mediumorchid", col_vector_tmp[10], 
                col_vector_tmp[11], "magenta", "yellow3")



pdf("QQplot_allpairs_withtissuecategory.pdf",width=7, height=7)


pvals_tmp <- as.data.frame(result_binddat$Pval)
pvals <- na.omit(pvals_tmp)
colnames(pvals) <- "Pval"
distribution_tmp=seq(0,1,length=nrow(pvals)+2)
distribution=distribution_tmp[-(length(distribution_tmp))]
distribution=distribution[-1]


#Adipose
filter= grep("ASC",rownames(result_binddat))
pvals_p_tmp <- as.data.frame(result_binddat$Pval[filter])
pvals_p <- na.omit(pvals_p_tmp)
dotforplot2_tmp <- c(0:round(nrow(pvals_p)/5),seq(round(nrow(pvals_p)/5),round(nrow(pvals_p)/4),5),seq(round(nrow(pvals_p)/4),round(nrow(pvals_p)/3),10),seq(round(nrow(pvals_p)/3),round(nrow(pvals_p)/2),20),seq(round(nrow(pvals_p)/2),nrow(pvals_p),40))
dotforplot2=nrow(pvals_p)-dotforplot2_tmp
distribution_tmp=seq(0,1,length=nrow(pvals_p)+2)
distribution=distribution_tmp[-(length(distribution_tmp))]
distribution=distribution[-1]

plot(sort(-log10(distribution))[dotforplot2],sort(-log10(pvals_p$`result_binddat$Pval[filter]`))[dotforplot2],
     xlab="", ylab="",
     xlim=(c(0,7)),ylim=(c(0,7)),
     cex.main=0.9,
     col=col_vector[1],pch=19,cex=0.7,
     xaxt="n",yaxt="n")
axis(side=1, at=c(0,2,4,6))
axis(side=2, at=c(0,2,4,6))


#Breast
filter= grep("HME1|MCF10A|BT20|HS578T|MCF7|MDAMB231|SKBR3",rownames(result_binddat))
pvals_p_tmp <- as.data.frame(result_binddat$Pval[filter])
pvals_p <- na.omit(pvals_p_tmp)
dotforplot2_tmp <- c(0:round(nrow(pvals_p)/5),seq(round(nrow(pvals_p)/5),round(nrow(pvals_p)/4),5),seq(round(nrow(pvals_p)/4),round(nrow(pvals_p)/3),10),seq(round(nrow(pvals_p)/3),round(nrow(pvals_p)/2),20),seq(round(nrow(pvals_p)/2),nrow(pvals_p),40))
dotforplot2=nrow(pvals_p)-dotforplot2_tmp
distribution_tmp=seq(0,1,length=nrow(pvals_p)+2)
distribution=distribution_tmp[-(length(distribution_tmp))]
distribution=distribution[-1]
par(new=T)
plot(sort(-log10(distribution))[dotforplot2],sort(-log10(pvals_p$`result_binddat$Pval[filter]`))[dotforplot2],
     xlab="", ylab="",
     xlim=(c(0,7)),ylim=(c(0,7)),
     cex.main=0.9,
     col=col_vector[2],pch=19,cex=0.7,
     xaxt="n",yaxt="n")
axis(side=1, at=c(0,2,4,6))
axis(side=2, at=c(0,2,4,6))


#CNS
filter= grep("NEU|NPC",rownames(result_binddat))
pvals_p_tmp <- as.data.frame(result_binddat$Pval[filter])
pvals_p <- na.omit(pvals_p_tmp)
dotforplot2_tmp <- c(0:round(nrow(pvals_p)/5),seq(round(nrow(pvals_p)/5),round(nrow(pvals_p)/4),5),seq(round(nrow(pvals_p)/4),round(nrow(pvals_p)/3),10),seq(round(nrow(pvals_p)/3),round(nrow(pvals_p)/2),20),seq(round(nrow(pvals_p)/2),nrow(pvals_p),40))
dotforplot2=nrow(pvals_p)-dotforplot2_tmp
distribution_tmp=seq(0,1,length=nrow(pvals_p)+2)
distribution=distribution_tmp[-(length(distribution_tmp))]
distribution=distribution[-1]
par(new=T)
plot(sort(-log10(distribution))[dotforplot2],sort(-log10(pvals_p$`result_binddat$Pval[filter]`))[dotforplot2],
     xlab="", ylab="",
     xlim=(c(0,7)),ylim=(c(0,7)),
     cex.main=0.9,
     col=col_vector[3],pch=19,cex=0.7,
     xaxt="n",yaxt="n")
axis(side=1, at=c(0,2,4,6))
axis(side=2, at=c(0,2,4,6))


#Digestive
filter= grep("CL34|HCT116|HT115|HT29|LOVO|MDST8|NCIH508|HCIH716|RKO|SNU|SW480|SW620|SW948",rownames(result_binddat))
pvals_p_tmp <- as.data.frame(result_binddat$Pval[filter])
pvals_p <- na.omit(pvals_p_tmp)
dotforplot2_tmp <- c(0:round(nrow(pvals_p)/5),seq(round(nrow(pvals_p)/5),round(nrow(pvals_p)/4),5),seq(round(nrow(pvals_p)/4),round(nrow(pvals_p)/3),10),seq(round(nrow(pvals_p)/3),round(nrow(pvals_p)/2),20),seq(round(nrow(pvals_p)/2),nrow(pvals_p),40))
dotforplot2=nrow(pvals_p)-dotforplot2_tmp
distribution_tmp=seq(0,1,length=nrow(pvals_p)+2)
distribution=distribution_tmp[-(length(distribution_tmp))]
distribution=distribution[-1]
par(new=T)
plot(sort(-log10(distribution))[dotforplot2],sort(-log10(pvals_p$`result_binddat$Pval[filter]`))[dotforplot2],
     xlab="", ylab="",
     xlim=(c(0,7)),ylim=(c(0,7)),
     cex.main=0.9,
     col=col_vector[4],pch=19,cex=0.7,
     xaxt="n",yaxt="n")
axis(side=1, at=c(0,2,4,6))
axis(side=2, at=c(0,2,4,6))




#Haematopoietic
filter= grep("U266|HL60|JURKAT|NOMO1|PL21|SKM1|THP1|U937|WSUDLCL2",rownames(result_binddat))
pvals_p_tmp <- as.data.frame(result_binddat$Pval[filter])
pvals_p <- na.omit(pvals_p_tmp)
dotforplot2_tmp <- c(0:round(nrow(pvals_p)/5),seq(round(nrow(pvals_p)/5),round(nrow(pvals_p)/4),5),seq(round(nrow(pvals_p)/4),round(nrow(pvals_p)/3),10),seq(round(nrow(pvals_p)/3),round(nrow(pvals_p)/2),20),seq(round(nrow(pvals_p)/2),nrow(pvals_p),40))
dotforplot2=nrow(pvals_p)-dotforplot2_tmp
distribution_tmp=seq(0,1,length=nrow(pvals_p)+2)
distribution=distribution_tmp[-(length(distribution_tmp))]
distribution=distribution[-1]
par(new=T)
plot(sort(-log10(distribution))[dotforplot2],sort(-log10(pvals_p$`result_binddat$Pval[filter]`))[dotforplot2],
     xlab="", ylab="",
     xlim=(c(0,7)),ylim=(c(0,7)),
     cex.main=0.9,
     col=col_vector[5],pch=19,cex=0.7,
     xaxt="n",yaxt="n")
axis(side=1, at=c(0,2,4,6))
axis(side=2, at=c(0,2,4,6))



#Liver
filter= grep("HEPG2|HUH7|PHH",rownames(result_binddat))
pvals_p_tmp <- as.data.frame(result_binddat$Pval[filter])
pvals_p <- na.omit(pvals_p_tmp)
dotforplot2_tmp <- c(0:round(nrow(pvals_p)/5),seq(round(nrow(pvals_p)/5),round(nrow(pvals_p)/4),5),seq(round(nrow(pvals_p)/4),round(nrow(pvals_p)/3),10),seq(round(nrow(pvals_p)/3),round(nrow(pvals_p)/2),20),seq(round(nrow(pvals_p)/2),nrow(pvals_p),40))
dotforplot2=nrow(pvals_p)-dotforplot2_tmp
distribution_tmp=seq(0,1,length=nrow(pvals_p)+2)
distribution=distribution_tmp[-(length(distribution_tmp))]
distribution=distribution[-1]
par(new=T)
plot(sort(-log10(distribution))[dotforplot2],sort(-log10(pvals_p$`result_binddat$Pval[filter]`))[dotforplot2],
     xlab="", ylab="",
     xlim=(c(0,7)),ylim=(c(0,7)),
     cex.main=0.9,
     col=col_vector[6],pch=19,cex=0.7,
     xaxt="n",yaxt="n")
axis(side=1, at=c(0,2,4,6))
axis(side=2, at=c(0,2,4,6))


#Lung
filter= grep("A549|CORL23|DV90|H1299|HCC|NCIH1694|NCIH1836|NCIH2073|NCIH596|SKLU1|T3M10",rownames(result_binddat))
pvals_p_tmp <- as.data.frame(result_binddat$Pval[filter])
pvals_p <- na.omit(pvals_p_tmp)
dotforplot2_tmp <- c(0:round(nrow(pvals_p)/5),seq(round(nrow(pvals_p)/5),round(nrow(pvals_p)/4),5),seq(round(nrow(pvals_p)/4),round(nrow(pvals_p)/3),10),seq(round(nrow(pvals_p)/3),round(nrow(pvals_p)/2),20),seq(round(nrow(pvals_p)/2),nrow(pvals_p),40))
dotforplot2=nrow(pvals_p)-dotforplot2_tmp
distribution_tmp=seq(0,1,length=nrow(pvals_p)+2)
distribution=distribution_tmp[-(length(distribution_tmp))]
distribution=distribution[-1]
par(new=T)
plot(sort(-log10(distribution))[dotforplot2],sort(-log10(pvals_p$`result_binddat$Pval[filter]`))[dotforplot2],
     xlab="", ylab="",
     xlim=(c(0,7)),ylim=(c(0,7)),
     cex.main=0.9,
     col=col_vector[7],pch=19,cex=0.7,
     xaxt="n",yaxt="n")
axis(side=1, at=c(0,2,4,6))
axis(side=2, at=c(0,2,4,6))




#Muscloskeletal
filter= grep("SKB_|SKL_|SKL.C",rownames(result_binddat))
pvals_p_tmp <- as.data.frame(result_binddat$Pval[filter])
pvals_p <- na.omit(pvals_p_tmp)
dotforplot2_tmp <- c(0:round(nrow(pvals_p)/5),seq(round(nrow(pvals_p)/5),round(nrow(pvals_p)/4),5),seq(round(nrow(pvals_p)/4),round(nrow(pvals_p)/3),10),seq(round(nrow(pvals_p)/3),round(nrow(pvals_p)/2),20),seq(round(nrow(pvals_p)/2),nrow(pvals_p),40))
dotforplot2=nrow(pvals_p)-dotforplot2_tmp
distribution_tmp=seq(0,1,length=nrow(pvals_p)+2)
distribution=distribution_tmp[-(length(distribution_tmp))]
distribution=distribution[-1]
par(new=T)
plot(sort(-log10(distribution))[dotforplot2],sort(-log10(pvals_p$`result_binddat$Pval[filter]`))[dotforplot2],
     xlab="", ylab="",
     xlim=(c(0,7)),ylim=(c(0,7)),
     cex.main=0.9,
     col=col_vector[8],pch=19,cex=0.7,
     xaxt="n",yaxt="n")
axis(side=1, at=c(0,2,4,6))
axis(side=2, at=c(0,2,4,6))



#Ovary
filter= grep("COV644|EFO27|OV7|RMGI|RMUGS|TYKNU",rownames(result_binddat))
pvals_p_tmp <- as.data.frame(result_binddat$Pval[filter])
pvals_p <- na.omit(pvals_p_tmp)
dotforplot2_tmp <- c(0:round(nrow(pvals_p)/5),seq(round(nrow(pvals_p)/5),round(nrow(pvals_p)/4),5),seq(round(nrow(pvals_p)/4),round(nrow(pvals_p)/3),10),seq(round(nrow(pvals_p)/3),round(nrow(pvals_p)/2),20),seq(round(nrow(pvals_p)/2),nrow(pvals_p),40))
dotforplot2=nrow(pvals_p)-dotforplot2_tmp
distribution_tmp=seq(0,1,length=nrow(pvals_p)+2)
distribution=distribution_tmp[-(length(distribution_tmp))]
distribution=distribution[-1]
par(new=T)
plot(sort(-log10(distribution))[dotforplot2],sort(-log10(pvals_p$`result_binddat$Pval[filter]`))[dotforplot2],
     xlab="", ylab="",
     xlim=(c(0,7)),ylim=(c(0,7)),
     cex.main=0.9,
     col=col_vector[9],pch=19,cex=0.7,
     xaxt="n",yaxt="n")
axis(side=1, at=c(0,2,4,6))
axis(side=2, at=c(0,2,4,6))


#Pancreas
filter= grep("YAPC",rownames(result_binddat))
pvals_p_tmp <- as.data.frame(result_binddat$Pval[filter])
pvals_p <- na.omit(pvals_p_tmp)
dotforplot2_tmp <- c(0:round(nrow(pvals_p)/5),seq(round(nrow(pvals_p)/5),round(nrow(pvals_p)/4),5),seq(round(nrow(pvals_p)/4),round(nrow(pvals_p)/3),10),seq(round(nrow(pvals_p)/3),round(nrow(pvals_p)/2),20),seq(round(nrow(pvals_p)/2),nrow(pvals_p),40))
dotforplot2=nrow(pvals_p)-dotforplot2_tmp
distribution_tmp=seq(0,1,length=nrow(pvals_p)+2)
distribution=distribution_tmp[-(length(distribution_tmp))]
distribution=distribution[-1]
par(new=T)
plot(sort(-log10(distribution))[dotforplot2],sort(-log10(pvals_p$`result_binddat$Pval[filter]`))[dotforplot2],
     xlab="", ylab="",
     xlim=(c(0,7)),ylim=(c(0,7)),
     cex.main=0.9,
     col=col_vector[10],pch=19,cex=0.7,
     xaxt="n",yaxt="n")
axis(side=1, at=c(0,2,4,6))
axis(side=2, at=c(0,2,4,6))


#Prostate
filter= grep("PC3|VCAP|LNCAP",rownames(result_binddat))
pvals_p_tmp <- as.data.frame(result_binddat$Pval[filter])
pvals_p <- na.omit(pvals_p_tmp)
dotforplot2_tmp <- c(0:round(nrow(pvals_p)/5),seq(round(nrow(pvals_p)/5),round(nrow(pvals_p)/4),5),seq(round(nrow(pvals_p)/4),round(nrow(pvals_p)/3),10),seq(round(nrow(pvals_p)/3),round(nrow(pvals_p)/2),20),seq(round(nrow(pvals_p)/2),nrow(pvals_p),40))
dotforplot2=nrow(pvals_p)-dotforplot2_tmp
distribution_tmp=seq(0,1,length=nrow(pvals_p)+2)
distribution=distribution_tmp[-(length(distribution_tmp))]
distribution=distribution[-1]
par(new=T)
plot(sort(-log10(distribution))[dotforplot2],sort(-log10(pvals_p$`result_binddat$Pval[filter]`))[dotforplot2],
     xlab="", ylab="",
     xlim=(c(0,7)),ylim=(c(0,7)),
     cex.main=0.9,
     col=col_vector[11],pch=19,cex=0.7,
     xaxt="n",yaxt="n")
axis(side=1, at=c(0,2,4,6))
axis(side=2, at=c(0,2,4,6))


#Skin
filter= grep("A375|FIBRNPC|MCH58|SKMEL",rownames(result_binddat))
pvals_p_tmp <- as.data.frame(result_binddat$Pval[filter])
pvals_p <- na.omit(pvals_p_tmp)
dotforplot2_tmp <- c(0:round(nrow(pvals_p)/5),seq(round(nrow(pvals_p)/5),round(nrow(pvals_p)/4),5),seq(round(nrow(pvals_p)/4),round(nrow(pvals_p)/3),10),seq(round(nrow(pvals_p)/3),round(nrow(pvals_p)/2),20),seq(round(nrow(pvals_p)/2),nrow(pvals_p),40))
dotforplot2=nrow(pvals_p)-dotforplot2_tmp
distribution_tmp=seq(0,1,length=nrow(pvals_p)+2)
distribution=distribution_tmp[-(length(distribution_tmp))]
distribution=distribution[-1]
par(new=T)
plot(sort(-log10(distribution))[dotforplot2],sort(-log10(pvals_p$`result_binddat$Pval[filter]`))[dotforplot2],
     xlab="", ylab="",
     xlim=(c(0,7)),ylim=(c(0,7)),
     cex.main=0.9,
     col=col_vector[12],pch=19,cex=0.7,
     xaxt="n",yaxt="n")
axis(side=1, at=c(0,2,4,6))
axis(side=2, at=c(0,2,4,6))




#Uterus
filter= grep("HEC108|HELA|JHUEM2|SNGM",rownames(result_binddat))
pvals_p_tmp <- as.data.frame(result_binddat$Pval[filter])
pvals_p <- na.omit(pvals_p_tmp)
dotforplot2_tmp <- c(0:round(nrow(pvals_p)/5),seq(round(nrow(pvals_p)/5),round(nrow(pvals_p)/4),5),seq(round(nrow(pvals_p)/4),round(nrow(pvals_p)/3),10),seq(round(nrow(pvals_p)/3),round(nrow(pvals_p)/2),20),seq(round(nrow(pvals_p)/2),nrow(pvals_p),40))
dotforplot2=nrow(pvals_p)-dotforplot2_tmp
distribution_tmp=seq(0,1,length=nrow(pvals_p)+2)
distribution=distribution_tmp[-(length(distribution_tmp))]
distribution=distribution[-1]
par(new=T)
plot(sort(-log10(distribution))[dotforplot2],sort(-log10(pvals_p$`result_binddat$Pval[filter]`))[dotforplot2],
     xlab="", ylab="",
     xlim=(c(0,7)),ylim=(c(0,7)),
     cex.main=0.9,
     col=col_vector[13],pch=19,cex=0.7,
     xaxt="n",yaxt="n")
axis(side=1, at=c(0,2,4,6))
axis(side=2, at=c(0,2,4,6))



#y=x line
par(new=T)
abline(0,1,lty=3)


cols=c(col_vector)
labels=c("Adipose","Breast","CNS","Digestive","Haematopoietic","Liver",
         "Lung","Musculoskeletal","Ovary","Pancreas","Prostate","Skin","Uterus")

pchs=c(rep(19,13))
legend("bottomright", legend=labels, col=cols,pch=pchs,cex=0.75)

dev.off()



