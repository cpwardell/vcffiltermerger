## Reads in VCFs, filter results and annotations and outputs a single file

## Load packages
library(VariantAnnotation)

## RUN IT LIKE THIS: /c/Program\ Files/R/R-3.3.2/bin/Rscript.exe example.R vcfpath filterpath annopath
## /c/Program\ Files/R/R-3.3.2/bin/Rscript.exe /g/Genomics\ Lab/MGP/vcffiltermerger/vcf_filter_merger.R "Z://chris_working_dir/project_erasmusmc/erasmusmc_strelka/_localcopy/T05208Tumor.all.somatic.snvs.vcf" "Z://chris_working_dir/project_erasmusmc/erasmusmc_strelka/filteredsnvs/T05208Tumor" "Z://cody_working_dir/erasmus_annotations/T05208Tumor.all.somatic.snvs.vcf.tsv"

## Command line options; specify paths to VCF, filters and annotations
args = commandArgs(trailingOnly=TRUE)

vcfpath=args[1]
filterpath=args[2]
annopath=args[3]

#vcfpath="Z://chris_working_dir/project_erasmusmc/erasmusmc_strelka/_localcopy/EMN_4095_diagnosis.all.somatic.indels.vcf"
#filterpath="Z://chris_working_dir/project_erasmusmc/erasmusmc_strelka/filteredindels/EMN_4095_diagnosis"
#annopath="Z://cody_working_dir/erasmus_annotations/EMN_4095_diagnosis.all.somatic.indels.vcf.tsv"

#vcfpath="Z://chris_working_dir/project_erasmusmc/erasmusmc_strelka/_localcopy/T05208Tumor.all.somatic.snvs.vcf"
#filterpath="Z://chris_working_dir/project_erasmusmc/erasmusmc_strelka/filteredsnvs/T05208Tumor"
#annopath="Z://cody_working_dir/erasmus_annotations/T05208Tumor.all.somatic.snvs.vcf.tsv"

## Indel files with duplicate rownames
#vcfpath="Z://cody_working_dir/mgp_oncotator/UK_data/Strelka_1470762457/_EGAR00001321238_EGAS00001001147_C273MACXX_2_378.seqvar/all.somatic.indels.vcf"
#filterpath="Z://chris_working_dir/project_mgp/SNEF/strelkaindels/retry/C273MACXX_2_378"
#annopath="Z://cody_working_dir/mgp_oncotator/UK_data/Strelka_1470762457/_EGAR00001321238_EGAS00001001147_C273MACXX_2_378.seqvar/all.somatic.indels.tsv"

#################
## Read in VCF ##
#################
vcf=readVcf(vcfpath,genome="GRCh37")
########################
## End of read in VCF ##
########################


#########################################################
## Read in filters and merge them into a single object ##
#########################################################
filterfiles=dir(filterpath,pattern="txt")
boundfilters=as.data.frame(matrix(nrow=0,ncol=0))

## Discover tumor name and normal name from filenames
## Tumor name will occur more often
samplenames=sort(table(sapply(strsplit(filterfiles,"\\."),"[",1)))
normalname=names(samplenames)[1]
tumorname=names(samplenames)[2]

## Reorder filterfiles to that tumor is always first
filterfiles=c(sort(grep(tumorname,filterfiles,value=TRUE)),sort(grep(normalname,filterfiles,value=TRUE)))

for(i in 1:length(filterfiles)){
  tumornormal=ifelse(grepl(tumorname,filterfiles[i]),"tumor.","normal.")
  #thisfilter=read.delim(paste(filterpath,filterfiles[i],sep="/"),header=FALSE,stringsAsFactors=FALSE,row.names=1)
  thisfilter=read.delim(paste(filterpath,filterfiles[i],sep="/"),header=FALSE,stringsAsFactors=FALSE)#,row.names=1)
  #thisfilter=thisfilter[,2:ncol(thisfilter)]
  thisfilter=as.data.frame(thisfilter[,2:ncol(thisfilter)])
  #rownames(thisfilter)=1:nrow(thisfilter)
  filtername=sapply(strsplit(filterfiles[i],"\\."),"[",2)
  filtername=ifelse(grepl(tumorname,filterfiles[i]),paste0(tumornormal,filtername),paste0(tumornormal,filtername))
  colnames(thisfilter)=ifelse(c(FALSE,rep(TRUE,ncol(thisfilter)-1)),paste(filtername,1:ncol(thisfilter),sep="."),filtername)
  
  ## Different behaviour for 1st filter
  if(i==1){
    boundfilters=rbind(boundfilters,thisfilter)
    }else{
    boundfilters=cbind(boundfilters,thisfilter)
  }
}
################################################################
## End of read in filters and merge them into a single object ##
################################################################

#########################################
## Calculate pass/fail for all filters ##
#########################################

## Calculate filter results and append these (i.e. TRUE/FALSE TABLE)
## 1st step; source() filter functions from a separate file
setwd("G:/Genomics Lab/MGP/vcffiltermerger")
source("filter_functions.R")
filterresults=filterset(boundfilters)

################################################
## End of calculate pass/fail for all filters ##
################################################


###################################
## Read in Oncotator annotations ##
###################################
anno=read.delim(annopath,comment.char = "#",stringsAsFactors=FALSE)
anno=anno[anno$"Matched_Norm_Sample_Barcode"%in%"TUMOR",]

## Remove duplicate records: CAUTION, WE MAY BE REMOVING THE INCORRECT RECORDS
coords=paste(anno$Chromosome,anno$Start_position,sep=":")
anno=anno[!(duplicated(coords)),]

##########################################
## End of read in Oncotator annotations ##
##########################################


##############################
## Merge all data and write ##
##############################

## Files to merge:
#vcf
#boundfilters
#filterresults
#anno

#CHROM   POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  TUMOR   NORMAL

## We dump any variant where the ALT allele wasn't recorded.  This appears to be a weird
## behaviour by Strelka; very low VAF alleles are not recorded... I think
nones=boundfilters$normal.vaf.2%in%"None"
vcf=vcf[!nones,]
boundfilters=boundfilters[!nones,]
filterresults=filterresults[!nones,]
anno=anno[!nones,]

## Make a fake vcf... this is easier than querying the original object
CHROM=sapply(strsplit(rownames(boundfilters),":"),"[",1)
POS=sapply(strsplit(rownames(boundfilters),":"),"[",2)
REF=boundfilters$normal.vaf
ALT=boundfilters$normal.vaf.2

final=cbind(CHROM,POS,REF,ALT,tumorname)
final=cbind(final,boundfilters,filterresults,anno)

write.table(final,file=paste0(tumorname,".txt"),quote=FALSE,sep="\t",row.names=FALSE,col.names=TRUE)

