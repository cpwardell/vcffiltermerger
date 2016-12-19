## Description - functions to process filter data

## Load packages (none presently)

## Make a call based on the current filter values
########################## Alignability - COMPLETE
aligntest=function(input){
  input==1
}
#aligntest(dfcifilters$tumor.alignability)
#table(aligntest(dfcifilters[dfcifilters$caller%in%"mutect2","tumor.alignability"]))
#table(aligntest(dfcifilters[dfcifilters$caller%in%"strelkaindels","tumor.alignability"]))
#table(aligntest(dfcifilters[dfcifilters$caller%in%"strelkaSNVs","tumor.alignability"]))

########################## On target - COMPLETE
ontargettest=function(input){
  input==1
}
#table(ontargettest(dfcifilters[dfcifilters$caller%in%"mutect2","tumor.ontarget"]))
#table(ontargettest(dfcifilters[dfcifilters$caller%in%"strelkaindels","tumor.ontarget"]))
#table(ontargettest(dfcifilters[dfcifilters$caller%in%"strelkaSNVs","tumor.ontarget"]))

########################## Base quality - COMPLETE
basequaltest=function(input){
  input>=25
}
#table(basequaltest(dfcifilters[dfcifilters$caller%in%"mutect2","tumor.base_quality"]))
#table(basequaltest(dfcifilters[dfcifilters$caller%in%"strelkaindels","tumor.base_quality"]))
#table(basequaltest(dfcifilters[dfcifilters$caller%in%"strelkaSNVs","tumor.base_quality"]))

########################## Mapping quality - COMPLETE
mapqualtest=function(input){
  input>=50
}
#table(mapqualtest(dfcifilters[dfcifilters$caller%in%"mutect2","tumor.mapping_quality"]))
#table(mapqualtest(dfcifilters[dfcifilters$caller%in%"strelkaindels","tumor.mapping_quality"]))
#table(mapqualtest(dfcifilters[dfcifilters$caller%in%"strelkaSNVs","tumor.mapping_quality"]))

########################## Read direction - COMPLETE
readdirectiontest=function(input){
  apply(input,1,function(x){fisher.test(matrix(x,nrow=2,byrow=TRUE))$p.value>0.05})
}
#x=readdirectiontest(dfcifilters[,c("tumor.read_direction","tumor.read_direction.2","tumor.read_direction.3","tumor.read_direction.4")])
#table(x)

########################## OxoG - COMPLETE
oxogtest=function(oxog,tvaf){ 
  
  ## Define error function for tryCatch
  #returnNA=function(anything){NA}
  
  x=-0.5+(1)*oxog<tvaf
  
  #x=tryCatch(-0.5+(1)*oxog<tvaf,error=returnNA)
  
  #if(length(x)==1 & is.na(x)[1]){
  #  x=rep(NA,length(oxog))
  #}
  x[is.na(x)]=TRUE
  x
}
#table(oxogtest(dfcifilters[dfcifilters$caller%in%"mutect2","tumor.oxog"],dfcifilters[dfcifilters$caller%in%"mutect2","tumor.vaf.6"]))
#table(oxogtest(dfcifilters[dfcifilters$caller%in%"strelkaindels","tumor.oxog"],dfcifilters[dfcifilters$caller%in%"strelkaindels","tumor.vaf.6"]))
#table(oxogtest(dfcifilters[dfcifilters$caller%in%"strelkaSNVs","tumor.oxog"],dfcifilters[dfcifilters$caller%in%"strelkaSNVs","tumor.vaf.6"]))
########################## Min Tdepth - COMPLETE
tdepthtest=function(input){
  input>=10
}
#table(tdepthtest(dfcifilters[dfcifilters$caller%in%"mutect2","tumor.depth"]))
#table(tdepthtest(dfcifilters[dfcifilters$caller%in%"strelkaindels","tumor.depth"]))
#table(tdepthtest(dfcifilters[dfcifilters$caller%in%"strelkaSNVs","tumor.depth"]))

########################## Min Ndepth - COMPLETE
ndepthtest=function(input){
  input>=10
}
#table(ndepthtest(dfcifilters[dfcifilters$caller%in%"mutect2","normal.depth"]))
#table(ndepthtest(dfcifilters[dfcifilters$caller%in%"strelkaindels","normal.depth"]))
#table(ndepthtest(dfcifilters[dfcifilters$caller%in%"strelkaSNVs","normal.depth"]))

########################## Max Tdepth - COMPLETE
tmaxdepthtest=function(inputdata){ 
  ## Read in mean depths; text file in same directory as this script
  ## These values are calculated using hsmetrics in Picard
  ## Format is three columns
  ## sample  MEAN_TARGET_COVERAGE    MEDIAN_TARGET_COVERAGE
  ## CTD_24_diagnosis        59.081339       51
  depths=read.table("depths.tsv",header=TRUE,stringsAsFactors=FALSE,row.names=1)
  ## Is depth less than 3x mean depth?
  results=inputdata$tumor.depth<3*depths[tumorname,"MEAN_TARGET_COVERAGE"]
  results[is.na(results)]=TRUE ## if there is no data, just accept all results
  results
}
#table(tmaxdepthtest(boundfilters))

########################## min TVAF - COMPLETE
tvaftest=function(input){ # could be 20% for indels?
  input>=0.05
}
#table(tvaftest(dfcifilters[dfcifilters$caller%in%"mutect2","tumor.vaf.6"]))
#table(tvaftest(dfcifilters[dfcifilters$caller%in%"strelkaindels","tumor.vaf.6"]))
#table(tvaftest(dfcifilters[dfcifilters$caller%in%"strelkaSNVs","tumor.vaf.6"]))
########################## max NVAF - COMPLETE
nvaftest=function(input){ 
  input<0.05
}
#table(nvaftest(dfcifilters[dfcifilters$caller%in%"mutect2","normal.vaf.6"]))
#table(nvaftest(dfcifilters[dfcifilters$caller%in%"strelkaindels","normal.vaf.6"]))
#table(nvaftest(dfcifilters[dfcifilters$caller%in%"strelkaSNVs","normal.vaf.6"]))

tlevtest=function(input){ # Max 2.5, 2.2 seems reasonable
  input<=2.2
}
#table(tlevtest(dfcifilters[dfcifilters$caller%in%"mutect2","tumor.altlevendist.7"]))
#table(tlevtest(dfcifilters[dfcifilters$caller%in%"strelkaindels","tumor.altlevendist.7"]))
#table(tlevtest(dfcifilters[dfcifilters$caller%in%"strelkaSNVs","tumor.altlevendist.7"]))

minalttest=function(input){ # Minimum 3 ALT reads required
  input>=3
}
#table(minalttest(dfcifilters[dfcifilters$caller%in%"mutect2","tumor.vaf.4"]))

## Calculate logical pass/fail for each filter
filterset=function(input){
  message("Processing alignresult")
  alignresult=aligntest(input$tumor.alignability)
  message("Processing ontargetresult")
  ontargetresult=ontargettest(input$tumor.ontarget)
  message("Processing basequalresult")
  basequalresult=basequaltest(input$tumor.base_quality)
  message("Processing mapqualtest")
  mapqualtest=mapqualtest(input$tumor.mapping_quality)
  #message("Processing readdirectionresult")
  #readdirectionresult=readdirectiontest(input[,c("tumor.read_direction","tumor.read_direction.2","tumor.read_direction.3","tumor.read_direction.4")])
  message("Processing oxogtest")
  oxogtest=oxogtest(input$tumor.oxog,input$tumor.vaf.6)
  message("Processing tdepthresult")
  tdepthresult=tdepthtest(input$normal.depth)
  message("Processing ndepthresult")
  ndepthresult=ndepthtest(input$normal.depth)
  message("Processing tmaxdepthresult")
  tmaxdepthresult=tmaxdepthtest(input)
  message("Processing tvafresult")
  tvafresult=tvaftest(input$tumor.vaf.6)
  message("Processing nvafresult")
  nvafresult=nvaftest(input$normal.vaf.6)
  message("Processing tlevresult")
  tlevresult=tlevtest(input$tumor.altlevendist.6)
  message("Processing minaltresult")
  minaltresult=minalttest(input$tumor.vaf.4)
  message("Concatenating results")
# results=cbind(alignresult,ontargetresult,basequalresult,mapqualtest,
#                readdirectionresult,oxogtest,tdepthresult,ndepthresult,
#                tmaxdepthresult,tvafresult,nvafresult)
  results=cbind(alignresult,ontargetresult,basequalresult,mapqualtest,
                oxogtest,tdepthresult,ndepthresult,
                tmaxdepthresult,tvafresult,nvafresult,tlevresult,minaltresult)
}
#dfciresults=filterset(dfcifilters)
#uamsresults=filterset(uamsfilters)
#mmrfresults=filterset(mmrffilters)

## Output final results file:
#write.table(dfcifilters,file="dfcifilters.20161101.txt",quote=FALSE,sep="\t",row.names=FALSE,col.names=TRUE)
#write.table(uamsfilters,file="uamsfilters.20161101.txt",quote=FALSE,sep="\t",row.names=FALSE,col.names=TRUE)
#write.table(mmrffilters,file="mmrffilters.20161101.txt",quote=FALSE,sep="\t",row.names=FALSE,col.names=TRUE)
