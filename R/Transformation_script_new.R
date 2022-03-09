#Function for data transformation 
require(DESeq2)
Transform <-function (data_es,method="log"){
  if(missing(method)){
    data_es_log<-log10(data_es)
    return(data_es_log)
  }
  else if(method=="vst"){
    data_es_vst<-vst(data_es)
    return(data_es_vst)
  }
   else if(method=="sqrt"){
     data_es_sqrt<-sqrt(data_es)
     return(data_es_sqrt)
   }
  else{
    data_es_log<-log10(data_es)
    return(data_es_log)
  }
}


