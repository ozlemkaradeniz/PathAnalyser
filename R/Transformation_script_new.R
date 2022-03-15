#Function for data transformation 
require(edgeR)
Transform <-function (data_es){
  
    data_es_log_cpm<-cpm(data_es, log=TRUE)
    return(data_es_log_cpm)
  
}


