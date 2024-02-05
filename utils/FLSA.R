FLSA = function(X=NULL,dtm=NULL,fcm=NULL,words=NULL,noc=NULL,maxW=30,lambda=0.6,fc_validity=TRUE,compute_stats=TRUE){
  #X: local/global weighted dtm
  #dtm: raw dtm
  #words: dictionary
  #noc: number of topics
  #maxW: maximum number of words per topic
  
  M = NCOL(X); N = NROW(X) 
  pw = matrix(apply(X,2,sum)/sum(X),M,1)
  #pd = matrix(1/N,N,1)
  pd = matrix(apply(X,1,sum)/sum(X),N,1)
  PW_D = diag(1/apply(X,1,sum))%*%X       #probability of word w in document d
  dtm_dimred = irlba::irlba(A = X, nv = 2) #via SVD as in the standard LSA technique
  #out = fclust::FKM.med(X = dtm_dimred$u,k = noc)
  out = fclust::FKM(X = dtm_dimred$u,k = noc)
  fc_meas=NULL;if(fc_validity==TRUE){fc_meas = fclust::Fclust.index(out,index = c("PE","MPC","SIL.F"));names(fc_meas)=c("PE","MPC","SIL.F")}
  PT_D = as.matrix(out$U)
  PTD = PT_D*(pd%*%matrix(1,1,noc))
  PD_T = PTD/(t(apply(PTD,2,sum)%*%matrix(1,1,N))); colnames(PD_T)=paste0("topic",1:noc);
  PW_T = t(PW_D) %*% PD_T; colnames(PW_T)=paste0("topic",1:noc);                                         
  
  PD_T_sorted = t(mapply(function(i)order(PD_T[i,],decreasing = TRUE),1:N)) #most likely topic for each document
  PW_T_sorted = mapply(function(j)names(PW_T[order(PW_T[,j],decreasing = TRUE),j][1:maxW]),1:noc) #most frequent terms for each topic
  colnames(PD_T_sorted)=paste0("topic",1:noc); colnames(PW_T_sorted)=paste0("topic",1:noc);
  
  PW_T_frex = stm::calcfrex(logbeta = log(t(PW_T)),wordcounts = apply(dtm,2,sum)); rownames(PW_T_frex) = rownames(PW_T)
  Frex_sorted = mapply(function(j)names(PW_T_frex[order(PW_T_frex[,j],decreasing = TRUE),j][1:maxW]),1:noc) #most frequent terms for each topic
  colnames(Frex_sorted)=paste0("topic",1:noc);
  
  Word_Relevance = lambda*log(PW_T) + (1-lambda)*log(PW_T/(pw%*%matrix(1,1,noc)))
  Word_Relevance_sorted = mapply(function(j)names(Word_Relevance[order(Word_Relevance[,j],decreasing = TRUE),j][1:maxW]),1:noc) 
  colnames(Word_Relevance_sorted)=paste0("topic",1:noc);
  
  doc_lengths = apply(dtm,1,sum)
  PT_unno = apply(PD_T*doc_lengths%*%matrix(1,1,noc),2,sum); 
  PT = PT_unno/sum(PT_unno) #marginal topic distribution P(T)
  pw2 = PW_T%*%PT
  
  FLSA=list(PW=list(fromData=pw,fromModel=pw2),PD=pd,PT=PT,PT_D=PT_D,PTD=PTD,PD_T=PD_T,PW_T=PW_T,PD_T_sort=PD_T_sorted,PW_T_sort=PW_T_sorted,Frex_W_T=PW_T_frex,Frex_W_T_sort=Frex_sorted,Relevance_W_T=Word_Relevance,Relevance_W_T_sort=Word_Relevance_sorted,dimred=dtm_dimred,fc_validity=fc_meas)
  
  if(compute_stats==TRUE){
    #coher = mapply(function(j)coherence(slam::as.simple_triplet_matrix(dtm),PW_T_sorted[1:top_n_tokens,j]),1:noc); names(coher)=paste0("topic",1:noc)
    #Coher = mapply(function(j)coherence(dtm,PW_T_sorted[,j]),1:noc); colnames(Coher)=paste0("topic",1:noc)
    Coher = coherence(dtm = dtm,fcm = fcm,topic_words_dist = t(PW_T),top_terms_topics = PW_T_sorted)
    topicSize = topic_size(PW_T)
    meanTopicLenght = mean_token_length(PW_T_sorted); names(meanTopicLenght)=paste0("topic",1:noc)
    corpusDistance = dist_from_corpus(PW_T,dtm); names(corpusDistance)=paste0("topic",1:noc)
    ftDfDistance = tf_df_dist(PW_T_sorted,dtm); names(ftDfDistance)=paste0("topic",1:noc)
    
    documentsProminence = doc_prominence(PD_T); names(documentsProminence)=paste0("topic",1:noc)
    exclusivity = topic_exclusivity(PW_T,PW_T_sorted,rownames(PW_T)); names(exclusivity)=paste0("topic",1:noc)
    
    ma10 = metric_arun2010(t(PW_T),PD_T,doc_lengths)
    mcj09 = metric_caoJuan2009(t(PW_T))
    
    #Coher_matrix = matrix(unlist(Coher),4,noc); rownames(Coher_matrix)=c("cos_sim","jacc","pmi","mimno2011"); colnames(Coher_matrix)=paste0("topic",1:noc)
    stats=list(topic_coherence=Coher,topic_size=topicSize,topic_length=meanTopicLenght,topic_exclusivity=exclusivity,dist_topics_corpus=corpusDistance,
               dist_ftDf=ftDfDistance,document_prominence=documentsProminence,metric_arun2010=ma10,metric_caoJuan09=mcj09)
  }else{
    stats=NULL
  }
  
  out = list(FLSA=FLSA,stats=stats)
  return(out)
}

