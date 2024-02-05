
#### User-defined functions ####

add_legend = function(...) {
  #From: https://stackoverflow.com/questions/3932038/plot-a-legend-outside-of-the-plotting-area-in-base-graphics
  opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), 
              mar=c(0, 0, 0, 0), new=TRUE)
  on.exit(par(opar))
  plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
  legend(...)
}


mean_inf.rm = function(x=NULL){
  x[is.infinite(x)]=NA
  return(mean(x,na.rm=TRUE))
}

compute_perplexity = function(dtm=NULL,PWT=NULL,PDT=NULL){
  # Note: in case of K-fold CV, dtm is the test dtm
  # PWT is obtained using the train model, PDT is obtained by running
  # another model on the test dtm
  
  PDT = PDT/apply(PDT,1,sum) #re-normalize PD_T so it sums 1 row-wise
  prpx = text2vec::perplexity(topic_word_distribution = t(PWT),
                              doc_topic_distribution = PDT,
                              X = Matrix::Matrix(dtm))
  return(prpx)
}

word_overlap = function(FLSA_out=NULL,words,maxW=30){
  #Compute overlapping among topics considering the first maxW words
  noc=NCOL(FLSA_out[[1]]$PW_T)
  U=diag(noc)
  Z=gtools::combinations(n = noc,r = 2,v = 1:noc,repeats.allowed = FALSE)
  U[Z]=mapply(function(i)sum(words[order(FLSA_out[[1]]$PW_T[,Z[i,1]],decreasing = TRUE)][1:maxW] %in% words[order(FLSA_out[[1]]$PW_T[,Z[i,2]],decreasing = TRUE)][1:maxW])/maxW,1:NROW(Z))
  U[lower.tri(U)]=U[upper.tri(U)]
  rownames(U)=colnames(U)=paste0("topic",1:noc)
  return(U)
}

exclusive_term_ratio = function(FLSA_out=NULL,maxW=30){
  noc=NCOL(FLSA_out[[1]]$PW_T)
  etr=rep(NA,noc)
  for(j in 1:noc){
    maxWords = names(sort(FLSA_out[[1]]$PW_T[,j],decreasing = TRUE))[1:maxW]
    iid=which(names(FLSA_out[[1]]$PW_T[,j])%in%maxWords)
    if(is.null(dim(FLSA_out[[1]]$PW_T))){
      z=sum(FLSA_out[[1]]$PW_T[iid])
    }else{
      z=FLSA_out[[1]]$PW_T[iid,-j]
    }
    etr[j]=mean(FLSA_out[[1]]$PW_T[iid,j]/z)
  }
  names(etr)=paste0("topic",1:noc)
  return(list(raw=etr,normalized=etr/sum(etr)))
}

intertopic_distanceMap = function(FLSA_out=NULL,doc_lengths=NULL,plotx=TRUE,log_adj=2,lbls_cex=1.5,cols=NULL,cols_gps=NULL,pch_topics=20,k_rect=2,lbls=NULL,gps=NULL,coh_measure=NULL,new_window=TRUE,...){
  noc=NCOL(FLSA_out[[1]]$PW_T)
  if(is.null(lbls)){lbls=1:noc}
  if(is.null(gps)){xgps=1:noc*0}else{xgps=gps}
  if(is.null(coh_measure)){coh_measure="mean_logratio"}
  if(is.null(cols)){cols=paletteer::paletteer_c("grDevices::Heat 2",n=25,direction = -1)}
  x_cols=seq(from=min(FLSA_out[[2]]$topic_coherence[,coh_measure]),to=max(FLSA_out[[2]]$topic_coherence[,coh_measure]),length=length(cols)) 
  iid = sapply(FLSA_out[[2]]$topic_coherence[,coh_measure],function(x)which.min(abs(x-x_cols)))
  #if(is.null(cols)){cols=paletteer::paletteer_c("grDevices::Heat",n=25,direction = -1)}
  Cs = text2vec::sim2(x = t(FLSA_out[[1]]$PW_T),method="cosine")
  ACs = diag(noc)*0; ACs[lower.tri(ACs)]=acos(Cs[lower.tri(Cs)]); ACs[upper.tri(ACs)]=ACs[lower.tri(ACs)] #the angular distance represents the angle between the vectors in the high-dimensional space
  mds_res = stats::cmdscale(ACs,k = 2)
  ds_df = data.frame(mds_res, topics = seq_len(noc), Freq = FLSA_out[[1]]$PT*100,stringsAsFactors = FALSE,names=lbls,groups=xgps)
  
  if(!is.null(gps)){
    ugps=unique(ds_df$groups)
    if(is.null(cols_gps)){cgps=paletteer::paletteer_d("ggsci::planetexpress_futurama",n = length(ugps))}else{cgps=cols_gps}
    #uiid=sapply(ugps,function(x)chull(x = ds_df[ds_df$groups==x,1],y = ds_df[ds_df$groups==x,2]))
    }
  
  if(plotx==TRUE){
    if(new_window==TRUE){x11()}
    xmax=max(ds_df[,1:2])+(max(ds_df[,1:2])*2); xmin=min(ds_df[,1:2])-(max(ds_df[,1:2])*2)
    plot(ds_df[,1],ds_df[,2],cex=log(ds_df$Freq)*log_adj,bty="n",ylim=c(xmin,xmax),xlim=c(xmin,xmax),pch=pch_topics,col=cols[iid],xlab="dim 1",ylab="dim 2",...); 
    abline(h = 0,v = 0,lty=2); text(ds_df[,1],ds_df[,2],labels=lbls,cex=lbls_cex,font=2)
    #plot(ds_df[,1],ds_df[,2],cex=log(ds_df$Freq)*log_adj,bty="n",ylim=c(xmin,xmax),xlim=c(xmin,xmax),pch=pch_topics,col=cols[iid],xlab="dim 1",ylab="dim 2"); 
    
    if(is.null(gps)==FALSE){
      for(g in 1:length(ugps)){
        #Xps = ds_df[ds_df$groups==ugps[1],1:2][uiid[[1]],]
        Xps = ds_df[ds_df$groups==ugps[g],1:2]
        #polygon(x = Xps[,1],y = Xps[,2],border="blue",col="lightblue",lwd=2)
        #k1=abs(mean(Xps[,1]))/k_rect; k2=abs(mean(Xps[,2]))/k_rect
        k1=k2=abs(mean(apply(Xps[,1:2],2,mean)))/k_rect
        rect(xleft = min(Xps[,1])-k1,xright = max(Xps[,1])+k1,ybottom = min(Xps[,2])-k2,ytop = max(Xps[,2])+k2,border=cgps[g],lwd=3,lty=2,mex=20)
      }
    }
    
  }else{
    return(list(Cosine_sim=Cs,Angular_dist=ACs,ds_df=ds_df))  
  }
} 

metric_caoJuan2009 = function(topic_word_distrib){
  #http://doi.org/10.1016/j.neucom.2008.06.011
  cos_sim = 1-text2vec::sim2(topic_word_distrib, method='cosine') #direction: minimize (based on cosine similarity)
  cos_sim[upper.tri(cos_sim)]=0
  return(round(cos_sim,5))
}

metric_arun2010 = function(topic_word_distrib,doc_topic_distrib,doc_lengths){
  #http://doi.org/10.1007/978-3-642-13657-3_43
  # Note:: It will fail when num. of words in the vocabulary is less then the num. of topics (which is very unusual)
  
  # CM1 – sing. value decomp. of topic-word distrib.
  cm1 = svd(topic_word_distrib)$d #take singular values only
  cm1 = cm1/sum(cm1)     # normalize
  
  # CM2 – topics scaled by document lengths
  cm2 = doc_lengths %*% doc_topic_distrib
  cm2 = -sort(-cm2)   # sort in desc. order (just like cm1 is already sorted in desc. order)
  cm2 = cm2/sum(cm2)     # normalize
  
  # symmetric Kullback-Leibler divergence KL(cm1||cm2) + KL(cm2||cm1)
  # note: using log(x/y) instead of log(x) - log(y) here because values in cm vectors are not small
  kl = sum(cm1 * (log(cm1 / cm2))) + sum(cm2 * (log(cm2 / cm1))) #direction: minimize
  
  return(kl)
}

coherence = function(dtm_data=NULL,fcm_data=NULL,topic_words_dist=NULL,top_terms_topics=NULL,M=10,smth=1e-9,average=FALSE,metrics_text2vec=c("mean_logratio","mean_pmi","mean_npmi","mean_difference","mean_npmi_cosim","mean_npmi_cosim2"),probcoh=TRUE){
  # Note:
  ## (textmineR)
  ## Probabilistic coherence: it measures how associated words are in a topic, controlling for statistical independence.
  ## Averaging the probabilistic coherence across all topics is similar to using the silhouette coefficient to select the number of clusters when clustering.
  ## (text2vec)
  ## 'logratio', 'pmi', 'nmpi' usually opt for smaller numbers; 'mean_difference' 'mean_npmi_cosim' 'mean_npmi_cosim2' tend to propose higher numbers of topics.
  ## (1) mean_logratio: similar to the UMass metric, however, with a smaller smoothing constant by default and using the mean for aggregation instead of the sum
  ## This metric is more intrinsic in nature. It attempts to confirm that the models learned data known to be in the corpus
  ## (2) mean_pmi: similar to the UCI metric, however, with a smaller smoothing constant by default and using the mean for aggregation instead of the sum
  ## This metric is more extrinsic in nature (it needs an external tcm).
  ## (3) mean_npmi: it may perform better than the simpler pmi metric
  ## (4) mean_difference: intrinsic metric in the standard use case
  ## (5) mean_nmpi_cosim: extrinsic metric in the standard use case
  ## (6) mean_nmpi_cosim2: extrinsic metric in the standard use case
  
  probcoh = if(probcoh==TRUE){textmineR::CalcProbCoherence(phi = topic_words_dist,dtm = Matrix::Matrix(dtm_data),M)}else{NA}
  cohs=matrix(NA,nrow = NCOL(dtm_data),ncol = 1); if(!is.null(metrics_text2vec)){
    cohs=text2vec::coherence(x = top_terms_topics,tcm = Matrix::Matrix(fcm_data),n_doc_tcm = NROW(dtm_data),smooth=smth,metrics=metrics_text2vec)
    }
  out=cbind(cohs,probcoh)
  if(average==TRUE){out=apply(out,2,mean_inf.rm)}
  return(out)  
}

# coherence = function(dtm,top_terms,a_smooth=1){
#   dtm_m = as.matrix(dtm[,top_terms]); N = NROW(dtm)
#   maxW=length(top_terms)
#   Z = gtools::combinations(n = maxW,r = 2,v = 1:maxW); nZ = NROW(Z)
#   
#   # Coherence based on Cosine similarity computed on TF-IDF dtm matrix
#   X = as.matrix(tm::weightTfIdf(dtm)[,top_terms])
#   x = mapply(function(i)sum(apply(X[,Z[i,]],1,prod))/(sqrt(sum(X[,Z[i,1]]^2))*sqrt(sum(X[,Z[i,2]]^2))),1:nZ); x[is.na(x)]=0
#   ch1 = (2/(maxW^2-maxW))*sum(x)
#   
#   # Coherence based on Jaccard similarity computed on Boolean dtm matrix
#   X = dtm_m[,top_terms]>0
#   jacc_wi_wj = rep(NA,nZ)
#   for(i in 1:nZ){
#     int_wi_wj = sum(X[,Z[i,1]]==TRUE & X[,Z[i,2]]==TRUE)
#     union_wi_wj = sum(X[,Z[i,1]]==TRUE)+sum(X[,Z[i,2]]==TRUE)-int_wi_wj
#     jacc_wi_wj[i] = int_wi_wj/union_wi_wj
#   }
#   jacc_wi_wj[is.na(jacc_wi_wj)]=0
#   ch2 = (2/(maxW^2-maxW))*sum(jacc_wi_wj)
#   
#   # Coherence based on PMI computed on fcm matrix
#   X = as.matrix(quanteda::fcm(quanteda::as.dfm(dtm_m)))
#   ch3 = (2/(maxW^2-maxW))*log2((a_smooth+X[Z[i,1],Z[i,2]]*N)/prod(apply(a_smooth+dtm_m[,Z[i,]],2,sum)))
#   
#   # Coherence based on Mimno's et al formula (see topicdoc package)
#   X = slam::tcrossprod_simple_triplet_matrix(slam::as.simple_triplet_matrix(t(dtm_m>0)))
#   ch4 = sum(mapply(function(i)log((X[Z[i,1],Z[i,2]]+a_smooth)/(X[Z[i,1],Z[i,1]]+a_smooth)),1:nZ))
#   
#   return(list(cos_sim=ch1,jacc=ch2,pmi=ch3,mimno2011=ch4))
# }

topic_coherence_SR <- function(top_words, #from the SpeedReader library
                            document_term_matrix,
                            vocabulary = NULL,
                            numeric_top_words = FALSE,
                            K = length(top_words)){
  
  # make sure the data is the right format
  vocabulary <- as.character(vocabulary)
  
  # perform some basic checks and throw errors if we see something weird.
  if(is.null(vocabulary) & !numeric_top_words){
    stop("You must provide a vocabulary vector!")
  }
  if(K > length(top_words)){
    K <- length(top_words)
    warning(paste("You must select a value for K that is less than length(top_words). K has automatically been set to :",K,sep = " "))
  }
  if(length(vocabulary) != ncol(document_term_matrix)){
    stop("The vocaublary vector must have the same number of entries as the number of columns in the document_term_matrix, and the word indicated by entries in the i'th column of document_term_matrix must correspond to the i'th entry in vocabulary.")
  }
  
  #if we are only using the K top words then reduce our top words vector
  top_words <- top_words[1:K]
  
  # binarize the document term matrix
  document_term_matrix <- matrix(as.numeric(document_term_matrix > 0),
                                 nrow = nrow(document_term_matrix),
                                 ncol = ncol(document_term_matrix))
  coherence_score <- 0
  for(i in 2:length(top_words)){
    for(j in 1:(i-1)){
      # we can either look up against vocab or just use indexes
      if(numeric_top_words){
        jindex <- top_words[j]
        iindex <- top_words[i]
      }else{
        jindex <- which(vocabulary == top_words[j])
        iindex <- which(vocabulary == top_words[i])
      }
      
      document_frequency <- sum(document_term_matrix[,jindex])
      j_positive <- which(document_term_matrix[,jindex] > 0)
      i_positive <- which(document_term_matrix[,iindex] > 0)
      co_document_frequency <- sum(i_positive %in% j_positive)
      
      coherence_score <- coherence_score + log((co_document_frequency + 1)/document_frequency)
      
    }
  }
  if(is.infinite(coherence_score)){
    coherence_score <- NA
    warning("The coherence score was not finite. Make sure that all words in your vocabulary appear atleast once.")
  }
  return(coherence_score)
}

calcCalinskiHarabasz = function(data, belongmatrix, centers){
  # using this formula : https://www.geeksforgeeks.org/calinski-harabasz-index-cluster-validity-indices-set-3/
  data <- as.matrix(data)
  c <- apply(data,2,mean)
  k <- ncol(belongmatrix)
  nk <- colSums(belongmatrix)
  p1 <- (sum(rowSums((centers - c) **2) * nk)) / (k-1)
  p2 <- sum(apply(centers, 1, function(ck){
    sum((sweep(data, 2, ck, `-`))**2)
  })) / (nrow(data)-k)
  return(p1/p2)
}

compute_communality = function(dfmx = NULL){
  m = NROW(dfmx)
  Y=as.matrix(dfmx); 
  Y[Y>1] = 1
  x = mapply(function(j)sum(Y[,j]==rep(1,m)),1:NCOL(Y))
  comun_overall = sum(x==m)/length(x)
  
  X = matrix(NA,m,m); rownames(X) = colnames(X) = paste0("t",rep(1:m))
  for(i in 1:m){
    for(j in 1:m){
      x = mapply(function(k)sum(Y[c(i,j),k]==rep(1,2)),1:NCOL(Y))
      X[i,j] = sum(x==2)/length(x[x>0])
    }
  }
  
  return(list(X=X,overall=comun_overall))
}




#### Functions from 'topicdoc' library (https://github.com/cran/topicdoc/) ####
# Note: functions have been adapted to work with fuzzythe -LSA implementation #
require(slam); require(topicmodels)

coherence2 <- function(dtm_data, top_terms, smoothing_beta=1){
  # Get the relevant entries of the document-term matrix
  rel_dtm <- dtm_data[,top_terms]

  # Turn it into a logical representing co-occurences
  df_dtm <- rel_dtm > 0

  # Calculate document frequencies for each term and all of its co-occurences
  cooc_mat <- slam::tcrossprod_simple_triplet_matrix(t(df_dtm))

  # Quickly get the number of top terms for the for-loop below
  top_n_tokens <- length(top_terms)

  # Using the syntax from the paper, calculate coherence
  c_l <- 0
  for (m in 2:top_n_tokens) {
    for (l in 1:(top_n_tokens - 1)) {
      df_ml <- cooc_mat[m,l]
      df_l <- cooc_mat[l,l]
      #c_l <- c_l + log((df_ml + smoothing_beta) / (df_l + smoothing_beta*top_n_tokens)) #with Laplacian smoothing
      c_l <- c_l + log((df_ml + smoothing_beta) / df_l)
    }
  }
  c_l
}

topic_size <- function(PW_T){
  # Obtain the beta matrix from the topicmodel object
  #beta_mat <- exp(topic_model@beta)
  beta_mat <- t(PW_T)
  # Normalize the beta values within each topic
  # SO link for reference - https://stats.stackexchange.com/a/51750
  beta_normed <- beta_mat %*% diag(1/colSums(beta_mat))
  
  # Sum the partial tokens per topic
  rowSums(beta_normed, na.rm = TRUE)
}

mean_token_length <- function(top_terms, top_n_tokens = 10){
  # Obtain the top terms from the topicmodel object
  #top_terms <- terms(topic_model, top_n_tokens)
  
  # Calculate the number of characters per token in each topic
  nchar_mat <- apply(top_terms, 2, nchar)
  
  # Calculate the averages for each topic
  unname(colMeans(nchar_mat))
}

dist_from_corpus <- function(PW_T, dtm_data){
  # Obtain the beta matrix from the topicmodel object
  #beta_mat <- exp(topic_model@beta)
  beta_mat <- t(PW_T)
  
  # Coerce dtm to slam format
  dtm_data <- slam::as.simple_triplet_matrix(dtm_data)
  
  # Calculate token frequency across all documents
  global_tf_counts <- slam::col_sums(dtm_data, na.rm = TRUE)
  
  # Get corpus-level probability of each token's occurence
  corpus_dist <- global_tf_counts/sum(global_tf_counts)
  
  # Using the Hellinger distance, calculate the distance
  # of each topic's token distribution from the corpus distribution
  topicmodels::distHellinger(beta_mat, matrix(corpus_dist, nrow = 1))[,1]
}

tf_df_dist <- function(top_terms, dtm_data){
  # Obtain the top terms from the topicmodel object
  #top_terms <- terms(topic_model, top_n_tokens)
  
  # Coerce dtm to slam format
  dtm_data <- slam::as.simple_triplet_matrix(dtm_data)
  
  # Calculate distance b/w token frequency/document frequency for each topic
  unname(apply(top_terms, 2, tf_df_dist_diff, dtm_data = dtm_data))
}

tf_df_dist_diff <- function(dtm_data, top_terms){
  # Obtain the indicies of the top terms in the dtm
  # and select only those
  top_terms_inds <- which(colnames(dtm_data) %in% top_terms)
  rel_dtm <- dtm_data[,top_terms_inds]
  
  # Calculate the token frequencies and document frequencies
  # using slam to keep things sparse
  tf_counts <- slam::col_sums(rel_dtm)
  df_counts <- slam::col_sums(rel_dtm > 0)
  
  # Using the Hellinger distance, calculate the distance
  # of each topic's token frequencies from its document frequencies
  topicmodels::distHellinger(matrix(tf_counts, nrow = 1),
                matrix(df_counts, nrow = 1))
}

doc_prominence <- function(PD_T, method = c("largest_gamma","gamma_threshold"),
                                      gamma_threshold = 0.2){
  # Ensure the user passed a valid method argument
  method <- match.arg(method)
  
  # Obtain the gamma matrix from the topicmodel object
  #gamma_mat <- topic_model@gamma
  gamma_mat <- PD_T
  
  if (method == "gamma_threshold") {
    # Count the number of documents per topic that exceed the gamma threshold
    colSums(gamma_mat > gamma_threshold)
  } else {
    # Find the topic with the largest gamma per document
    row_maxs <- max.col(gamma_mat, ties.method = "first")
    
    # Sum up the results
    as.vector(table(row_maxs))
  }
}

topic_exclusivity <- function(PW_T, top_terms,terms,excl_weight = 0.5){
  # excl_weight: a numeric between 0 and 1 indicating the weight to place on exclusivity
  # versus frequency in the calculation (0.5 is the default)
   
  # Obtain the beta matrix from the topicmodel object
  #beta_mat <- exp(topic_model@beta)
  beta_mat <- t(PW_T)
  
  # Normalize the beta values within each topic
  # SO link for reference - https://stats.stackexchange.com/a/51750
  beta_normed <- beta_mat %*% diag(1/colSums(beta_mat))
  
  # Calculate exclusivity (using approx ECDF)
  excls <- apply(beta_normed, 1, rank)/ncol(beta_normed)
  
  # Calculate frequency (using approx ECDF)
  freqs <- apply(beta_mat, 1, rank)/ncol(beta_mat)
  
  # Obtain the indicies of the top terms in the model
  # and select only those from the exclusivity and frequency matrices
  #top_terms <- terms(topic_model, top_n_tokens)
  
  # Create an empty vector for the frex scores
  frex <- vector(length = ncol(freqs))
  
  # Loop through each topic's terms and calculate its total frex score
  for (i in 1:ncol(freqs)) {
    # Identifying and selecting rows for this topic's terms
    term_inds <- which(terms %in% top_terms[,i])
    this_excls <- excls[term_inds, i]
    this_freqs <- freqs[term_inds, i]
    
    # Calculate frex score using the provided exclusivity weight
    excl_term <- excl_weight / this_excls
    freq_term <- (1 - excl_weight) / this_freqs
    frex[i] <- sum(1 / (excl_term + freq_term))
  }
  frex
}




