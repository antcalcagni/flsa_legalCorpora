#setwd("/home/antonio.calcagni/storage/FLSA_sds24")
rm(list=ls()); options(warn = -1)
source("utils/utils.R"); source("utils/FLSA.R")
library(doParallel)
set.seed(111);

ncores=44
cl = parallel::makeCluster(ncores)
doParallel::registerDoParallel(cl,cores=ncores)
combAbind = abind::abind

cv_evalTopics = function(K=5,cvs=NULL,Xw=NULL,dtm_data=NULL,words_dic=NULL,noc=NULL){
  Cvout = matrix(nrow=12,ncol=K)
  for(k in 1:K){
    cat("\n CV..",k)
    train_iid = cvs$subsets[cvs$which==k]; test_iid = cvs$subsets[cvs$which!=k]
    train_model = FLSA(X=Xw[train_iid,],dtm=dtm_data[train_iid,],words=words_dic[train_iid],noc=noc,compute_stats=FALSE)
    fcm_test = Matrix::crossprod(Matrix::Matrix(dtm_data[test_iid,]),Matrix::Matrix(dtm_data[test_iid,]))
    pdT_test = FLSA(X=Xw[test_iid,],dtm=dtm_data[test_iid,],words=words_dic[test_iid],noc=noc,fc_validity=FALSE,compute_stats=FALSE)[[1]]$PD_T
    Cvout[,k] = c(coherence(dtm_data=dtm_data[test_iid,],fcm_data=fcm_test,topic_words_dist=t(train_model[[1]]$PW_T),top_terms_topics=train_model[[1]]$PW_T_sort,average=TRUE,smth=1),
                  mean_inf.rm(mapply(function(j)dist_from_corpus(PW_T=train_model[[1]]$PW_T[,j],dtm_data=dtm_data[test_iid,]),1:noc)),
                  train_model[[1]]$fc_validity[3],
                  metric_arun2010(t(train_model[[1]]$PW_T),pdT_test,apply(dtm_data[test_iid,],1,sum)),
                  compute_perplexity(dtm_data[test_iid,],train_model[[1]]$PW_T,pdT_test),
                  mean_inf.rm(mapply(function(j)topic_coherence_SR(train_model[[1]]$PW_T_sort[,j],dtm_data[test_iid,],colnames(dtm_data[test_iid,]),K=NROW(train_model[[1]]$PW_T_sort)),1:noc))
    ) #coherences, hellinger dist, silh, arun index, perplex, mimno2011
  }
  return(Cvout)
}

cv_evalTopics_lda = function(K=5,cvs=NULL,Xw=NULL,dtm_data=NULL,words_dic=NULL,noc=NULL){
  Cvout = matrix(nrow=12,ncol=K)
  for(k in 1:K){
    cat("\n CV..",k)
    train_iid = cvs$subsets[cvs$which==k]; test_iid = cvs$subsets[cvs$which!=k]
    fcm_test = Matrix::crossprod(Matrix::Matrix(dtm_data[test_iid,]),Matrix::Matrix(dtm_data[test_iid,]))
    train_model = topicmodels::LDA(slam::as.simple_triplet_matrix(dtm_data[train_iid,]),noc,method="Gibbs",control=list(burnin=500,iter=2500))
    pdT_test = topicmodels::LDA(slam::as.simple_triplet_matrix(dtm_data[test_iid,]),noc,method="Gibbs",control=list(burnin=500,iter=2500))@gamma
    Cvout[,k] = c(coherence(dtm_data=dtm_data[test_iid,],fcm_data=fcm_test,topic_words_dist=exp(train_model@beta),top_terms_topics=topicmodels::terms(train_model,30),average=TRUE,smth=1,probcoh = FALSE),
                  mean_inf.rm(mapply(function(j)dist_from_corpus(PW_T=t(exp(train_model@beta)),dtm_data=dtm_data[test_iid,]),1:noc)),
                  NA,
                  metric_arun2010(exp(train_model@beta),pdT_test,apply(dtm_data[test_iid,],1,sum)),
                  compute_perplexity(dtm_data[test_iid,],exp(t(train_model@beta)),pdT_test),
                  mean_inf.rm(mapply(function(j)topic_coherence_SR(topicmodels::terms(train_model,30),dtm_data[test_iid,],colnames(dtm_data[test_iid,]),K=30),1:noc))
    ) #coherences, hellinger dist, silh, arun index, perplex, mimno2011
  }
  return(Cvout)
}


cat("\n Running..")
load("casestudy_data.RData")

ntopics = 2:45
B = length(ntopics)
chunk.size = floor(B/ncores)
n_resid = rep(0,ncores)
if((ncores*chunk.size)<B){nresid = B - ncores*chunk.size;
n_resid = rmultinom(1, size = nresid, prob = rep(1/ncores,ncores))}

KK=25 #maximum steps of the outer loop (repetitions)
for(kk in 1:KK){ #repeated-Kfold-CV
cat("\n |------- KK = ",kk," -------|")
  
K=10; cvs = cvTools::cvFolds(NROW(dtm_matrix_filtered),K,type="random")[c(4,5)]

cat("\n FLSA..")
out <- foreach(h=1:ncores, .combine = "combAbind", .export = c(NULL,NULL)) %dopar% {
  res = array(NA,dim=c(12,K,chunk.size+n_resid[h]))
  iid_core = ((h-1)*chunk.size+1):((h*chunk.size)+n_resid[h])
  for(u in iid_core){
    res[,,u-(h-1)*chunk.size] = cv_evalTopics(K,cvs,X,dtm_matrix_filtered,words,ntopics[u])
  }
  res
}
names(out)=ntopics
save(out,file = paste0("casestudy_cv_results_flsa_",kk,".RData"))
rm(out)

cat("\n LDA..")
out <- foreach(h=1:ncores, .combine = "combAbind", .export = c(NULL,NULL)) %dopar% {
  res = array(NA,dim=c(12,K,chunk.size+n_resid[h]))
  iid_core = ((h-1)*chunk.size+1):((h*chunk.size)+n_resid[h])
  for(u in iid_core){
    res[,,u-(h-1)*chunk.size] = cv_evalTopics_lda(K,cvs,X,dtm_matrix_filtered,words,ntopics[u])
  }
  res
}
names(out)=ntopics
save(out,file = paste0("casestudy_cv_results_lda_",kk,".RData"))
rm(out)

cat("\n |------- -- - -- -------| \n")
}

doParallel::stopImplicitCluster(); parallel::stopCluster(cl)
cat("\n Done")
