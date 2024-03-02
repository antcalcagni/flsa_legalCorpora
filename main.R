
# Initial settings --------------------------------------------------------
rm(list=ls()); graphics.off(); options(warn = -1)
source("utils/utils.R"); source("utils/FLSA.R")
library(tm); library(quanteda); library(stringr)


# Data --------------------------------------------------------------------
txt_paths = tm::DirSource(directory = "data/out_corpora/",mode = "text")
corpus = tm::VCorpus(txt_paths); N=length(corpus)
corpus = paste(mapply(function(i)corpus[[i]]$content,1:N),collapse = " ")

# Dtm and Vocab with raw data
dtm = tm::DocumentTermMatrix(corpus); dtm_matrix = as.matrix(dtm)
out = apply(dtm_matrix,2,sum); 
Vocab = data.frame(term=names(out),freq=as.numeric(out),cums=cumsum(as.numeric(out))); Vocab = Vocab[order(Vocab$freq,decreasing = FALSE),]; head(Vocab)
Vocab$freq_rel = round(Vocab$freq/sum(Vocab$freq),4); Vocab$cums_rel = round(Vocab$cums/sum(Vocab$freq),4)
Vocab = Vocab[order(Vocab$freq,decreasing = TRUE),]; head(Vocab)
M_raw = NROW(Vocab)

# Save and manually mark words that need to be removed 
Vocab$keep = 1
write.csv(x = Vocab,file = "data/dictionary.csv",dec = ".")
X=read.csv(file = "data/dictionary.csv")[c(1,7)]

#x11(); plot(1:M,Vocab$freq,type="h",bty="n")
summary(Vocab$freq)
hpx = sum(Vocab$freq==1)/M_raw*100

# Removing words appearing rarely
length(Vocab$term[Vocab$freq<9])/M_raw
1-length(Vocab$term[Vocab$freq<9])/M_raw

# Filter-out the dtm matrix
jjd = c(X[!X[,2],1], #words marked manually
        as.numeric(which(apply(dtm_matrix,2,sum)<8)), #remove rare words
        grep(x = colnames(dtm_matrix),pattern = "per_cento"),
        grep(x = colnames(dtm_matrix),pattern = "anno_"),
        grep(x = colnames(dtm_matrix),pattern = "unità"),
        grep(x = colnames(dtm_matrix),pattern = "comma|commi|lett"),
        grep(x = colnames(dtm_matrix),pattern = "[0-9]{1,3}[a-z]{1}"), #posizioni economiche
        grep(x = colnames(dtm_matrix),pattern = "\\b[0-9]{4}-[0-9]{4}\\b") #years
)
length(jjd)/M_raw

dtm_matrix_filtered = dtm_matrix[,setdiff(1:M_raw,jjd)] #direction: keep
dtm_matrix_filtered = dtm_matrix_filtered[which(apply(dtm_matrix_filtered,1,sum)>0),] #keep documents with at least 1 word!
out = apply(dtm_matrix_filtered,2,sum); 
Vocab_filtered = data.frame(term=names(out),freq=as.numeric(out)); Vocab_filtered = Vocab_filtered[order(Vocab_filtered$freq,decreasing = TRUE),]; head(Vocab_filtered)
M = NROW(Vocab_filtered); 
summary(Vocab_filtered$freq)

corpus = c()
for(i in 1:NROW(dtm_filtered)){
  x=dtm_matrix_filtered[1,dtm_matrix_filtered[1,]>0]
  corpus=c(corpus,rep(names(x),x))
}
corpus = paste(corpus,collapse = " ")
dtm = tm::DocumentTermMatrix(corpus); dtm_matrix = as.matrix(dtm)

dtm_filtered = as.DocumentTermMatrix(dtm_matrix_filtered,weighting = tm::weightTf)
fcm_matrix_filtered = t(dtm_matrix_filtered)%*%dtm_matrix_filtered
X = lsa::lw_tf(m = dtm_matrix_filtered)*lsa::gw_idf(m = dtm_matrix_filtered) 
words = colnames(X)
save(X,dtm_matrix_filtered,dtm_filtered,fcm_matrix_filtered,words,file = "data/casestudy_data.RData")


# Repeated 10-fold CV -----------------------------------------------------
# The procedure has been run on a HPC remote system (see: 'casestudy_parallel_cv_all.R').
# The 25 x 10-fold CVs have been averaged over the outer 25 repetitions.
# The results are: 
#     - results/casestudy_cv_results_flsa.RData
#     - results/casestudy_cv_results_lda.RData

load("results/casestudy_cv_results_flsa.RData")
nocs = 2:45
K=dim(OUT_flsa)[2]
ntopics = dim(OUT_flsa)[3]+1
cvres_flsa = OUT_flsa

load("results/casestudy_cv_results_lda.RData")
cvres_lda = OUT_lda


# Running models on selected n topics -------------------------------------
# These results refer to the number of topics selected via cross-validation

# FLSA
noc=12
out_flsa = FLSA(X=X,dtm = dtm_matrix_filtered,fcm = fcm_matrix_filtered,words = words,noc = noc,fc_validity = FALSE,compute_stats = TRUE)
apply(out_flsa[[2]]$topic_coherence,2,mean)
write.csv(out_flsa[[1]]$Frex_W_T_sort,file = "results/fLSA_Frexs_K12.csv")

# LDA
noc=14; out_lda = vector("list",3)
out_lda[[3]] = topicmodels::LDA(slam::as.simple_triplet_matrix(dtm_matrix_filtered),noc,method="Gibbs",control=list(burnin=500,iter=2500))
out_lda[[1]]$PW_T = t(exp(out_lda[[3]]@beta)); out_lda[[1]]$PD_T = out_lda[[3]]@gamma
rownames(out_lda[[1]]$PW_T) = words
out_lda[[2]]$topic_coherence = coherence(dtm = dtm_matrix_filtered,fcm = fcm_matrix_filtered,topic_words_dist = t(out_lda[[1]]$PW_T),top_terms_topics = topicmodels::terms(out_lda[[3]],30))
doc_lengths = apply(dtm_matrix_filtered,1,sum); 
PT_unno = apply(out_lda[[1]]$PD_T*doc_lengths%*%matrix(1,1,noc),2,sum); out_lda[[1]]$PT = PT_unno/sum(PT_unno)
apply(out_lda[[2]]$topic_coherence,2,mean)

PW_T_frex = stm::calcfrex(logbeta = log(t(out_lda[[1]]$PW_T)),wordcounts = apply(dtm_matrix_filtered,2,sum)); rownames(PW_T_frex) = rownames(out_lda[[1]]$PW_T)
Frex_sorted = mapply(function(j)names(PW_T_frex[order(PW_T_frex[,j],decreasing = TRUE),j][1:30]),1:noc) #most frequent terms for each topic
colnames(Frex_sorted)=paste0("topic",1:noc); out_lda[[1]]$Frex_W_T_sort = Frex_sorted
write.csv(out_lda[[1]]$Frex_W_T_sort,file = "results/LDA_Frexs_K14.csv")



# Paper materials ---------------------------------------------------------

cls = c("#36648B", #steelblue4 
        "#CD4F39") #tomato3 (red-like)


## Figure 1
tikzDevice::tikz(file='../contribution/fig1.tex',width=7.5,height=2)
par(mfrow=c(1,4))
#Panels A-B
m=1; lwdx=2
Z2 = apply(OUT_lda[m,,,],c(2),function(x)c(quantile(x,c(0.25,0.75),na.rm=TRUE),mean(x,na.rm=TRUE)))
Z1 = apply(OUT_flsa[m,,,],c(2),function(x)c(quantile(x,c(0.25,0.75)),mean(x)))
#A
plot(nocs,Z1[3,],bty="n",type="l",lwd=lwdx,ylim = c(min(Z1,Z2),max(Z1,Z2)),main="UMass-like",xlab="",ylab="",col=cls[1]); 
polygon(c(nocs, rev(nocs)), c(Z1[2,], rev(Z1[1,])),col = "#EEE9E9",border = FALSE); lines(nocs,Z1[3,],lwd=lwdx,col=cls[1])
px1=kneePoint(x = nocs,y = Z1[3,],plot = FALSE,df = 1,bty="n",sign = -1,xQuery = nocs)
abline(v = c(px1),col=cls[1],lty=2,lwd=1.25); text(c(px1)+2,c(Z1[3,1]),labels = c(px1),cex = 1.25,col=cls[1])
title("(A)", line = 1.5,adj=0); 
#B
plot(nocs,Z2[3,],bty="n",type="l",lwd=lwdx,ylim = c(min(Z1,Z2),max(Z1,Z2)),main="UMass-like",xlab="",ylab="",col=cls[2]); 
polygon(c(nocs, rev(nocs)), c(Z2[2,], rev(Z2[1,])),col = "#EEE9E9",border = FALSE); lines(nocs,Z2[3,],lwd=lwdx,col=cls[2])
px1=kneePoint(x = nocs,y = Z2[3,],plot = FALSE,df = 1,bty="n",sign = -1,xQuery = nocs)
abline(v = c(px1),col=cls[2],lty=2,lwd=1.25); text(c(px1+2),c(Z2[3,1]),labels = c(1+px1),cex = 1.25,col=cls[2])
title("(B)", line = 1.5,adj=0); 

m=2; lwdx=2
Z2 = apply(OUT_lda[m,,,],c(2),function(x)c(quantile(x,c(0.25,0.75),na.rm=TRUE),mean(x,na.rm=TRUE)))
Z1 = apply(OUT_flsa[m,,,],c(2),function(x)c(quantile(x,c(0.25,0.75)),mean(x)))
#C
plot(nocs,Z1[3,],bty="n",type="l",lwd=lwdx,ylim = c(min(Z1,Z2),max(Z1,Z2)),main="UCI-like",xlab="",ylab="",col=cls[1]); 
polygon(c(nocs, rev(nocs)), c(Z1[2,], rev(Z1[1,])),col = "#EEE9E9",border = FALSE); lines(nocs,Z1[3,],lwd=lwdx,col=cls[1])
px1=kneePoint(x = nocs,y = Z1[3,],plot = FALSE,df = 1,bty="n",sign = 1,xQuery = nocs)
abline(v = c(px1),col=cls[1],lty=2,lwd=1.25); text(c(px1+2),c(Z1[3,1]),labels = c(px1-2),cex = 1.25,col=cls[1])
title("(C)", line = 1.5,adj=0); 
#D
plot(nocs,Z2[3,],bty="n",type="l",lwd=lwdx,ylim = c(min(Z1,Z2),max(Z1,Z2)),main="UCI-like",xlab="",ylab="",col=cls[2]); 
polygon(c(nocs, rev(nocs)), c(Z2[2,], rev(Z2[1,])),col = "#EEE9E9",border = FALSE); lines(nocs,Z2[3,],lwd=lwdx,col=cls[2])
px1=kneePoint(x = nocs,y = Z2[3,],plot = FALSE,df = 1,bty="n",sign = 1,xQuery = nocs)
abline(v = c(px1),col=cls[2],lty=2,lwd=1.25); text(c(px1+2),c(Z2[3,1]),labels = c(2+px1),cex = 1.25,col=cls[2])
title("(D)", line = 1.5,adj=0); 

add_legend("bottom",fill = cls,legend = c("flsa","lda"),border = FALSE,bty = "n",ncol = 2,cex=1.5)
dev.off()


## Figure 2
tikzDevice::tikz(file='fig2.tex',width=5.5,height=2.5)
cls2=paletteer::paletteer_c("grDevices::Heat",n=25,direction = -1) #https://r-charts.com/color-palettes/#dynamic
x_cls2=seq(from=0,to=1,length=length(cls2)) 
par(mfrow=c(1,2))
#Panel A
noc=12
pt=out_flsa$FLSA$PT/max(out_flsa$FLSA$PT)
exc=round(exclusive_term_ratio(out_flsa)$normalized,2); exc=exc/max(exc); iid = sapply(exc,function(x)which.min(abs(x-x_cls2)))
plot(pt,ylab="",xlab="",bty="n",type="h",lwd=5,col=cls2[iid],axes = FALSE,xlim=c(1,noc),ylim=c(0,1)); 
axis(side = 1,at = 1:noc,labels = 1:noc);axis(side = 2,at = round(seq(from=0,to=1,length=noc),2));
title("(A) fLSA", line = 1,adj=0)
#Panel B
noc=14
pt=out_lda[[1]]$PT
exc=round(exclusive_term_ratio(out_lda)$normalized,2); exc=exc/max(exc); iid = sapply(exc,function(x)which.min(abs(x-x_cls2)))
plot(1:noc,pt/max(pt),col=cls2[iid],bty="n",ylab="",xlab="",type="h",lwd=5,axes = FALSE,ylim=c(0,1)); 
title("(B) LDA", line = 1,adj=0)
axis(side = 1,at = 1:noc,labels = 1:noc); axis(side = 2,at = round(seq(from=0,to=1,length=noc),2));
add_legend("bottom",legend = as.numeric(summary(x_cls2)[-4]),fill = cls2[as.numeric(sapply(summary(x_cls2)[-4],function(x)which(x==x_cls2)))],border = FALSE,bty = "n",cex=1.25,ncol=5)
dev.off()


## Figure 3 
tikzDevice::tikz(file='fig3.tex',width=8.5,height=4.5)
par(mfrow=c(1,2))
#Panel A
distx = text2vec::dist2(x = t(out_flsa$FLSA$PW_T),method="cosine") #higher similarity scores indicate greater redundancy between topics
hc=hclust(d = as.dist(distx),method = "ward.D2")
ngx=which.max(sapply(2:NCOL(distx),function(x)clusterSim::index.S(d = as.dist(distx),cl = cutree(hc,k = x))))+1; gx=cutree(hc,ngx)
intertopic_distanceMap(FLSA_out = out_flsa,doc_lengths = doc_lengths,plotx = TRUE,lbls_cex = 1.5,log_adj = 2.5,new_window = FALSE,gps=cutree(hc,ngx),k_rect=1.5,
                       coh_measure = "mean_pmi",cols = paletteer::paletteer_c("ggthemes::Orange Light", 10,1),
                       cols_gps = rep("#133955",length(gx))); title("(A) fLSA", line = 1,adj=0)
#Panel B
intertopic_distanceMap(FLSA_out = out_lda,doc_lengths = doc_lengths,plotx = TRUE,lbls_cex = 1.5,log_adj = 2.5,new_window = FALSE,gps=NULL,
                       coh_measure = "mean_pmi",cols = paletteer::paletteer_c("ggthemes::Orange Light", 10,1),
                       cols_gps = rep("#133955",length(gx))); title("(B) LDA", line = 1,adj=0)
dev.off()


## Now, run FLSA on a smaller number of topics (n=3)
## This is needed to create Figure 4
noc=3
out_flsa = FLSA(X=X,dtm = dtm_matrix_filtered,fcm = fcm_matrix_filtered,words = words,noc = noc,fc_validity = FALSE,compute_stats = TRUE)

out_flsa[[1]]$Frex_W_T_sort
out_flsa[[1]]$PT
exclusive_term_ratio(out_flsa)
apply(out_flsa[[2]]$topic_coherence,2,mean)

# Frexs and topic proportion
Pwts = out_flsa[[1]]$PW_T; Pwts_n = t(apply(Pwts,1,function(x)x/sum(x)))
Wds = out_flsa[[1]]$Frex_W_T_sort
Iid = matrix(as.numeric(sapply(as.vector(Wds),function(x)which(x==rownames(Pwts)))),ncol = NCOL(Wds))
rownames(Pwts_n)[Iid[10,2]]="politiche_agric_alim" #shorten a particular long word for the sake of graph representation
write.csv(out_flsa[[1]]$Frex_W_T_sort,file = "results/fLSA_Frexs_K3.csv")

## Figure 4
tikzDevice::tikz(file='fig4.tex',width=7,height=3,sanitize = TRUE)
par(mar=c(5,8,5,2)+0.1,mfrow=c(1,3)) # Doubles left margin.
j=1;barplot(t(round(Pwts_n[Iid[1:10,j],],3)),border = FALSE,beside = FALSE,horiz = TRUE,las=2,col = c("#EEB422","#00B2EE","#B4EEB4"),cex.names = 1.35); title(paste0("Topic ",j), line = 1,adj=0,cex=1.35)
j=2;barplot(t(round(Pwts_n[Iid[1:10,j],],3)),border = FALSE,beside = FALSE,horiz = TRUE,las=2,col = c("#EEB422","#00B2EE","#B4EEB4"),cex.names = 1.35); title(paste0("Topic ",j), line = 1,adj=0,cex=1.35)
j=3;barplot(t(round(Pwts_n[Iid[1:10,j],],3)),border = FALSE,beside = FALSE,horiz = TRUE,las=2,col = c("#EEB422","#00B2EE","#B4EEB4"),cex.names = 1.35); title(paste0("Topic ",j), line = 1,adj=0,cex=1.35)
add_legend("bottom",legend = c("Topic 1","Topic 2","Topic 3"),fill = c("#EEB422","#00B2EE","#B4EEB4"),border = FALSE,bty = "n",cex=1.5,ncol=3)
dev.off()



# Supplementary Materials -------------------------------------------------

## Figure 1-sup
tikzDevice::tikz(file='fig1_supp.tex',width=6.5,height=4)
par(mfcol=c(2,3))
lwdx=2; mains = c("Mean-diff","Hellinger dist", "Arun 2010")
for(m in 3:5){
  
  Z2 = apply(OUT_lda[m,,,],c(2),function(x)c(quantile(x,c(0.25,0.75),na.rm=TRUE),mean(x,na.rm=TRUE)))
  Z1 = apply(OUT_flsa[m,,,],c(2),function(x)c(quantile(x,c(0.25,0.75)),mean(x)))
  #A
  plot(nocs,Z1[3,],bty="n",type="l",lwd=lwdx,ylim = c(min(Z1,Z2),max(Z1,Z2)),main=mains[m-2],xlab="",ylab="",col=cls[1]); 
  polygon(c(nocs, rev(nocs)), c(Z1[2,], rev(Z1[1,])),col = "#EEE9E9",border = FALSE); lines(nocs,Z1[3,],lwd=lwdx,col=cls[1])
  if(m<5){px1=kneePoint(x = nocs,y = Z1[3,],plot = FALSE,df = 1,bty="n",sign = 1,xQuery = nocs)
  abline(v = c(px1),col=cls[1],lty=2,lwd=1.25); text(c(px1)+2,c(Z1[3,1]),labels = c(px1),cex = 1.25,col=cls[1])}
  title(paste0("(",LETTERS[m-2],")"), line = 1.5,adj=0); 
  #B
  plot(nocs,Z2[3,],bty="n",type="l",lwd=lwdx,ylim = c(min(Z1,Z2),max(Z1,Z2)),main=mains[m-2],xlab="",ylab="",col=cls[2]); 
  polygon(c(nocs, rev(nocs)), c(Z2[2,], rev(Z2[1,])),col = "#EEE9E9",border = FALSE); lines(nocs,Z2[3,],lwd=lwdx,col=cls[2])
  if(m<5){px1=kneePoint(x = nocs,y = Z2[3,],plot = FALSE,df = 1,bty="n",sign = 1,xQuery = nocs)
  abline(v = c(px1),col=cls[2],lty=2,lwd=1.25); text(c(px1+2),c(Z2[3,1]),labels = c(1+px1),cex = 1.25,col=cls[2])}
  title(paste0("(",LETTERS[m-2],")"), line = 1.5,adj=0); 
}
add_legend("bottom",fill = cls,legend = c("flsa","lda"),border = FALSE,bty = "n",ncol = 2,cex=1.5)
dev.off()


## Figure 2-sup
tikzDevice::tikz(file='fig2_supp.tex',width=4.5,height=6)
par(mfrow=c(3,2))

mains = c("UMass-like","UCI-like","Mean-diff","Hellinger dist", "Arun 2010")
for(m in 1:5){
  Z1 = t(apply(OUT_flsa[m,,,],c(1,3),mean,na.rm=TRUE))
  Z2 = t(apply(OUT_lda[m,,,],c(1,3),mean,na.rm=TRUE))
  boxplot(Z1,at=seq(1,20,by=2)-1,frame=FALSE,col=cls[1],ylim=c(min(Z1,Z2),max(Z1,Z2)),xlim=c(0,22),names=rep("",10),xlab="",main=mains[m])
  boxplot(Z2,at=seq(1,20,by=2),frame=FALSE,add=TRUE,col=cls[2]); 
  title(paste0("(",LETTERS[m],")"), line = 1.5,adj=0); 
}
add_legend("bottom",fill = cls,legend = c("flsa","lda"),border = FALSE,bty = "n",ncol = 2,cex=1.5)
dev.off()




