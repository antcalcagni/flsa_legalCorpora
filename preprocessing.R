
# Initial settings --------------------------------------------------------
rm(list=ls()); graphics.off()
source("utils/utils.R")
library(tm); library(quanteda); library(stringr)

options(warn = -1)

## Function used to handle with {legge, decreto legge,..}
extract_replace_abbjus = function(txt=NULL,current_abbjus=NULL,new_abbjus=NULL,regex_pattern=NULL){
  # search for plain abbjus (no abbreviation)
  x = str_extract_all(string = txt,pattern = paste0("\\b",str_squish(current_abbjus),regex_pattern))[[1]] 
  if(length(x)>0){
    for(i in 1:length(x)){
      z = gsub(x = gsub(x = x[i],pattern = str_squish(current_abbjus),replacement = str_squish(new_abbjus)),pattern = " ",replacement = "_")
      txt = str_replace(string = txt,pattern = x[i],replacement = z)
    }
  }
  # search for abbreviated abbjus (if any)
  x = str_extract_all(string = txt,pattern = paste0("\\b",str_squish(new_abbjus),"\\b"))[[1]] 
  if(length(x)>0){
    for(i in 1:length(x)){
      z = gsub(x = gsub(x = x[i],pattern = str_squish(new_abbjus),replacement = str_squish(new_abbjus)),pattern = " ",replacement = "_")
      txt = str_replace(string = txt,pattern = x[i],replacement = z)
    }
  }
  return(txt)
}

## Regex expressions useful for detecting {legge,decreto legge,..}
p1 = "\\s[0-9]{1,}\\s[0-9]{1,}\\s[0-9]{4}\\sn\\s[0-9]{1,}" #case 1 (eg, dpcm 23 10 2020 n 101)
p2 = "\\s[0-9]{1,}\\s[0-9]{1,}\\s[0-9]{4}" #case 2 (eg, dpcm 23 10 2020)
p3 = "(\\sn\\s)*[0-9]{1,}\\sdel\\s[0-9]{4}" #case 3 (eg, n 101 del 2020, 101 del 2020)


# Data --------------------------------------------------------------------
month_names = c("gennaio","febbraio","marzo","aprile","maggio","giugno","luglio","agosto","settembre","ottobre","novembre","dicembre")
out_name = "data/out_corpora"; if(!dir.exists(out_name)){dir.create(out_name)} #output folder
txts = read.csv(file = "data/originale_testo legge 178_2020_bis.csv",header = FALSE,encoding = "UTF-8") #read raw corpus
abbjus = read.csv(file = "data/abbreviazioni_giuridiche.csv",header = FALSE); #read juridical terms
stpwds = read.delim(file = "data/stopwords-it_integrate.txt",header = FALSE,sep = "\n",row.names = NULL,stringsAsFactors = FALSE,fileEncoding = "ISO-8859-1");
#polirems = read.delim(file = "data/polirematiche_it.txt",header = FALSE,sep = "\n",row.names = NULL,stringsAsFactors = FALSE);
polirems = read.delim(file = "data/parole_tecniche-it.txt",header = FALSE,sep = "\n",row.names = NULL,stringsAsFactors = FALSE);


# Pre-processing 1 --------------------------------------------------------
library(doParallel)
source("utils/preprocessFun_1.R")

ncores=8
cl = parallel::makeCluster(ncores)
doParallel::registerDoParallel(cl,cores=ncores)
B = NROW(txts)
chunk.size = floor(B/ncores)
if((ncores*chunk.size)<B){nresid = B - ncores*chunk.size}
n_resid = rep(0,ncores)
set.seed(111);n_resid = rmultinom(1, size = nresid, prob = rep(1/ncores,ncores))

foreach(h=1:ncores, .export = c(NULL,NULL),.packages = c("stringr")) %dopar% {
  iid_core = ((h-1)*chunk.size+1):((h*chunk.size)+n_resid[h])
  for(u in iid_core){
    preprocessFun_1(n=u,txt=str_squish(tolower(txts[u,1])),month_names,abbjus,stpwds,polirems)
  }
}



# Pre-processing 2 --------------------------------------------------------
rm(list=ls()); graphics.off()
out_name = "data/out_ngrams"; if(!dir.exists(out_name)){dir.create(out_name)} #output folder
library(ngram)

score_ngram = function(Vocab=NULL,ngram=NULL,f_ngram=NULL){
  unigram = str_split(string = str_squish(ngram),pattern = " ")[[1]]
  ni = mapply(function(j)Vocab$freq[Vocab$term==unigram[j]],1:length(unigram))
  jjd = unlist(lapply(ni,length))
  ni = unlist(ni)*jjd
  fi = (ni+1e-02)/sum(Vocab$freq)

  pmi = log2(f_ngram/prod(fi))
  npmi = pmi/-sum(log2(fi))

  return(c(pmi,npmi,round(fi,4)))
}
#Note: The larger the PMI score, the larger P(w1,..,wn) if compared to P(w1) x..x P(wn)
#      If P(w1,..,wn) is higher then P(w1),..,P(wn) this means that the n-gram occurrs more frequently then the single words upon which is based

txt_paths = tm::DirSource(directory = "data/out_corpora/",mode = "text")
corpus = tm::VCorpus(txt_paths); N=length(corpus) #no need to load the corpus K times (see line 125 of this code!)

#### Looking for K-grams ####
## Note:  1) K-grams need to be found by varying (manually) the variable K (eg, K in {5,4,3,2})
##        2) Avoid creating ngrams containing verbs!
K=2

dtm = tm::DocumentTermMatrix(corpus); dtm_matrix = as.matrix(dtm)
out = apply(dtm_matrix,2,sum); 
Vocab = data.frame(term=names(out),freq=as.numeric(out)); Vocab = Vocab[order(Vocab$freq,decreasing = TRUE),]; head(Vocab)
corpus_all = paste(mapply(function(i)corpus[[i]]$content,1:N),collapse = " ")

ngs = ngram(corpus_all,n=K)
ng_out = get.phrasetable(ngs)
summary(ng_out$freq)
100*sum(ng_out$freq>9)/NROW(ng_out)
ng_out = ng_out[ng_out$freq>9,] #discard those ngrams occurring rarely
out = t(mapply(function(i){score_ngram(Vocab = Vocab,ngram = ng_out[i,1],f_ngram = ng_out[i,3])},i=1:NROW(ng_out)))
ng_out$pmi = out[,1]; ng_out$npmi = out[,2]; ng_out$ni = out[,3:NCOL(out)];

Xh=ng_out[order(ng_out$pmi,decreasing = TRUE),]; rownames(Xh)=1:NROW(Xh); Xh$keep = 0
#tail(ng_out[order(ng_out$pmi,decreasing = TRUE),],20) 
write.csv(x = Xh,file = paste0(out_name,"/output_",K,"grams.csv"),row.names = FALSE,na = "")
    ## Note: Now K-grams can be selected by adding 1/0 on the column 'keep' in the output csv file ##
Xh=read.csv(file = paste0(out_name,"/output_",K,"grams.csv"),header = TRUE)

jjd = which(Xh$keep==1)
ngrams_kept = str_squish(Xh[jjd,1]); print(ngrams_kept)
ngrams_kept_iid = mapply(function(j)which(unlist(mapply(function(i)str_detect(pattern = ngrams_kept[j],string = corpus[[i]]$content)[[1]],1:N))),1:length(jjd))
sum(unlist(lapply(ngrams_kept_iid,function(x)length(x)>0)))==length(ngrams_kept) #check for empty docs

ll = unlist(lapply(ngrams_kept_iid,length)); ll_max = max(ll)
X = matrix(data = NA,nrow = length(ngrams_kept),ncol = ll_max+1); 
for(j in 1:length(jjd)){X[j,] = c(ngrams_kept[j],c(ngrams_kept_iid[[j]],rep(NA,ll_max-ll[j])))}
write.csv(x = X,file = paste0(out_name,"/matrix_",K,"grams_docs.csv"),row.names = FALSE,na = "")

for(j in 1:length(ngrams_kept_iid)){
  for(n in 1:length(ngrams_kept_iid[[j]])){
    txt=str_replace_all(string = corpus[[ngrams_kept_iid[[j]][n]]]$content,pattern = ngrams_kept[j],replacement = gsub(ngrams_kept[j],pattern=" ",replacement="_"))
    corpus[[ngrams_kept_iid[[j]][n]]]$content = txt 
    write(x = str_squish(txt),file = paste0("data/out_corpora/txt_",ngrams_kept_iid[[j]][n],".txt"))
  }
}



# Pre-processing 3 --------------------------------------------------------
#library(doParallel)
rm(list=ls()); graphics.off()
source("utils/verbs2list.R")

replace_text = function(n=NULL,txt=NULL,list_words=NULL,M=NULL){
  for(j in 1:M){txt = str_replace_all(string = txt,pattern = paste0("\\s",as.character(list_words[j]),"\\s"),replacement = " ")}
  write(x = str_squish(txt),file = paste0("data/out_corpora/txt_",n,".txt"))    
}

vbs = c("essere","avere","dovere","potere","fare")
vbs_list = unlist(lapply(vbs,verbs2list))

txt_paths = tm::DirSource(directory = "data/out_corpora/",mode = "text")
corpus = tm::VCorpus(txt_paths); N=length(corpus); N = length(corpus)

for(u in 1:N){
  replace_text(n = u,txt = corpus[[u]]$content,list_words = vbs_list,M = length(vbs_list))  
}





