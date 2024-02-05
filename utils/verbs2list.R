verbs2list = function(verb="essere"){
  vbs_tree = list(c("presente","imperfetto","passato remoto","futuro semplice","passato prossimo","trapassato prossimo","trapassato remoto","futuro anteriore"),
                  c("presente","passato"),c("presente","passato","imperfetto","trapassato"),c("presente"),c("presente","passato"),c("presente","passato"),c("presente","passato")); names(vbs_tree) = c("indicativo","condizionale","congiuntivo","imperativo","infinito","participio","gerundio")
  tmps = str_unique(as.vector(unlist(vbs_tree)))                
                  
  vbs=c(verb,paste(verb,"femminile",sep = "-"))
  out_list=c()
  for(i in 1:2){
    rm(html,html_raw,out)
    html_raw = tolower(readLines(paste0("https://www.coniugazione.it/verbo/",vbs[i],".php")))
    html = html_raw[grep(x = html_raw,pattern = '<div id=\"contenu\"'):length(html_raw)]
    iid_mode = grep(x = html,pattern = '<h2 class=\"mode\">')[1:7] #Note: length(iid_mode) = length(vbs_tree)
    
    for(m in 2:(length(iid_mode))){  #from m=2..m=8
      x_html = gsub(x = html[iid_mode[m-1]:(iid_mode[m])],pattern = "<b>|</b>",replacement = "") 
      x_html = gsub(x = x_html,pattern = "<.*?>",replacement = " ") #current mode
      x_html = str_squish(x_html); x_html = x_html[!is.na(x_html)]
      for(j in 1:7){x_html=gsub(x = x_html,pattern = names(vbs_tree)[j],replacement = " ")} #x_html = gsub(x = x_html,pattern = names(vbs_tree)[m-1],replacement = " ")
      x_html = gsub(x = x_html,pattern = "-|;|,",replacement = " ")
      x_html = gsub(x = x_html,pattern = "&agrave",replacement = "à")
      x_html = gsub(x = x_html,pattern = "&ograve",replacement = "ò")
      x_html = gsub(x = x_html,pattern = "&igrave",replacement = "ì")
      x_html = gsub(x = x_html,pattern = "&egrave",replacement = "è")
      x_html = paste(x_html,collapse = " ")
      x_html = str_squish(x_html)
      
      out = x_html; for(t in 1:length(tmps)){out=str_replace_all(string = out,pattern = tmps[t],replacement = "XXA")} #out = x_html; for(t in 1:length(vbs_tree[[m-1]])){out=str_replace_all(string = out,pattern = vbs_tree[[m-1]][t],replacement = "XXA")}
      out = unlist(str_split(string = out,pattern = "XXA"))
      out = str_replace_all(string = out,pattern = "\\bio\\s|\\btu\\s|\\blui\\s|\\blei\\s|\\bnoi\\s|\\bvoi\\s|\\bloro\\s|\\bessi\\s|\\besse\\s|\\besso\\s",replacement = " XXA ")
      out = str_replace_all(string = out,pattern = "\\bche\\b",replacement = " ")
      out = str_squish(unlist(str_split(string = out,pattern = "XXA")))
      out = out[str_length(out)>0]
      out_list = c(out_list,out)
    }
  }
  
  iid_trm = which(str_detect(string = out_list,pattern = "/"))
  if(length(iid_trm)>0){
    for(j in 1:length(iid_trm)){ 
      x=unlist(str_split(string = out_list[iid_trm[j]],pattern = "/"))
      out_list[iid_trm[j]]=x[1]
      #out_list = c(out_list,x[-1])
    }
  }
  
  iid_trm = which(lapply(str_split(string = out_list,pattern = " "),length)>2)
  if(length(iid_trm)>0){
    x = unlist(str_split(string = out_list[iid_trm]," "))
    out_list = out_list[-iid_trm]
    out_list = c(out_list,x)
  }
  
  return(str_unique(out_list))
  #return(out_list)
}
