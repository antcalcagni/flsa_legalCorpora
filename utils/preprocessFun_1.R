preprocessFun_1 = function(n=NULL,txt=NULL,month_names=NULL,abbjus=NULL,stpwds=NULL,polirems=NULL){

  for(i in 1:length(month_names)){
    x=str_extract_all(string = txt,pattern = paste0("[0-9]{1,}\\s",month_names[i]))[[1]]  
    if(length(x)>0){
      y=gsub(x = x,pattern = month_names[i],replacement = as.character(i))
      for(j in 1:length(x)){txt=gsub(x = txt,pattern = x[j],replacement = y[j])}
    }
  }
  
  txt = gsub(x = txt,pattern = "[,.:;°]",replacement = "")  #remove basic punctuaction
  txt = gsub(x = txt,pattern = "[«»'/]",replacement = " ") #remove other symbols 
  txt = str_replace_all(string = txt,pattern = "\\([a-z]{1,}\\)",replacement = "") #remove strings in parentheses
  txt = str_replace_all(string = txt,pattern = "\\([0-9]{1,}\\)",replacement = "") #remove numbers in parentheses --ie, comma
  
  for(b in 1:NROW(abbjus)){
    txt = extract_replace_abbjus(txt,abbjus[b,1],abbjus[b,2],p1) #case 1 (no abbr)
    txt = extract_replace_abbjus(txt,abbjus[b,1],abbjus[b,2],p2) #case 2
    txt = extract_replace_abbjus(txt,abbjus[b,1],abbjus[b,2],p3) #case 3
    txt = str_replace_all(string = txt,pattern = str_squish(abbjus[b,1]),replacement = str_squish(abbjus[b,2]))
  }
  txt = str_squish(txt);
  
  # euro, mila euro, milioni di euro
  x=str_extract_all(string = txt,pattern = "(euro|mila euro|milioni di euro|milione di euro)")[[1]]
  if(length(x)>0){for(i in 1:length(x)){txt = str_replace(string = txt,pattern = x[i],replacement = gsub(x = x[i],pattern = " ",replacement = "_"))}}
  
  # anno + yyyy
  x=str_extract_all(string = txt,pattern = "anno\\s[0-9]{4,}")[[1]]
  if(length(x)>0){for(i in 1:length(x)){txt = str_replace(string = txt,pattern = x[i],replacement = gsub(x = x[i],pattern = " ",replacement = "_"))}}
  
  # num + 'per cento'
  x=str_extract_all(string = txt,pattern = "[0-9]{1,}\\sper\\scento")[[1]]
  if(length(x)>0){for(i in 1:length(x)){txt = str_replace(string = txt,pattern = x[i],replacement = gsub(x = x[i],pattern = " ",replacement = "_"))}}
  
  # num + unità
  x=str_extract_all(string = txt,pattern = "[0-9]{1,}\\sunità")[[1]]
  if(length(x)>0){for(i in 1:length(x)){txt = str_replace(string = txt,pattern = x[i],replacement = gsub(x = x[i],pattern = " ",replacement = "_"))}}
  
  # posizione_economica + codice (eg, posizione economica f1)
  x=str_extract_all(string = txt,pattern = "(posizione_economica)\\s[a-z]{1,}")[[1]]
  if(length(x)>0){for(i in 1:length(x)){txt = str_replace(string = txt,pattern = x[i],replacement = gsub(x = x[i],pattern = " ",replacement = "_"))}}
  
  # Parole tecniche/polirematiche
  for(p in 1:NROW(polirems)){
    x=str_squish(str_extract_all(string = txt,pattern = paste0("\\s",tolower(as.character(polirems[p,1])),"\\s"))[[1]])
    if(length(x)>0){for(i in 1:length(x)){txt = str_replace(string = txt,pattern = x[i],replacement = gsub(x = x[i],pattern = " ",replacement = "_"))}}}
  
  # Common Italian stopwords
  for(i in 1:NROW(stpwds)){txt = str_replace_all(string = txt,pattern = paste0("\\s",as.character(stpwds[i,1]),"\\s"),replacement = " ")}
  
  # Remove links to commas in the same budget law 
  txt=str_replace_all(string = txt,pattern = "(comma|art|artt|articolo)\\s[0-9]{1,}",replacement = "") 
  txt=str_replace_all(string = txt,pattern = "-[a-z]{3,}",replacement = "") # remove 1-bis, 2-ter, 3-quinquies
  txt=str_replace_all(string = txt,pattern = "(lettera|lett)\\s[a-z]{1}\\)",replacement = "")
  
  # Remove floating symbols (if any)
  txt=str_replace_all(string = txt,pattern = '\\"',replacement = "")
  txt=str_replace_all(string = txt,pattern = '\\(',replacement = "")
  txt=str_replace_all(string = txt,pattern = '\\)',replacement = "")
  
  # Remove floating numbers and words with length less K characters
  x=str_split(string = txt,pattern = " ")[[1]][-1] #-1 remove the number in front of the text
  txt=paste(x[nchar(x)>3 & is.na(as.numeric(x))],collapse = " ")
  
  #txts_processed[n] = str_squish(txt)
  write(x = str_squish(txt),file = paste0(out_name,"/txt_",n,".txt"))
}
