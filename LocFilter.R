# LOC filter function
library(dplyr)
LocFilter <- function(list.txt, out.txt){
  list <- read.delim(list.txt, header = FALSE)
  vector <- unlist(list)
  condition <- !grepl('LOC', vector)
  print(condition)
  out <- vector[condition]
  print(out)
  write.table(out, file = out.txt, row.names = FALSE, col.names = FALSE, quote = FALSE)
  return(out)
}


inFile <- '~/ConGen/DAVID/DAVID_horse_rhino/H_R_logFC_P_cut/Newtop100.txt'
noLocHRtop <- LocFilter(inFile, 'LOCFILTEROUT.txt')
