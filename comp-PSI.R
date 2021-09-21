

args = commandArgs(trailingOnly=TRUE)

case <- read.delim(file=args[1],stringsAsFactors=F)
con <- read.delim(file=args[2],stringsAsFactors=F)

comp <- data.frame(stringsAsFactors=F, "exon_ID"=case$exon_ID,"case.PSI"=case$ave, "con.PSI"=con$ave)
comp$exon <- as.character(sapply(comp$exon_ID, function(x) unlist(strsplit(split=";", x))[1]))
comp$transcript <- as.character(sapply(comp$exon_ID, function(x) unlist(strsplit(split=";", x))[2]))

comp$p.value.uncorrected <- as.numeric(sapply(1:nrow(comp), function(x) tryCatch({t.test(case[x,2:(ncol(case)-1)],con[x,2:(ncol(con)-1)])$p.value},error = {function(e) return(NA)})))

write.table(comp, file=args[3], col.names=T, row.names=F, sep="\t", quote=F)

