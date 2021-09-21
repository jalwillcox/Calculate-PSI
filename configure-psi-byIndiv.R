args = commandArgs(trailingOnly=TRUE)

psi <- read.delim(file=args[1],stringsAsFactors=F)

psi$ave <- rowMeans(psi[,-1],na.rm=T)

write.table(psi, file=args[1], col.names=T, row.names=F, sep="\t", quote=F)

