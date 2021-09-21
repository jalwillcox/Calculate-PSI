library(ggplot2)

args = commandArgs(trailingOnly=TRUE)

# Read in PSI data

psi <- read.delim(file=args[1], stringsAsFactors=F)

psi$transcript <- as.character(sapply(psi$exon_ID, function(x) unlist(strsplit(split=";", x))[2]))

# Only use the longest transcript
lt <- tail(names(sort(table(psi$transcript))),1)
psi <- psi[psi$transcript == lt,]

psi$exon_num <- 1:nrow(psi)

# Format for plotting

pdf(args[2],width=14,height=7)
ggplot(psi) +
  geom_area(aes(exon_num,ave),fill="deepskyblue") +
  xlab("Exon") +
  ylab("PSI") +
  scale_y_continuous(breaks=seq(0,1,0.1),limits=c(0,1.1)) +
  theme(
    panel.background = element_rect(fill = "white", color = "black"),
    panel.grid.major = element_line(colour = "light gray"),
    panel.grid.minor = element_line(colour = "light gray"),
  )
dev.off()


