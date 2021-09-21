library(ggplot2)

args = commandArgs(trailingOnly=TRUE)

# Read in PSI data

psi <- read.delim(file=args[1], stringsAsFactors=F)

# Only use the longest transcript
lt <- tail(names(sort(table(psi$transcript))),1)
psi <- psi[psi$transcript == lt,]

psi$exon_num <- 1:nrow(psi)


# Format for plotting

psi_plot <- psi[,c("exon_num","con.PSI","p.value.uncorrected")]
colnames(psi_plot)=c("Exon","PSI","p.value") 
psi_plot$case.con="control"

psi_plot <- rbind(psi_plot,data.frame(stringsAsFactors=F, "Exon"=psi$exon_num, "PSI"=psi$case.PSI, "p.value"=psi$p.value.uncorrected, "case.con"="case"))

psi_plot$case.con <- factor(psi_plot$case.con, levels=c("control","case"))
psi_plot$log.p.value <- -log10(psi_plot$p.value)

psi_plot$p.mark <- 1.04

cutoff=-log10(0.05/nrow(psi))
cols=c("case"="black", "control"="deepskyblue")

pdf(args[2],width=14,height=7)
ggplot(psi_plot) +
  geom_area(data=psi_plot[psi_plot$case.con=="control",],aes(Exon,PSI,fill=case.con)) +
  scale_fill_manual(values=cols, guide=F) +
  geom_line(aes(Exon,PSI,color=case.con)) +
  scale_color_manual(values=cols) +
  geom_point(data=psi_plot[psi_plot$case.con=="case" & psi_plot$log.p.value > cutoff,],
             aes(x=Exon,y=p.mark,size=log.p.value),shape=25,fill="goldenrod1") +
  geom_label(data=psi_plot[psi_plot$case.con=="case" & psi_plot$log.p.value > cutoff,],
            aes(Exon,p.mark,label=Exon),hjust=0.5, vjust=-0.5) +
  guides(size = guide_legend(expression(-log[10](p-value))),color = guide_legend("")) +
  scale_y_continuous(breaks=seq(0,1,0.1),limits=c(0,1.1)) +
  theme(
    panel.background = element_rect(fill = "white", color = "black"),
    panel.grid.major = element_line(colour = "light gray"),
    panel.grid.minor = element_line(colour = "light gray"),
  )
dev.off()


