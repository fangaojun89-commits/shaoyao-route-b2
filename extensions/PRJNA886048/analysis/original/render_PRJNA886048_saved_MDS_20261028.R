# Presentation-only rendering from the saved gene model; no gene models refitted.
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==2L)
.libPaths(c(file.path(args[1],"environment/R-library/4.6"),"C:/R/Library",.libPaths()))
out <- args[2]
fit <- readRDS(file.path(out,"genes/gene_model_objects.rds"))
samples <- read.delim(file.path(out,"samples.tsv"))
mds <- limma::plotMDS(fit$dge,plot=FALSE)
coords <- data.frame(run=samples$run,treatment=samples$treatment,x=mds$x,y=mds$y)
write.table(coords,file.path(out,"complete_cohort_MDS_coordinates.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
draw <- function() {
  par(mar=c(5,5,3,2))
  expand <- function(x) range(x)+c(-1,1)*diff(range(x))*.18
  colors <- ifelse(samples$treatment=="Baclofen","#1F77B4","#C55A11")
  plot(mds$x,mds$y,type="n",xlim=expand(mds$x),ylim=expand(mds$y),
       xlab=sprintf("Leading logFC dimension 1 (%.0f%%)",100*mds$var.explained[1]),
       ylab=sprintf("Leading logFC dimension 2 (%.0f%%)",100*mds$var.explained[2]),
       main="Injured spinal tissue: all 12 libraries")
  points(mds$x,mds$y,pch=19,col=colors)
  text(mds$x,mds$y,labels=samples$run,pos=3,cex=.8,col=colors)
  legend("bottomleft",legend=c("Baclofen (n = 6)","Saline (n = 6)"),pch=19,
         col=c("#1F77B4","#C55A11"),bty="n",cex=.8)
}
pdf(file.path(out,"complete_cohort_MDS_reviewed.pdf"),width=10,height=7)
draw();dev.off()
png(file.path(out,"complete_cohort_MDS_reviewed.png"),width=1500,height=1050,res=150)
draw();dev.off()
