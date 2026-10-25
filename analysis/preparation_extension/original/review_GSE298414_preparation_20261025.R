# Read saved new-source fits only: no model estimation or data acquisition.
base <- 'R:/Shaoyao_storage/data/GSE298414_preparation/analysis_20261025_03'
stopifnot(file.exists(file.path(base,'receipt.json')))
review <- file.path(base,'review');stopifnot(!dir.exists(review));dir.create(review)
load_tables <- function(label) {
  a <- read.delim(gzfile(file.path(base,label,'edgeR_all_genes.tsv.gz')))
  b <- read.delim(gzfile(file.path(base,label,'voom_all_genes.tsv.gz')))
  stopifnot(!anyDuplicated(paste(a$contrast_id,a$gene_id)),identical(a$gene_id,b$gene_id),
    identical(a$contrast_id,b$contrast_id),length(unique(a$contrast_id))==12L,
    isTRUE(all.equal(a$FDR_all_12_contrasts,p.adjust(a$PValue,'BH'),tolerance=1e-12)),
    isTRUE(all.equal(b$FDR_all_12_contrasts,p.adjust(b$P.Value,'BH'),tolerance=1e-12)))
  obj <- readRDS(file.path(base,label,'model_objects.rds'))
  stopifnot(qr(obj$design)$rank==24L,identical(rownames(obj$design),obj$samples$sample_id),
    identical(rownames(obj$contrasts),colnames(obj$design)),
    all(colSums(obj$contrasts)==0),all(colSums(abs(obj$contrasts))==2))
  for(id in unique(a$contrast_id)) {
    aa <- a[a$contrast_id==id,];bb <- b[b$contrast_id==id,]
    stopifnot(identical(aa$gene_id,rownames(obj$voom_fit$coefficients)),
      isTRUE(all.equal(bb$logFC,unname(obj$voom_fit$coefficients[,id]),tolerance=1e-12)))
  }
  list(a=a,b=b,s=obj$samples)
}
full <- load_tables('full_72');omit <- load_tables('omit_disputed_71')
stopifnot(nrow(full$s)==72L,nrow(omit$s)==71L,
          identical(setdiff(full$s$sample_id,omit$s$sample_id),'GSM9013914'))
affected <- 'JGE02_L1_minus_JGE06_L1'
x <- full$a[full$a$contrast_id==affected,];y <- omit$a[omit$a$contrast_id==affected,]
stopifnot(identical(x$gene_id,y$gene_id))
big <- abs(x$logFC)>=1 & x$FDR_all_12_contrasts<=.05
both <- x$FDR_all_12_contrasts<=.05 & y$FDR_all_12_contrasts<=.05
pc <- read.delim(file.path(base,'full_72/PCA_coordinates.tsv'))
pv <- read.delim(file.path(base,'full_72/PCA_variance.tsv'))
qc <- read.delim(file.path(base,'full_72/sample_qc.tsv'))
stopifnot(identical(pc$sample_id,qc$sample_id))
within <- read.delim(file.path(base,'full_72/logCPM_correlation.tsv'),check.names=FALSE)
cm <- as.matrix(within[,-1]);rownames(cm)<-within$sample_id
rep_cor <- do.call(rbind,lapply(unique(qc$group),function(g) {
  ids <- qc$sample_id[qc$group==g];z<-cm[ids,ids]
  data.frame(group=g,min_pair_correlation=min(z[lower.tri(z)]),mean_pair_correlation=mean(z[lower.tri(z)]))
}))
write.table(rep_cor,file.path(review,'within_condition_correlations.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
plot_panels <- function() {
  par(mfrow=c(1,2),mar=c(4.5,4.5,3,1),cex=.9)
  cols <- ifelse(qc$solvent=='water','#1971A4','#B76314')
  shapes <- ifelse(qc$extraction=='CEM',16,17)
  plot(pc$PC1,pc$PC2,col=cols,pch=shapes,cex=1.1,
    xlab=sprintf('PC1 (%.1f%%)',100*pv$variance_fraction[1]),
    ylab=sprintf('PC2 (%.1f%%)',100*pv$variance_fraction[2]),main='72 PC12 libraries')
  j <- which(pc$sample_id=='GSM9013914');points(pc$PC1[j],pc$PC2[j],pch=1,cex=1.8,lwd=1.6)
  legend('topright',c('Water CEM','Water IEM','Ethanol CEM','Ethanol IEM','Disputed label'),
    col=c('#1971A4','#1971A4','#B76314','#B76314','black'),pch=c(16,17,16,17,1),bty='n',cex=.65)
  plot(x$logFC,y$logFC,pch=16,cex=.35,col=adjustcolor('#1971A4',alpha.f=.18),
    xlab='CEM - IEM log2 FC: full 72',ylab='CEM - IEM log2 FC: disputed library omitted',
    main='Ethanol, PR:GR 2:1, low level')
  abline(0,1,lty=2,col='gray40')
}
png(file.path(review,'preparation_QC_and_label_influence.png'),width=1800,height=850,res=160)
plot_panels();dev.off()
pdf(file.path(review,'preparation_QC_and_label_influence.pdf'),width=11.25,height=5.3)
plot_panels();dev.off()
result <- list(status='SAVED_RESULT_CHECKS_PASSED',gene_contrast_rows_per_method=nrow(full$a),
  gene_filter_identical=identical(unique(full$a$gene_id),unique(omit$a$gene_id)),
  p_adjustment_checked='BH across complete gene x 12 family, both methods and both sample sets',
  saved_voom_coefficients_match=TRUE,contrast_direction_checked=TRUE,
  full_libraries=nrow(full$s),omission_libraries=nrow(omit$s),
  within_condition_correlation_min=min(rep_cor$min_pair_correlation),
  within_condition_correlation_median=median(rep_cor$mean_pair_correlation),
  TMM_range=range(qc$TMM),PC1_PC2_variance_fraction=pv$variance_fraction[1:2],
  affected_contrast=affected,affected_effect_spearman=cor(x$logFC,y$logFC,method='spearman'),
  affected_full_global_discoveries=sum(x$FDR_all_12_contrasts<=.05),
  affected_omission_global_discoveries=sum(y$FDR_all_12_contrasts<=.05),
  affected_shared_global_discoveries=sum(both),affected_shared_sign_conflicts=sum(both & sign(x$logFC)!=sign(y$logFC)),
  affected_full_large_discoveries=sum(big),affected_large_direction_retained=sum(big & sign(x$logFC)==sign(y$logFC)),
  model_refits=0)
jsonlite::write_json(result,file.path(review,'review_receipt.json'),auto_unbox=TRUE,pretty=TRUE)
print(result)
