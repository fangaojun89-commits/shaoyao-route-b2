# Post-fit QC-triggered influence analysis. Keep every primary sample/result.
# Omit each of all six libraries in the affected comparison once, not only the
# visually unusual library. Reuse imported counts and the full retained universe.
suppressPackageStartupMessages({library(edgeR);library(limma)})
base <- 'R:/Shaoyao_storage/data/GSE298414_preparation/analysis_20261025_03'
out <- file.path(base,'high_ethanol_six_library_influence')
stopifnot(!dir.exists(out));dir.create(out)
obj <- readRDS(file.path(base,'full_72/model_objects.rds'))
input <- readRDS(file.path(base,'imported_counts.rds'))
id <- 'JGE02_L3_minus_JGE06_L3'
selected <- obj$samples$sample_id[obj$samples$group %in% c('JGE02_L3','JGE06_L3')]
stopifnot(length(selected)==6L, 'GSM9013891' %in% selected)
primary <- read.delim(gzfile(file.path(base,'full_72/edgeR_all_genes.tsv.gz')))
primary <- primary[primary$contrast_id==id,]
summary <- list()
for(omit in selected) {
  keep_sample <- obj$samples$sample_id!=omit
  d <- DGEList(input$counts[obj$keep,keep_sample])
  design <- obj$design[keep_sample,,drop=FALSE]
  stopifnot(qr(design)$rank==24L,nrow(design)==71L)
  d <- calcNormFactors(d,method='TMM')
  d <- estimateDisp(d,design,robust=TRUE)
  ql <- glmQLFit(d,design,robust=TRUE)
  v <- voom(d,design,plot=FALSE)
  vf <- eBayes(contrasts.fit(lmFit(v,design),obj$contrasts),robust=TRUE)
  ep <- ef <- matrix(NA_real_,nrow(d),ncol(obj$contrasts),dimnames=list(rownames(d),colnames(obj$contrasts)))
  for(j in seq_len(ncol(obj$contrasts))) {
    tab <- glmQLFTest(ql,contrast=obj$contrasts[,j])$table
    ep[,j] <- tab$PValue;ef[,j] <- tab$logFC
  }
  eq <- matrix(p.adjust(ep,'BH'),nrow(ep),dimnames=dimnames(ep))
  vq <- matrix(p.adjust(vf$p.value,'BH'),nrow(ep),dimnames=dimnames(ep))
  tab <- data.frame(gene_id=rownames(ep),edgeR_logFC=ef[,id],edgeR_P=ep[,id],
    edgeR_FDR_within=p.adjust(ep[,id],'BH'),edgeR_FDR_all12=eq[,id],
    voom_logFC=vf$coefficients[,id],voom_P=vf$p.value[,id],
    voom_FDR_within=p.adjust(vf$p.value[,id],'BH'),voom_FDR_all12=vq[,id])
  stopifnot(identical(tab$gene_id,primary$gene_id),all(is.finite(ep)))
  con <- gzfile(file.path(out,paste0(omit,'_affected_contrast.tsv.gz')),'wt')
  write.table(tab,con,sep='\t',row.names=FALSE,quote=FALSE);close(con)
  saveRDS(list(omitted=omit,design=design,contrasts=obj$contrasts,
    edgeR_p=ep,edgeR_logFC=ef,voom_p=vf$p.value,voom_logFC=vf$coefficients),
    file.path(out,paste0(omit,'_complete_family_statistics.rds')))
  both <- tab$edgeR_FDR_all12<=.05 & tab$voom_FDR_all12<=.05
  primary_sig <- primary$FDR_all_12_contrasts<=.05
  summary[[omit]] <- data.frame(omitted=omit,source_group=obj$samples$group[match(omit,obj$samples$sample_id)],
    omitted_PCA_outlier=omit=='GSM9013891',edgeR_global_FDR05=sum(tab$edgeR_FDR_all12<=.05),
    voom_global_FDR05=sum(tab$voom_FDR_all12<=.05),both_methods_global=sum(both),
    both_methods_same_sign=sum(both & sign(tab$edgeR_logFC)==sign(tab$voom_logFC)),
    edgeR_within_FDR05=sum(tab$edgeR_FDR_within<=.05),voom_within_FDR05=sum(tab$voom_FDR_within<=.05),
    edgeR_effect_rho_vs_full=cor(tab$edgeR_logFC,primary$logFC,method='spearman'),
    retained_primary_discoveries=sum(primary_sig & tab$edgeR_FDR_all12<=.05),
    retained_primary_same_direction=sum(primary_sig & tab$edgeR_FDR_all12<=.05 & sign(tab$edgeR_logFC)==sign(primary$logFC)),
    median_absolute_edgeR_logFC=median(abs(tab$edgeR_logFC)))
  write.table(do.call(rbind,summary),file.path(out,'summary.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
  print(summary[[omit]])
  rm(d,ql,v,vf,ep,ef,eq,vq,tab);gc()
}
jsonlite::write_json(list(status='SIX_FIXED_OMISSIONS_COMPLETE_REVIEW_PENDING',comparison=id,
  libraries_omitted_once=selected,primary_unchanged=TRUE,retained_genes=sum(obj$keep),
  rationale='Visible PCA outlier and low within-condition correlation; all six comparison libraries checked',
  accepted_analysis_reruns=0),file.path(out,'receipt.json'),pretty=TRUE,auto_unbox=TRUE)
