# Validate saved complete-family statistics; no refit.
b <- 'R:/Shaoyao_storage/data/GSE298414_preparation/analysis_20261025_03/high_ethanol_six_library_influence'
stopifnot(file.exists(file.path(b,'receipt.json')))
s <- read.delim(file.path(b,'summary.tsv'))
stopifnot(nrow(s)==6L,!anyDuplicated(s$omitted),sum(s$omitted_PCA_outlier)==1L)
for(i in seq_len(nrow(s))) {
  r <- readRDS(file.path(b,paste0(s$omitted[i],'_complete_family_statistics.rds')))
  t <- read.delim(gzfile(file.path(b,paste0(s$omitted[i],'_affected_contrast.tsv.gz'))))
  id <- 'JGE02_L3_minus_JGE06_L3'
  stopifnot(nrow(r$edgeR_p)==13761L,ncol(r$edgeR_p)==12L,
    identical(dimnames(r$edgeR_p),dimnames(r$voom_p)),all(is.finite(r$voom_p)),
    identical(t$gene_id,rownames(r$edgeR_p)),nrow(r$design)==71L,qr(r$design)$rank==24L)
  eq <- matrix(p.adjust(r$edgeR_p,'BH'),13761,12,dimnames=dimnames(r$edgeR_p))
  vq <- matrix(p.adjust(r$voom_p,'BH'),13761,12,dimnames=dimnames(r$voom_p))
  stopifnot(isTRUE(all.equal(t$edgeR_FDR_all12,unname(eq[,id]),tolerance=1e-12)),
    isTRUE(all.equal(t$voom_FDR_all12,unname(vq[,id]),tolerance=1e-12)),
    s$edgeR_global_FDR05[i]==sum(eq[,id]<=.05),s$voom_global_FDR05[i]==sum(vq[,id]<=.05),
    s$both_methods_global[i]==sum(eq[,id]<=.05 & vq[,id]<=.05))
}
jsonlite::write_json(list(status='ALL_SIX_SAVED_OMISSIONS_VERIFIED',model_refits=0,
  across_12_contrast_edgeR_discovery_range=range(s$edgeR_global_FDR05),
  across_12_contrast_voom_discovery_range=range(s$voom_global_FDR05),
  effect_spearman_range=range(s$edgeR_effect_rho_vs_full),
  outlier_omission=s[s$omitted_PCA_outlier,]),file.path(b,'review_receipt.json'),pretty=TRUE,auto_unbox=TRUE)
print(s)
