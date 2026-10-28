# Portable saved-table readback. Base R only; no model fitting or input writes.
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==1L)
root <- normalizePath(args[1],winslash='/',mustWork=TRUE)
read <- function(name) read.delim(file.path(root,name),check.names=FALSE,stringsAsFactors=FALSE)
eq <- function(x,y) stopifnot(isTRUE(all.equal(x,y,tolerance=1e-9,check.attributes=FALSE)))
samples <- read('data/derived/samples.tsv')
qc <- read('data/output/complete_library_qc.tsv')
ena <- read('data/raw_ledger/ENA_quantification_manifest.tsv')
stopifnot(nrow(samples)==12L,!anyDuplicated(samples$run),!anyDuplicated(samples$biosample),
          setequal(samples$run,paste0('SRR217758',46:57)),all(table(samples$treatment)==6L),
          setequal(ena$run_accession,samples$run),identical(qc$run,samples$run),
          all(qc$input_pairs>=qc$retained_pairs),all(qc$retained_pairs>=qc$pseudoaligned_pairs))
eq(qc$pseudoaligned_fraction,qc$pseudoaligned_pairs/qc$retained_pairs)
primary <- read('data/output/genes/primary_edgeR.tsv')
voom <- read('data/output/genes/sensitivity_voom.tsv')
filter <- read('data/output/genes/gene_filter.tsv')
stopifnot(nrow(primary)==17013L,!anyDuplicated(primary$gene_id),identical(primary$gene_id,voom$gene_id),
          setequal(primary$gene_id,filter$gene_id[filter$retained]),
          all(is.finite(primary$logFC)),all(is.finite(voom$logFC)))
eq(primary$FDR,p.adjust(primary$PValue,'BH'));eq(voom$adj.P.Val,p.adjust(voom$P.Value,'BH'))
stopifnot(sum(primary$FDR<.05)==0L,sum(voom$adj.P.Val<.05)==0L,
          all(voom$CI.L<=voom$logFC & voom$logFC<=voom$CI.R))
cam <- read('data/output/camera_complete_results.tsv')
act <- read('data/output/signed_activity_contrasts.tsv')
stopifnot(!anyDuplicated(cam[c('collection','term_id')]),!anyDuplicated(act[c('collection','term_id')]))
expected <- list(Reactome=c(1588L,65L),DoRothEA=c(91L,4L),PROGENy=c(14L,0L))
for(cl in names(expected)) {
  z <- cam[cam$collection==cl,]
  stopifnot(nrow(z)==expected[[cl]][1],sum(z$FDR<.05)==expected[[cl]][2])
  eq(z$FDR,p.adjust(z$PValue,'BH'))
}
for(cl in c('PROGENy','DoRothEA')) {
  z <- act[act$collection==cl,]
  scores <- read(paste0('data/output/',cl,'_per_library_scores.tsv'))
  stopifnot(nrow(z)==if(cl=='PROGENy')14L else 91L,all(z$df==10),all(z$status=='TESTED'),
            sum(z$FDR<.05)==0L,identical(scores$term_id,z$term_id),
            identical(names(scores)[-1],samples$run))
  m <- as.matrix(scores[,-1]);storage.mode(m)<-'double'
  b <- m[,samples$treatment=='Baclofen',drop=FALSE]
  v <- m[,samples$treatment=='Saline solution',drop=FALSE]
  vb <- apply(b,1,var);vv <- apply(v,1,var)
  se <- sqrt(((5*vb+5*vv)/10)*(1/6+1/6))
  eq(z$mean_baclofen,rowMeans(b));eq(z$mean_vehicle,rowMeans(v))
  eq(z$estimate,rowMeans(b)-rowMeans(v));eq(z$SE,se);eq(z$t,z$estimate/se)
  eq(z$PValue,2*pt(-abs(z$t),z$df));eq(z$FDR,p.adjust(z$PValue,'BH'))
  eq(z$ci_low,z$estimate-qt(.975,z$df)*z$SE);eq(z$ci_high,z$estimate+qt(.975,z$df)*z$SE)
}
pairs <- read('data/output/comparison/complete_program_comparisons.tsv')
summary <- read('data/output/comparison/descriptive_comparison_summary.tsv')
basis <- read('data/derived/retained_PROGENy_comparator_basis.tsv')
stopifnot(nrow(pairs)==266L,nrow(summary)==19L,nrow(basis)==266L,
          !anyDuplicated(pairs[c('comparator_id','program_id')]),
          identical(pairs$program_id,basis$program_id),identical(pairs$comparator_id,basis$comparator_id))
for(key in setdiff(names(basis),'SCI_effect')) eq(pairs[[key]],basis[[key]])
pa <- act[act$collection=='PROGENy',]
for(id in summary$comparator_id) {
  p <- pairs[pairs$comparator_id==id,];s <- summary[summary$comparator_id==id,]
  a <- pa[match(p$program_id,pa$term_id),]
  stopifnot(nrow(p)==14L,all(is.finite(p$effect)),all(is.finite(p$SCI_effect)))
  for(map in list(c('SCI_effect','estimate'),c('SCI_ci_low','ci_low'),c('SCI_ci_high','ci_high'),c('SCI_PValue','PValue'),c('SCI_FDR','FDR'))) eq(p[[map[1]]],a[[map[2]]])
  same <- sign(p$effect)==sign(p$SCI_effect)
  stopifnot(all(p$effect!=0 & p$SCI_effect!=0),s$finite_pairs==14L,s$missing_pairs==0L,
            s$exact_zero_pairs==0L,s$same_sign==sum(same),s$opposite_sign==sum(!same),
            identical(p$direction_relation,ifelse(same,'SAME_SIGN','OPPOSITE_SIGN')))
  eq(s$descriptive_spearman,cor(p$effect,p$SCI_effect,method='spearman'))
}
stopifnot(sum(summary$panel_role=='core_context')==9L)
cat('PASS_PORTABLE_TABLE_READBACK\n12 libraries; 17013 genes; 1693 CAMERA rows; 105 signed-score rows; 266 pairs; 19 comparator summaries\nModels fitted: 0. Raw inputs modified: 0.\n')
