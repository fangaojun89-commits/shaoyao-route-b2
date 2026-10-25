# New within-source preparation analysis; no accepted analysis is read or refitted.
suppressPackageStartupMessages({library(edgeR); library(limma); library(tximport)})
base <- 'R:/Shaoyao_storage/data/GSE298414_preparation'
root <- 'C:/Users/fanga/Desktop/Shaoyao/route_b2'
out <- file.path(base, 'analysis_20261025_03')
stopifnot(!dir.exists(out))
dir.create(out)
write_tsv <- function(x, path) write.table(x, path, sep='\t', quote=FALSE, row.names=FALSE, na='NA')
write_gz <- function(x, path) { con <- gzfile(path, 'wt'); on.exit(close(con)); write_tsv(x, con) }
samples <- read.delim(file.path(base,'selected_rsem_72/samples.tsv'),check.names=FALSE)
inventory <- read.delim(file.path(base,'inputs_20261025/cem_iem_contrasts_supplement_verified.tsv'),check.names=FALSE)
stopifnot(nrow(samples)==72L, !anyDuplicated(samples$sample_id), !anyDuplicated(samples$biosample),
          length(unique(samples$group))==24L, all(table(samples$group)==3L), all(samples$batch=='A'),
          nrow(inventory)==12L, !anyDuplicated(inventory$contrast_id))
setup <- function(s) {
  groups <- factor(s$group, levels=sort(unique(samples$group)))
  design <- model.matrix(~0+groups)
  colnames(design) <- levels(groups); rownames(design) <- s$sample_id
  contrast <- matrix(0, ncol(design), nrow(inventory), dimnames=list(colnames(design),inventory$contrast_id))
  for(j in seq_len(nrow(inventory))) {
    contrast[inventory$positive_group[j],j] <- 1
    contrast[inventory$negative_group[j],j] <- -1
    p <- s[s$group==inventory$positive_group[j],]
    n <- s[s$group==inventory$negative_group[j],]
    stopifnot(all(p$extraction=='CEM'),all(n$extraction=='IEM'),
      length(unique(c(p$ratio_PR_GR,n$ratio_PR_GR)))==1L,
      length(unique(c(p$solvent,n$solvent)))==1L,
      length(unique(c(p$ordinal_level,n$ordinal_level)))==1L)
  }
  stopifnot(qr(design)$rank==24L, all(colSums(contrast)==0),all(colSums(abs(contrast))==2))
  # Independently recover each contrast from synthetic group means, testing coding/direction.
  mu <- seq_len(ncol(design)); y <- as.vector(design %*% mu)
  beta <- qr.solve(design,y)
  expected <- vapply(seq_len(nrow(inventory)),function(j)
    mean(y[s$group==inventory$positive_group[j]])-mean(y[s$group==inventory$negative_group[j]]),numeric(1))
  stopifnot(isTRUE(all.equal(as.vector(t(contrast)%*%beta),expected,tolerance=1e-10)))
  list(design=design,contrast=contrast,groups=groups)
}
full_setup <- setup(samples)
paths <- setNames(file.path(base,'selected_rsem_72',samples$member),samples$sample_id)
importer <- function(x) read.delim(gzfile(x),check.names=FALSE)
raw <- tximport(paths,type='rsem',txIn=FALSE,txOut=FALSE,countsFromAbundance='no',importer=importer)
# RSEM gene import deliberately returns unadjusted counts even if an abundance
# conversion argument is supplied. Convert explicitly using deposited gene TPM
# and sample-specific effective gene lengths; verify the implemented formula.
stopifnot(identical(raw$countsFromAbundance,'no'),all(is.finite(raw$length)),all(raw$length>=0),
          all(raw$abundance[raw$length==0]==0))
# Tiny assignments to genes with zero effective length and TPM are retained in
# the source record but not treated as valid expression measurements.
length_valid <- rowSums(raw$length<=0)==0
bad <- which(raw$length==0 & raw$counts>0,arr.ind=TRUE)
write_tsv(data.frame(gene_id=rownames(raw$counts)[bad[,1]],sample_id=colnames(raw$counts)[bad[,2]],
  expected_count=raw$counts[bad],TPM=raw$abundance[bad]),file.path(out,'zero_effective_length_assignments.tsv'))
write_tsv(data.frame(gene_id=rownames(raw$counts),all_selected_effective_lengths_positive=length_valid),
          file.path(out,'effective_length_eligibility.tsv'))
adjusted <- raw
adjusted$counts <- tximport::makeCountsFromAbundance(raw$counts,raw$abundance,raw$length,'lengthScaledTPM')
adjusted$countsFromAbundance <- 'lengthScaledTPM'
manual <- raw$abundance * rowMeans(raw$length)
manual <- sweep(manual,2,colSums(manual),'/')
manual <- sweep(manual,2,colSums(raw$counts),'*')
stopifnot(isTRUE(all.equal(manual,adjusted$counts,tolerance=1e-12)),
          max(abs(adjusted$counts-raw$counts))>1)
rm(manual)
stopifnot(identical(dim(raw$counts),c(32883L,72L)),identical(dimnames(raw$counts),dimnames(adjusted$counts)),
          identical(colnames(adjusted$counts),samples$sample_id),
          all(is.finite(adjusted$counts)),all(adjusted$counts>=0),
          isTRUE(all.equal(colSums(raw$counts),colSums(adjusted$counts),tolerance=1e-8)))
annotation <- read.delim(file.path(root,'data/processed/jgw03_three_cell_lines/GSE289929/gene_annotation.tsv'))
stopifnot(identical(annotation$raw_gene_id,rownames(raw$counts)))
saveRDS(list(counts=adjusted$counts,samples=samples,length_valid=length_valid,countsFromAbundance='lengthScaledTPM'),file.path(out,'imported_counts.rds'))
write_tsv(data.frame(sample_id=samples$sample_id,raw_expected_total=colSums(raw$counts),
  adjusted_total=colSums(adjusted$counts)),file.path(out,'import_totals.tsv'))
rm(raw,adjusted); gc()
imported <- readRDS(file.path(out,'imported_counts.rds'))
fit_set <- function(counts,s,label,keep=NULL) {
  target <- file.path(out,label);dir.create(target)
  st <- setup(s)
  stopifnot(identical(colnames(counts),s$sample_id))
  d <- calcNormFactors(DGEList(counts=counts),method='TMM')
  if(is.null(keep)) keep <- filterByExpr(d,design=st$design) & imported$length_valid
  stopifnot(sum(keep)>100, length(keep)==nrow(counts))
  write_tsv(data.frame(gene_id=rownames(counts),retained=keep,length_eligible=imported$length_valid),file.path(target,'gene_filter.tsv'))
  d <- calcNormFactors(d[keep,,keep.lib.sizes=FALSE],method='TMM')
  d <- estimateDisp(d,st$design,robust=TRUE)
  ql <- glmQLFit(d,st$design,robust=TRUE)
  v <- voom(d,st$design,plot=FALSE)
  vf <- eBayes(contrasts.fit(lmFit(v,st$design),st$contrast),robust=TRUE)
  er <- lv <- vector('list',nrow(inventory))
  for(j in seq_len(nrow(inventory))) {
    a <- topTags(glmQLFTest(ql,contrast=st$contrast[,j]),n=Inf,sort.by='none')$table
    b <- topTable(vf,coef=j,number=Inf,sort.by='none',confint=TRUE)
    stopifnot(identical(rownames(a),rownames(b)))
    er[[j]] <- data.frame(contrast_id=inventory$contrast_id[j],gene_id=rownames(a),a,row.names=NULL)
    lv[[j]] <- data.frame(contrast_id=inventory$contrast_id[j],gene_id=rownames(b),b,row.names=NULL)
  }
  er <- do.call(rbind,er); lv <- do.call(rbind,lv)
  er$FDR_all_12_contrasts <- p.adjust(er$PValue,'BH')
  lv$FDR_all_12_contrasts <- p.adjust(lv$P.Value,'BH')
  stopifnot(all(is.finite(er$PValue)),all(er$PValue>=0 & er$PValue<=1),
            all(is.finite(lv$P.Value)),all(lv$CI.L<=lv$logFC & lv$logFC<=lv$CI.R))
  for(j in seq_len(nrow(inventory))) {
    ix <- er$contrast_id==inventory$contrast_id[j]
    stopifnot(isTRUE(all.equal(er$FDR[ix],p.adjust(er$PValue[ix],'BH'))),
      isTRUE(all.equal(lv$adj.P.Val[ix],p.adjust(lv$P.Value[ix],'BH'))))
  }
  er$gene_symbol <- annotation$gene_symbol[match(er$gene_id,annotation$raw_gene_id)]
  lv$gene_symbol <- annotation$gene_symbol[match(lv$gene_id,annotation$raw_gene_id)]
  write_gz(er,file.path(target,'edgeR_all_genes.tsv.gz'));write_gz(lv,file.path(target,'voom_all_genes.tsv.gz'))
  logcpm <- cpm(d,log=TRUE,prior.count=2)
  corr <- cor(logcpm)
  pc <- prcomp(t(logcpm),center=TRUE,scale.=FALSE,rank.=5)
  qc <- data.frame(s,adjusted_filtered_library=d$samples$lib.size,TMM=d$samples$norm.factors,
                   median_correlation=sapply(seq_len(ncol(corr)),function(i)median(corr[-i,i])))
  write_tsv(qc,file.path(target,'sample_qc.tsv'))
  write_tsv(data.frame(sample_id=s$sample_id,pc$x),file.path(target,'PCA_coordinates.tsv'))
  write_tsv(data.frame(component=seq_along(pc$sdev),variance_fraction=pc$sdev^2/sum(pc$sdev^2)),file.path(target,'PCA_variance.tsv'))
  write_tsv(data.frame(sample_id=colnames(corr),corr,check.names=FALSE),file.path(target,'logCPM_correlation.tsv'))
  saveRDS(list(dge=d,voom=v,voom_fit=vf,ql_fit=ql,design=st$design,contrasts=st$contrast,
               samples=s,keep=keep),file.path(target,'model_objects.rds'))
  summaries <- lapply(seq_len(nrow(inventory)),function(j) {
    a <- er[er$contrast_id==inventory$contrast_id[j],]; b <- lv[lv$contrast_id==inventory$contrast_id[j],]
    stopifnot(identical(a$gene_id,b$gene_id))
    data.frame(analysis=label,inventory[j,c('contrast_id','ratio_PR_GR','solvent','ordinal_level')],
      tested_genes=nrow(a),positive_libraries=sum(s$group==inventory$positive_group[j]),
      negative_libraries=sum(s$group==inventory$negative_group[j]),
      edgeR_FDR05=sum(a$FDR<=.05),edgeR_global_FDR05=sum(a$FDR_all_12_contrasts<=.05),
      edgeR_global_FDR05_absLogFC1=sum(a$FDR_all_12_contrasts<=.05 & abs(a$logFC)>=1),
      voom_FDR05=sum(b$adj.P.Val<=.05),voom_global_FDR05=sum(b$FDR_all_12_contrasts<=.05),
      effect_spearman=cor(a$logFC,b$logFC,method='spearman'),
      shared_global_significant=sum(a$FDR_all_12_contrasts<=.05 & b$FDR_all_12_contrasts<=.05),
      median_abs_edgeR_logFC=median(abs(a$logFC)),nominal_residual_df=nrow(s)-qr(st$design)$rank)
  })
  summary <- do.call(rbind,summaries);write_tsv(summary,file.path(target,'contrast_summary.tsv'))
  print(summary[,c('analysis','contrast_id','tested_genes','edgeR_global_FDR05','voom_global_FDR05')])
  list(keep=keep,summary=summary,edgeR=er,voom=lv)
}
primary <- fit_set(imported$counts,samples,'full_72')
retain <- samples$sample_id!='GSM9013914'
influence <- fit_set(imported$counts[,retain],samples[retain,],'omit_disputed_71',primary$keep)
comparison <- lapply(inventory$contrast_id,function(id) {
  a <- primary$edgeR[primary$edgeR$contrast_id==id,]
  b <- influence$edgeR[influence$edgeR$contrast_id==id,]
  stopifnot(identical(a$gene_id,b$gene_id))
  da <- a$FDR_all_12_contrasts<=.05;db <- b$FDR_all_12_contrasts<=.05
  data.frame(contrast_id=id,effect_spearman=cor(a$logFC,b$logFC,method='spearman'),
    all_gene_sign_agreement=mean(sign(a$logFC)==sign(b$logFC)),
    full_global_FDR05=sum(da),omission_global_FDR05=sum(db),shared=sum(da&db),
    shared_same_direction=sum(da&db&sign(a$logFC)==sign(b$logFC)),
    max_absolute_logFC_change=max(abs(a$logFC-b$logFC)))
})
write_tsv(do.call(rbind,comparison),file.path(out,'disputed_library_influence.tsv'))
write_tsv(rbind(primary$summary,influence$summary),file.path(out,'all_contrast_summaries.tsv'))
writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
sources <- c(file.path(base,'GSE298414_RAW.tar'),file.path(base,'inputs_20261025/sample_ledger_supplement_verified.tsv'),
  file.path(base,'inputs_20261025/cem_iem_contrasts_supplement_verified.tsv'),
  file.path(root,'code/R/analyze_GSE298414_preparation_20261025.R'),
  file.path(root,'protocol/GSE298414_WITHIN_SOURCE_PREPARATION_ANALYSIS_20261025.md'))
jsonlite::write_json(list(status='NEW_SOURCE_GENE_COMPARISONS_FITTED_REVIEW_PENDING',
  input_libraries=72,conditions=24,contrasts=12,retained_genes=sum(primary$keep),
  primary_method='edgeR robust QL; TMM; explicit gene TPM x across-sample mean effective length rescaled to original count totals',secondary_method='limma-voom robust eBayes',
  influence_omission='GSM9013914 only; fixed full-analysis gene filter and imported length scaling',
  absolute_dose_inference=FALSE,candidate_gate_eligible=FALSE,accepted_analysis_reruns=0,
  inputs=data.frame(path=sources,sha256=vapply(sources,digest::digest,character(1),algo='sha256',file=TRUE))),
  file.path(out,'receipt.json'),pretty=TRUE,auto_unbox=TRUE)
cat('NEW SOURCE ANALYSIS FINISHED; SCIENTIFIC REVIEW PENDING\n')
