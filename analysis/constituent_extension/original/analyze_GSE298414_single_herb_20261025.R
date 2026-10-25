# New mixture-minus-single-herb estimands; prior 72-library outputs are not modified.
suppressPackageStartupMessages({library(edgeR);library(limma);library(tximport);library(jsonlite)})
root <- 'C:/Users/fanga/Desktop/Shaoyao/route_b2'
out <- 'R:/Shaoyao_storage/data/GSE298414_preparation/single_herb_20261025_01'
wtsv <- function(x,p) write.table(x,p,sep='\t',quote=FALSE,row.names=FALSE,na='NA')
wgz <- function(x,p) {z<-gzfile(p,'wt');on.exit(close(z));wtsv(x,z)}
samples<-read.delim(file.path(out,'samples.tsv'),check.names=FALSE)
inventory<-read.delim(file.path(out,'contrast_inventory.tsv'),check.names=FALSE)
stopifnot(nrow(samples)==126L,length(unique(samples$group))==42L,all(table(samples$group)==3L),
          !anyDuplicated(samples$sample_id),!anyDuplicated(samples$biosample),all(samples$batch=='A'),
          nrow(inventory)==60L,all(table(inventory$contrast_class)==30L))
setup<-function(s) {
  g<-factor(s$group,levels=sort(unique(samples$group)))
  des<-model.matrix(~0+g);colnames(des)<-levels(g);rownames(des)<-s$sample_id
  cmat<-matrix(0,ncol(des),nrow(inventory),dimnames=list(colnames(des),inventory$contrast_id))
  for(j in seq_len(nrow(inventory))) {
    cmat[inventory$positive_group[j],j]<-1;cmat[inventory$negative_group[j],j]<- -1
    a<-s[s$group==inventory$positive_group[j],];b<-s[s$group==inventory$negative_group[j],]
    stopifnot(nrow(a)>=2L,nrow(b)>=2L,all(a$preparation_class%in%c('CEM','IEM')),
      all(b$preparation_class==if(inventory$contrast_class[j]=='mixture_minus_PR') 'single_herb_PR' else 'single_herb_GR'),
      length(unique(c(a$solvent,b$solvent)))==1L,length(unique(c(a$ordinal_level,b$ordinal_level)))==1L)
  }
  stopifnot(qr(des)$rank==42L,all(colSums(cmat)==0),all(colSums(abs(cmat))==2))
  mu<-seq_len(ncol(des));y<-as.vector(des%*%mu);b<-qr.solve(des,y)
  direct<-vapply(seq_len(nrow(inventory)),function(j)mean(y[s$group==inventory$positive_group[j]])-mean(y[s$group==inventory$negative_group[j]]),numeric(1))
  stopifnot(isTRUE(all.equal(as.vector(t(cmat)%*%b),direct,tolerance=1e-10)))
  list(design=des,contrasts=cmat)
}
invisible(setup(samples))
import_path<-file.path(out,'imported_counts.rds')
if(!file.exists(import_path)) {
  paths<-setNames(samples$rsem_path,samples$sample_id)
  raw<-tximport(paths,type='rsem',txIn=FALSE,txOut=FALSE,countsFromAbundance='no',
               importer=function(p)read.delim(gzfile(p),check.names=FALSE))
  stopifnot(identical(dim(raw$counts),c(32883L,126L)),identical(colnames(raw$counts),samples$sample_id),
    identical(raw$countsFromAbundance,'no'),all(is.finite(raw$length)),all(raw$length>=0),all(raw$abundance[raw$length==0]==0))
  counts<-tximport::makeCountsFromAbundance(raw$counts,raw$abundance,raw$length,'lengthScaledTPM')
  manual<-raw$abundance*rowMeans(raw$length);manual<-sweep(manual,2,colSums(manual),'/');manual<-sweep(manual,2,colSums(raw$counts),'*')
  stopifnot(isTRUE(all.equal(counts,manual,tolerance=1e-12)),all(is.finite(counts)),all(counts>=0),
    isTRUE(all.equal(colSums(counts),colSums(raw$counts),tolerance=1e-8)))
  valid<-rowSums(raw$length<=0)==0
  wtsv(data.frame(sample_id=samples$sample_id,raw_expected_total=colSums(raw$counts),adjusted_total=colSums(counts)),file.path(out,'import_totals.tsv'))
  wtsv(data.frame(gene_id=rownames(counts),all_effective_lengths_positive=valid),file.path(out,'effective_length_eligibility.tsv'))
  saveRDS(list(counts=counts,samples=samples,length_valid=valid,countsFromAbundance='lengthScaledTPM'),import_path)
  rm(raw,counts,manual,valid);gc()
}
input<-readRDS(import_path)
stopifnot(identical(input$samples,samples))
annotation<-read.delim(file.path(root,'data/processed/jgw03_three_cell_lines/GSE289929/gene_annotation.tsv'))
stopifnot(identical(annotation$raw_gene_id,rownames(input$counts)))
adjust_class<-function(x,p) {
 x$FDR_class<-NA_real_
 for(cl in unique(x$contrast_class)) {i<-x$contrast_class==cl;x$FDR_class[i]<-p.adjust(x[[p]][i],'BH')}
 x
}
fit_set<-function(s,label,keep=NULL) {
 target<-file.path(out,label)
 marker<-file.path(target,'COMPLETE.json')
 if(file.exists(marker)) {cat('Reusing completed phase',label,'\n');return(readRDS(file.path(target,'model_objects.rds'))$keep)}
 stopifnot(!dir.exists(target));dir.create(target)
 st<-setup(s);counts<-input$counts[,s$sample_id,drop=FALSE]
 d<-calcNormFactors(DGEList(counts),method='TMM')
 if(is.null(keep)) keep<-filterByExpr(d,design=st$design)&input$length_valid
 wtsv(data.frame(gene_id=rownames(counts),retained=keep,length_eligible=input$length_valid),file.path(target,'gene_filter.tsv'))
 d<-calcNormFactors(d[keep,,keep.lib.sizes=FALSE],method='TMM')
 d<-estimateDisp(d,st$design,robust=TRUE);ql<-glmQLFit(d,st$design,robust=TRUE)
 v<-voom(d,st$design,plot=FALSE);base_fit<-lmFit(v,st$design)
 vf<-eBayes(contrasts.fit(base_fit,st$contrasts),robust=TRUE)
 saveRDS(list(dge=d,voom=v,voom_base_fit=base_fit,voom_fit=vf,ql_fit=ql,design=st$design,contrasts=st$contrasts,samples=s,keep=keep),file.path(target,'model_objects.rds'))
 er<-lv<-vector('list',60)
 for(j in 1:60) {
  a<-topTags(glmQLFTest(ql,contrast=st$contrasts[,j]),n=Inf,sort.by='none')$table
  b<-topTable(vf,coef=j,number=Inf,sort.by='none',confint=TRUE)
  stopifnot(identical(rownames(a),rownames(b)))
  er[[j]]<-data.frame(inventory[j,c('contrast_id','contrast_class')],gene_id=rownames(a),a,row.names=NULL)
  lv[[j]]<-data.frame(inventory[j,c('contrast_id','contrast_class')],gene_id=rownames(b),b,row.names=NULL)
 }
 er<-adjust_class(do.call(rbind,er),'PValue');lv<-adjust_class(do.call(rbind,lv),'P.Value')
 stopifnot(all(is.finite(er$PValue)),all(is.finite(lv$P.Value)),all(lv$CI.L<=lv$logFC & lv$logFC<=lv$CI.R))
 er$gene_symbol<-annotation$gene_symbol[match(er$gene_id,annotation$raw_gene_id)]
 lv$gene_symbol<-annotation$gene_symbol[match(lv$gene_id,annotation$raw_gene_id)]
 wgz(er,file.path(target,'edgeR_all_genes.tsv.gz'));wgz(lv,file.path(target,'voom_all_genes.tsv.gz'))
 summary<-do.call(rbind,lapply(inventory$contrast_id,function(id){
  a<-er[er$contrast_id==id,];b<-lv[lv$contrast_id==id,]
  data.frame(analysis=label,contrast_id=id,contrast_class=a$contrast_class[1],tested_genes=nrow(a),
    edgeR_within05=sum(a$FDR<=.05),edgeR_class05=sum(a$FDR_class<=.05),
    voom_within05=sum(b$adj.P.Val<=.05),voom_class05=sum(b$FDR_class<=.05),
    shared_class05=sum(a$FDR_class<=.05 & b$FDR_class<=.05),effect_spearman=cor(a$logFC,b$logFC,method='spearman'),
    median_abs_edgeR_logFC=median(abs(a$logFC)))
 }))
 wtsv(summary,file.path(target,'contrast_summary.tsv'));rm(er,lv);gc()
 expr<-cpm(d,log=TRUE,prior.count=2);pc<-prcomp(t(expr),rank.=5);corr<-cor(expr)
 wtsv(data.frame(s,adjusted_filtered_library=d$samples$lib.size,TMM=d$samples$norm.factors,
  median_correlation=sapply(seq_len(ncol(corr)),function(i)median(corr[-i,i]))),file.path(target,'sample_qc.tsv'))
 wtsv(data.frame(sample_id=s$sample_id,pc$x),file.path(target,'PCA_coordinates.tsv'))
 wtsv(data.frame(component=seq_along(pc$sdev),variance_fraction=pc$sdev^2/sum(pc$sdev^2)),file.path(target,'PCA_variance.tsv'))
 write_json(list(status='NEW_GENE_CONTRASTS_COMPLETE',libraries=nrow(s),groups=42,contrasts=60,genes=nrow(d),residual_df=nrow(s)-42),marker,pretty=TRUE,auto_unbox=TRUE)
 cat(label,'completed',nrow(d),'genes, 60 new contrasts\n')
 keep
}
keep<-fit_set(samples,'full_126');gc()
invisible(fit_set(samples[samples$sample_id!='GSM9013891',],'omit_GSM9013891_125',keep));gc()
writeLines(capture.output(sessionInfo()),file.path(out,'gene_sessionInfo.txt'))
write_json(list(status='GENE_PRIMARY_AND_PLANNED_OMISSION_COMPLETE_REVIEW_PENDING',libraries=126,conditions=42,new_contrasts=60,
 primary_genes=sum(keep),planned_omission='GSM9013891',prior_36_contrasts_rerun=FALSE,candidate_gate_changes=0),
 file.path(out,'gene_receipt.json'),pretty=TRUE,auto_unbox=TRUE)
