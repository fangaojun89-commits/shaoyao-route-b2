# New contrasts and program outcomes from saved 72-library inputs/fits.
suppressPackageStartupMessages({library(edgeR);library(limma);library(decoupleR)})
root <- 'C:/Users/fanga/Desktop/Shaoyao/route_b2'
base <- 'R:/Shaoyao_storage/data/GSE298414_preparation/analysis_20261025_03'
out <- file.path(base,'program_and_stratified_extension_01')
stopifnot(!dir.exists(out) || length(list.files(out,all.files=TRUE,no..=TRUE))==0L)
if(!dir.exists(out)) dir.create(out)
wtsv <- function(x,p) write.table(x,p,sep='\t',row.names=FALSE,quote=FALSE,na='NA')
wgz <- function(x,p) {z<-gzfile(p,'wt');on.exit(close(z));wtsv(x,z)}
obj <- readRDS(file.path(base,'full_72/model_objects.rds'))
s <- obj$samples;design <- obj$design
stopifnot(nrow(s)==72L,ncol(design)==24L,qr(design)$rank==24L)
old <- read.delim('R:/Shaoyao_storage/data/GSE298414_preparation/inputs_20261025/cem_iem_contrasts_supplement_verified.tsv')
inventory <- data.frame(contrast_id=old$contrast_id,contrast_class='extraction',
                       positive_group=old$positive_group,negative_group=old$negative_group)
for(code in c('02','04','06','08')) for(level in 1:3) {
  pos<-paste0('JGE',code,'_L',level);neg<-paste0('JGW',code,'_L',level)
  inventory<-rbind(inventory,data.frame(contrast_id=paste0(pos,'_minus_',neg),contrast_class='solvent',positive_group=pos,negative_group=neg))
}
for(solvent in c('W','E')) for(pair in list(c('02','04'),c('06','08'))) for(level in 1:3) {
  pos<-paste0('JG',solvent,pair[1],'_L',level);neg<-paste0('JG',solvent,pair[2],'_L',level)
  inventory<-rbind(inventory,data.frame(contrast_id=paste0(pos,'_minus_',neg),contrast_class='ratio',positive_group=pos,negative_group=neg))
}
stopifnot(nrow(inventory)==36L,!anyDuplicated(inventory$contrast_id),all(table(inventory$contrast_class)==12L))
contrasts<-matrix(0,24,36,dimnames=list(colnames(design),inventory$contrast_id))
for(j in 1:36) {
  contrasts[inventory$positive_group[j],j]<-1;contrasts[inventory$negative_group[j],j]<- -1
  p<-s[s$group==inventory$positive_group[j],];n<-s[s$group==inventory$negative_group[j],]
  stopifnot(nrow(p)==3L,nrow(n)==3L,length(unique(c(p$ordinal_level,n$ordinal_level)))==1L)
  cls<-inventory$contrast_class[j]
  if(cls=='solvent') stopifnot(all(p$solvent=='70_percent_ethanol'),all(n$solvent=='water'),
    length(unique(c(p$extraction,n$extraction)))==1L,length(unique(c(p$ratio_PR_GR,n$ratio_PR_GR)))==1L)
  if(cls=='ratio') stopifnot(all(p$ratio_PR_GR=='2:1'),all(n$ratio_PR_GR=='1:2'),
    length(unique(c(p$extraction,n$extraction)))==1L,length(unique(c(p$solvent,n$solvent)))==1L)
}
stopifnot(identical(contrasts[,1:12],obj$contrasts),all(colSums(contrasts)==0))
wtsv(inventory,file.path(out,'contrast_inventory.tsv'))
adjust_class <- function(tab,p='PValue',collection=NULL) {
  tab$FDR_class <- NA_real_
  group<-if(is.null(collection)) tab$contrast_class else interaction(tab[[collection]],tab$contrast_class,drop=TRUE)
  for(g in unique(group)) {ix<-which(group==g);tab$FDR_class[ix]<-p.adjust(tab[[p]][ix],'BH')}
  tab
}

# Saved QL fit needs only new contrasts. The base voom coefficients/SEs are
# reconstructed once with the already saved expression and precision weights.
lm_base<-lmFit(obj$voom,design)
vf<-eBayes(contrasts.fit(lm_base,contrasts),robust=TRUE)
stopifnot(isTRUE(all.equal(vf$coefficients[,1:12],obj$voom_fit$coefficients,tolerance=1e-12)),
          isTRUE(all.equal(vf$p.value[,1:12],obj$voom_fit$p.value,tolerance=1e-12)))
saveRDS(lm_base,file.path(out,'saved_voom_base_linear_fit.rds'))
er<-lv<-list()
for(j in 13:36) {
  a<-topTags(glmQLFTest(obj$ql_fit,contrast=contrasts[,j]),n=Inf,sort.by='none')$table
  b<-topTable(vf,coef=j,number=Inf,sort.by='none',confint=TRUE)
  stopifnot(identical(rownames(a),rownames(b)))
  er[[j-12]]<-data.frame(inventory[j,c('contrast_id','contrast_class')],gene_id=rownames(a),a,row.names=NULL)
  lv[[j-12]]<-data.frame(inventory[j,c('contrast_id','contrast_class')],gene_id=rownames(b),b,row.names=NULL)
}
er<-adjust_class(do.call(rbind,er));lv<-adjust_class(do.call(rbind,lv),'P.Value')
wgz(er,file.path(out,'new_24_edgeR_genes.tsv.gz'));wgz(lv,file.path(out,'new_24_voom_genes.tsv.gz'))
gene_summary<-do.call(rbind,lapply(inventory$contrast_id[13:36],function(id) {
  a<-er[er$contrast_id==id,];b<-lv[lv$contrast_id==id,]
  data.frame(contrast_id=id,contrast_class=a$contrast_class[1],tested_genes=nrow(a),
    edgeR_within05=sum(a$FDR<=.05),edgeR_class05=sum(a$FDR_class<=.05),
    voom_within05=sum(b$adj.P.Val<=.05),voom_class05=sum(b$FDR_class<=.05),
    shared_class05=sum(a$FDR_class<=.05 & b$FDR_class<=.05),
    effect_spearman=cor(a$logFC,b$logFC,method='spearman'),median_abs_edgeR_logFC=median(abs(a$logFC)))
}))
wtsv(gene_summary,file.path(out,'new_24_gene_summary.tsv'));rm(er,lv);gc()
cat('New gene contrasts saved; old 12 voom contrasts bridged exactly.\n')

# Canonical label-independent symbol selection; all mapping decisions preserved.
ann<-read.delim(file.path(root,'data/processed/jgw03_three_cell_lines/GSE289929/gene_annotation.tsv'))
ann<-ann[match(rownames(obj$dge),ann$raw_gene_id),]
stopifnot(identical(ann$raw_gene_id,rownames(obj$dge)))
expr<-cpm(obj$dge,log=TRUE,prior.count=2)
collapse<-function(column) {
  symbol<-ann[[column]];idx<-which(!is.na(symbol)&nzchar(symbol)&symbol!='NA')
  idx<-idx[order(symbol[idx],-rowMeans(expr[idx,,drop=FALSE]),ann$raw_gene_id[idx])]
  ix<-idx[!duplicated(symbol[idx])]
  mat<-expr[ix,,drop=FALSE];rownames(mat)<-symbol[ix]
  list(index=ix,expression=mat,mapping=data.frame(raw_gene_id=ann$raw_gene_id[ix],symbol=symbol[ix],mean_logCPM=rowMeans(mat)))
}
human<-collapse('human_symbol');native<-collapse('gene_symbol')
wtsv(human$mapping,file.path(out,'human_symbol_selection.tsv'));wtsv(native$mapping,file.path(out,'rat_symbol_selection.tsv'))
resource<-file.path(root,'data/processed/gene_sets')
reactome<-read.delim(file.path(resource,'reactome97_human_long.tsv'))
syngo<-read.delim(file.path(resource,'quickgo_syngo_mammalian_direct_human.tsv'))
collections<-list(Reactome97=unique(reactome[,c('term_id','term_name','gene_symbol')]),
                 SynGO_sensitivity=unique(syngo[,c('term_id','term_name','gene_symbol')]))
vh<-obj$voom[human$index,];rownames(vh$E)<-rownames(human$expression)
if(!is.null(vh$weights)) rownames(vh$weights)<-rownames(vh$E)
cam_all<-list();coverage<-list()
for(collection in names(collections)) {
  r<-collections[[collection]]
  sets<-lapply(split(r$gene_symbol,r$term_id),unique)
  index<-lapply(sets,function(g)which(rownames(vh$E)%in%g))
  eligible<-lengths(index)>=10L & lengths(index)<=500L
  coverage[[collection]]<-data.frame(collection=collection,term_id=names(sets),
    reference_members=lengths(sets),mapped_expressed_members=lengths(index),eligible=eligible)
  idx<-index[eligible];stopifnot(length(idx)>0)
  lab<-setNames(r$term_name[!duplicated(r$term_id)],r$term_id[!duplicated(r$term_id)])
  for(j in 1:36) {
    a<-camera(vh,index=idx,design=design,contrast=contrasts[,j],inter.gene.cor=.01,use.ranks=FALSE,sort=FALSE)
    stopifnot(identical(rownames(a),names(idx)),all(is.finite(a$PValue)))
    effect<-vf$coefficients[human$index,j]
    cam_all[[paste(collection,j)]]<-data.frame(collection=collection,inventory[j,c('contrast_id','contrast_class')],
      term_id=rownames(a),term_name=unname(lab[rownames(a)]),a,
      mean_member_logFC=vapply(idx,function(k)mean(effect[k]),numeric(1)),
      median_member_logFC=vapply(idx,function(k)median(effect[k]),numeric(1)),background_genes=nrow(vh$E),row.names=NULL)
  }
  cat(collection,'CAMERA comparisons complete\n')
}
cam<-adjust_class(do.call(rbind,cam_all),collection='collection')
wtsv(do.call(rbind,coverage),file.path(out,'unsigned_program_coverage.tsv'))
wgz(cam,file.path(out,'camera_all_programs.tsv.gz'))
saveRDS(list(human=human,native=native,contrasts=contrasts,inventory=inventory),file.path(out,'fixed_program_mapping.rds'))
rm(cam_all,vh,reactome,syngo);gc()

# Footprints: ULM scores, followed by library-level treatment inference.
score_fit<-function(scores,keep_sample=rep(TRUE,nrow(s)),analysis='full_72',collection) {
  des<-design[keep_sample,,drop=FALSE]
  stopifnot(qr(des)$rank==24L)
  fit<-eBayes(contrasts.fit(lmFit(scores[,keep_sample,drop=FALSE],des),contrasts),robust=TRUE,trend=FALSE)
  ans<-do.call(rbind,lapply(1:36,function(j) {
    a<-topTable(fit,coef=j,number=Inf,sort.by='none',confint=TRUE)
    data.frame(analysis=analysis,collection=collection,inventory[j,c('contrast_id','contrast_class')],term_id=rownames(a),
      estimate=a$logFC,CI_low=a$CI.L,CI_high=a$CI.R,t=a$t,PValue=a$P.Value,FDR=a$adj.P.Val,
      residual_df=nrow(des)-24L,row.names=NULL)
  }))
  stopifnot(all(is.finite(ans$PValue)),all(ans$CI_low<=ans$estimate & ans$estimate<=ans$CI_high))
  adjust_class(ans,collection='collection')
}
pro<-read.delim(file.path(resource,'progeny_human_top500.tsv'))
dor<-read.delim(file.path(resource,'dorothea_rat_AB_one_to_one.tsv'))
nets<-list(PROGENy=data.frame(source=pro$pathway,target=pro$gene_symbol,mor=pro$weight),
           DoRothEA=data.frame(source=dor$tf,target=dor$target,mor=dor$mor))
activities<-fits<-act_cov<-influence<-list()
omit_ids<-c('GSM9013914',s$sample_id[s$group%in%c('JGE02_L3','JGE06_L3')])
stopifnot(length(omit_ids)==7L,!anyDuplicated(omit_ids))
for(collection in names(nets)) {
  expression<-if(collection=='PROGENy')human$expression else native$expression
  threshold<-if(collection=='PROGENy')50L else 5L
  ref<-unique(nets[[collection]])
  stopifnot(all(is.finite(ref$mor)),!anyDuplicated(ref[c('source','target')]))
  terms<-sort(unique(ref$source));net<-ref[ref$target%in%rownames(expression)&ref$mor!=0,]
  n<-table(factor(net$source,levels=terms));eligible<-n>=threshold
  act_cov[[collection]]<-data.frame(collection=collection,term_id=terms,expressed_nonzero_targets=as.integer(n),
    required_minimum=threshold,eligible=as.logical(eligible))
  net<-net[net$source%in%terms[eligible],]
  a<-decoupleR::run_ulm(expression,net,.source='source',.target='target',.mor='mor',minsize=threshold,center=FALSE)
  scores<-matrix(NA_real_,sum(eligible),ncol(expression),dimnames=list(terms[eligible],colnames(expression)))
  stopifnot(!anyDuplicated(a[c('source','condition')]))
  scores[cbind(match(a$source,rownames(scores)),match(a$condition,colnames(scores)))]<-a$score
  stopifnot(all(is.finite(scores)))
  activities[[collection]]<-scores
  fits[[collection]]<-score_fit(scores,collection=collection)
  for(omit in omit_ids) influence[[paste(collection,omit)]]<-score_fit(scores,s$sample_id!=omit,paste0('omit_',omit),collection)
  # Independently verify one ULM score against a direct ordinary regression.
  weights<-numeric(nrow(expression));names(weights)<-rownames(expression)
  first<-net[net$source==rownames(scores)[1],];weights[first$target]<-first$mor
  direct<-summary(lm(expression[,1]~weights))$coefficients[2,'t value']
  stopifnot(isTRUE(all.equal(unname(scores[1,1]),unname(direct),tolerance=1e-8)))
  cat(collection,'signed scores and all seven fixed-score omissions complete\n')
}
signed<-do.call(rbind,fits);infl<-do.call(rbind,influence)
wtsv(do.call(rbind,act_cov),file.path(out,'signed_program_coverage.tsv'))
wtsv(signed,file.path(out,'signed_activity_all_programs.tsv'))
wgz(infl,file.path(out,'signed_activity_seven_omissions.tsv.gz'))
saveRDS(activities,file.path(out,'signed_activity_scores.rds'))
summarize<-function(tab,method) do.call(rbind,lapply(split(tab,interaction(tab$collection,tab$contrast_id,drop=TRUE)),function(a)
  data.frame(method=method,collection=a$collection[1],contrast_id=a$contrast_id[1],contrast_class=a$contrast_class[1],
             tested=nrow(a),within_FDR05=sum(a$FDR<=.05),class_FDR05=sum(a$FDR_class<=.05))))
summary<-rbind(summarize(cam,'CAMERA_membership'),summarize(signed,'ULM_score_contrast'))
wtsv(summary,file.path(out,'program_family_summary.tsv'))
writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
files<-c(file.path(base,'full_72/model_objects.rds'),file.path(root,'code/R/extend_GSE298414_preparation_programs_20261025.R'),
  file.path(root,'protocol/GSE298414_PROGRAM_AND_STRATIFIED_EXTENSION_20261025.md'),
  file.path(resource,c('reactome97_human_long.tsv','quickgo_syngo_mammalian_direct_human.tsv',
                      'progeny_human_top500.tsv','dorothea_rat_AB_one_to_one.tsv')))
jsonlite::write_json(list(status='NEW_STRATIFIED_CONTRASTS_AND_PROGRAMS_COMPLETE_REVIEW_PENDING',
  libraries=72,conditions=24,total_contrasts=36,new_gene_contrasts=24,existing_gene_contrasts_preserved=12,
  gene_filter=13761,human_mapped_background=nrow(human$expression),rat_mapped_background=nrow(native$expression),
  signed_activity_omissions=omit_ids,omission_scope='Fixed full-sample scores; no CAMERA or full-pipeline omission claim',
  gene_normalization_dispersion_voom_weights_reused=TRUE,candidate_gate_changes=0,
  source_manifest=data.frame(path=files,sha256=vapply(files,digest::digest,character(1),algo='sha256',file=TRUE))),
  file.path(out,'receipt.json'),pretty=TRUE,auto_unbox=TRUE)
cat('EXTENSION COMPLETE; REVIEW PENDING\n')
