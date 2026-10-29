# Read-only comparison; original row order and fixed tolerance 1e-9.
# spec TSV: name, kind (table/RDS), baseline, rebuilt. Output must be new.
args <- commandArgs(TRUE)
stopifnot(length(args) == 2L, !file.exists(args[2L]))
spec <- read.delim(args[1L], stringsAsFactors=FALSE, check.names=FALSE)
stopifnot(!anyDuplicated(spec$name), all(spec$kind %in% c('table','RDS')))
dir.create(args[2L], recursive=TRUE)
records <- list()
add <- function(file, component, rule, pass, detail='', maximum=NA_real_) {
  records[[length(records)+1L]] <<- data.frame(file=file, component=component,
    rule=rule, pass=isTRUE(pass), detail=paste(detail,collapse=' | '),
    max_absolute_difference=maximum, stringsAsFactors=FALSE)
}
same <- function(x,y) isTRUE(all.equal(x,y,tolerance=0,check.attributes=TRUE))
numeric_check <- function(x,y,file,component,exact=FALSE) {
  missing <- identical(is.na(x),is.na(y)) && identical(is.nan(x),is.nan(y))
  add(file,component,'missingness',missing)
  eq <- all.equal(x,y,tolerance=1e-9,check.attributes=TRUE)
  finite <- is.finite(x)&is.finite(y)
  delta <- if(any(finite)) max(abs(x[finite]-y[finite])) else 0
  add(file,component,if(exact) 'exact_discrete_numeric' else 'R_all.equal_1e-9',
    if(exact) same(x,y) else isTRUE(eq),if(isTRUE(eq)) '' else eq,delta)
}
exact_names <- c('entrez_id','time_day','N','exercise_hour','sample_order','replicate',
 'n_samples','design_rank','design_columns','residual_df','input_features',
 'filterByExpr_features','cpm1_in_at_least_3_features','retained_intersection_features',
 'retained_annotated_symbols','library_size','detected_cpm1','tested_features',
 'FDR_0_05','FDR_0_05_up','FDR_0_05_down','programs_tested','within_contrast_q_0_05',
 'positive_q_0_05','negative_q_0_05','ranked_symbols','exact_ties_before','exact_ties_after',
 'tie_break_iterations','pathways_tested','seed','combined_direction','activity_sign',
 'camera_sign','mroast_sign','fgsea_sign','programs')
read_table <- function(p) read.delim(p,quote='"',comment.char='',check.names=FALSE,
  stringsAsFactors=FALSE,na.strings=c('NA','NaN',''),fileEncoding='UTF-8')
for(i in seq_len(nrow(spec))) {
  f <- spec$name[i]
  tryCatch({
    if(spec$kind[i]=='table') {
      a<-read_table(spec$baseline[i]);b<-read_table(spec$rebuilt[i])
      shape<-identical(names(a),names(b)) && identical(dim(a),dim(b)) && !anyDuplicated(names(a))
      add(f,'$','schema_and_row_count',shape)
      if(!shape) next
      for(n in names(a)) {
        x<-a[[n]];y<-b[[n]]
        if(is.numeric(x)&&is.numeric(y)) {
          discrete<-n %in% exact_names || grepl('_n$|_p_0_05$|_q_0_05$|_locked_q_0_05$',n) || grepl('design|contrast_matrix',f)
          numeric_check(x,y,f,n,discrete)
          # Preserve the exact row-level decision domain at the original 0.05 threshold.
          if(grepl('(^[pq]_|^PValue$|^FDR$|_[pq]$|_[pq]_global$|_q_all_locked$)',n))
            add(f,n,'exact_0.05_threshold_flags_in_original_row_order',identical(x<0.05,y<0.05))
          if(grepl('(^logFC$|effect|^ES$|^NES$|signed_sqrt_F|_sign$)',n))
            add(f,n,'exact_effect_direction_in_original_row_order',identical(sign(x),sign(y)))
        } else add(f,n,'exact_identity_category_boolean',identical(x,y))
      }
    } else {
      a<-readRDS(spec$baseline[i]);b<-readRDS(spec$rebuilt[i])
      count<-0L
      walk<-function(x,y,path='$',depth=0L) {
        count<<-count+1L
        stopifnot(count<=100000L,depth<=30L)
        opaque<-function(z) is.environment(z)||is.function(z)||typeof(z)%in%c('externalptr','weakref')||isS4(z)
        if(isS4(x)&&isS4(y)&&all(class(x)%in%c('DGEList','DGEGLM','EList'))&&identical(class(x),class(y))) {
          stopifnot(requireNamespace('edgeR',quietly=TRUE),requireNamespace('limma',quietly=TRUE))
          sx<-methods::slotNames(x);sy<-methods::slotNames(y)
          add(f,path,'S4_slot_names',identical(sx,sy))
          if(identical(sx,sy)) for(sn in sx) walk(methods::slot(x,sn),methods::slot(y,sn),paste0(path,'@',sn),depth+1L)
          return()
        }
        if(opaque(x)||opaque(y)) {add(f,path,'OPAQUE_REQUIRES_REVIEW',FALSE,paste(class(x),collapse=';'));return()}
        structural<-identical(typeof(x),typeof(y))&&identical(class(x),class(y))
        add(f,path,'storage_type_and_class',structural)
        if(!structural) return()
        # Inspect attributes as well as every underlying element; no model fitting or method dispatch.
        atx<-attributes(x);aty<-attributes(y)
        add(f,paste0(path,'@attributes'),'exact_attributes',identical(atx,aty))
        if(is.object(x)) {x<-unclass(x);y<-unclass(y)}
        attributes(x)<-NULL;attributes(y)<-NULL
        add(f,path,'storage_length',length(x)==length(y))
        if(length(x)!=length(y)) return()
        if(is.list(x)) {
          for(j in seq_along(x)) walk(x[[j]],y[[j]],paste0(path,'[[',j,']]'),depth+1L)
        } else if(is.numeric(x)) numeric_check(x,y,f,path,typeof(x)=='integer')
        else add(f,path,'exact_atomic_storage',identical(x,y))
      }
      walk(a,b)
    }
  },error=function(e) add(f,'$','READ_OR_COMPARE_ERROR',FALSE,conditionMessage(e)))
}
result<-do.call(rbind,records)
write.table(result,file.path(args[2L],'COMPONENT_COMPARISON.tsv'),sep='\t',row.names=FALSE,quote=TRUE,na='NA')
write.table(result[!result$pass,,drop=FALSE],file.path(args[2L],'DIFFERENCES.tsv'),sep='\t',row.names=FALSE,quote=TRUE,na='NA')
capture.output(sessionInfo(),file=file.path(args[2L],'sessionInfo.txt'))
cat('FILES',nrow(spec),'CHECKS',nrow(result),'FAILED_CHECKS',sum(!result$pass),'\n')
# A returned comparator is not full acceptance: file identities, projections,
# graphics, provenance and any exact numeric role ambiguities remain separate.
