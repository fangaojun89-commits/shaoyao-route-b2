args<-commandArgs(TRUE);stopifnot(length(args)==2L,!dir.exists(args[2L]))
dir.create(args[2L],recursive=TRUE);root<-args[2L]
a<-data.frame(id=c('a','b'),effect=c(1,-2),q_global=c(.01,.1),N=c(3L,4L))
variants<-list(equal=a,within=a,outside=a,count=a,missing=a,order=a[2:1,],threshold=a)
variants$within$effect[1]<-1+1e-12
variants$outside$effect[1]<-1+1e-4
variants$count$N[1]<-4L
variants$missing$effect[1]<-NA_real_
variants$threshold$q_global[1]<-.06
spec<-list()
for(n in names(variants)) {
  old<-file.path(root,paste0(n,'_old.tsv'));new<-file.path(root,paste0(n,'_new.tsv'))
  write.table(a,old,sep='\t',row.names=FALSE);write.table(variants[[n]],new,sep='\t',row.names=FALSE)
  spec[[length(spec)+1]]<-data.frame(name=n,kind='table',baseline=old,rebuilt=new)
}
for(n in c('rds_equal','rds_changed')) {
  old<-file.path(root,paste0(n,'_old.rds'));new<-file.path(root,paste0(n,'_new.rds'))
  x<-list(matrix=matrix(1:4,2),frame=a);y<-x
  if(n=='rds_changed') y$matrix[1,1]<-8L
  saveRDS(x,old);saveRDS(y,new)
  spec[[length(spec)+1]]<-data.frame(name=n,kind='RDS',baseline=old,rebuilt=new)
}
sp<-file.path(root,'spec.tsv');write.table(do.call(rbind,spec),sp,sep='\t',row.names=FALSE)
stopifnot(system2(file.path(R.home('bin'),'Rscript'),c('--vanilla',shQuote(args[1L]),shQuote(sp),shQuote(file.path(root,'out'))))==0)
r<-read.delim(file.path(root,'out','COMPONENT_COMPARISON.tsv'))
passed<-vapply(split(r$pass,r$file),all,logical(1))
stopifnot(all(passed[c('equal','within','rds_equal')]),!any(passed[c('outside','count','missing','order','threshold','rds_changed')]))
write.table(data.frame(case=names(passed),all_checks_pass=passed),file.path(root,'CASE_RESULTS.tsv'),sep='\t',row.names=FALSE)
cat('NINE_SYNTHETIC_COMPARISON_CASES_CONFIRMED_NOT_REAL_G277_ACCEPTANCE\n')
