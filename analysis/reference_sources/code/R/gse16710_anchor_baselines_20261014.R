#!/usr/bin/env Rscript
# Explicit cached-score comparison; never sources the original analyzer.
options(stringsAsFactors = FALSE, digits = 17, scipen = 15, warn = 1)

assert <- function(ok, message) if (!isTRUE(ok)) stop(message, call. = FALSE)
read_tsv <- function(path) read.delim(path, check.names = FALSE, quote = "",
                                    stringsAsFactors = FALSE, na.strings = "NA")
write_tsv <- function(x, path) write.table(x, path, sep = "\t", row.names = FALSE,
                                        col.names = TRUE, quote = FALSE, na = "NA")
write_json <- function(x, path) jsonlite::write_json(x, path, auto_unbox = TRUE,
                                                   pretty = TRUE, digits = 17,
                                                   null = "null", na = "null")
sha <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)
utc <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
contrast_name <- "spasticity_progression_interaction"
families <- c("Reactome", "PROGENy", "DoRothEA")
expected_counts <- c(1656L, 14L, 101L)
score_names <- c("reactome_mean_z_activity_primary.rds", "progeny_ulm_activity_primary.rds",
                 "dorothea_ulm_activity_primary.rds")
native_names <- c("reactome_mean_z_activity_all_contrasts.tsv",
                  "progeny_ulm_activity_all_contrasts.tsv", "dorothea_ulm_activity_all_contrasts.tsv")
group_levels <- c("Sham21", "Sham60", "Spi21", "Spi60")

check_schema <- function(tab, required, numeric_fields=character(), probability_fields=character()) {
  assert(all(required %in% names(tab)),paste("Missing source fields:",paste(setdiff(required,names(tab)),collapse=",")))
  for(field in numeric_fields) {
    assert(length(tab[[field]])==nrow(tab) &&
             (is.numeric(tab[[field]]) || (is.logical(tab[[field]]) && all(is.na(tab[[field]])))),
           paste("Non-numeric source field",field))
    assert(all(is.na(tab[[field]])|is.finite(tab[[field]])),paste("Infinite source field",field))
  }
  for(field in probability_fields) assert(all(is.na(tab[[field]]) |
    (tab[[field]]>=0 & tab[[field]]<=1)),paste("Invalid source probability",field))
  invisible(TRUE)
}
safe_max <- function(x) if(any(is.finite(x))) max(x[is.finite(x)]) else NA_real_
safe_quantiles <- function(x) {
  probs<-c(0,.5,.9,1)
  if(any(is.finite(x))) return(as.list(quantile(x[is.finite(x)],probs)))
  as.list(setNames(rep(NA_real_,4),c("0%","50%","90%","100%")))
}

source_paths <- function() {
  c(setNames(file.path("data/processed/GSE16710", score_names), paste0("score_", families)),
    setNames(file.path("results/GSE16710/program_level", native_names), paste0("native_", families)),
    ledger = "data/processed/GSE16710/sample_ledger_verified.tsv",
    qc = "data/processed/GSE16710/qc_inclusion_decision.tsv",
    design = "results/GSE16710/diagnostics/design_matrix_primary.tsv",
    camera = "results/GSE16710/program_level/reactome_camera_primary_all_contrasts.tsv",
    observed_t = "results/GSE16710/diagnostics/program_Freedman_Lane_10000.tsv",
    s1_gate = "results/GSE16710/program_level/S1_gate_all_collections.tsv",
    master_gate = "results/master_gates/master_gate_summary.tsv",
    sap = "protocol/statistical_analysis_plan.md", yaml = "config/contrasts.yaml",
    original_analyzer = "code/R/analyze_GSE16710_spasticity_anchor.R",
    worker = "code/R/gse16710_anchor_baselines_20261014.R",
    definition = "protocol/GSE16710_ANCHOR_BASELINE_COMPARISON_20261014.md",
    method_review = paste0("results/diagnostics/GSE16710_anchor_baseline_method_review_20261014/",
                           "INDEPENDENT_METHOD_REVIEW.md"),
    source_manifest = "data/processed/GSE16710/source_file_manifest.tsv",
    preprocessing_versions = "data/processed/GSE16710/preprocessing_package_versions.tsv")
}

inventory <- function(root) {
  paths <- source_paths()
  absolute <- file.path(root, paths)
  assert(all(file.exists(absolute)), "Missing fixed input")
  data.frame(role = names(paths), relative_path = unname(paths),
             bytes = unname(file.info(absolute)$size),
             sha256 = vapply(absolute, sha, character(1)), row.names = NULL)
}

same_inventory <- function(a,b) identical(a$role,b$role) &&
  identical(a$relative_path,b$relative_path) && identical(a$sha256,b$sha256) &&
  isTRUE(all(a$bytes==b$bytes))

same_design <- function(x,old,samples) identical(old$sample_id,samples) &&
  identical(colnames(x),names(old)[-1L]) && identical(dim(x),dim(as.matrix(old[-1L]))) &&
  isTRUE(all(x==as.matrix(old[-1L])))

load_inputs <- function(root) {
  path <- function(role) file.path(root, source_paths()[[role]])
  ledger <- read_tsv(path("ledger")); qc <- read_tsv(path("qc"))
  assert(nrow(ledger) == 27L && !anyDuplicated(ledger$sample_id), "Ledger must be 27 unique arrays")
  assert(nrow(qc) == 27L && !anyDuplicated(qc$sample_id) &&
           setequal(qc$sample_id, ledger$sample_id), "QC and ledger keys differ")
  ledger <- ledger[order(as.integer(ledger$sample_order)), , drop = FALSE]
  q <- qc[match(ledger$sample_id, qc$sample_id), , drop = FALSE]
  keep <- q$primary_analysis_decision == "retain"
  assert(!anyNA(keep) && sum(keep) == 26L &&
           identical(ledger$sample_id[!keep], "GSM418682"), "Unexpected QC selection")
  meta <- ledger[keep, , drop = FALSE]
  meta$injury <- factor(meta$experimental_arm, levels = c("Sham", "Spi"))
  meta$time_factor <- factor(as.character(meta$time_day), levels = c("21", "60"))
  assert(!anyNA(meta$injury) && !anyNA(meta$time_factor), "Unknown primary group")
  groups <- paste0(meta$injury, meta$time_factor)
  assert(identical(as.integer(table(factor(groups, levels = group_levels))), c(6L,4L,8L,8L)),
         "Wrong four-cell sample counts")
  x <- model.matrix(~ injury * time_factor, meta); rownames(x) <- meta$sample_id
  old <- read_tsv(path("design"))
  assert(same_design(x,old,meta$sample_id), "Saved design identity mismatch")
  assert(qr(x)$rank == 4L && qr(x[,1:3])$rank == 3L, "Design rank failure")
  conf <- yaml::read_yaml(path("yaml"))
  assert(conf$datasets$GSE16710$expected_n == 27L, "Original pre-QC expected_n changed")
  mats <- lapply(families, function(fam) readRDS(path(paste0("score_", fam))))
  names(mats) <- families
  natives <- lapply(seq_along(families), function(i) {
    mat <- mats[[i]]
    assert(is.matrix(mat) && is.numeric(mat) && nrow(mat) == expected_counts[i] &&
             ncol(mat) == 26L && !is.null(rownames(mat)) && !anyDuplicated(rownames(mat)) &&
             !anyNA(rownames(mat)) && all(nzchar(rownames(mat))) &&
             identical(colnames(mat), meta$sample_id), paste("Score identity failure", families[i]))
    native <- read_tsv(path(paste0("native_", families[i])))
    native_numeric<-c("estimate","standard_error","ci_low","ci_high","moderated_t","df_total","p_value","q_value_BH")
    check_schema(native,c("feature_id","contrast",native_numeric),native_numeric,c("p_value","q_value_BH"))
    assert(all(is.na(native$standard_error)|native$standard_error>=0),"Negative retained standard error")
    native <- native[native$contrast == contrast_name, , drop = FALSE]
    assert(nrow(native) == nrow(mat) && !anyDuplicated(native$feature_id) &&
             setequal(native$feature_id, rownames(mat)), "Native result key mismatch")
    native[match(rownames(mat), native$feature_id), , drop = FALSE]
  })
  names(natives) <- families
  registry <- do.call(rbind, lapply(seq_along(families), function(i)
    data.frame(collection = families[i], feature_id = rownames(mats[[i]]))))
  rownames(registry) <- NULL
  assert(nrow(registry) == 1771L && !anyDuplicated(paste(registry$collection, registry$feature_id)),
         "Feature registry failure")
  list(meta = meta, x = x, groups = groups, mats = mats, natives = natives, features = registry)
}

fit_one <- function(y, x, groups) {
  assert(nrow(x)==26L && ncol(x)==4L && length(y)==26L && length(groups)==26L &&
           qr(x)$rank==4L && qr(x[,1:3,drop=FALSE])$rank==3L,
         "Unestimable or wrong-dimension baseline design")
  empty <- list(status = "nonfinite", estimate = NA_real_, se = NA_real_, t = NA_real_,
                p = NA_real_, ci_low = NA_real_, ci_high = NA_real_, mse = NA_real_,
                sse_full = NA_real_, sse_reduced = NA_real_, nested_f = NA_real_,
                negative_sse_roundoff = FALSE, group_effect_check = NA_real_,
                se_formula_check = NA_real_, rank_effect = NA_real_, rank_unique = NA_integer_,
                rank_tied_observations = NA_integer_, max_abs_studentized = NA_real_,
                max_cook = NA_real_, max_cook_sample = NA_character_,
                group_variance_ratio = NA_real_)
  for (g in group_levels) {
    empty[[paste0("mean_",g)]] <- NA_real_
    empty[[paste0("reduced_mean_",g)]] <- NA_real_
    empty[[paste0("variance_",g)]] <- NA_real_
    empty[[paste0("rank_mean_",g)]] <- NA_real_
  }
  if (!all(is.finite(y))) return(empty)
  full <- lm.fit(x, y); reduced <- lm.fit(x[,1:3,drop=FALSE], y)
  sse <- sum(full$residuals^2); sse0 <- sum(reduced$residuals^2)
  total <- sum((y-mean(y))^2)
  degenerate <- length(unique(y)) == 1L || sse <= 64*.Machine$double.eps*total
  out <- empty
  out$status <- if (degenerate) "numerically_degenerate" else "ok_nominal_working_model"
  out$estimate <- unname(full$coefficients[4]); out$mse <- sse/22
  out$sse_full <- sse; out$sse_reduced <- sse0
  assert(sse-sse0 <= 1e-10*(1+sse0), "Full SSE exceeds reduced SSE beyond roundoff")
  out$negative_sse_roundoff <- sse0 < sse
  r <- rank(y, ties.method = "average")
  out$rank_unique <- length(unique(r))
  out$rank_tied_observations <- sum(duplicated(y) | duplicated(y, fromLast = TRUE))
  means <- vapply(group_levels, function(g) mean(y[groups == g]), numeric(1))
  rmeans <- vapply(group_levels, function(g) mean(r[groups == g]), numeric(1))
  variances <- vapply(group_levels, function(g) var(y[groups == g]), numeric(1))
  direct <- means["Spi60"]-means["Spi21"]-means["Sham60"]+means["Sham21"]
  out$group_effect_check <- unname(out$estimate-direct)
  out$rank_effect <- unname(rmeans["Spi60"]-rmeans["Spi21"]-rmeans["Sham60"]+rmeans["Sham21"])
  out$group_variance_ratio <- if (min(variances)>0) max(variances)/min(variances) else NA_real_
  for (g in group_levels) {
    out[[paste0("mean_",g)]] <- means[[g]]
    out[[paste0("rank_mean_",g)]] <- rmeans[[g]]
    out[[paste0("variance_",g)]] <- variances[[g]]
    out[[paste0("reduced_mean_",g)]] <- mean(reduced$fitted.values[groups == g])
  }
  if (!degenerate) {
    v <- solve(crossprod(x)); h <- rowSums((x %*% v)*x)
    out$se <- sqrt(out$mse*v[4,4]); out$t <- out$estimate/out$se
    out$p <- 2*pt(-abs(out$t), df=22)
    out$ci_low <- out$estimate-qt(.975,22)*out$se
    out$ci_high <- out$estimate+qt(.975,22)*out$se
    out$nested_f <- max(0,sse0-sse)/out$mse
    out$se_formula_check <- out$se-sqrt(out$mse*sum(1/table(groups)))
    student <- full$residuals/sqrt(out$mse*(1-h))
    cook <- student^2*h/(4*(1-h))
    out$max_abs_studentized <- max(abs(student))
    out$max_cook <- max(cook); out$max_cook_sample <- rownames(x)[which.max(cook)]
  }
  out
}

close_enough <- function(a,b,tolerance=1e-8) length(a)>0L && length(a)==length(b) &&
  all(is.finite(a) & is.finite(b) & abs(a-b) <= tolerance*(1+abs(b)))
self_tests <- function() {
  grp <- rep(group_levels,c(6,4,8,8))
  meta <- data.frame(injury=factor(sub("[0-9]+$","",grp),levels=c("Sham","Spi")),
                     time_factor=factor(sub("^[A-Za-z]+","",grp),levels=c("21","60")))
  x <- model.matrix(~injury*time_factor,meta); rownames(x) <- paste0("s",seq_len(26))
  # Center noise within cells so expected interaction is exactly 1.25.
  noise <- seq_len(26) %% 5
  noise <- noise-ave(noise,grp,FUN=mean)
  y <- as.vector(x %*% c(2,-.7,.3,1.25))+noise/3
  z <- fit_one(y,x,grp)
  assert(close_enough(z$estimate,1.25) && abs(z$se_formula_check)<1e-12 &&
           abs(z$group_effect_check)<1e-12 && close_enough(z$nested_f,z$t^2), "Synthetic OLS failure")
  beta <- solve(crossprod(x),crossprod(x,y))
  assert(close_enough(beta[4],z$estimate),"QR versus normal equations failure")
  additive <- fit_one(y-1.25*x[,4],x,grp)
  assert(abs(additive$estimate)<1e-12, "Native additive fixture failure")
  for (bad in list(rep(4,26),as.vector(x %*% c(1,2,3,4)))) {
    b <- fit_one(bad,x,grp)
    assert(b$status=="numerically_degenerate" && is.na(b$p),"Degeneracy handling failure")
  }
  missing <- y; missing[1] <- NA_real_
  assert(fit_one(missing,x,grp)$status=="nonfinite","Missing-score retention failure")
  infinite <- y; infinite[1] <- Inf
  assert(fit_one(infinite,x,grp)$status=="nonfinite","Nonfinite retention failure")
  tied <- c(1,1,3,3,3,7,7,9,rep(11,18))
  assert(identical(rank(tied)[1:5],c(1.5,1.5,4,4,4)),"Average tie rank failure")
  a <- fit_one(tied,x,grp); b <- fit_one(3*tied+7,x,grp); c <- fit_one(-tied,x,grp)
  assert(close_enough(a$rank_effect,b$rank_effect) && close_enough(a$rank_effect,-c$rank_effect),
         "Rank affine/reversal failure")
  perm <- c(26:1)
  assert(close_enough(fit_one(y[perm],x[perm,,drop=FALSE],grp[perm])$estimate,z$estimate),
         "Joint order invariance failure")
  # Unequal cells with an exactly additive raw score, but nonlinear rank interaction.
  nonlin <- as.vector(x %*% c(0,10,1,0))
  nonlin_fit <- fit_one(nonlin,x,grp)
  assert(abs(nonlin_fit$estimate)<1e-12 && abs(nonlin_fit$rank_effect)>0,
         "Distinct native/rank estimand fixture failure")
  fake_old<-data.frame(sample_id=rownames(x),x,check.names=FALSE)
  assert(same_design(x,fake_old,rownames(x)),"Saved design fixture identity failure")
  fake_old$sample_id<-rev(fake_old$sample_id)
  assert(!same_design(x,fake_old,rownames(x)),"Sample-order rejection failure")
  broken<-x; broken[,4]<-0
  assert(inherits(tryCatch(fit_one(y,broken,grp),error=function(e)e),"error"),
         "Rank-deficient design rejection failure")
  assert(!close_enough(numeric(),numeric()) && is.na(safe_max(numeric())) &&
           all(is.na(unlist(safe_quantiles(numeric())))),"Empty comparisons/summary failure")
  assert(inherits(tryCatch(check_schema(data.frame(a=1),c("a","p"),"p","p"),error=function(e)e),"error"),
         "Missing source column rejection failure")
  fixture_names<-c("OLS_effect_SE_F", "QR_normal_equations", "native_additive", "constant_degenerate",
    "perfect_cell_degenerate", "missing_retained", "infinite_retained", "average_ties",
    "rank_affine_reversal", "joint_order_invariance", "distinct_rank_estimand",
    "saved_design_identity", "sample_order_rejection", "rank_deficiency_rejection",
    "empty_not_estimable", "missing_column_rejection")
  list(status="PASS", fixtures=length(fixture_names), fixture_names=fixture_names,
       method="synthetic only; no real score used", time_utc=utc())
}

safe_cor <- function(a,b) {
  k <- is.finite(a)&is.finite(b)
  if(sum(k)<2 || length(unique(a[k]))<2 || length(unique(b[k]))<2) return(NA_real_)
  cor(a[k],b[k],method="spearman")
}
overlap_summary <- function(p1,p2) {
  k <- is.finite(p1)&is.finite(p2); n <- sum(k)
  if(n==0) return(list(paired_n=0L,rho=NA_real_,top_k=0L,top_left=0L,top_right=0L,top_overlap=0L))
  kk <- min(20L,n); a <- p1[k]<=sort(p1[k])[kk]; b <- p2[k]<=sort(p2[k])[kk]
  list(paired_n=n,rho=safe_cor(p1,p2),top_k=kk,top_left=sum(a),top_right=sum(b),top_overlap=sum(a&b))
}

run_baselines <- function(root,out,input) {
  result_path <- file.path(out,"BASELINE_RESULTS.tsv")
  run_products<-c("RUN_ATTEMPT.json","BASELINE_RESULTS.tsv","SCORE_EXPORT_FOR_VERIFICATION.tsv",
                  "SAMPLE_DESIGN.tsv","SYNTHETIC_TESTS.json","INPUT_MANIFEST_AFTER.tsv","RUN_SUMMARY.json")
  assert(!any(file.exists(file.path(out,run_products))),"Refusing to overwrite any existing run attempt")
  write_json(list(started_utc=utc(),status="started_attempt_preserve_even_if_incomplete"),file.path(out,"RUN_ATTEMPT.json"))
  before <- inventory(root); frozen <- read_tsv(file.path(out,"INPUT_MANIFEST.tsv"))
  assert(same_inventory(before,frozen),"Input pin drift before execution")
  assert(identical(input$features,read_tsv(file.path(out,"FEATURE_REGISTRY.tsv"))),"Feature pin drift")
  tests <- self_tests(); write_json(tests,file.path(out,"SYNTHETIC_TESTS.json"))
  fl <- read_tsv(file.path(root,source_paths()[["observed_t"]]))
  check_schema(fl,c("family","feature_id_original","ordinary_t_observed"),"ordinary_t_observed")
  camera <- read_tsv(file.path(root,source_paths()[["camera"]]))
  check_schema(camera,c("feature_id","contrast","direction","p_value","q_value_BH"),
               c("p_value","q_value_BH"),c("p_value","q_value_BH"))
  assert(all(is.na(camera$direction)|camera$direction %in% c("Up","Down")),"Unknown CAMERA direction")
  camera <- camera[camera$contrast==contrast_name,,drop=FALSE]
  assert(!anyDuplicated(camera$feature_id),"CAMERA primary keys not unique")
  all <- list(); summaries <- list(); exports <- list()
  for (i in seq_along(families)) {
    fam <- families[i]; mat <- input$mats[[i]]; native <- input$natives[[i]]
    rows <- lapply(seq_len(nrow(mat)),function(j)
      as.data.frame(fit_one(mat[j,],input$x,input$groups),stringsAsFactors=FALSE))
    rows <- do.call(rbind,rows); rownames(rows)<-NULL
    rows <- cbind(data.frame(collection=fam,feature_id=rownames(mat)),rows)
    rows$q_value_BH_nominal <- p.adjust(rows$p,method="BH",n=nrow(rows))
    native_fields <- c("estimate","standard_error","ci_low","ci_high","moderated_t","df_total","p_value","q_value_BH")
    for (field in native_fields) rows[[paste0("retained_moderated_",field)]] <- native[[field]]
    ff <- fl[fl$family==fam,,drop=FALSE]
    assert(nrow(ff)==nrow(rows) && !anyDuplicated(ff$feature_id_original) &&
             setequal(ff$feature_id_original,rows$feature_id),"Observed-t key mismatch")
    rows$retained_FL_observed_t <- ff$ordinary_t_observed[match(rows$feature_id,ff$feature_id_original)]
    rows$native_effect_difference <- rows$estimate-native$estimate
    ok <- rows$status=="ok_nominal_working_model"
    finite <- is.finite(rows$estimate)&is.finite(native$estimate)
    fl_ok<-ok&is.finite(rows$retained_FL_observed_t)
    if(any(finite)) assert(close_enough(rows$estimate[finite],native$estimate[finite]),"Stored native effect mismatch")
    if(any(fl_ok)) assert(close_enough(rows$t[fl_ok],rows$retained_FL_observed_t[fl_ok],1e-7),"Observed ordinary-t mismatch")
    if(any(ok)) assert(close_enough(rows$nested_f[ok],rows$t[ok]^2,1e-7),"F=t squared algebra failure")
    rows$retained_effect_check_status<-ifelse(finite,"numeric_compared","unavailable_reference_or_score")
    assert(all(abs(rows$group_effect_check[finite])<=1e-8*(1+abs(rows$estimate[finite]))),
           "Group means contrast mismatch")
    assert(all(abs(rows$se_formula_check[ok])<=1e-8*(1+abs(rows$se[ok]))),"SE formula mismatch")
    rows$no_interaction_fixed_effect <- 0
    rows$rank_inference_status <- "descriptive_only_no_P_q_CI"
    rows$nominal_inference_status <- "conditional_homoskedastic_Gaussian_array_independence_not_verified"
    rows$camera_direction <- NA_character_; rows$camera_p_value <- NA_real_; rows$camera_q_value_BH <- NA_real_
    cm <- NULL
    if(fam=="Reactome") {
      cc <- camera[match(rows$feature_id,camera$feature_id),,drop=FALSE]
      rows$camera_direction<-cc$direction; rows$camera_p_value<-cc$p_value; rows$camera_q_value_BH<-cc$q_value_BH
      camera_sign <- ifelse(cc$direction=="Up",1,ifelse(cc$direction=="Down",-1,NA_real_))
      dir_ok <- is.finite(camera_sign)&is.finite(rows$estimate)&rows$estimate!=0
      cm <- c(overlap_summary(rows$p,cc$p_value),list(direction_paired_n=sum(dir_ok),
                                                   direction_agree_n=sum(camera_sign[dir_ok]==sign(rows$estimate[dir_ok])),
                                                   missing_keys=sum(!rows$feature_id %in% camera$feature_id),
                                                   extra_keys=sum(!camera$feature_id %in% rows$feature_id),
                                                   null="competitive, not same score-level null"))
    }
    directions <- is.finite(rows$estimate)&is.finite(rows$rank_effect)&rows$estimate!=0&rows$rank_effect!=0
    summaries[[fam]] <- list(n=nrow(rows),status_counts=as.list(table(rows$status)),
      native_effect_compared_n=sum(finite),ordinary_t_compared_n=sum(fl_ok),
      native_effect_max_abs_difference=safe_max(abs(rows$native_effect_difference[finite])),
      ordinary_t_max_abs_difference=safe_max(abs(rows$t[fl_ok]-rows$retained_FL_observed_t[fl_ok])),
      f_t_squared_max_abs_difference=safe_max(abs(rows$nested_f[ok]-rows$t[ok]^2)),
      moderated_vs_OLS=overlap_summary(rows$p,native$p_value),
      OLS_q_eligible_n=sum(is.finite(rows$q_value_BH_nominal)),
      retained_q_eligible_n=sum(is.finite(native$q_value_BH)),
      q_common_eligible_n=sum(is.finite(rows$q_value_BH_nominal)&is.finite(native$q_value_BH)),
      OLS_nominal_q05_n=sum(rows$q_value_BH_nominal<=.05,na.rm=TRUE),
      retained_moderated_q05_n=sum(native$q_value_BH<=.05,na.rm=TRUE),
      q05_overlap_n=sum(rows$q_value_BH_nominal<=.05 & native$q_value_BH<=.05,na.rm=TRUE),
      rank_effect_paired_n=sum(directions),rank_effect_direction_agree_n=sum(sign(rows$estimate[directions])==sign(rows$rank_effect[directions])),
      rank_effect_zero_n=sum(is.finite(rows$estimate)&is.finite(rows$rank_effect)&(rows$estimate==0|rows$rank_effect==0)),
      rank_effect_missing_n=sum(!is.finite(rows$estimate)|!is.finite(rows$rank_effect)),
      rank_effect_zero_or_missing_n=sum(!directions),absolute_native_rank_effect_rho=safe_cor(abs(rows$estimate),abs(rows$rank_effect)),
      variance_ratio_quantiles=safe_quantiles(rows$group_variance_ratio),
      max_cook_quantiles=safe_quantiles(rows$max_cook),
      tied_score_feature_n=sum(rows$rank_tied_observations>0,na.rm=TRUE),
      negative_sse_roundoff_n=sum(rows$negative_sse_roundoff),camera_vs_OLS=cm)
    all[[fam]]<-rows
    exports[[fam]]<-cbind(data.frame(collection=fam,feature_id=rownames(mat)),as.data.frame(mat,check.names=FALSE))
  }
  result<-do.call(rbind,all); rownames(result)<-NULL
  write_tsv(result,result_path)
  write_tsv(do.call(rbind,exports),file.path(out,"SCORE_EXPORT_FOR_VERIFICATION.tsv"))
  write_tsv(data.frame(sample_id=input$meta$sample_id,group=input$groups,
    leverage=rowSums((input$x %*% solve(crossprod(input$x)))*input$x)),file.path(out,"SAMPLE_DESIGN.tsv"))
  after <- inventory(root); write_tsv(after,file.path(out,"INPUT_MANIFEST_AFTER.tsv"))
  assert(identical(before,after),"Source mutated during baseline execution")
  output_files<-c("RUN_ATTEMPT.json","BASELINE_RESULTS.tsv","SCORE_EXPORT_FOR_VERIFICATION.tsv","SAMPLE_DESIGN.tsv",
                   "SYNTHETIC_TESTS.json","INPUT_MANIFEST_AFTER.tsv")
  write_json(list(status=if(all(result$status=="ok_nominal_working_model")) "COMPUTED_NUMERIC_CHECKS_PASS" else "PARTIAL_WITH_INVALID_ROWS",
    started_from_registry=sha(file.path(out,"REGISTRY.json")),completed_utc=utc(),
    sample_n=26L,family_n=length(families),features=nrow(result),rank_full=4L,rank_reduced=3L,residual_df=22L,
    sources_unchanged=identical(before,after),synthetic_tests=tests,collections=summaries,
    boundaries=c("Post-outcome detail; no change to primary gates", "Nominal OLS assumptions unresolved",
      "Rank descriptive only, no rank inference", "No independent replication or out-of-sample calibration",
      "Anchor score-level object only; entire SAP and project remain active"),
    outputs=lapply(output_files,function(f) list(file=f,bytes=file.info(file.path(out,f))$size,sha256=sha(file.path(out,f))))),
    file.path(out,"RUN_SUMMARY.json"))
  cat("Computed",nrow(result),"programs; inputs unchanged; see RUN_SUMMARY.json\n")
}

main <- function() {
  for(pkg in c("digest","jsonlite","yaml")) assert(requireNamespace(pkg,quietly=TRUE),paste("Missing existing package",pkg))
  args<-commandArgs(trailingOnly=TRUE)
  assert(length(args)==3L && args[1]%in%c("inspect","run"),"Usage: script inspect|run PROJECT_ROOT OUTPUT_DIR")
  root<-normalizePath(args[2],winslash="/",mustWork=TRUE); out<-args[3]
  if(args[1]=="inspect") {
    assert(!dir.exists(out),"Inspection output directory must be fresh")
    dir.create(out,recursive=TRUE)
    tests<-self_tests(); write_json(tests,file.path(out,"SYNTHETIC_TESTS_PREFIT.json"))
    before<-inventory(root); input<-load_inputs(root)
    write_tsv(before,file.path(out,"INPUT_MANIFEST.tsv"))
    write_tsv(input$features,file.path(out,"FEATURE_REGISTRY.tsv"))
    write_tsv(data.frame(sample_id=input$meta$sample_id,group=input$groups),file.path(out,"SAMPLE_REGISTRY.tsv"))
    write_json(list(schema="anchor_baseline_v1",frozen_utc=utc(),stage="before_new_baseline_fits",
      expected_features=as.list(setNames(expected_counts,families)),samples=26L,
      finite_rows=as.list(vapply(input$mats,function(m)sum(rowSums(!is.finite(m))==0L),integer(1))),
      source_manifest_sha256=sha(file.path(out,"INPUT_MANIFEST.tsv")),
      feature_registry_sha256=sha(file.path(out,"FEATURE_REGISTRY.tsv")),
      sample_registry_sha256=sha(file.path(out,"SAMPLE_REGISTRY.tsv")),
      R_version=R.version.string,packages=as.list(vapply(c("digest","jsonlite","yaml"),function(p)as.character(utils::packageVersion(p)),character(1))),
      synthetic_tests=tests),file.path(out,"REGISTRY.json"))
    assert(identical(before,inventory(root)),"Inspection input drift")
    cat("Input registry fixed; no real baseline effects computed\n")
  } else {
    assert(file.exists(file.path(out,"REGISTRY.json")),"Missing pre-fit registry")
    reg<-jsonlite::read_json(file.path(out,"REGISTRY.json"),simplifyVector=TRUE)
    assert(sha(file.path(out,"INPUT_MANIFEST.tsv"))==reg$source_manifest_sha256 &&
             sha(file.path(out,"FEATURE_REGISTRY.tsv"))==reg$feature_registry_sha256 &&
             sha(file.path(out,"SAMPLE_REGISTRY.tsv"))==reg$sample_registry_sha256,"Registry file drift")
    run_baselines(root,out,load_inputs(root))
  }
}
if(sys.nframe()==0L) main()
