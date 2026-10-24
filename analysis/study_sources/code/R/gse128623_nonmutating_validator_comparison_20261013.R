# Compatibility only: keep frozen statistical bodies, isolate comparison sorting.
# Caller authenticates the original kernel and this helper separately. This is
# not an execution authorization or an alternate statistical method.
g128_nonmutating_validator <- function(original_validator) {
  stopifnot(is.environment(original_validator),
            !identical(original_validator, globalenv()),
            exists("v2v_compare_table", original_validator, inherits = FALSE))
  original_compare <- get("v2v_compare_table", original_validator, inherits = FALSE)
  stopifnot(is.function(original_compare),
            identical(names(formals(original_compare)),
                      c("path", "expected", "keys", "numeric", "exact", "label")))

  # Rebind every closure: branch functions must resolve the isolated comparator
  # rather than retaining their original enclosing namespace. Private assignments
  # below never modify the original namespace or the global environment.
  compatible_validator <- new.env(parent = parent.env(original_validator))
  bindings <- ls(original_validator, all.names = TRUE)
  for (name in bindings) {
    value <- get(name, original_validator, inherits = FALSE)
    if (is.function(value)) {
      cloned <- as.function(c(as.list(formals(value)), list(body(value))),
                            envir = compatible_validator)
      attributes(cloned) <- attributes(value)
      assign(name, cloned, envir = compatible_validator)
    } else {
      assign(name, data.table::copy(value), envir = compatible_validator)
    }
  }

  # Original checks, tolerance and missingness semantics are unchanged. Only the
  # object exposed to by-reference sorting is a private copy. In particular, the
  # caller's edge table keeps its gene order for positional symbol-score pairing.
  compatible_validator$v2v_compare_table <- function(path, expected, keys,
                                                     numeric, exact, label) {
    original_compare(path, data.table::copy(expected), keys, numeric, exact, label)
  }
  lockEnvironment(compatible_validator, bindings = TRUE)
  compatible_validator
}
