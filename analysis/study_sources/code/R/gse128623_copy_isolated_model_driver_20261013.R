# Compatibility factory only; does not execute or authorize a model.
# All v4 code is preserved except its single validator namespace assignment.
g128_copy_isolated_model_driver <- function(driver_raw, helper_raw) {
  need <- function(ok, code) {
    if (!isTRUE(ok)) stop(paste0("COPY_ISOLATED_DRIVER:", code), call. = FALSE)
  }
  pins <- list(
    driver = list(bytes = 13656L,
      sha256 = "b4151ad1c11b07d54d805eb8f94fe336250fba529f619127d90114e5e8cf4181"),
    helper = list(bytes = 2127L,
      sha256 = "a84f8f660d5bf2e3633732eb7f88a3fd992172990eff942edd6a2080a2d48205")
  )
  authenticate <- function(value, pin) {
    need(is.raw(value) && is.null(attributes(value)) && length(value) == pin$bytes &&
           identical(tolower(as.vector(as.character(openssl::sha256(value)), "character")),
                     pin$sha256), "SOURCE_IDENTITY")
  }
  authenticate(driver_raw, pins$driver)
  authenticate(helper_raw, pins$helper)
  driver_nodes <- parse(text = rawToChar(driver_raw), keep.source = FALSE)
  helper_nodes <- parse(text = rawToChar(helper_raw), keep.source = FALSE)
  assignment <- function(nodes, name) {
    length(nodes) == 1L && is.call(nodes[[1L]]) && length(nodes[[1L]]) == 3L &&
      identical(nodes[[1L]][[1L]], as.name("<-")) &&
      identical(nodes[[1L]][[2L]], as.name(name)) &&
      is.call(nodes[[1L]][[3L]]) &&
      identical(nodes[[1L]][[3L]][[1L]], as.name("function"))
  }
  need(assignment(driver_nodes, "g128_model_execution_driver") &&
         assignment(helper_nodes, "g128_nonmutating_validator"), "SOURCE_AST")
  original_module <- new.env(parent = baseenv())
  eval(driver_nodes, original_module)
  original <- original_module$g128_model_execution_driver
  lockEnvironment(original_module, bindings = TRUE)
  target <- quote(validator <- kernels$validator)
  replacement <- quote(validator <- g128_nonmutating_validator(kernels$validator))
  replace_assignment <- function(node, from, to) {
    if (identical(node, from)) return(list(node = to, matches = 1L))
    if (!is.call(node)) return(list(node = node, matches = 0L))
    children <- lapply(as.list(node), replace_assignment, from = from, to = to)
    list(node = as.call(lapply(children, `[[`, "node")),
         matches = sum(vapply(children, `[[`, 0L, "matches")))
  }
  changed <- replace_assignment(body(original), target, replacement)
  restored <- replace_assignment(changed$node, replacement, target)
  need(changed$matches == 1L && restored$matches == 1L &&
         identical(restored$node, body(original)), "EXACTLY_ONE_ASSIGNMENT_CHANGE")
  compatible_module <- new.env(parent = baseenv())
  eval(helper_nodes, compatible_module)
  compatible <- as.function(c(as.list(formals(original)), list(changed$node)),
                            envir = compatible_module)
  need(identical(formals(compatible), formals(original)), "FORMALS_UNCHANGED")
  lockEnvironment(compatible_module, bindings = TRUE)
  function(analyzer_raw, validator_raw, loader_raw, member_guard_raw,
           comparison_guard_raw, paths) {
    result <- compatible(analyzer_raw, validator_raw, loader_raw, member_guard_raw,
                         comparison_guard_raw, paths)
    # The checks and statistical functions returned through a compatibility
    # execution. It is not a success of the unmodified original driver.
    result$status <- "COPY_ISOLATED_MODEL_SEQUENCE_RETURNED_CONTENT_VALIDATION_ONLY"
    result$compatibility_schema <- "gse128623-copy-isolated-v4-driver-v1"
    result$original_unmodified_driver_returned <- FALSE
    result$validator_assignment_changes <- 1L
    result$original_driver_source <- pins$driver
    result$comparison_helper_source <- pins$helper
    result$corrected_verification_lineage <-
      "Authenticated v4 driver cloned in RAM; only validator namespace assignment invokes the authenticated copy-isolation helper. Statistical bodies and all original driver gates are retained."
    result
  }
}
