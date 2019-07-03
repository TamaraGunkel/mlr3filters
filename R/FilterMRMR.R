#' @title Minimum redundancy maximal relevancy filter
#'
#' @aliases mlr_filters_MRMR
#' @format [R6::R6Class] inheriting from [Filter].
#' @include Filter.R
#'
#' @description Minimum redundancy maximal relevancy filter. Calls
#'   [praznik::MRMR()].
#'
#' @details This filter supports partial scoring via hyperparameter `k`. To use
#'   it, set `k` during construction via `param_vals.` By default all filter
#'   scores are calculated and the default of `k = 3` in the ParamSet does not
#'   apply.
#'
#' @family Filter
#' @export
#' @examples
#' task = mlr3::mlr_tasks$get("iris")
#' filter = FilterMRMR$new()
#' filter$calculate(task)
#' head(as.data.table(filter), 3)
FilterMRMR = R6Class("FilterMRMR", inherit = Filter,
  public = list(
    initialize = function(id = "MRMR", param_vals = list()) {
      super$initialize(
        id = id,
        packages = "praznik",
        feature_types = c("numeric", "factor", "integer", "character", "logical"),
        task_type = c("classif", "regr"),
        param_set = ParamSet$new(list(
          ParamInt$new("k", lower = 1L, default = 3L, tags = "filter"),
          ParamInt$new("threads", lower = 0L, default = 0L, tags = "filter"),
          ParamInt$new("nfeat", lower = 1, tags = "generic"),
          ParamDbl$new("frac", lower = 0, upper = 1, tags = "generic"),
          ParamDbl$new("cutoff", tags = "generic")
        )),
        param_vals = param_vals
      )
    }
  ),

  private = list(
    .calculate = function(task, n = NULL) {

      # setting params
      k = self$param_set$values$k
      threads = self$param_set$values$threads

      if (is.null(k)) {
        # by default we calculate all scores
        # partial scoring need to be specifically requested by setting k during construction
        k = length(task$feature_names)
      }
      if (is.null(threads)) {
        threads = self$param_set$default$threads
      }

      # n overwrites any param_vals
      if (!is.null(n)) {
        if (!is.null(k)) {
          warningf("Overwriting hyperparameter 'k' with the value given in `$filter_*().")
        }
        k = n
      }

      X = task$data(cols = task$feature_names)
      Y = task$truth()
      praznik::MRMR(X = X, Y = Y, k = k, threads = threads)$score
    }
  )
)