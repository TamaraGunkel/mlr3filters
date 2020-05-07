#' @title Boruta Feature Selection Wrapper
#'
#' @name mlr_filters_jmim
#'
#' @description Boruta Feature Selection Wrapper calling
#' [Boruta::Boruta()] in package \CRANpkg{Boruta}.
#' Results: 2 - Confirmed, 1 - Tentative, 0 - Rejected
#'
#' @family Filter
#' @template seealso_filter
#' @export
#' @examples
#' task = mlr3::tsk("iris")
#' filter = flt("boruta")
#' filter$calculate(task, nfeat = 2)
#' as.data.table(filter)
FilterBoruta = R6Class("FilterBoruta", inherit = Filter,

                     public = list(

                       #' @description Create a FilterJMIM object.
                       #' @param id (`character(1)`)\cr
                       #'   Identifier for the filter.
                       #' @param task_type (`character()`)\cr
                       #'   Types of the task the filter can operator on. E.g., `"classif"` or
                       #'   `"regr"`.
                       #' @param param_set ([paradox::ParamSet])\cr
                       #'   Set of hyperparameters.
                       #' @param feature_types (`character()`)\cr
                       #'   Feature types the filter operates on.
                       #'   Must be a subset of
                       #'   [`mlr_reflections$task_feature_types`][mlr3::mlr_reflections].
                       #' @param packages (`character()`)\cr
                       #'   Set of required packages.
                       #'   Note that these packages will be loaded via [requireNamespace()], and
                       #'   are not attached.
                       initialize = function(id = "boruta",
                                             task_type = c("classif", "regr"),
                                             param_set = ParamSet$new(list(
                                               ParamDbl$new("pValue", default = 0.01),
                                               ParamLgl$new("mcAdj", default = TRUE),
                                               ParamInt$new("maxRuns", default = 100),
                                               ParamInt$new("doTrace", default = 0),
                                               ParamLgl$new("holdHistory", default = TRUE)
                                             )),
                                             packages = "Boruta",
                                             feature_types = c("integer", "numeric", "factor", "ordered")) {
                         super$initialize(
                           id = id,
                           task_type = task_type,
                           param_set = param_set,
                           feature_types = feature_types,
                           packages = packages,
                           man = "mlr3filters::mlr_filters_boruta"
                         )
                       }
                     ),

                     private = list(

                       .calculate = function(task, nfeat) {
                         pValue = self$param_set$values$pValue %??% 0.01
                         mcAdj = self$param_set$values$mcAdj %??% TRUE
                         maxRuns = self$param_set$values$maxRuns %??% 100
                         doTrace = self$param_set$values$doTrace %??% 0
                         holdHistory = self$param_set$values$holdHistory %??% TRUE
                         X = task$data(cols = task$feature_names)
                         Y = task$truth()
                         results = Boruta::Boruta(x = X, y = Y, pValue = pValue, mcAdj = mcAdj,
                                                  maxRuns = maxRuns, doTrace = doTrace, holdHistory = holdHistory)$finalDecision
                         transformed_results = as.numeric(results)
                         names(transformed_results) = names(results)
                         transformed_results
                       }
                     )
)

#' @include mlr_filters.R
mlr_filters$add("boruta", FilterBoruta)
