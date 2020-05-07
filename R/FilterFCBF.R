entropy <- function(x, base = exp(1)) {
  if (!is.factor(x)) {
    stop("For calculating the entropy, the vector must be a factor")
  }
  t <- table(x)
  probabily_of_t <- t / sum(t)
  if (any(t == 0)) {
    probabily_of_t <- probabily_of_t[-which(t == 0)]
  }
  ent <- -1 * sum(probabily_of_t * log(probabily_of_t) / log(base))
  if (is.na(ent)) {
    ent <- 0
  }
  ent
}

# H(X,Y) - joint entropy
entropy.joint <- function(x, y, base = exp(1)) {
  if (!is.factor(x) || !is.factor(y)) {
    stop("For calculating the joint entropy, the vector x & y must be factors")
  }
  t <- table(x, y)
  probabily_of_t <- as.numeric(t / sum(t))
  if (any(probabily_of_t == 0)) {
    probabily_of_t <- probabily_of_t[-which(probabily_of_t == 0)]
  }
  ent <- -1 * sum(probabily_of_t * log(probabily_of_t) / log(base))
  if (is.na(ent)) {
    ent <- 0
  }
  ent
}

SU <- function(x, y, base = exp(1)) {
  if (is.character(x)) {
    x <- as.factor(x)
  }
  if (!is.factor(x) || !is.factor(y)) {
    stop(
      "For calculating the symmetrical uncertainty, the vectors x & y must be factors.
      Using a continuous(numeric) feature set leads to this error."
    )
  }
  Ht <- entropy.joint(x, y, base)
  Hx <- entropy(x, base)
  Hy <- entropy(y, base)
  #cat(Ht,' ',Hx,' ',Hy,'\n')

  # Returns the symmetrical uncertainty value for the vector pair
  2 * (Hy + Hx - Ht) / (Hx + Hy)

}

.get.next.elem <- function(s, first_prime) {
  index <- which(s == first_prime)
  if (index == length(s)) {
    NA
  } else {
    s[index + 1]
  }
}

fcbf <-
  function(x,
           y,
           thresh = 0.25,
           n_genes = NULL,
           verbose = FALSE,
           samples_in_rows = FALSE,
           balance_classes = FALSE) {
    if (!samples_in_rows) {
      x <- t(x)
    }
    if (!balance_classes) {
      x <- data.frame(x)
      nvar <- ncol(x)
      if (verbose) {
        message("Calculating symmetrical uncertainties")
      }
      su_ic <- apply(x, 2, function(xx, yy) {
        SU(xx, yy)
      }, y)

      if (length(n_genes)) {
        thresh <- sort(su_ic, decreasing = TRUE)[n_genes - 1]
      }

      s_prime <-
        data.frame(f = (seq_len(nvar))[which(su_ic >= thresh)], su = su_ic[which(su_ic >= thresh)])

      s_prime <- s_prime[sort.list(s_prime$su, decreasing = TRUE), ]

      # s_prime is the list of selected features ranked by su_ic
      s_prime <- s_prime[, 1]

      if (length(s_prime) == 1) {
        s_prime
      } else if (length(s_prime) == 0) {
        stop("No prospective features for this threshold level. Threshold: ",
             thresh)
      }

      print(paste('Number of prospective features = ', length(s_prime)))


      first_prime  <- s_prime[1]
      cnt <- 1
      while (TRUE) {
        if (verbose) {
          cat("Round ")
          cat(cnt, "\n")
          cnt <- cnt + 1
          print(paste(
            'first_prime  round ( |s_prime| =',
            length(s_prime),
            ')' ,
            sep =
              ' '
          ))
        }

        next_prime <- .get.next.elem(s_prime, first_prime)
        if (!is.na(next_prime)) {
          while (TRUE) {
            prime_to_be_compared <- next_prime
            su1 = SU(x[, first_prime], x[, next_prime])
            su2 = SU(x[, next_prime], y)
            if (su1 > su2) {
              next_prime <- .get.next.elem(s_prime, next_prime)
              s_prime <-
                s_prime[-which(s_prime == prime_to_be_compared)]
              if (verbose) {
                cat("  ",
                    su1,
                    " ",
                    su2,
                    " ",
                    "Removed feature ",
                    prime_to_be_compared,
                    "\n")
              }
            }
            else {
              next_prime <- .get.next.elem(s_prime, next_prime)
            }
            if (is.na(next_prime)) {
              break
            }
          }
        }
        first_prime  <- .get.next.elem(s_prime, first_prime)

        if (is.na(first_prime)) {
          break
        }
      }
      if (length(s_prime) > 1) {
        suvalues <- apply(x[, s_prime], 2, function(xx, yy) {
          SU(xx, yy)
        }, y)
        data.frame(index = s_prime, SU = suvalues)
      } else {
        data.frame(index = s_prime, SU = SU(x[, s_prime], y))
      }
    }
    else{
      instances_in_minor_class <- min(table(y))
      n_x <- as.data.frame(cbind(x, y))
      final_x <- data.frame(n_x[1, ])
      final_x <- final_x[-1, ]
      for (i in levels(as.factor(n_x[, ncol(n_x)]))) {
        n_x_i <- n_x[n_x[, ncol(n_x)] == i, ]
        n_x_i <-
          n_x_i[sample(seq_len(length.out = nrow(n_x_i)), instances_in_minor_class), ]
        final_x <- rbind(n_x_i, final_x)
      }
      final_x$y <- NULL
      fcbf(t(final_x),
           y,
           thresh,
           verbose,
           samples_in_rows,
           balance_classes = FALSE)

    }
  }


#' @title Fast Correlation Based Filter
#'
#' @name mlr_filters_fcbf
#'
#' @description Fast Correlation Based Filter  (embedded implementation taken from FCBF (package on Bioconductor)).
#'
#' @family Filter
#' @template seealso_filter
#' @export
#' @examples
#' task = mlr3::tsk("iris")
#' filter = flt("fcbf")
#' filter$calculate(task, nfeat = 2)
#' as.data.table(filter)
FilterFCBF = R6Class("FilterFCBF", inherit = Filter,

                     public = list(

                       #' @description Create a FilterFCBF object.
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
                       initialize = function(id = "fcbf",
                                             task_type = c("classif"),
                                             param_set = ParamSet$new(list(
                                               ParamLgl$new("verbose", default = FALSE),
                                               ParamLgl$new("balance_classes", default = FALSE)
                                             )),
                                             packages = "FSelectorRcpp",
                                             feature_types = c("integer", "numeric", "factor", "ordered")) {
                         super$initialize(
                           id = id,
                           task_type = task_type,
                           param_set = param_set,
                           feature_types = feature_types,
                           packages = packages,
                           man = "mlr3filters::mlr_filters_fcbf"
                         )
                       }
                     ),

                     private = list(

                       .calculate = function(task, nfeat) {
                         verbose = self$param_set$values$verbose %??% FALSE
                         balance_classes = self$param_set$values$balance_classes %??% FALSE
                         X = task$data(cols = task$feature_names)
                         Y = task$truth()
                         disc_data =  FSelectorRcpp::discretize(X, Y)
                         Y = disc_data$Y
                         disc_data$Y = NULL
                         result = fcbf(disc_data, Y, thresh = 0, n_genes = NULL, verbose = verbose,
                              samples_in_rows = TRUE, balance_classes)
                         new_result = result$SU
                         names(new_result) = colnames(disc_data)[result$index]
                         new_result
                       }
                     )
)

#' @include mlr_filters.R
mlr_filters$add("fcbf", FilterFCBF)
