#' Disaggregate binned hazards to single year of age and single calendar year
#'
#' @description
#' Expand a period-by-agegroup hazard array from \code{\link{predict_apc}} (or
#' \code{\link{predict_strata}}) onto a single-year-of-age by single-calendar-year grid, as
#' required by an individual-level microsimulation such as IMPACTncd. bamp models age
#' \emph{groups} and (possibly multi-year) \emph{periods}; a microsimulation that advances one
#' year of age per one-year cycle needs a hazard for every single year of age and calendar year.
#'
#' @param x a hazard array whose \strong{first dimension is period} and \strong{second dimension is
#'   age group}, optionally with a third dimension (e.g. MCMC samples) that is carried through
#'   unchanged. This matches \code{predict_apc(..., hazard=TRUE)$samples$hazard} and the per-stratum
#'   \code{$samples$hazard} from \code{\link{predict_strata}}. A plain \code{[period, agegroup]}
#'   matrix is also accepted. (For the quantile array \code{[quantile, period, agegroup]}, permute
#'   it to put period first, e.g. \code{aperm(h, c(2, 3, 1))}.)
#' @param agegroup_width single positive integer, or an integer vector of length equal to the number
#'   of age groups, giving how many single years of age each group spans (e.g. 5 for five-year age
#'   groups).
#' @param period_width single positive integer: how many calendar years each period spans (often the
#'   same number you passed as \code{period_length} to \code{predict_apc}). Default 1.
#' @param method disaggregation method. Currently only \code{"constant"} (piecewise constant): the
#'   group's per-person-year hazard is assigned to each single year of age it covers and each
#'   calendar year of the period.
#' @param start_age,start_period optional integers labelling the first single year of age and the
#'   first calendar year, used only to set \code{dimnames} on the result.
#'
#' @return an array with the period dimension expanded to \code{period_width} times as many single
#'   calendar years, the age-group dimension expanded to \code{sum(agegroup_width)} single years of
#'   age, and any further dimensions unchanged. Dimension order is
#'   \code{[year, age, ...]}.
#'
#' @details
#' Piecewise-constant disaggregation treats the per-person-year hazard as constant within each age
#' group and within each period. This is the standard, conservative default and is exact when the
#' underlying hazard is flat across the bin; where the age gradient is steep a smoother
#' interpolation may be preferable (not yet implemented). The hazards must already be expressed per
#' the desired time unit -- i.e. produced with \code{predict_apc(..., hazard=TRUE, period_length=)}
#' where \code{period_length} is the period width in years -- so that replicating a value across the
#' single years of a period is correct rather than double-counting.
#'
#' @seealso \code{\link{predict_apc}}, \code{\link{predict_strata}}
#' @examples
#' \dontrun{
#' data(apc)
#' m <- bamp(cases, population, age="rw1", period="rw1", cohort="rw1", periods_per_agegroup=5)
#' pred <- predict_apc(m, periods=2, hazard=TRUE, period_length=5)
#' # five-year age groups, five-year periods -> single year of age and calendar year
#' h_yr <- disaggregate_hazard(pred$samples$hazard, agegroup_width=5, period_width=5)
#' dim(h_yr)   # [years, single-ages, draws]
#' }
#' @export
disaggregate_hazard <- function(x, agegroup_width, period_width = 1L,
                                method = c("constant"),
                                start_age = NULL, start_period = NULL)
{
  method <- match.arg(method)
  d <- dim(x)
  if (is.null(d) || length(d) < 2L || length(d) > 3L)
    stop("'x' must be a [period, agegroup] matrix or a [period, agegroup, sample] array.")
  P <- d[1]; A <- d[2]

  if (length(agegroup_width) == 1L) agegroup_width <- rep(agegroup_width, A)
  agw <- as.integer(round(agegroup_width))
  if (length(agw) != A || any(agw < 1L))
    stop("'agegroup_width' must be one positive integer, or one per age group (", A, ").")
  pw <- as.integer(round(period_width))
  if (length(pw) != 1L || pw < 1L)
    stop("'period_width' must be a single positive integer.")

  pidx <- rep(seq_len(P), each = pw)       # calendar year -> period
  aidx <- rep(seq_len(A), times = agw)     # single year of age -> age group

  out <- if (length(d) == 2L) x[pidx, aidx, drop = FALSE]
         else                 x[pidx, aidx, , drop = FALSE]

  yr_names  <- if (!is.null(start_period)) as.character(seq.int(start_period, length.out = length(pidx)))
  age_names <- if (!is.null(start_age))    as.character(seq.int(start_age,    length.out = length(aidx)))
  if (!is.null(yr_names) || !is.null(age_names)) {
    dn <- vector("list", length(d)); dn[[1]] <- yr_names; dn[[2]] <- age_names
    dimnames(out) <- dn
  }
  out
}
