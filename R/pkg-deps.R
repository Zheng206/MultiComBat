# R/pkg-deps.R

#' Package imports
#'
#' Centralized namespace imports for this package.
#'
#' @keywords internal
#'
#' @importFrom mgcv gam anova.gam
#' @importFrom lme4 lmer
#' @importFrom MASS ginv mvrnorm
#' @importFrom stats lm median model.matrix predict update var qnorm manova anova rt
#'    fligner.test bartlett.test kruskal.test aov as.formula coef p.adjust resid residuals
#'    rnorm rbinom runif
#' @importFrom dplyr bind_rows mutate select arrange filter group_by group_split pull
#'    distinct left_join across everything
#' @importFrom magrittr %>%
#' @importFrom utils globalVariables
#' @importFrom purrr map_dfr map
#' @importFrom tidyr pivot_wider extract pivot_longer
#' @importFrom broom tidy
#' @importFrom car leveneTest
#' @importFrom ggplot2 ggplot aes geom_tile geom_text scale_fill_gradientn geom_density facet_wrap
#'    labs theme_minimal theme element_text element_blank geom_histogram scale_fill_manual after_stat
#'
#' @importFrom MCMCpack rinvgamma riwish
#'
#'
#' @noRd
NULL
