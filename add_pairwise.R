
library(gt)
library(gtsummary)
library(tidyverse)

#### from github issue (from author of package)
# works, but only if all variables are continous (or could edit the fns value in add_stat if all categorical)
# define a pairwise test in the gtsummary format to calculate pairwise p-values
gts_pairwise.t.test <- function(data, variable, by, ...) {
  pairwise.t.test(data[[variable]], data[[by]]) %>% 
    broom::tidy() %>%
    mutate(label = glue::glue("**{group2} vs. {group1}**")) %>%
    select(label, p.value) %>%
    spread(label, p.value)
}

# test function
gts_pairwise.t.test(trial, "age", "grade")
#> # A tibble: 1 x 3
#>   `**I vs. II**` `**I vs. III**` `**II vs. III**`
#>            <dbl>           <dbl>            <dbl>
#> 1              1               1                1

trial %>%
  select(age, marker, grade) %>%
  # create summary table split by Grade
  tbl_summary(
    by = grade,
    missing = "no",
    statistic = everything() ~ "{mean} ({sd})"
  ) %>%
  # add pairwise comparisons
  add_stat(fns = everything() ~ gts_pairwise.t.test) %>%
  modify_fmt_fun(c('**I vs. II**', '**I vs. III**', '**II vs. III**') ~ style_pvalue)




### new
add_pairwise <- function(
    x,
    include = dplyr::everything(),
    p.adjust.method = c("none", "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr"),
    continuous_method = c("t.test", "wilcox.test"),
    categorical_method = c("chisq", "fisher"),   
    header_fmt = "**{level1} vs. {level2}**",
    pvalue_fun = gtsummary::style_pvalue
) {
  p.adjust.method <- match.arg(p.adjust.method)
  continuous_method <- match.arg(continuous_method)
  categorical_method <- match.arg(categorical_method)
  
  if (!inherits(x, "tbl_summary")) {
    stop("`add_pairwise()` expects a gtsummary::tbl_summary() object.")
  }
  by_var <- x$inputs$by
  if (is.null(by_var) || !nzchar(by_var)) {
    stop("`add_pairwise()` requires tbl_summary(by = ...) to be set.")
  }
  
  dat <- x$inputs$data
  g <- droplevels(as.factor(dat[[by_var]]))
  lvls <- levels(g)
  
  if (length(lvls) < 2) return(x)
  
  pairs <- utils::combn(lvls, 2, simplify = FALSE)
  
  # Canonical column names (syntactic) + display headers (pretty)
  col_names <- vapply(
    pairs,
    function(p) make.names(glue::glue("p_{p[1]}_vs_{p[2]}")),
    character(1)
  )
  header_map <- rlang::set_names(
    as.list(vapply(pairs, function(p) glue::glue(header_fmt, level1 = p[1], level2 = p[2]), character(1))),
    col_names
  )
  
  # ---- helper: compute pairwise p-values for ONE variable, returning a 1-row tibble
  gts_pairwise_generic <- function(data, variable, by, ...) {
    y <- data[[variable]]
    g <- data[[by]]
    
    ok <- stats::complete.cases(y, g)
    y <- y[ok]
    g <- droplevels(as.factor(g[ok]))
    
    # Always return 1 row with the full set of columns
    out <- as.list(rep(NA_real_, length(pairs)))
    names(out) <- col_names
    
    if (nlevels(g) < 2) return(dplyr::as_tibble(out))
    
    is_cont <- is.numeric(y) || is.integer(y)
    
    pvals <- rep(NA_real_, length(pairs))
    
    for (i in seq_along(pairs)) {
      a <- pairs[[i]][1]
      b <- pairs[[i]][2]
      
      keep <- g %in% c(a, b)
      y2 <- y[keep]
      g2 <- droplevels(g[keep])
      
      if (nlevels(g2) < 2) {
        pvals[i] <- NA_real_
        next
      }
      
      if (is_cont) {
        pvals[i] <- tryCatch({
          if (continuous_method == "t.test") {
            stats::t.test(y2 ~ g2)$p.value
          } else {
            stats::wilcox.test(y2 ~ g2)$p.value
          }
        }, error = function(e) NA_real_)
      } else {
        tab <- table(y2, g2)
        # Need at least 2 columns and at least 2 rows in the outcome to test
        if (ncol(tab) < 2 || nrow(tab) < 2) {
          pvals[i] <- NA_real_
          next
        }
        
        pvals[i] <- tryCatch({
          if (categorical_method == "fisher") {
            stats::fisher.test(tab)$p.value
          } else {
            suppressWarnings(stats::chisq.test(tab, correct = FALSE)$p.value)
          }
        }, error = function(e) NA_real_)
      }
    }
    
    # adjust within-variable across its pairs if requested
    if (p.adjust.method != "none") {
      pvals <- stats::p.adjust(pvals, method = p.adjust.method)
    }
    
    out <- as.list(pvals)
    names(out) <- col_names
    dplyr::as_tibble(out)
  }
  
  # Add all pairwise columns (like in Daniel Sjoberg’s github example)
  x2 <- gtsummary::add_stat(
    x = x,
    fns = rlang::new_formula(rlang::enquo(include), rlang::expr(gts_pairwise_generic)),
    location = gtsummary::all_stat_cols() ~ "label"
  )
  
  # headers + p-value formatting
  x2 <-
    x2 %>%
    #gtsummary::modify_header(rlang::!!!header_map) %>%
    gtsummary::modify_fmt_fun(dplyr::all_of(col_names) ~ pvalue_fun)
  
  # ---- build footnote text ----
  cont_label <- if (continuous_method == "t.test") {
    "t-test"
  } else {
    "Wilcoxon rank-sum test"
  }
  
  cat_label <- if (categorical_method == "fisher") {
    "Fisher's exact test"
  } else {
    "Chi-squared test"
  }
  
  footnote_text <- paste0(
    "Pairwise comparisons performed using ",
    cont_label,
    " for continuous variables and ",
    cat_label,
    " for categorical variables",
    if (p.adjust.method != "none") {
      paste0("; p-values adjusted using ", p.adjust.method, " method")
    } else {
      ""
    },
    "."
  )
  
  # ---- attach footnote to all pairwise columns ----
  x2 <- x2 %>%
    gtsummary::modify_footnote(
      dplyr::all_of(col_names) ~ footnote_text
    )
  
  x2
}





## test new
trial %>%
  select(age, marker, grade, stage) %>%
  # create summary table split by Grade
  tbl_summary(
    by = stage,
    missing = "no",
    statistic = list(
      gtsummary::all_continuous() ~ "{mean} ({sd})",
      gtsummary::all_categorical() ~ "{n} ({p}%)"
    )
  ) %>% add_pairwise(p.adjust.method = "fdr")





















