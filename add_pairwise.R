
library(gt)
library(gtsummary)
library(tidyverse)
library(emmeans)


## remove prop... keep wilcox and fisher
## additional options are anova and chi square

### new version
add_pairwise <- function(
    x,
    include = dplyr::everything(),
    p.adjust.method = c("holm", "none", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr"),
    continuous_method = c("wilcox.test", "lm_emmeans"),
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
  # col_names <- vapply(
  #   pairs,
  #   function(p) make.names(glue::glue("p_{p[1]}_vs_{p[2]}")),
  #   character(1)
  # )
  # header_map <- rlang::set_names(
  #   as.list(vapply(pairs, function(p) glue::glue(header_fmt, level1 = p[1], level2 = p[2]), character(1))),
  #   col_names
  # )
  
  # --- build pair list first ---
  pairs <- utils::combn(lvls, 2, simplify = FALSE)
  
  # --- internal (syntactic) column names ---
  col_names <- vapply(
    pairs,
    function(p) make.names(glue::glue("p_{p[1]}_vs_{p[2]}")),
    character(1)
  )
  
  # --- display headers (pretty) shown in the table ---
  header_map <- rlang::set_names(
    as.list(vapply(pairs, function(p) glue::glue("{p[1]} vs. {p[2]}"), character(1))),
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
    
    # ---- Optional precompute for lm + emmeans (continuous only) ----
    # We'll compute all pairwise Tukey-adjusted p-values once per variable, then look up.
    pw_lookup <- NULL
    if (is_cont && continuous_method == "lm_emmeans") {
      
      if (!requireNamespace("emmeans", quietly = TRUE)) {
        stop("continuous_method = 'lm_emmeans' requires the emmeans package.")
      }
      
      fit <- stats::lm(y ~ g)
      
      # pw <- emmeans::emmeans(fit, ~ g) |>
      #   emmeans::pairs(adjust = "tukey") |>
      #   as.data.frame()
      
      pw <- emmeans::contrast(emmeans::emmeans(fit, ~ g), method = "pairwise", adjust = "tukey") |>
        as.data.frame()
      
      # Build lookup
      parts <- strsplit(as.character(pw$contrast), " - ", fixed = TRUE)
      keys  <- vapply(parts, function(z) paste(sort(z), collapse = "||"), character(1))
      pw_lookup <- stats::setNames(as.numeric(pw$p.value), keys)
    }
    
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
          if (continuous_method == "wilcox.test") {
            stats::wilcox.test(y2 ~ g2)$p.value
            
          } else if (continuous_method == "lm_emmeans") {
            if (is.null(pw_lookup)) return(NA_real_)
            key <- paste(sort(c(a, b)), collapse = "||")
            unname(pw_lookup[[key]])
            
          } else {
            NA_real_
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
    
    # Adjust within-variable across its pairs if requested:
    # - If continuous uses lm_emmeans, Tukey is already applied, so DO NOT re-adjust.
    # - Otherwise (wilcox or categorical), apply p.adjust.method if requested.
    if (p.adjust.method != "none") {
      if (!(is_cont && continuous_method == "lm_emmeans")) {
        pvals <- stats::p.adjust(pvals, method = p.adjust.method)
      }
    }
    
    out <- as.list(pvals)
    names(out) <- col_names
    out_tbl <- dplyr::as_tibble(out)
    
    # debug:
    message("pvals for variable = ", variable, ": ",
            paste(sprintf("%s=%s", names(out), unlist(out)), collapse = ", "))
    
    out_tbl
  }
  
  # Add all pairwise columns
  x2 <- gtsummary::add_stat(
    x = x,
    fns = rlang::new_formula(rlang::enquo(include), rlang::expr(gts_pairwise_generic)),
    location = gtsummary::all_stat_cols() ~ "label"
  )
  
  # headers + p-value formatting
  x2 <- do.call(
    gtsummary::modify_header,
    c(list(x2), header_map)
  ) %>%
    # spanning header over the pairwise columns
    gtsummary::modify_spanning_header(
      dplyr::all_of(col_names) ~ "**Pairwise comparisons (p-value)**"
    ) %>%
    gtsummary::modify_fmt_fun(dplyr::all_of(col_names) ~ pvalue_fun)
  
  # ---- build footnote text ----
  cont_label <- if (continuous_method == "wilcox.test") {
    "Wilcoxon rank-sum test"
  } else {
    "linear model with Tukey-adjusted pairwise comparisons"
  }
  
  cat_label <- if (categorical_method == "fisher") {
    "Fisher's exact test"
  } else {
    "Chi-squared test"
  }
  
  # If lm_emmeans is used, don't apply p.adjust.method to continuous p-values.
  adj_clause <- ""
  if (p.adjust.method != "none") {
    if (continuous_method == "lm_emmeans") {
      adj_clause <- paste0("; ", p.adjust.method, " p-value adjustment applied for categorical variables")
    } else {
      adj_clause <- paste0("; p-values adjusted using ", p.adjust.method, " method")
    }
  }
  
  footnote_text <- paste0(
    "Pairwise comparisons performed using ",
    cont_label,
    " for continuous variables and ",
    cat_label,
    " for categorical variables",
    adj_clause,
    "."
  )
  
  # ---- attach footnote to all pairwise columns ----
  x2 <- x2 %>%
    gtsummary::modify_footnote(
      dplyr::all_of(col_names) ~ footnote_text
    )
  
  x2
}


## test
trial$death <- factor(trial$death)

trial %>%
  select(age, marker, grade, stage, death) %>%
  # create summary table split by Grade
  tbl_summary(
    by = stage,
    missing = "no",
    statistic = list(
      gtsummary::all_continuous() ~ "{mean} ({sd})",
      gtsummary::all_categorical() ~ "{n} ({p}%)"
    )
  ) %>% add_pairwise(continuous_method = "lm_emmeans",
                     categorical_method = "fisher",
                     p.adjust.method = "holm")



