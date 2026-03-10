#' prplot
#'
#' Build a partial residual / component-plus-residual plot for a fitted model
#' using raw ggplot2.
#'
#' For numeric predictors, the plotted y-values are:
#'   residuals(mdl) + term_contribution(xvariable)
#' where term_contribution is taken from predict(type = "terms").
#'
#' For simple linear-model main effects, this is the standard
#' component-plus-residual construction.
#'
#' @param mdl fitted model object supporting model.frame(), residuals(),
#'   predict(type = "terms"), coef(), vcov().
#' @param xvariable character scalar naming the predictor to plot.
#' @param byvariable optional character scalar for faceting.
#' @param titlestring plot title.
#' @param ystring y-axis label.
#' @param addpoints numeric point size; if <= 0, omit points.
#' @param palette kept for backward compatibility; not used directly.
#' @param colorvar optional character scalar naming a grouping variable for point color.
#' @param extradata optional data.frame aligned row-wise with model.frame(mdl).
#' @param fit_line one of c("model", "global_lm", "group_lm", "none").
#' @param show_annotation logical; if TRUE, annotate effect summary and p-value.
#' @param annotation_effect one of c("standardized_beta", "raw", "partial_r2", "none").
#' @param annotation_digits integer digits for displayed statistics.
#' @param annotation_x numeric in \code{[0, 1]}, fractional x placement within data range.
#' @param annotation_y numeric in \code{[0, 1]}, fractional y placement within data range.
#'
#' @return ggplot object.
#' @export
prplot <- function(
  mdl,
  xvariable,
  byvariable = NULL,
  titlestring = "",
  ystring = "Partial residual",
  addpoints = 0,
  palette = "npg",
  colorvar = NULL,
  extradata = NULL,
  fit_line = c("model", "global_lm", "group_lm", "none"),
  show_annotation = TRUE,
  annotation_effect = c("standardized_beta", "raw", "partial_r2", "none"),
  annotation_digits = 3,
  annotation_x = 0.02,
  annotation_y = 0.98
) {
  fit_line <- match.arg(fit_line)
  annotation_effect <- match.arg(annotation_effect)

  .prplot_validate_inputs(
    mdl = mdl,
    xvariable = xvariable,
    byvariable = byvariable,
    titlestring = titlestring,
    ystring = ystring,
    addpoints = addpoints,
    palette = palette,
    colorvar = colorvar,
    extradata = extradata,
    fit_line = fit_line,
    show_annotation = show_annotation,
    annotation_effect = annotation_effect,
    annotation_digits = annotation_digits,
    annotation_x = annotation_x,
    annotation_y = annotation_y
  )

  model_data <- stats::model.frame(mdl)
  .prplot_validate_model_variables(
    model_data = model_data,
    xvariable = xvariable,
    byvariable = byvariable,
    colorvar = colorvar,
    extradata = extradata
  )

  plot_data <- .prplot_build_plot_data(
    mdl = mdl,
    model_data = model_data,
    xvariable = xvariable,
    byvariable = byvariable,
    colorvar = colorvar,
    extradata = extradata
  )

  color_name <- .prplot_resolve_color_mapping(model_data, colorvar, extradata)
  annotation_text <- if (show_annotation) {
    .prplot_term_annotation(
      mdl = mdl,
      model_data = model_data,
      xvariable = xvariable,
      effect_metric = annotation_effect,
      digits = annotation_digits
    )
  } else {
    NULL
  }

  if (.prplot_is_discrete(plot_data[[xvariable]])) {
    return(
      .prplot_build_discrete_plot(
        plot_data = plot_data,
        xvariable = xvariable,
        byvariable = byvariable,
        titlestring = titlestring,
        ystring = ystring,
        addpoints = addpoints,
        colorvar = color_name,
        annotation_text = annotation_text,
        annotation_x = annotation_x,
        annotation_y = annotation_y
      )
    )
  }

  .prplot_build_continuous_plot(
    mdl = mdl,
    plot_data = plot_data,
    xvariable = xvariable,
    byvariable = byvariable,
    titlestring = titlestring,
    ystring = ystring,
    addpoints = addpoints,
    colorvar = color_name,
    fit_line = fit_line,
    annotation_text = annotation_text,
    annotation_x = annotation_x,
    annotation_y = annotation_y
  )
}

.prplot_validate_inputs <- function(
  mdl,
  xvariable,
  byvariable,
  titlestring,
  ystring,
  addpoints,
  palette,
  colorvar,
  extradata,
  fit_line,
  show_annotation,
  annotation_effect,
  annotation_digits,
  annotation_x,
  annotation_y
) {
  if (missing(mdl) || is.null(mdl)) stop("`mdl` must be a fitted model.", call. = FALSE)
  if (missing(xvariable) || !is.character(xvariable) || length(xvariable) != 1L) {
    stop("`xvariable` must be a single character string.", call. = FALSE)
  }
  if (!is.null(byvariable) && (!is.character(byvariable) || length(byvariable) != 1L)) {
    stop("`byvariable` must be NULL or a single character string.", call. = FALSE)
  }
  if (!is.character(titlestring) || length(titlestring) != 1L) {
    stop("`titlestring` must be a single character string.", call. = FALSE)
  }
  if (!is.character(ystring) || length(ystring) != 1L) {
    stop("`ystring` must be a single character string.", call. = FALSE)
  }
  if (!is.numeric(addpoints) || length(addpoints) != 1L || is.na(addpoints) || addpoints < 0) {
    stop("`addpoints` must be a single numeric value >= 0.", call. = FALSE)
  }
  if (!is.character(palette) || length(palette) != 1L) {
    stop("`palette` must be a single character string.", call. = FALSE)
  }
  if (!is.null(colorvar) && (!is.character(colorvar) || length(colorvar) != 1L)) {
    stop("`colorvar` must be NULL or a single character string.", call. = FALSE)
  }
  if (!is.null(extradata) && !is.data.frame(extradata)) {
    stop("`extradata` must be a data.frame when provided.", call. = FALSE)
  }
  if (!fit_line %in% c("model", "global_lm", "group_lm", "none")) {
    stop("`fit_line` must be one of 'model', 'global_lm', 'group_lm', 'none'.", call. = FALSE)
  }
  if (!is.logical(show_annotation) || length(show_annotation) != 1L || is.na(show_annotation)) {
    stop("`show_annotation` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!annotation_effect %in% c("standardized_beta", "raw", "partial_r2", "none")) {
    stop("`annotation_effect` must be one of 'standardized_beta', 'raw', 'partial_r2', 'none'.", call. = FALSE)
  }
  if (!is.numeric(annotation_digits) || length(annotation_digits) != 1L || is.na(annotation_digits)) {
    stop("`annotation_digits` must be a numeric scalar.", call. = FALSE)
  }
  if (!is.numeric(annotation_x) || length(annotation_x) != 1L || is.na(annotation_x)) {
    stop("`annotation_x` must be a numeric scalar.", call. = FALSE)
  }
  if (!is.numeric(annotation_y) || length(annotation_y) != 1L || is.na(annotation_y)) {
    stop("`annotation_y` must be a numeric scalar.", call. = FALSE)
  }
}

.prplot_validate_model_variables <- function(model_data, xvariable, byvariable, colorvar, extradata) {
  model_names <- names(model_data)

  if (!xvariable %in% model_names) {
    stop(sprintf("`xvariable` '%s' not found in model frame.", xvariable), call. = FALSE)
  }
  if (!is.null(byvariable) && !byvariable %in% model_names) {
    stop(sprintf("`byvariable` '%s' not found in model frame.", byvariable), call. = FALSE)
  }
  if (!is.null(extradata) && nrow(extradata) != nrow(model_data)) {
    stop("`extradata` must have the same number of rows as `model.frame(mdl)`.", call. = FALSE)
  }
  if (!is.null(colorvar)) {
    color_in_model <- colorvar %in% model_names
    color_in_extra <- !is.null(extradata) && colorvar %in% names(extradata)
    if (!color_in_model && !color_in_extra) {
      warning(sprintf("`colorvar` '%s' not found in model frame or extradata; ignoring color mapping.", colorvar), call. = FALSE)
    }
  }
}

.prplot_build_plot_data <- function(mdl, model_data, xvariable, byvariable, colorvar, extradata) {
  term_info <- .prplot_extract_training_term(mdl, model_data, xvariable)
  partial_resid <- stats::residuals(mdl) + term_info$fit

  out <- data.frame(
    partial_residual = as.numeric(partial_resid),
    term_fit = as.numeric(term_info$fit),
    stringsAsFactors = FALSE
  )
  out[[xvariable]] <- model_data[[xvariable]]

  if (!is.null(byvariable)) {
    out[[byvariable]] <- model_data[[byvariable]]
  }

  color_name <- .prplot_resolve_color_mapping(model_data, colorvar, extradata)
  if (!is.null(color_name)) {
    if (color_name %in% names(model_data)) {
      out[[color_name]] <- model_data[[color_name]]
    } else if (!is.null(extradata) && color_name %in% names(extradata)) {
      out[[color_name]] <- extradata[[color_name]]
    }
  }

  out <- .prplot_restore_factor_levels(out, model_data, xvariable)
  if (!is.null(byvariable)) out <- .prplot_restore_factor_levels(out, model_data, byvariable)
  if (!is.null(color_name) && color_name %in% names(model_data)) {
    out <- .prplot_restore_factor_levels(out, model_data, color_name)
  }

  out
}

.prplot_extract_training_term <- function(mdl, model_data, xvariable) {
  pred_terms <- tryCatch(
    stats::predict(mdl, type = "terms", se.fit = FALSE),
    error = function(e) NULL
  )
  if (is.null(pred_terms)) {
    stop("Could not compute predict(mdl, type = 'terms').", call. = FALSE)
  }

  term_names <- colnames(pred_terms)
  term_name <- .prplot_find_term_column(term_names, xvariable)
  if (is.null(term_name)) {
    stop(
      paste0(
        "Could not identify a unique simple fitted term for '",
        xvariable,
        "'. This plot currently supports simple main-effect terms."
      ),
      call. = FALSE
    )
  }

  fit <- pred_terms[, term_name]
  list(term_name = term_name, fit = as.numeric(fit))
}

.prplot_find_term_column <- function(term_names, xvariable) {
  if (is.null(term_names)) return(NULL)

  exact_hit <- which(term_names == xvariable)
  if (length(exact_hit) == 1L) {
    return(term_names[exact_hit])
  }

  prefix_hits <- grep(paste0("^", xvariable), term_names)
  if (length(prefix_hits) == 1L) {
    return(term_names[prefix_hits])
  }

  NULL
}

.prplot_restore_factor_levels <- function(out_data, source_data, variable) {
  if (!variable %in% names(out_data) || !variable %in% names(source_data)) {
    return(out_data)
  }
  if (is.factor(source_data[[variable]])) {
    out_data[[variable]] <- factor(out_data[[variable]], levels = levels(source_data[[variable]]))
  }
  out_data
}

.prplot_resolve_color_mapping <- function(model_data, colorvar, extradata) {
  if (is.null(colorvar) || !nzchar(colorvar)) return(NULL)
  if (colorvar %in% names(model_data)) return(colorvar)
  if (!is.null(extradata) && colorvar %in% names(extradata)) return(colorvar)
  NULL
}

.prplot_is_discrete <- function(x) {
  is.factor(x) || is.character(x) || is.logical(x)
}

.prplot_common_theme <- function() {
  ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      legend.position = "top",
      legend.title = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold"),
      panel.grid.minor = ggplot2::element_blank()
    )
}

.prplot_apply_faceting <- function(plot_object, byvariable) {
  if (is.null(byvariable)) return(plot_object)
  plot_object + ggplot2::facet_wrap(stats::as.formula(paste("~", byvariable)))
}

.prplot_term_annotation <- function(mdl, model_data, xvariable, effect_metric = "standardized_beta", digits = 3) {
  cf <- tryCatch(stats::coef(summary(mdl)), error = function(e) NULL)
  if (is.null(cf) || is.null(rownames(cf))) {
    return(NULL)
  }

  rn <- rownames(cf)
  idx <- which(rn == xvariable)
  if (length(idx) != 1L) {
    prefix_idx <- grep(paste0("^", xvariable), rn)
    if (length(prefix_idx) == 1L) {
      idx <- prefix_idx
    } else {
      return(NULL)
    }
  }

  pcol <- grep("Pr\\(>.*\\)", colnames(cf))
  pval <- if (length(pcol) >= 1L) cf[idx, pcol[1]] else NA_real_
  effect_line <- .prplot_effect_label(mdl, model_data, xvariable, idx, effect_metric, digits)
  p_line <- paste0("p = ", .prplot_format_p(pval, digits = digits))

  if (is.null(effect_line) || !nzchar(effect_line)) {
    return(p_line)
  }
  paste(effect_line, p_line, sep = "\n")
}

.prplot_effect_label <- function(mdl, model_data, xvariable, coef_idx, effect_metric, digits) {
  if (identical(effect_metric, "none")) {
    return(NULL)
  }

  cf <- tryCatch(stats::coef(summary(mdl)), error = function(e) NULL)
  if (is.null(cf)) return(NULL)
  est <- unname(cf[coef_idx, 1])

  if (identical(effect_metric, "raw")) {
    return(paste0("Effect = ", formatC(est, digits = digits, format = "fg")))
  }

  if (identical(effect_metric, "standardized_beta")) {
    if (!xvariable %in% names(model_data)) return(NULL)
    x <- model_data[[xvariable]]
    y <- stats::model.response(model_data)
    if (!is.numeric(x) || !is.numeric(y)) return(NULL)
    sx <- stats::sd(x, na.rm = TRUE)
    sy <- stats::sd(y, na.rm = TRUE)
    if (is.na(sx) || is.na(sy) || sx == 0 || sy == 0) return(NULL)
    std_beta <- est * sx / sy
    return(paste0("Std. beta = ", formatC(std_beta, digits = digits, format = "fg")))
  }

  if (identical(effect_metric, "partial_r2")) {
    tcol <- grep("t value", colnames(cf), fixed = TRUE)
    if (length(tcol) != 1L) return(NULL)
    tval <- unname(cf[coef_idx, tcol])
    dfres <- tryCatch(stats::df.residual(mdl), error = function(e) NA_real_)
    if (is.na(tval) || is.na(dfres)) return(NULL)
    pr2 <- (tval^2) / ((tval^2) + dfres)
    return(paste0("Partial R2 = ", formatC(pr2, digits = digits, format = "fg")))
  }

  NULL
}

.prplot_format_p <- function(p, digits = 3) {
  if (is.null(p) || length(p) != 1L || is.na(p)) return("NA")
  threshold <- 10^(-digits)
  if (p < threshold) {
    return(paste0("< ", format(threshold, scientific = FALSE, trim = TRUE)))
  }
  formatC(p, digits = digits, format = "fg")
}

.prplot_annotation_coords <- function(plot_data, xvariable, annotation_x, annotation_y) {
  x <- plot_data[[xvariable]]
  y <- plot_data$partial_residual

  xr <- range(x, na.rm = TRUE)
  yr <- range(y, na.rm = TRUE)

  if (!all(is.finite(xr)) || !all(is.finite(yr))) return(NULL)

  xpad <- if ((xr[2] - xr[1]) > 0) 0 else 1
  ypad <- if ((yr[2] - yr[1]) > 0) 0 else 1

  list(
    x = xr[1] + annotation_x * ((xr[2] - xr[1]) + xpad),
    y = yr[1] + annotation_y * ((yr[2] - yr[1]) + ypad)
  )
}

.prplot_build_annotation_layer <- function(plot_data, xvariable, annotation_text, annotation_x, annotation_y) {
  if (is.null(annotation_text) || !nzchar(annotation_text)) return(NULL)
  coords <- .prplot_annotation_coords(plot_data, xvariable, annotation_x, annotation_y)
  if (is.null(coords)) return(NULL)

  ggplot2::annotate(
    "text",
    x = coords$x,
    y = coords$y,
    label = annotation_text,
    hjust = 0,
    vjust = 1,
    size = 3.8,
    fontface = "plain",
    inherit.aes = FALSE
  )
}

.prplot_model_line_data <- function(mdl, plot_data, xvariable) {
  x <- plot_data[[xvariable]]
  if (!is.numeric(x)) return(NULL)

  model_data <- stats::model.frame(mdl)
  term_info <- .prplot_extract_training_term(mdl, model_data, xvariable)
  xgrid <- seq(min(x, na.rm = TRUE), max(x, na.rm = TRUE), length.out = 200)

  newdata <- model_data[rep(1L, length(xgrid)), , drop = FALSE]
  for (nm in names(newdata)) {
    if (identical(nm, xvariable)) {
      newdata[[nm]] <- xgrid
    } else {
      col <- model_data[[nm]]
      if (is.numeric(col)) {
        newdata[[nm]] <- rep(stats::median(col, na.rm = TRUE), length(xgrid))
      } else if (is.factor(col)) {
        ref <- levels(col)[1]
        newdata[[nm]] <- factor(rep(ref, length(xgrid)), levels = levels(col))
      } else if (is.character(col)) {
        tab <- sort(table(col), decreasing = TRUE)
        ref <- names(tab)[1]
        newdata[[nm]] <- rep(ref, length(xgrid))
      } else if (is.logical(col)) {
        tab <- sort(table(col), decreasing = TRUE)
        ref <- as.logical(names(tab)[1])
        newdata[[nm]] <- rep(ref, length(xgrid))
      } else {
        newdata[[nm]] <- rep(col[1], length(xgrid))
      }
    }
  }

  pred <- tryCatch(
    stats::predict(mdl, newdata = newdata, type = "terms", se.fit = TRUE),
    error = function(e) NULL
  )
  if (is.null(pred) || is.null(pred$fit)) return(NULL)

  term_name <- term_info$term_name
  term_cols <- colnames(pred$fit)
  if (is.null(term_cols) || !term_name %in% term_cols) return(NULL)

  fit <- as.numeric(pred$fit[, term_name])
  se_fit <- if (!is.null(pred$se.fit)) as.numeric(pred$se.fit[, term_name]) else rep(NA_real_, length(fit))
  crit <- stats::qnorm(0.975)

  data.frame(
    x = xgrid,
    fit = fit,
    lwr = fit - crit * se_fit,
    upr = fit + crit * se_fit,
    stringsAsFactors = FALSE
  )
}

.prplot_build_discrete_plot <- function(
  plot_data,
  xvariable,
  byvariable,
  titlestring,
  ystring,
  addpoints,
  colorvar,
  annotation_text,
  annotation_x,
  annotation_y
) {
  aes_base <- ggplot2::aes_string(x = xvariable, y = "partial_residual")
  p <- ggplot2::ggplot(plot_data, aes_base)

  if (!is.null(colorvar)) {
    p <- p + ggplot2::geom_boxplot(
      ggplot2::aes_string(color = colorvar),
      outlier.shape = NA,
      alpha = 0.4
    )
  } else {
    p <- p + ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.4)
  }

  if (addpoints > 0) {
    if (!is.null(colorvar)) {
      p <- p + ggplot2::geom_jitter(
        ggplot2::aes_string(color = colorvar),
        width = 0.15,
        height = 0,
        size = addpoints,
        alpha = 0.8
      )
    } else {
      p <- p + ggplot2::geom_jitter(
        width = 0.15,
        height = 0,
        size = addpoints,
        alpha = 0.8
      )
    }
  }

  ann <- .prplot_build_annotation_layer(plot_data, xvariable, annotation_text, annotation_x, annotation_y)
  if (!is.null(ann)) {
    p <- p + ann + ggplot2::coord_cartesian(clip = "off")
  }

  p <- .prplot_apply_faceting(p, byvariable)
  p +
    ggplot2::labs(title = titlestring, x = xvariable, y = ystring) +
    .prplot_common_theme()
}

.prplot_build_continuous_plot <- function(
  mdl,
  plot_data,
  xvariable,
  byvariable,
  titlestring,
  ystring,
  addpoints,
  colorvar,
  fit_line,
  annotation_text,
  annotation_x,
  annotation_y
) {
  p <- ggplot2::ggplot(plot_data, ggplot2::aes_string(x = xvariable, y = "partial_residual"))

  if (addpoints > 0) {
    if (!is.null(colorvar)) {
      p <- p + ggplot2::geom_point(
        ggplot2::aes_string(color = colorvar),
        size = addpoints,
        alpha = 0.8
      )
    } else {
      p <- p + ggplot2::geom_point(size = addpoints, alpha = 0.8)
    }
  }

  if (fit_line == "model") {
    line_data <- .prplot_model_line_data(mdl, plot_data, xvariable)
    if (!is.null(line_data)) {
      p <- p +
        ggplot2::geom_ribbon(
          data = line_data,
          ggplot2::aes(x = .data$x, ymin = .data$lwr, ymax = .data$upr),
          inherit.aes = FALSE,
          alpha = 0.15,
          fill = "grey50"
        ) +
        ggplot2::geom_line(
          data = line_data,
          ggplot2::aes(x = .data$x, y = .data$fit),
          inherit.aes = FALSE,
          linewidth = 1
        )
    }
  } else if (fit_line == "global_lm") {
    p <- p + ggplot2::geom_smooth(method = "lm", se = TRUE, formula = y ~ x)
  } else if (fit_line == "group_lm" && !is.null(colorvar)) {
    p <- p + ggplot2::geom_smooth(
      ggplot2::aes_string(color = colorvar),
      method = "lm",
      se = TRUE,
      formula = y ~ x
    )
  }

  ann <- .prplot_build_annotation_layer(plot_data, xvariable, annotation_text, annotation_x, annotation_y)
  if (!is.null(ann)) {
    p <- p + ann + ggplot2::coord_cartesian(clip = "off")
  }

  p <- .prplot_apply_faceting(p, byvariable)
  p +
    ggplot2::labs(title = titlestring, x = xvariable, y = ystring) +
    .prplot_common_theme()
}
