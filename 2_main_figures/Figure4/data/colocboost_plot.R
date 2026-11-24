#' @rdname colocboost_plot
#'
#' @title Plot visualization plot from a ColocBoost output.
#'
#' @description `colocboost_plot` generates visualization plots for colocalization events from a ColocBoost analysis.
#'
#' @param cb_output Output object from `colocboost` analysis
#' @param y Specifies the y-axis values, default is "log10p" for -log10 transformed marginal association p-values.
#' @param grange Optional plotting range of x-axis to zoom in to a specific region.
#' @param plot_cos_idx Optional indices of CoS to plot
#' @param outcome_idx Optional indices of outcomes to include in the plot. \code{outcome_idx=NULL} to plot only the outcomes having colocalization.
#' @param plot_all_outcome Optional to plot all outcome in the same figure.
#' @param plot_focal_only Logical, if TRUE only plots colocalization with focal outcome, default is FALSE.
#' @param plot_focal_cos_outcome_only Logical, if TRUE only plots colocalization including at least on colocalized outcome with focal outcome, default is FALSE.
#' @param points_color Background color for non-colocalized variables, default is "grey80".
#' @param cos_color Optional custom colors for CoS.
#' @param add_highlight Logical, if TRUE adds vertical lines at specified positions, default is FALSE
#' @param add_highlight_style Highlight marker style; supports vertical lines, hollow circles, coral-filled circles, or a coral lead plus hollow coding variant.
#' @param add_highlight_idx Optional indices for add_highlight variables.
#' @param add_highlight_outcome_idx Optional list, same length as add_highlight_idx; each
#'   element is a vector of outcome indices for which that highlight position is actually
#'   drawn (e.g. only outcomes truly in the corresponding CoS). NULL (default) draws every
#'   highlight position on every outcome, as before.
#' @param outcome_names Optional vector of outcomes names for the subtitle of each figure. \code{outcome_names=NULL} for the outcome name shown in \code{data_info}.
#' @param plot_cols Number of columns in the plot grid, default is 2. If you have many colocalization. please consider increasing this.
#' @param variant_coord Logical, if TRUE uses variant coordinates on x-axis, default is FALSE. This is required the variable names including position information.
#' @param show_top_variables Logical, if TRUE shows top variables for each CoS, default is FALSE
#' @param show_cos_to_uncoloc Logical, if TRUE shows colocalization to uncolocalized outcomes to diagnose, default is FALSE
#' @param show_cos_to_uncoloc_idx Optional indices for showing CoS to all uncolocalized outcomes
#' @param show_cos_to_uncoloc_outcome Optional outcomes for showing CoS to uncolocalized outcomes
#' @param plot_ucos Logical, if TRUE plots also trait-specific (uncolocalized) sets , default is FALSE
#' @param plot_ucos_idx Optional indices of trait-specific (uncolocalized) sets to plot when included
#' @param title_specific Optional specific title to display in plot title
#' @param ylim_each Logical, if TRUE uses separate y-axis limits for each plot, default is TRUE
#' @param outcome_legend_pos Position for outcome legend, default is "top"
#' @param outcome_legend_size Size for outcome legend text, default is 1.2
#' @param cos_legend_pos Proportion of the legend from (left edge, bottom edge) following by the distance for multiple labels (horizontal-x space, vertical-y space), default as (0.05, 0.4, 0.1, 0.5) at the left - median position
#' @param show_variable Logical, if TRUE displays variant IDs, default is FALSE
#' @param lab_style Vector of two numbers for label style (size, boldness), default is c(2, 1)
#' @param axis_style Vector of two numbers for axis style (size, boldness), default is c(2, 1)
#' @param title_style Vector of two numbers for title style (size, boldness), default is c(2.5, 2)
#' @param ... Additional parameters passed to `plot` functions
#'
#' @return Visualization plot for each colocalization event.
#'
#' @examples
#' # colocboost example
#' set.seed(1)
#' N <- 1000
#' P <- 100
#' # Generate X with LD structure
#' sigma <- 0.9^abs(outer(1:P, 1:P, "-"))
#' X <- MASS::mvrnorm(N, rep(0, P), sigma)
#' colnames(X) <- paste0("SNP", 1:P)
#' L <- 3
#' true_beta <- matrix(0, P, L)
#' true_beta[10, 1] <- 0.5 # SNP10 affects trait 1
#' true_beta[10, 2] <- 0.4 # SNP10 also affects trait 2 (colocalized)
#' true_beta[50, 2] <- 0.3 # SNP50 only affects trait 2
#' true_beta[80, 3] <- 0.6 # SNP80 only affects trait 3
#' Y <- matrix(0, N, L)
#' for (l in 1:L) {
#'   Y[, l] <- X %*% true_beta[, l] + rnorm(N, 0, 1)
#' }
#' res <- colocboost(X = X, Y = Y)
#' colocboost_plot(res, plot_cols = 1)
#' colocboost_plot(res, plot_cols = 1, outcome_idx = 1:3)
#'
#' @importFrom utils head tail
#' @importFrom graphics abline axis legend mtext par points text
#' @importFrom grDevices adjustcolor
#' 
#' @source See detailed instructions in our tutorial portal: \url{https://statfungen.github.io/colocboost/articles/Visualization_ColocBoost_Output.html}
#'
#' @family colocboost_plot
#' @export
colocboost_plot <- function(cb_output, y = "log10p",
                            grange = NULL,
                            plot_cos_idx = NULL,
                            outcome_idx = NULL,
                            plot_all_outcome = FALSE,
                            plot_focal_only = FALSE,
                            plot_focal_cos_outcome_only = FALSE,
                            points_color = "grey80",
                            cos_color = NULL,
                            add_highlight = FALSE,
                            add_highlight_idx = NULL,
                            add_highlight_style = "vertical_lines",
                            add_highlight_outcome_idx = NULL,
                            add_highlight_uncoloc_idx = NULL,
                            outcome_annotations = NULL,
                            annotation_scale = 1,
                            bottom_track = NULL,
                            compact_x_axes = FALSE,
                            track_height = 0.35,
                            panel_top_margin = 0.8,
                            outer_top_margin = 1,
                            panel_axis_size = NULL,
                            panel_y_ticks = NULL,
                            panel_y_limits = NULL,
                            panel_left_margin = 5,
                            panel_right_margin = 1,
                            panel_inner_right_margin = NULL,
                            panel_column_widths = NULL,
                            panel_left_margin_by_col = NULL,
                            hide_ylab_cols = NULL,
                            track_top_margin = 0.2,
                            outcome_names = NULL,
                            plot_cols = 2,
                            variant_coord = FALSE,
                            genomic_ticks = NULL,
                            show_top_variables = FALSE,
                            show_cos_to_uncoloc = FALSE,
                            show_cos_to_uncoloc_idx = NULL,
                            show_cos_to_uncoloc_outcome = NULL,
                            plot_ucos = FALSE,
                            plot_ucos_idx = NULL,
                            title_specific = NULL,
                            ylim_each = TRUE,
                            outcome_legend_pos = "top",
                            outcome_legend_size = 1.8,
                            cos_legend_pos = c(0.05, 0.4, 0.1, 0.5),
                            show_variable = FALSE,
                            lab_style = c(2, 1),
                            axis_style = c(2, 1),
                            title_style = c(2.5, 2),
                            ...) {
  if (!inherits(cb_output, "colocboost")) {
    stop("Input of colocboost_plot must be a 'colocboost' object!")
  }
  
  # get cb_plot_input data from colocboost results
  cb_plot_input <- get_input_plot(cb_output,
                                  plot_cos_idx = plot_cos_idx,
                                  variant_coord = variant_coord,
                                  outcome_names = outcome_names,
                                  plot_focal_only = plot_focal_only,
                                  plot_focal_cos_outcome_only = plot_focal_cos_outcome_only,
                                  show_cos_to_uncoloc = show_cos_to_uncoloc,
                                  show_cos_to_uncoloc_idx = show_cos_to_uncoloc_idx,
                                  show_cos_to_uncoloc_outcome = show_cos_to_uncoloc_outcome,
                                  plot_ucos = plot_ucos, plot_ucos_idx = plot_ucos_idx
  )
  # get initial set up of plot
  cb_plot_init <- plot_initial(cb_plot_input,
                               y = y, points_color = points_color, cos_color = cos_color,
                               ylim_each = ylim_each, title_specific = title_specific,
                               outcome_legend_pos = outcome_legend_pos, outcome_legend_size = outcome_legend_size,
                               cos_legend_pos = cos_legend_pos, plot_ucos = plot_ucos,
                               show_variable = show_variable, lab_style = lab_style, axis_style = axis_style,
                               title_style = title_style, ...
  )
  
  colocboost_plot_basic <- function(cb_plot_input, cb_plot_init,
                                    outcome_idx = NULL, plot_all_outcome = FALSE,
                                    grange = NULL, plot_cols = 2,
                                    add_highlight = FALSE, add_highlight_idx = NULL,
                                    add_highlight_style = "vertical_lines",
                                    add_highlight_outcome_idx = NULL,
                                    add_highlight_uncoloc_idx = NULL,
                                    outcome_annotations = NULL,
                                    annotation_scale = 1,
                                    bottom_track = NULL,
                                    compact_x_axes = FALSE,
                                    track_height = 0.35,
                                    panel_top_margin = 0.8,
                                    outer_top_margin = 1,
                                    panel_axis_size = NULL,
                                    panel_y_ticks = NULL,
                                    panel_y_limits = NULL,
                                    panel_left_margin = 5,
                                    panel_right_margin = 1,
                                    panel_inner_right_margin = NULL,
                                    panel_column_widths = NULL,
                                    panel_left_margin_by_col = NULL,
                                    hide_ylab_cols = NULL,
                                    track_top_margin = 0.2,
                                    genomic_ticks = NULL,
                                    show_top_variables = TRUE,
                                    ...) {
    args <- list(...)
    args$xlab <- cb_plot_init$xlab
    args$ylab <- cb_plot_init$ylab
    args$col <- cb_plot_init$bg
    if (is.null(grange)) {
      args$x <- cb_plot_init$x
      y <- cb_plot_init$y
    } else {
      args$x <- cb_plot_init$x[grange]
      y <- lapply(cb_plot_init$y, function(yy) yy[grange])
    }
    args$pch <- cb_plot_init$pch
    args$cex.axis <- if (is.null(panel_axis_size)) cb_plot_init$axis_size else panel_axis_size
    args$cex.lab <- cb_plot_init$lab_size
    args$font.lab <- cb_plot_init$lab_face
    if (is.null(panel_inner_right_margin)) {
      panel_inner_right_margin <- panel_right_margin
    }

    draw_outcome_annotation <- function(annotation) {
      if (is.null(annotation)) return(invisible())

      usr <- par("usr")
      x_span <- usr[2] - usr[1]
      y_span <- usr[4] - usr[3]
      x_left <- usr[1] + 0.035 * x_span
      y_top <- usr[4] - 0.035 * y_span

      halo_dx <- 0.0012 * x_span
      halo_dy <- 0.0012 * y_span
      halo_offsets <- expand.grid(
        dx = c(-halo_dx, 0, halo_dx),
        dy = c(-halo_dy, 0, halo_dy)
      )
      halo_offsets <- halo_offsets[halo_offsets$dx != 0 | halo_offsets$dy != 0, ]
      draw_halo_text <- function(x, y, labels, adj = c(0, 1), cex = 1, font = 1) {
        for (jj in seq_len(nrow(halo_offsets))) {
          text(
            x + halo_offsets$dx[jj], y + halo_offsets$dy[jj], labels = labels,
            adj = adj, cex = cex, font = font, col = "white", xpd = FALSE
          )
        }
        text(
          x, y, labels = labels, adj = adj, cex = cex, font = font,
          col = "black", xpd = FALSE
        )
      }
      draw_halo_text(
        x_left, y_top, labels = annotation$line1,
        adj = c(0, 1), cex = 1.38 * annotation_scale, font = 2
      )
      line_gap_fraction <- if (!is.null(annotation$line_gap)) {
        annotation$line_gap
      } else {
        0.095
      }
      if (!is.null(annotation$line2)) {
        draw_halo_text(
          x_left, y_top - line_gap_fraction * y_span, labels = annotation$line2,
          adj = c(0, 1), cex = 1.28 * annotation_scale
        )
      }

      if (!is.null(annotation$glycan)) {
        glycan_y_offset_fraction <- if (!is.null(annotation$glycan_y_offset)) {
          annotation$glycan_y_offset
        } else {
          line_gap_fraction
        }
        glycan_y <- y_top - glycan_y_offset_fraction * y_span
        component_gap_fraction <- if (!is.null(annotation$component_gap)) {
          annotation$component_gap
        } else {
          0.092
        }
        symbol_offset_fraction <- if (!is.null(annotation$symbol_offset)) {
          annotation$symbol_offset
        } else {
          0.048
        }
        component_gap <- component_gap_fraction * x_span
        for (ii in seq_along(annotation$glycan)) {
          component <- annotation$glycan[[ii]]
          component_x <- x_left + (ii - 1) * component_gap
          draw_halo_text(
            component_x, glycan_y, labels = paste0(component$count, "\u00d7"),
            adj = c(0, 0.5), cex = 1.28 * annotation_scale
          )
          points(
            component_x + symbol_offset_fraction * x_span, glycan_y,
            pch = component$pch, bg = "white", col = "white",
            cex = 2.02 * annotation_scale, lwd = 1.2, xpd = FALSE
          )
          points(
            component_x + symbol_offset_fraction * x_span, glycan_y,
            pch = component$pch, bg = component$fill,
            col = component$border, cex = 1.82 * annotation_scale, lwd = 1.0, xpd = FALSE
          )
        }
      }
    }
    # - change position
    cb_plot_init$outcome_legend_pos <- switch(cb_plot_init$outcome_legend_pos,
                                              "right" = 4,
                                              "left" = 2,
                                              "top" = 3,
                                              "bottom" = 1
    )
    
    # - begin plotting
    coloc_cos <- cb_plot_input$cos
    outcomes <- cb_plot_input$outcomes
    if (length(y)==1) outcome_idx <- 1
    if (plot_all_outcome){
      outcome_idx <- 1:length(y)
      if (!is.null(cb_plot_input$focal_outcome)) {
        p_focal <- grep(cb_plot_input$focal_outcome, outcomes)
        outcome_idx <- c(p_focal, setdiff(outcome_idx, p_focal))
      }
    } else {
      if (is.null(outcome_idx)) {
        if (is.null(coloc_cos)) {
          # - no colocalized effects, draw all outcomes in this region
          if (length(cb_plot_input$outcomes) == 1) {
            message("There is no fine-mapped causal effect in this region!. Showing margianl for this outcome!")
          } else {
            if (length(y) == 1){
              message("There is no colocalization in this region!. Showing VCP = 0!")
            } else {
              message("There is no colocalization in this region!. Showing margianl for all outcomes!")
            }
          }
          outcome_idx <- 1:length(y)
        } else {
          n.coloc <- length(coloc_cos)
          coloc_index <- cb_plot_input$coloc_index
          outcome_idx <- Reduce(union, coloc_index)
        }
        if (!is.null(cb_plot_input$focal_outcome)) {
          p_focal <- grep(cb_plot_input$focal_outcome, outcomes)
          include_focal <- sapply(cb_plot_input$coloc_index, function(ci) {
            p_focal %in% ci
          })
          if (any(include_focal)) {
            coloc_index <- cb_plot_input$coloc_index[order(include_focal == "FALSE")]
            coloc_index <- Reduce(union, coloc_index)
            outcome_idx <- c(p_focal, setdiff(coloc_index, p_focal))
          }
        }
      }
    }
    if (length(outcome_idx) == 1) {
      plot_cols <- 1
    }
    nrow <- ceiling(length(outcome_idx) / plot_cols)
    if (!is.null(cb_plot_init$xtext)) {
      bottom <- 6
    } else {
      bottom <- 2
    }
    # - restore users' options
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar))
    if (!is.null(bottom_track)) {
      layout(
        matrix(seq_len((nrow + 1) * plot_cols), nrow = nrow + 1, byrow = TRUE),
        heights = c(rep(1, nrow), track_height),
        widths = if (is.null(panel_column_widths)) rep(1, plot_cols) else panel_column_widths
      )
      par(oma = if (!is.null(cb_plot_init$title)) c(0, 0, 3, 0) else c(0, 0, outer_top_margin, 0))
    } else if (!is.null(cb_plot_init$title)) {
      par(mfrow = c(nrow, plot_cols), mar = c(bottom, 5, 2, 1), oma = c(0, 0, 3, 0))
    } else {
      par(mfrow = c(nrow, plot_cols), mar = c(bottom, 5, 2, 1), oma = c(0, 0, 1, 0))
    }
    ylab_all <- args$ylab
    left_margin_of <- function(col) if (is.null(panel_left_margin_by_col)) panel_left_margin else panel_left_margin_by_col[col]
    shared_x_range <- range(args$x, finite = TRUE)
    if (is.null(genomic_ticks)) {
      shared_tick_step <- if (diff(shared_x_range) <= 2e5) 5e4 else 1e5
      shared_tick_at <- seq(
        ceiling(shared_x_range[1] / shared_tick_step) * shared_tick_step,
        floor(shared_x_range[2] / shared_tick_step) * shared_tick_step,
        by = shared_tick_step
      )
    } else {
      shared_tick_at <- genomic_ticks[
        genomic_ticks >= shared_x_range[1] & genomic_ticks <= shared_x_range[2]
      ]
    }
    for (plot_i in seq_along(outcome_idx)) {
      iy <- outcome_idx[plot_i]
      plot_row <- ceiling(plot_i / plot_cols)
      plot_col <- (plot_i - 1) %% plot_cols + 1
      current_right_margin <- if (plot_col < plot_cols) {
        panel_inner_right_margin
      } else {
        panel_right_margin
      }
      custom_y_ticks <- NULL
      if (!is.null(panel_y_ticks)) {
        tick_key <- as.character(plot_i)
        if (!is.null(names(panel_y_ticks)) && tick_key %in% names(panel_y_ticks)) {
          custom_y_ticks <- panel_y_ticks[[tick_key]]
        } else if (is.null(names(panel_y_ticks)) && length(panel_y_ticks) >= plot_i) {
          custom_y_ticks <- panel_y_ticks[[plot_i]]
        }
      }
      custom_y_limits <- NULL
      if (!is.null(panel_y_limits)) {
        limit_key <- as.character(plot_i)
        if (!is.null(names(panel_y_limits)) && limit_key %in% names(panel_y_limits)) {
          custom_y_limits <- panel_y_limits[[limit_key]]
        } else if (is.null(names(panel_y_limits)) && length(panel_y_limits) >= plot_i) {
          custom_y_limits <- panel_y_limits[[plot_i]]
        }
      }
      args$ylab <- if (plot_col %in% hide_ylab_cols) "" else ylab_all
      args$yaxt <- if (is.null(custom_y_ticks)) NULL else "n"
      show_x_axis <- !compact_x_axes || (plot_row == nrow && is.null(bottom_track))
      if (!is.null(bottom_track)) {
        par(mar = if (show_x_axis) c(1.8, left_margin_of(plot_col), panel_top_margin, current_right_margin) else c(0.4, left_margin_of(plot_col), panel_top_margin, current_right_margin))
      }
      args$y <- y[[iy]]
      args$ylim <- if (is.null(custom_y_limits)) {
        c(cb_plot_init$ymin[iy], cb_plot_init$ymax[iy])
      } else {
        custom_y_limits
      }
      if (!is.null(cb_plot_init$xtext)) {
        args$xaxt <- "n"
        do.call(plot, args)
        axis(1, at = args$x, labels = FALSE)
        text(x = args$x, y = par("usr")[3] - 0.1, labels = cb_plot_init$xtext, srt = 45, adj = 1, xpd = TRUE)
      } else if (isTRUE(cb_plot_input$variant_coord)) {
        args$xaxt <- "n"
        do.call(plot, args)
        if (show_x_axis) {
          axis(
            1, at = shared_tick_at,
            labels = paste0(formatC(shared_tick_at / 1e6, format = "f", digits = 2), "M"),
            cex.axis = cb_plot_init$axis_size
          )
        }
      } else {
        do.call(plot, args)
      }
      if (!is.null(custom_y_ticks)) {
        axis(2, at = custom_y_ticks, cex.axis = args$cex.axis)
      }
      if (length(y) != 1 && is.null(outcome_annotations)){
        mtext(outcomes[iy],
              side = cb_plot_init$outcome_legend_pos, line = 0.2, adj = 0.5,
              cex = cb_plot_init$outcome_legend_size, font = 1
        )
      }
      
      # mark variables in CoS to colocalized outcomes
      if (!is.null(coloc_cos)) {
        n.coloc <- length(coloc_cos)
        coloc_index <- cb_plot_input$coloc_index
        if (length(y)==1){coloc_index = lapply(coloc_index, function(i) 1 )}
        legend_text <- list(col = vector())
        legend_text$col <- head(cb_plot_init$col, n.coloc)
        
        # check which coloc set for this outcome
        p.coloc <- sapply(coloc_index, function(idx) sum(idx == iy) != 0)
        p.coloc <- which(p.coloc)
        coloc_cos.idx <- coloc_cos[p.coloc]
        for (i.cs in p.coloc) {
          # add the points with specific color
          cs <- as.numeric(coloc_cos[[i.cs]])
          x0 <- intersect(args$x, cs)
          y1 <- args$y[match(x0, args$x)]
          points(x0, y1, pch = 21, bg = legend_text$col[i.cs], col = NA, cex = 2.5, lwd = 2.5)
          if (show_top_variables) {
            # add the hits points with "red"
            cs_hits <- as.numeric(cb_plot_input$cos_hits[[i.cs]])
            x_hits <- intersect(args$x, cs_hits)
            y_hits <- args$y[match(x_hits, args$x)]
            points(x_hits, y_hits, pch = 21, bg = legend_text$col[i.cs], col = "#E31A1C", cex = 3, lwd = 3)
          }
        }
        
        # mark variables in CoS to uncolocalized outcomes
        uncoloc <- cb_plot_input$uncoloc
        if (!is.null(uncoloc)) {
          p.uncoloc <- sapply(uncoloc$outcome_to_uncoloc, function(idx) sum(idx == iy) != 0)
          p.uncoloc <- which(p.uncoloc)
          shape_col <- texts_col <- c()
          for (i.uncoloc in p.uncoloc) {
            uncoloc_outcome <- uncoloc$outcome_to_uncoloc[[i.uncoloc]]
            if (iy %in% uncoloc_outcome) {
              # add the points with specific color
              cs <- as.numeric(uncoloc$cos_to_uncoloc[[i.uncoloc]])
              x0 <- intersect(args$x, cs)
              y1 <- args$y[match(x0, args$x)]
              points(x0, y1,
                     pch = 4, col = adjustcolor(legend_text$col[uncoloc$cos_idx_to_uncoloc[i.uncoloc]], alpha.f = 0.4),
                     cex = 1.5, lwd = 1.5
              )
              shape_col <- c(shape_col, adjustcolor(legend_text$col[uncoloc$cos_idx_to_uncoloc[i.uncoloc]], alpha.f = 1))
              texts_col <- c(texts_col, adjustcolor(legend_text$col[uncoloc$cos_idx_to_uncoloc[i.uncoloc]], alpha.f = 0.8))
            }
          }
          texts <- unlist(sapply(names(uncoloc$outcome_to_uncoloc), function(name) {
            positions <- which(uncoloc$outcome_to_uncoloc[[name]] == iy)
            if(length(positions) > 0) {
              return(uncoloc$cos_uncoloc_texts[[name]][positions])
            } else {
              return(NULL)
            }
          }))
          # NOTE: this used to be `if (length(texts) == 0) next`, which -- because `next`
          # always targets the nearest enclosing loop, i.e. the outer per-outcome `for (iy ...)`
          # -- skipped add_highlight and draw_outcome_annotation below for EVERY outcome that
          # simply had nothing to report here (the common case: any outcome colocalizing
          # cleanly with a CoS actually shown). That silently dropped panel titles/highlights
          # for all "normal" outcomes as soon as show_cos_to_uncoloc = TRUE and only drew them
          # for the (rare) outcomes uncolocalized from every shown CoS. Skip just the legend.
          if (length(texts) > 0) {
            # Get current plot area coordinates
            usr <- par("usr")
            x_pos <- usr[1] + cb_plot_init$cos_legend_pos[1] * (usr[2] - usr[1])  # 5% from left edge
            y_pos <- usr[3] + cb_plot_init$cos_legend_pos[2] * (usr[4] - usr[3])  # 50% from bottom edge
            legend(x = x_pos, y = y_pos, texts,
                   bty = "n", col = shape_col, text.col = texts_col,
                   cex = 1.5, pt.cex = 1.5, pch = 4, x.intersp = cos_legend_pos[3], y.intersp = cos_legend_pos[4]
            )
          }
        }
      }
      if (add_highlight) {
        highlight_colors <- c("#D73027", "#542788", "#1B9E77")

        for (iii in seq_along(add_highlight_idx)) {
          # Skip drawing this highlight position on outcomes it doesn't actually apply to
          # (e.g. a CoS's lead variant on an outcome that isn't a member of that CoS) --
          # show_cos_to_uncoloc's x-mark already covers "not colocalized here" on its own.
          if (!is.null(add_highlight_outcome_idx) && !is.null(add_highlight_outcome_idx[[iii]]) &&
              !(iy %in% add_highlight_outcome_idx[[iii]])) {
            next
          }
          if (!(add_highlight_style %in% c(
            "vertical_lines", "star", "coral_circle", "lead_coral_coding_hollow"
          ))) {
            message("Unknown add_highlight_style; defaulting to vertical_lines")
            add_highlight_style <- "vertical_lines"
          }

          this_col <- highlight_colors[(iii - 1) %% length(highlight_colors) + 1]  # cycle colors

          if (add_highlight_style == "vertical_lines") {
            abline(v = add_highlight_idx[iii], col = this_col, lwd = 1.5, lty = "dashed")
          } else if (add_highlight_style == "star") {
            points(add_highlight_idx[iii], args$y[match(add_highlight_idx[iii], args$x)], 
                   pch = 21, bg = NA, col = "white", cex = 3.2, lwd = 4)
            points(add_highlight_idx[iii], args$y[match(add_highlight_idx[iii], args$x)], 
                   pch = 21, bg = NA, col = this_col, cex = 2.8, lwd = 2.8)
          } else if (add_highlight_style == "coral_circle") {
            points(add_highlight_idx[iii], args$y[match(add_highlight_idx[iii], args$x)],
                   pch = 21, bg = "white", col = "white", cex = 3.25, lwd = 3.8)
            points(add_highlight_idx[iii], args$y[match(add_highlight_idx[iii], args$x)],
                   pch = 21, bg = "#F46D43", col = "#B2182B", cex = 2.85, lwd = 2.6)
          } else if (add_highlight_style == "lead_coral_coding_hollow" && iy %in% add_highlight_uncoloc_idx) {
            # Outcome not colocalized with the CoS: no filled marker; coloured ring with an x inside
            ring_col <- if (iii == 1) "#B2182B" else "#542788"
            hx <- add_highlight_idx[iii]; hy <- args$y[match(hx, args$x)]
            points(hx, hy, pch = 21, bg = NA, col = "white", cex = 3.25, lwd = 3.8)
            points(hx, hy, pch = 21, bg = NA, col = ring_col, cex = 2.85, lwd = 2.6)
            points(hx, hy, pch = 4, col = ring_col, cex = 1.6, lwd = 2)
          } else if (add_highlight_style == "lead_coral_coding_hollow") {
            points(add_highlight_idx[iii], args$y[match(add_highlight_idx[iii], args$x)],
                   pch = 21, bg = if (iii == 1) "white" else NA,
                   col = "white", cex = 3.25, lwd = 3.8)
            if (iii == 1) {
              points(add_highlight_idx[iii], args$y[match(add_highlight_idx[iii], args$x)],
                     pch = 21, bg = "#F46D43", col = "#B2182B", cex = 2.85, lwd = 2.6)
            } else {
              points(add_highlight_idx[iii], args$y[match(add_highlight_idx[iii], args$x)],
                     pch = 21, bg = NA, col = "#542788", cex = 2.85, lwd = 2.6)
            }
          }
        }
      }
      if (!is.null(outcome_annotations) && length(outcome_annotations) >= iy) {
        draw_outcome_annotation(outcome_annotations[[iy]])
      }
    }
    if (!is.null(bottom_track)) {
      empty_panels <- nrow * plot_cols - length(outcome_idx)
      if (empty_panels > 0) {
        for (empty_i in seq_len(empty_panels)) {
          plot.new()
        }
      }
      for (track_col in seq_len(plot_cols)) {
        current_right_margin <- if (track_col < plot_cols) {
          panel_inner_right_margin
        } else {
          panel_right_margin
        }
        par(mar = c(3.8, left_margin_of(track_col), track_top_margin, current_right_margin))
        bottom_track(
          x_range = shared_x_range,
          tick_at = shared_tick_at,
          column = track_col,
          axis_size = cb_plot_init$axis_size
        )
      }
    }
    if (!is.null(cb_plot_init$title)) {
      mtext(cb_plot_init$title,
            side = 3, line = 0, outer = TRUE,
            cex = cb_plot_init$title_size, font = cb_plot_init$title_face
      )
    }
    return(invisible())
  }
  
  colocboost_plot_basic(cb_plot_input, cb_plot_init,
                        grange = grange,
                        outcome_idx = outcome_idx, 
                        plot_all_outcome = plot_all_outcome,
                        plot_cols = plot_cols,
                        add_highlight = add_highlight, add_highlight_idx = add_highlight_idx,
                        add_highlight_style = add_highlight_style,
                        add_highlight_outcome_idx = add_highlight_outcome_idx,
                        add_highlight_uncoloc_idx = add_highlight_uncoloc_idx,
                        outcome_annotations = outcome_annotations,
                        annotation_scale = annotation_scale,
                        bottom_track = bottom_track,
                        compact_x_axes = compact_x_axes,
                        track_height = track_height,
                        panel_top_margin = panel_top_margin,
                        outer_top_margin = outer_top_margin,
                        panel_axis_size = panel_axis_size,
                        panel_y_ticks = panel_y_ticks,
                        panel_y_limits = panel_y_limits,
                        panel_left_margin = panel_left_margin,
                        panel_right_margin = panel_right_margin,
                        panel_inner_right_margin = panel_inner_right_margin,
                        panel_column_widths = panel_column_widths,
                        panel_left_margin_by_col = panel_left_margin_by_col,
                        hide_ylab_cols = hide_ylab_cols,
                        track_top_margin = track_top_margin,
                        genomic_ticks = genomic_ticks,
                        show_top_variables = show_top_variables,
                        ...
  )
}



# get input data for cb_plot
get_input_plot <- function(cb_output, plot_cos_idx = NULL,
                           variant_coord = FALSE,
                           outcome_names = NULL,
                           plot_focal_only = FALSE,
                           plot_focal_cos_outcome_only = FALSE,
                           show_cos_to_uncoloc = FALSE,
                           show_cos_to_uncoloc_idx = NULL,
                           show_cos_to_uncoloc_outcome = NULL,
                           plot_ucos = FALSE,
                           plot_ucos_idx = NULL) {
  
  # check ucos exists
  if (plot_ucos && !"ucos_details" %in% names(cb_output)) {
    warning(
      "Since you want to plot trait-specific (uncolocalized) sets with plot_ucos = TRUE,",
      " but there is no output of ucos from colocboost.\n",
      " Please run colocboost model with output_level=2!",
      " Only show colocalization results if exists."
    )
  }
  # redefined outcome names
  if (!is.null(outcome_names)) {
    cb_output$data_info$outcome_info$outcome_names <- outcome_names
    names(cb_output$data_info$z) <- outcome_names
  }
  
  # extract results from colocboost
  analysis_outcome <- cb_output$data_info$outcome_info$outcome_names
  focal_outcome_idx <- which(cb_output$data_info$outcome_info$is_focal)
  if (length(focal_outcome_idx) != 0) {
    focal_outcome <- analysis_outcome[focal_outcome_idx]
  } else {
    focal_outcome <- NULL
  }
  # check if focal cos
  focal_cos <- cb_output$cos_summary$cos_id[cb_output$cos_summary$focal_outcome != FALSE]
  if_focal <- !is.na(match(names(cb_output$cos_details$cos$cos_variables), focal_cos))
  # extract z-scores
  variables <- cb_output$data_info$variables
  Z <- cb_output$data_info$z
  coef <- cb_output$data_info$coef
  vcp <- list(as.numeric(cb_output$vcp))
  
  # if finemapping
  if (cb_output$data_info$n_outcomes == 1) {
    cb_output$cos_details$cos$cos_variables <- cb_output$ucos_details$ucos$ucos_variables
    cb_output$cos_details$cos$cos_index <- cb_output$ucos_details$ucos$ucos_index
    cb_output$cos_details$cos_top_variables <- cb_output$ucos_details$ucos_top_variables
    cb_output$cos_details$cos_outcomes <- cb_output$ucos_details$ucos_outcomes
    if (!is.null(cb_output$cos_details$cos$cos_variables)) {
      cb_output$cos_details$cos_vcp <- cb_output$ucos_details$ucos_weight
      vcp <- list(as.numeric(cb_output$vpa))
    }
    if_focal <- rep(TRUE, length(cb_output$cos_details$cos$cos_variables))
  }
  # extract coloc_cos
  coloc_variables <- cb_output$cos_details$cos$cos_variables
  coloc_cos <- cb_output$cos_details$cos$cos_index
  coloc_index <- cb_output$cos_details$cos_outcomes$outcome_index
  # top_variables
  coloc_hits <- lapply(names(coloc_cos), function(cn) {
    p <- grep(cn, rownames(cb_output$cos_details$cos_top_variables))
    cb_output$cos_details$cos_top_variables$top_index[p]
  })
  names(coloc_hits) <- names(coloc_cos)
  if (!is.null(coloc_cos)) {
    # extract vcp each outcome
    cos_vcp <- lapply(1:length(analysis_outcome), function(iy) {
      pos <- which(sapply(coloc_index, function(idx) iy %in% idx))
      if (length(pos) != 0) {
        w <- do.call(cbind, cb_output$cos_details$cos_vcp[pos])
        return(1 - apply(1 - w, 1, prod))
      } else {
        return(rep(0, length(variables)))
      }
    })
    ncos <- length(cb_output$cos_details$cos$cos_index)
    
    select_cs <- 1:ncos
    if (!is.null(plot_cos_idx)) {
      if (length(setdiff(plot_cos_idx, c(1:ncos))) != 0) {
        stop("Please check plot_cos_idx!")
      }
      select_cs <- plot_cos_idx
    } else {
      if (plot_focal_only || plot_focal_cos_outcome_only) {
        if (sum(if_focal) == 0) {
          message("No focal CoS, draw all CoS.")
        } else if (plot_focal_only) {
          select_cs <- which(if_focal)
        } else {  # plot_focal_cos_outcome_only is true here
          # Get all outcomes colocalized with focal CoS
          focal_outcomes <- unique(unlist(coloc_index[if_focal]))
          # Find CoS that include at least one of these focal outcomes
          cos_with_focal_outcomes <- sapply(coloc_index, function(cos_outcomes) {
            length(intersect(cos_outcomes, focal_outcomes)) > 0
          })
          # Extract the indices of these CoS
          select_cs <- which(cos_with_focal_outcomes)
        }
      }
    }
    coloc_variables <- coloc_variables[select_cs]
    coloc_cos <- coloc_cos[select_cs]
    coloc_index <- coloc_index[select_cs]
    coloc_hits <- coloc_hits[select_cs]
  } else {
    if (cb_output$data_info$n_outcomes == 1) {
      warning("No fine-mapped causal effects in this region!")
    } else {
      warning("No colocalized effects in this region!")
    }
    ncos <- 0
    coloc_index <- select_cs <- NULL
    vcp <- list(rep(0, cb_output$data_info$n_variables))
    cos_vcp <- lapply(1:length(analysis_outcome), function(iy) rep(0, cb_output$data_info$n_variables) )
  }
  # extract x axis
  if (variant_coord) {
    variable.info <- do.call(rbind, lapply(variables, function(variable) strsplit(variable, ":")[[1]]))
    chrom <- as.numeric(sapply(variable.info[, 1], function(ss) strsplit(ss, "chr")[[1]][2]))
    x <- data.frame(
      "chrom" = chrom,
      "pos" = as.numeric(variable.info[, 2])
    )
    coloc_cos <- lapply(coloc_cos, function(cos) x$pos[cos])
    coloc_hits <- lapply(coloc_hits, function(cos) x$pos[cos])
  } else {
    x <- data.frame(
      "chrom" = NA,
      "pos" = c(1:length(variables))
    )
  }
  plot_input <- list(
    "outcomes" = analysis_outcome,
    "focal_outcome" = focal_outcome,
    "variables" = variables,
    "x" = x,
    "Zscores" = Z,
    "vcp" = vcp,
    "cos_vcp" = cos_vcp,
    "ucos_cos_int_weights" = cos_vcp,
    "coef" = coef,
    "cos" = coloc_cos,
    "cos_hits" = coloc_hits,
    "coloc_index" = coloc_index,
    "select_cos" = select_cs
  )
  
  # check plot uncolocalizaed confidence sets from ucos_details
  if (plot_ucos & cb_output$data_info$n_outcomes > 1){
    ucos_details <- cb_output$ucos_details
    if (!is.null(ucos_details)){
      ucos <- ucos_details$ucos$ucos_index
      ucos_outcome_index <- ucos_details$ucos_outcomes$outcome_index
      ucos_hits <- lapply(ucos, function(x) x[[1]])
      # check inclusion of other options
      select_ucos <- 1:length(ucos)
      if (!is.null(plot_ucos_idx)) {
        if (length(setdiff(plot_cos_idx, select_ucos)) != 0) {
          stop("Please check plot_ucos_idx!")
        }
        select_ucos <- plot_ucos_idx
      } else if (plot_focal_cos_outcome_only && sum(if_focal) != 0) {
        # Get all outcomes colocalized with focal CoS
        focal_outcomes <- unique(unlist(plot_input$coloc_index))
        # Find uCoS that include at least one of these focal outcomes
        ucos_with_focal_outcomes <- sapply(ucos_outcome_index, function(ucos_outcomes) {
          length(intersect(ucos_outcomes, focal_outcomes)) > 0
        })
        # Extract the indices of these uCoS
        select_ucos <- which(ucos_with_focal_outcomes)
      }
      plot_input$cos <- c(plot_input$cos, ucos[select_ucos])
      plot_input$cos_hits <- c(plot_input$cos_hits, ucos_hits[select_ucos])
      plot_input$coloc_index <- c(plot_input$coloc_index, ucos_outcome_index[select_ucos])
      plot_input$select_cos <- c(plot_input$select_cos, ncos + select_ucos)
      # add the recalibrated cos and ucos weights
      weights <- ucos_details$ucos_weight
      ucos_weights <- lapply(1:cb_output$data_info$n_outcomes, function(l){
        pp <- which(sapply(ucos_outcome_index, function(oo) oo == l))
        if (length(pp) == 0){
          return(rep(0, length(weights[[1]])))
        } else {
          w <- weights[pp]
          w <- do.call(cbind, w)
          return(1 - apply(1-w, 1, prod))
        }
      })
      cos_vcp <- plot_input$cos_vcp
      combined <- lapply(1:cb_output$data_info$n_outcomes, function(l){
        ucos_w <- ucos_weights[[l]]
        cos_w <- cos_vcp[[l]]
        1 - (1-ucos_w)*(1-cos_w)
      })
      names(combined) <- names(plot_input$cos_vcp)
      plot_input$ucos_cos_int_weights <- combined
    }
  }
  
  # check if plot cos to uncolocalized outcome
  # use the updated coloc_cos and related components from plot_input if available
  coloc_cos <- plot_input$cos
  coloc_index <- plot_input$coloc_index
  coloc_hits <- plot_input$cos_hits
  if (show_cos_to_uncoloc & !is.null(coloc_cos)) {
    if (is.null(show_cos_to_uncoloc_idx)) {
      cos_to_uncoloc <- coloc_cos
      cos_idx_to_uncoloc <- 1:length(coloc_index)
      if (is.null(show_cos_to_uncoloc_outcome)) {
        message("Show all CoSs to uncolocalized outcomes.")
        outcome_to_uncoloc <- lapply(coloc_index, function(cidx) {
          setdiff(1:length(analysis_outcome), cidx)
        })
      } else {
        message("Show all CoSs to uncolocalized outcomes ", paste(show_cos_to_uncoloc_outcome, collapse = ","))
        outcome_to_uncoloc <- lapply(coloc_index, function(cidx) {
          setdiff(show_cos_to_uncoloc_outcome, cidx)
        })
      }
    } else {
      if (show_cos_to_uncoloc_idx > length(coloc_cos)) {
        warning(
          "There are only ", length(coloc_cos), " CoS in this region. ",
          "Cannot show the ordered ", paste(show_cos_to_uncoloc_idx, collapse = ","),
          " CoS, which does not exist"
        )
        outcome_to_uncoloc <- cos_to_uncoloc <- cos_idx_to_uncoloc <- NULL
      } else {
        cos_idx_to_uncoloc <- show_cos_to_uncoloc_idx
        cos_to_uncoloc <- coloc_cos[show_cos_to_uncoloc_idx]
        if (is.null(show_cos_to_uncoloc_outcome)) {
          message(
            "Show the ordered ", paste(cos_idx_to_uncoloc, collapse = ","),
            " CoS for all uncolocalized outcomes."
          )
          outcome_to_uncoloc <- sapply(show_cos_to_uncoloc_idx, function(idx) {
            l <- list(setdiff(1:length(analysis_outcome), coloc_index[[idx]]))
            names(l) <- names(coloc_index[idx])
            return(l)
          })
        } else {
          message(
            "Show the ordered ", paste(cos_idx_to_uncoloc, collapse = ","),
            " CoS for outcomes ", paste(show_cos_to_uncoloc_outcome, collapse = ",")
          )
          outcome_to_uncoloc <- sapply(show_cos_to_uncoloc_idx, function(idx) {
            l <- list(setdiff(show_cos_to_uncoloc_outcome, coloc_index[[idx]]))
            names(l) <- names(coloc_index[idx])
            return(l)
          })
        }
      }
    }
    if (!is.null(outcome_to_uncoloc)) {
      cos_uncoloc_texts <- lapply(1:length(outcome_to_uncoloc), function(idx){
        cname <- names(outcome_to_uncoloc)[idx]
        if (grepl("^ucos", cname)){
          texts_tmp <- "Uncolocalized effect (uCoS)"
        } else {
          texts_tmp <- getOption("colocboost_plot.uncoloc_cos_label", "Uncolocalized effect (CoS)")   # legend wording can be overridden per figure; default unchanged
        }
        return(rep(texts_tmp, length(outcome_to_uncoloc[[idx]])))
      })
      names(cos_uncoloc_texts) <- names(outcome_to_uncoloc)
    } else {
      cos_uncoloc_texts <- NULL
    }
    draw_uncoloc <- sapply(outcome_to_uncoloc, length) != 0
    if (any(draw_uncoloc)) {
      pos <- which(draw_uncoloc)
      outcome_to_uncoloc <- outcome_to_uncoloc[pos]
      cos_to_uncoloc <- cos_to_uncoloc[pos]
      cos_idx_to_uncoloc <- cos_idx_to_uncoloc[pos]
      cos_uncoloc_texts <- cos_uncoloc_texts[pos]
    }
    uncoloc <- list(
      "outcome_to_uncoloc" = outcome_to_uncoloc,
      "cos_to_uncoloc" = cos_to_uncoloc,
      "cos_idx_to_uncoloc" = cos_idx_to_uncoloc,
      "cos_uncoloc_texts" = cos_uncoloc_texts
    )
    plot_input$uncoloc <- uncoloc
  }
  
  plot_input$variant_coord <- variant_coord
  class(plot_input) <- "colocboost"
  return(plot_input)
}


#' @importFrom stats pnorm
plot_initial <- function(cb_plot_input, y = "log10p",
                         points_color = "grey80", cos_color = NULL,
                         ylim_each = TRUE, title_specific = NULL,
                         outcome_legend_size = 1.5,
                         outcome_legend_pos = "right",
                         cos_legend_pos = "bottomleft",
                         plot_ucos = FALSE,
                         show_variable = FALSE,
                         lab_style = c(2, 1),
                         axis_style = c(1, 1),
                         title_style = c(2.5, 2),
                         ...) {
  args <- list(...)
  if (!exists("pch", args)) args$pch <- 16
  
  # - set background point color and cos color pools
  args$bg <- points_color
  if (is.null(cos_color)) {
    
    # cos_color <- c(
    #   "#377EB8", "#E69F00", "#33A02C", "#984EA3", "#F46D43",
    #   "#A65628", "#1F4E79", "#B2182B", "#D73027", "#F781BF",
    #   "#4DAF4A", "#E41A1C", "#FF7F00", "#6A3D9A", "#1B7837",
    #   "#E6AB02", "#542788", "#74ADD1", "#A50F15", "#01665E"
    # )
    cos_color <- c(
      # Original 20 colors (preserved exactly)
      "#377EB8", "#E69F00", "#33A02C", "#984EA3", "#F46D43",
      "#A65628", "#1F4E79", "#B2182B", "#D73027", "#F781BF",
      "#4DAF4A", "#E41A1C", "#FF7F00", "#6A3D9A", "#1B7837",
      "#E6AB02", "#542788", "#74ADD1", "#A50F15", "#01665E",
      
      # Additional 80 colors - carefully selected to complement the original palette
      # Blues and Teals
      "#5DADE2", "#85C1E9", "#2E86C1", "#1B4F72", "#154360",
      "#AED6F1", "#3498DB", "#21618C", "#1A5490", "#17A2B8",
      "#20B2AA", "#008B8B", "#4682B4", "#87CEEB", "#B0E0E6",
      
      # Greens
      "#58D68D", "#82E0AA", "#27AE60", "#1E8449", "#186A3B",
      "#ABEBC6", "#2ECC71", "#229954", "#196F3D", "#52C41A",
      "#90EE90", "#32CD32", "#228B22", "#006400", "#9ACD32",
      
      # Reds and Pinks
      "#EC7063", "#F1948A", "#CB4335", "#A93226", "#922B21",
      "#FADBD8", "#E74C3C", "#C0392B", "#943126", "#FF6B6B",
      "#FF69B4", "#FFB6C1", "#FF1493", "#DC143C", "#B22222",
      
      # Oranges and Yellows
      "#F8C471", "#FAD7A0", "#E67E22", "#CA6F1E", "#B7950B",
      "#FCF3CF", "#F39C12", "#D68910", "#B7950B", "#FFA500",
      "#FFFF00", "#FFD700", "#FFA500", "#FF8C00", "#DAA520",
      
      # Purples and Violets
      "#BB8FCE", "#D2B4DE", "#8E44AD", "#7D3C98", "#6C3483",
      "#E8DAEF", "#9B59B6", "#7E57C2", "#673AB7", "#9C27B0",
      "#8A2BE2", "#9370DB", "#BA55D3", "#DA70D6", "#DDA0DD",
      
      # Browns and Earth Tones
      "#CD853F", "#D2691E", "#A0522D", "#8B4513", "#654321",
      "#DEB887", "#F4A460", "#BC8F8F", "#D2B48C", "#F5DEB3",
      "#8B7355", "#A0522D", "#CD853F", "#DEB887", "#D2691E",
      
      # Grays and Neutrals
      "#85929E", "#AEB6BF", "#5D6D7E", "#34495E", "#2C3E50",
      "#D5D8DC", "#BDC3C7", "#95A5A6", "#7F8C8D", "#17202A",
      "#808080", "#A9A9A9", "#C0C0C0", "#696969", "#778899"
    )
    
    
    # cos_color <- c("#1F70A9", "#33A02C", "#CAB2D6", "#EA7827")
  }
  args$col <- cos_color[cb_plot_input$select_cos]
  
  # - set data and x-lab and y-lab
  if (y == "log10p") {
    plot_data <- lapply(cb_plot_input$Zscores, function(z) {
      -log10(2 * pnorm(-abs(z)))
    })
    ylab <- expression(-log[10](P))
  } else if (y == "z_original") {
    plot_data <- cb_plot_input$Zscores
    ylab <- "Z score"
  } else if (y == "cos_vcp") {
    if (plot_ucos){
      plot_data <- cb_plot_input$ucos_cos_int_weights
      ylab <- "Integrated VPA from CoS and uCoS"
    } else {
      plot_data <- cb_plot_input$cos_vcp
      ylab <- "CoS-specific VCP"
    }
    args$ylim <- c(0, 1)
  } else if (y == "vcp") {
    plot_data <- cb_plot_input$vcp
    ylab <- "VCP"
    if (length(cb_plot_input$outcomes) == 1) {
      ylab <- "VPA"
    }
    args$ylim <- c(0, 1)
  }else if (y == "coef") {
    plot_data <- cb_plot_input$coef
    ylab <- "Coefficients"
  } else {
    stop("Invalid y value! Choose from 'log10p', 'z_original', 'vcp', 'coef', or 'cos_vcp'!")
  }
  if (!exists("xlab", args)) args$xlab <- "variables"
  if (!exists("ylab", args)) args$ylab <- ylab
  args$lab_size <- as.numeric(lab_style[1])
  args$lab_face <- lab_style[2]
  
  # - set title format
  args$title <- title_specific
  args$title_size <- as.numeric(title_style[1])
  args$title_face <- title_style[2]
  
  # - set x-axis and y-axis
  args$x <- cb_plot_input$x$pos
  args$y <- plot_data
  if (show_variable) {
    args$xaxt <- "n"
    args$xtext <- cb_plot_input$variables
  }
  args$axis_size <- as.numeric(axis_style[1])
  args$axis_face <- axis_style[2]
  
  # - set ylim for each subfigure
  if (exists("ylim", args)) {
    ymax <- rep(args$ylim[2], length(args$y))
    ymin <- rep(args$ylim[1], length(args$y))
  } else {
    ymax <- NULL
    ymin <- rep(0, length(args$y))
    
    # Check if ylim_each is FALSE but no ylim is provided
    if (!ylim_each) {
      warning("Since you want to plot with the same y-axis range (ylim_each=FALSE), but ylim was not provided, ", 
              "still plotting with trait-specific y ranges (equivalent to ylim_each=TRUE).")
      ylim_each <- TRUE
    }
  }
  if (ylim_each & is.null(ymax)) {
    ymax <- sapply(plot_data, function(p) {
      valid_values <- p[is.finite(p)]
      if (length(valid_values) == 0) {
        ymax <- 330 # Default if no valid values
      } else {
        ymax <- max(valid_values) * 1.05
      }
      return(ymax)
    })
    if (y == "coef") {
      ymin <- sapply(plot_data, function(p) min(p) * 1.05)
    }
  }
  args$ymax <- ymax
  args$ymin <- ymin
  
  # - set legend text position and format
  args$outcome_legend_pos <- outcome_legend_pos
  args$outcome_legend_size <- outcome_legend_size
  if (outcome_legend_pos == "right") {
    args$outcome_legend_angle <- 90
  } else if (outcome_legend_pos == "left") {
    args$outcome_legend_angle <- 270
  } else {
    args$outcome_legend_angle <- 0
  }
  args$cos_legend_pos <- cos_legend_pos
  
  return(args)
}

