#' Phytoplankton Class Abundance Figure
#'
#' This function plots the class abundances as output by `simulated_annealing`.
#'
#' @param c_matrix C matrix of class abundance concentrations
#'
#' @return A stacked line plot with sample number on x axis, chl a 
#'         concentrations on y axis, and phytoplankton groups as colors
#' @examples
#' # ADD_EXAMPLES_HERE
phyto_figure <- function(c_matrix) {

  # NOTE: NULL assignment to stop NOTE during the package "R-CMD-check"
  #   error -  `no visible binding for global variable`
  vals    <- NULL
  row_num <- NULL

  
  if (length(names(c_matrix)) <= 11) {
    color_pal <- c(
      "#999999", "#E69F00", "#56B4E9", "#009E73",
      "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#009E73",
      "#001E73", "#013E73"
    )
  } else if (length(names(c_matrix)) <= 15) {
    # color palette from:
    # https://msk678.github.io/Farbenblind.html#large-palette-15-colors
    message("Using 15 color palette for `class abundance plot`.")
    color_pal <- c(
      "#004949", "#009292", "#FF6DB6", "#FFB6DB", "#490092",
      "#006DDB", "#B66DFF", "#6DB6FF", "#B6DBFF", "#920000",
      "#924900", "#DB6D00", "#24FF24", "#FFFF6D", "#999999"
    )
  } else {
    message("More than 15 pigments used, color palette set to NULL.")
    color_pal <- NULL
  }

  PLE <- tidyr::pivot_longer(
    data      = cbind(c_matrix, "row_num" = seq(nrow(c_matrix))),
    cols      = -row_num,
    names_to  = "names",
    values_to = "vals"
  )

  plt <-
    ggplot2::ggplot(PLE, ggplot2::aes(x = row_num, y = vals, fill = names)) +
    ggplot2::geom_area() +
    ggplot2::scale_y_continuous(
      limits = c(0, NA),
      expand = ggplot2::expansion(mult = c(0, 1.2))
    ) +
    ggplot2::labs(
      x = "Sample number",
      y = expression("Chl a concentrations" ~ (mg ~ m^-3)),
      fill = "Phytoplankton\nNames"
    ) +
    ggplot2::theme_bw()

  if (length(names(c_matrix)) <= 15) {
    plt <-
      plt +
      ggplot2::scale_color_manual(values = color_pal) +
      ggplot2::scale_fill_manual(values = color_pal)
  }
  return(plt)
}



#' Convergence Figure
#'
#' A figure to show the pigment ratios for each phytoplankton group for each
#' iteration.
#'
#' @param fm_iter A data.frame with columns of iter, phyto, pigment and ratio
#' @param niter Optional: the number of iterations on the x axis. 
#'              If `NULL`, will extract from the `iter` column of `fm_iter`.
#'
#' @return A figure with each pigment ratio per iteration per group
#' @examples
#' # ADD_EXAMPLES_HERE
convergence_figure <- function(fm_iter, niter = NULL) {
  
  # NOTE: NULL assignment to stop NOTE during the package "R-CMD-check"
  #   error -  `no visible binding for global variable`
  iter    <- NULL
  ratio   <- NULL
  pigment <- NULL

  
  if (length(unique(fm_iter$pigment)) <= 11) {
    color_pal <- c(
      "#999999", "#E69F00", "#56B4E9", "#009E73",
      "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#009E73",
      "#001E73", "#013E73"
    )
  } else if (length(unique(fm_iter$pigment)) <= 15) {
    # color palette from:
    # https://msk678.github.io/Farbenblind.html#large-palette-15-colors
    message("Using 15 color palette for `convergence plot`.")
    color_pal <- c(
      "#004949", "#009292", "#FF6DB6", "#FFB6DB", "#490092",
      "#006DDB", "#B66DFF", "#6DB6FF", "#B6DBFF", "#920000",
      "#924900", "#DB6D00", "#24FF24", "#FFFF6D", "#999999"
    )
  } else {
    message("More than 15 pigments used, color palette set to NULL.")
    color_pal <- NULL
  }

  if (is.null(niter)) niter <- max(fm_iter$iter)

  # add pretty breaks
  break_pts <- if (niter < 10) {
    seq(0, niter, 1)
  } else if (niter < 50) {
    seq(0, niter, 5)
  } else if (niter <= 100) {
    seq(0, niter, 10)
  } else if (niter <= 500) {
    seq(0, niter, 50)
  } else {
    pretty(x = c(0, niter), 5)
  }

  converge_plt <-
    ggplot2::ggplot(fm_iter, ggplot2::aes(x = iter, y = ratio, color = pigment)) +
    ggplot2::geom_line() +
    ggplot2::geom_point() +
    ggplot2::facet_wrap(~phyto) +
    ggplot2::scale_x_continuous(breaks = break_pts) +
    ggplot2::scale_y_continuous(
      limits = c(0, NA),
      expand = ggplot2::expansion(mult = c(0, 1.2))
    ) +
    ggplot2::labs(
      x = "Iteration Number",
      y = "Pigment Ratio",
      color = "Pigment"
    ) +
    ggplot2::theme_bw() +
    ggplot2::coord_cartesian(expand = FALSE)

  if (length(unique(fm_iter$pigment)) <= 15) {
    converge_plt <- 
      converge_plt +
      ggplot2::scale_color_manual(values = color_pal) 
  }
  
  converge <- list("F_mat_iter" = fm_iter, "converge_plot" = converge_plt)

  return(converge)
}
