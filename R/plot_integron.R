# Build and draw a map of ONE integron from a parsed .integrons data frame.
suppressPackageStartupMessages({library(dplyr); library(ggplot2)})

classify <- function(df) {
  df |> mutate(category = dplyr::case_when(
    type_elt == "attc" ~ "attC",
    grepl("inti", tolower(paste(annotation, model))) ~ "integrase",
    type_elt == "protein" ~ "cassette",
    TRUE ~ "other"))
}

prepare_integron <- function(df, integron_id = NULL, replicon = NULL) {
  d <- classify(df)
  if (!is.null(replicon))    d <- filter(d, ID_replicon == replicon)
  if (!is.null(integron_id)) d <- filter(d, ID_integron == integron_id)
  d <- filter(d, !is.na(pos_beg), !is.na(pos_end))
  if (nrow(d) == 0) stop("no features for that selection")
  a0 <- min(d$pos_beg); a1 <- max(d$pos_end)
  inti <- filter(d, category == "integrase")
  if (nrow(inti) > 0 && mean(inti$pos_beg) > (a0 + a1) / 2) {   # put intI on the left
    d <- d |> mutate(nb = a1 - pos_end + a0, ne = a1 - pos_beg + a0,
                     pos_beg = nb, pos_end = ne,
                     strand = ifelse(is.na(strand), strand, -strand)) |>
      select(-nb, -ne)
  }
  attr(d, "start") <- min(d$pos_beg); attr(d, "end") <- max(d$pos_end)
  d
}

plot_integron <- function(d, row_width = 25000, title = NULL) {
  s <- attr(d, "start"); e <- attr(d, "end")
  d <- d |> mutate(
    mid = (pos_beg + pos_end) / 2,
    row = floor((mid - s) / row_width),
    xb  = pmax(pos_beg - s - row * row_width, 0),
    xe  = pmin(pos_end - s - row * row_width, row_width),
    fwd = is.na(strand) | strand >= 0,
    x0  = ifelse(fwd, xb, xe), x1 = ifelse(fwd, xe, xb),
    y   = -row)
  genes <- filter(d, category %in% c("cassette", "integrase"))
  attc  <- filter(d, category == "attC")
  nrows <- max(d$row) + 1
  labs  <- data.frame(y = -(0:(nrows - 1)),
                      lab = paste0(round((s + (0:(nrows - 1)) * row_width) / 1000),
                                   " kb"))
  cols <- c(cassette = "#4C78A8", integrase = "#E45756")
  ggplot() +
    geom_segment(data = labs, aes(0, y, xend = row_width, yend = y),
                 colour = "grey88", linewidth = 0.4) +
    geom_segment(data = genes,
                 aes(x0, y, xend = x1, yend = y, colour = category),
                 linewidth = 3.5, lineend = "butt",
                 arrow = arrow(length = unit(0.12, "cm"), type = "closed")) +
    geom_point(data = attc, aes((xb + xe) / 2, y - 0.30),
               shape = 18, size = 2, colour = "#F58518") +
    geom_text(data = labs, aes(-1600, y, label = lab),
              hjust = 1, size = 2.6, colour = "grey45") +
    scale_colour_manual(values = cols, name = NULL,
                        breaks = c("integrase", "cassette"),
                        labels = c("intI integrase", "cassette CDS")) +
    scale_x_continuous(expand = expansion(mult = c(0.14, 0.02))) +
    coord_cartesian(clip = "off") +
    labs(title = title,
         subtitle = sprintf("%d cassette CDS  |  %d attC (orange)  |  %.1f kb",
                            sum(genes$category == "cassette"), nrow(attc),
                            (e - s) / 1000),
         x = "position within row (bp)", y = NULL) +
    theme_minimal(base_size = 11) +
    theme(axis.text.y = element_blank(), panel.grid = element_blank(),
          legend.position = "top", plot.subtitle = element_text(colour = "grey45"))
}

# ---- v0.2: layered multi-track region view (ggcoverage-style) ----
suppressPackageStartupMessages({library(patchwork)})

# Build 7-point arrow polygons for a set of genes (proper gene glyphs).
gene_arrows <- function(g, hb = 0.22, ht = 0.34, head_bp = 300) {
  do.call(rbind, lapply(seq_len(nrow(g)), function(i) {
    x0 <- g$gx0[i]; x1 <- g$gx1[i]; fwd <- isTRUE(g$fwd[i])
    hl <- min(head_bp, x1 - x0)
    if (fwd) { hs <- x1 - hl
      px <- c(x0, hs, hs, x1, hs, hs, x0)
    } else {  hs <- x0 + hl
      px <- c(x1, hs, hs, x0, hs, hs, x1) }
    py <- c(-hb, -hb, -ht, 0, ht, hb, hb)
    data.frame(id = g$element[i], x = px, y = py, fill = g$fill[i])
  }))
}

plot_region_tracks <- function(d, from = NULL, to = NULL, title = NULL) {
  s <- attr(d, "start"); e <- attr(d, "end")
  if (is.null(from)) from <- s
  if (is.null(to))   to   <- min(e, s + 20000)
  d <- dplyr::mutate(d, gx0 = pmin(pos_beg, pos_end),
                        gx1 = pmax(pos_beg, pos_end),
                        fwd = is.na(strand) | strand >= 0)
  genes <- d |> dplyr::filter(category %in% c("cassette","integrase"),
                              gx1 >= from, gx0 <= to) |>
    dplyr::mutate(fill = ifelse(category == "integrase", "intI",
                         ifelse(fwd, "forward", "reverse")))
  attc  <- d |> dplyr::filter(category == "attC",
                              pos_beg >= from, pos_beg <= to) |>
    dplyr::mutate(mid = (pos_beg + pos_end)/2)
  cds   <- genes |> dplyr::filter(category == "cassette") |>
    dplyr::mutate(mid = (gx0 + gx1)/2, len = gx1 - gx0)
  arr <- gene_arrows(genes)

  pal <- c(intI = "#E45756", forward = "#4C78A8", reverse = "#72B7B2")
  xs  <- ggplot2::scale_x_continuous(
           limits = c(from, to), expand = ggplot2::expansion(mult = 0.01),
           labels = function(v) round(v/1000, 1))
  blank_x <- ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                            axis.title.x = ggplot2::element_blank())

  t1 <- ggplot2::ggplot(arr, ggplot2::aes(x, y, group = id, fill = fill)) +
    ggplot2::geom_polygon(colour = "grey25", linewidth = 0.15) +
    ggplot2::scale_fill_manual(values = pal, name = NULL) + xs +
    ggplot2::coord_cartesian(ylim = c(-0.55, 0.55)) +
    ggplot2::labs(title = title, y = "cassettes") +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                   panel.grid = ggplot2::element_blank(),
                   legend.position = "top") + blank_x

  t2 <- ggplot2::ggplot(attc) +
    ggplot2::geom_segment(ggplot2::aes(mid, 0, xend = mid, yend = 1),
                          colour = "#F58518", linewidth = 0.4) +
    ggplot2::geom_point(ggplot2::aes(mid, 1), shape = 18, size = 2.2,
                        colour = "#F58518") + xs +
    ggplot2::coord_cartesian(ylim = c(0, 1.25)) +
    ggplot2::labs(y = "attC") + ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                   panel.grid = ggplot2::element_blank()) + blank_x

  t3 <- ggplot2::ggplot(cds) +
    ggplot2::geom_col(ggplot2::aes(mid, len, fill = fill), width = 280) +
    ggplot2::scale_fill_manual(values = pal, guide = "none") + xs +
    ggplot2::labs(y = "length (bp)", x = "genomic position (kb)") +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())

  t1 / t2 / t3 + patchwork::plot_layout(heights = c(3, 1, 2))
}
