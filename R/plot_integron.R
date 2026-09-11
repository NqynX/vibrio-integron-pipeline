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
