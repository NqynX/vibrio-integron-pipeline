source("R/read_integronfinder.R"); source("R/plot_integron.R")
df <- read_integronfinder("data/example/N16961.integrons")
big <- df |> dplyr::filter(type_elt == "attc") |>
  dplyr::count(ID_replicon, ID_integron, sort = TRUE) |> dplyr::slice(1)
d <- prepare_integron(df, big$ID_integron, big$ID_replicon)
s <- attr(d, "start")
p <- plot_region_tracks(d, from = s, to = s + 20000,
       title = "V. cholerae N16961 superintegron — first 20 kb")
dir.create("results/visualisation", showWarnings = FALSE, recursive = TRUE)
ggplot2::ggsave("results/visualisation/N16961_detail_v0.2.png", p,
                width = 12, height = 7, dpi = 150)
cat("wrote results/visualisation/N16961_detail_v0.2.png\n")
