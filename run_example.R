# First milestone: render the N16961 superintegron.
source("R/read_integronfinder.R"); source("R/plot_integron.R")
df <- read_integronfinder("data/example/N16961.integrons")
cat("rows:", nrow(df), " integrons:", length(unique(df$ID_integron)), "\n")
print(table(df$type_elt))
big <- df |> filter(type_elt == "attc") |>
  count(ID_replicon, ID_integron, sort = TRUE) |> slice(1)
cat("plotting", big$ID_integron, "on", big$ID_replicon, "-", big$n, "attC\n")
d <- prepare_integron(df, big$ID_integron, big$ID_replicon)
p <- plot_integron(d, 25000, "Vibrio cholerae N16961 superintegron")
dir.create("results/visualisation", showWarnings = FALSE, recursive = TRUE)
ggsave("results/visualisation/N16961_v0.1.png", p, width = 11, height = 7, dpi = 150)
cat("wrote results/visualisation/N16961_v0.1.png\n")
