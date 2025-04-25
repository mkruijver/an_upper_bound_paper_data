source("relationships.R")

relationship_names <- c("HS", "N", "GP")

x_plot <- seq(from = 0, to = 300, length = 500)

require(ibdsegments)
cdfs <- pdfs <- list()

pdfs <- sapply(relationship_names, function(relationship_name){
  kappa1 <- d_ibd(ibd = 1, pedigree = peds[[relationship_name]],
                  ids = persons[[relationship_name]])

  Vectorize(function(l)
    d_cibd(x = c(l,0), ibd = c(1,0),
           pedigree = peds[[relationship_name]],
           ids = persons[[relationship_name]]) /
      kappa1)
})

cdfs <- sapply(relationship_names, function(relationship_name){
  Vectorize(function(l) integrate(f = pdfs[[relationship_name]],
                                  lower = 0, upper = l)$val)
})

hazards <- sapply(relationship_names, function(relationship_name){
  function(x) pdfs[[relationship_name]](x) / (1 - cdfs[[relationship_name]](x))
})


plot_dfs <- list()
for (relationship_name in relationship_names){
  plot_dfs[[1 + length(plot_dfs)]] <-
    data.frame(x = x_plot,
               Relationship = relationship_name,
               Variable = "Density",
               y = pdfs[[relationship_name]](x_plot))
}

for (relationship_name in relationship_names){
  plot_dfs[[1 + length(plot_dfs)]] <-
    data.frame(x = x_plot,
               Relationship = relationship_name,
               Variable = "Hazard rate",
               y = hazards[[relationship_name]](x_plot))
}


plot_df <- do.call(rbind, plot_dfs)

require(ggplot2)

ggplot(plot_df) + aes(x=x, y=y, col=Relationship) + geom_line() +
  xlab("Segment length (cM)") + ylab("") +
  facet_wrap(~Variable) + theme_bw()

ggsave("fig4_pdf_and_hazard_rate.pdf", width = 7.5, height = 3.5)
