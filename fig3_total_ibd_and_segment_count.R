
source("relationships.R")

require(ibdsegments)

number_of_samples <- 1e5

set.seed(1)
sims <- list()
for (relationship_name in names(peds)){
  cat(relationship_name, " ")
  r <- r_cibd(number_of_samples,
              pedigree = peds[[relationship_name]],
              ids = persons[[relationship_name]], chromosome_length = L)$stats

  sims[[relationship_name]] <- data.frame(relationship = relationship_name,
                                          r, stringsAsFactors = FALSE)
}

sims_merged <- do.call(rbind, sims)

# readr::write_csv(sims_merged, "sims_1_1e5.csv")
#
# sims_merged <- readr::read_csv("sims_1_1e5.csv")

plot_df <- tidyr::pivot_longer(sims_merged, cols = c(total_length, segments),
                               names_to = "variable",
                               values_to = "value")

# reorder factor labels
o <- c("GP", "HS", "N", "GGP", "HN", "GN", "1C", "G^2GP", "H1C",
       "G^2N", "1C1R", "G^3GP",
       "H1C1R",
       "G^3N",
       "2C", "G^4GP", "H2C", "G^4N", "2C1R", "G^5GP", "H2C1R",
       "G^5N", "3C", "G^6GP", "H3C", "G^6N", "3C1R", "G^7GP",
       "H3C1R", "G^7N", "4C", "G^8GP", "H4C", "G^8N", "4C1R",
       "G^9GP", "H4C1R", "G^9N", "5C")


plot_df$relationship_ordered <- factor(plot_df$relationship, levels = rev(o))
plot_df$variable <- factor(plot_df$variable, levels = c("total_length", "segments"))

require(ggplot2)

ggplot(plot_df) + aes(y=relationship_ordered, x = value) +
  geom_boxplot(outlier.shape = NA) +
  theme_bw() +
  xlab("") +
  ylab("Relationship") + facet_wrap(~ variable, scales= "free",
                                    labeller = labeller(variable = c("total_length" = "Total IBD (cM)",
                                                                     "segments" = "Segment count"))) +
  scale_x_sqrt()

ggsave(filename = "fig3_sim_total_ibd_and_segments_by_relationship.pdf",
       width = 8)
