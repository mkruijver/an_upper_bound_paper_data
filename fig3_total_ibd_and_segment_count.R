source("relationships.R")
require(ibdsegments)

kappa_by_ped <- sapply(names(peds), function(relationship_name){

  paste0("1/",
         as.character(1/ibdsegments::d_ibd(ibd = 1,
                                           pedigree = peds[[relationship_name]],
                                           ids = persons[[relationship_name]],
                                           states = "kappa")))
})


# reorder factor labels
o <- rev(c("GP", "HS", "N", "GGP", "HN", "GN", "1C", "G^2GP", "H1C",
       "G^2N", "1C1R", "G^3GP",
       "H1C1R",
       "G^3N",
       "2C", "G^4GP", "H2C", "G^4N", "2C1R", "G^5GP", "H2C1R",
       "G^5N", "3C", "G^6GP", "H3C", "G^6N", "3C1R", "G^7GP",
       "H3C1R", "G^7N", "4C", "G^8GP", "H4C", "G^8N", "4C1R",
       "G^9GP", "H4C1R", "G^9N", "5C"))

o_split <- split(o, factor(as.character(kappa_by_ped[o]),
       levels = unique(as.character(kappa_by_ped[o]))))
o_split_with_gaps <- list()

for(i in 1:(length(o_split)-1)){
  o_split_with_gaps[[i]] <- c(o_split[[i]], paste0("gap", i))
}
o_split_with_gaps[[length(o_split)]] <- o_split[[length(o_split)]]

relationship_names_with_gaps <- unlist(o_split_with_gaps)

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

plot_df <- tidyr::pivot_longer(sims_merged, cols = c(total_length, segment_count),
                               names_to = "variable",
                               values_to = "value")

plot_df$relationship_ordered <- factor(plot_df$relationship,
                                       levels = relationship_names_with_gaps)

plot_df$variable <- factor(plot_df$variable, levels = c("total_length", "segment_count"))

require(ggplot2)


plot_df_gaps <- rbind(plot_df[-1],
      data.frame(relationship_ordered = rep(grep("gap", relationship_names_with_gaps, value = TRUE), each=2),
           variable = c("total_length", "segment_count"),
           value = NA))


axis_labels <- relationship_names_with_gaps
axis_labels[grepl("^gap", axis_labels)] <- ""


# add shading for clusters
df_cluster_positions <- do.call(rbind, lapply(seq_along(o_split), function(i) {
  cats <- o_split[[i]]
  start <- which(relationship_names_with_gaps == cats[1]) - 0.5
  end <- which(relationship_names_with_gaps == cats[length(cats)]) + 0.5
  data.frame(ymin = start, ymax = end, cluster = paste("Cluster", i))
}))
df_cluster_positions$fill <- c("gray95", "gray90")

sqrt2 <- scales::new_transform(
  "sqrt2",
  transform = \(x) sign(x) * sqrt(abs(x)),
  inverse = \(x) sign(x) * abs(x)^2
)

ggplot() +
  geom_rect(data = df_cluster_positions, aes(xmin = -Inf, xmax = Inf,
  ymin = ymin, ymax = ymax, fill = fill), alpha = 0.7) +
  scale_fill_identity() +  # use exact color codes without legend

  geom_boxplot(data = plot_df_gaps, aes(y=relationship_ordered, x = value),
               outlier.shape = NA, na.rm = TRUE) +
  scale_y_discrete(labels = ifelse(startsWith(relationship_names_with_gaps,
                                              prefix = "gap"),
                                   yes = "",
                                   no = relationship_names_with_gaps)) +

  theme_bw() + theme(axis.ticks.y=element_blank()) +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 10)) +
  xlab("") +
  ylab("Relationship") + facet_wrap(~ variable, scales= "free",
                                    labeller = labeller(variable = c("total_length" = "Total IBD (cM)",
                                                                     "segment_count" = "Segment count"))) +
  coord_trans(x = sqrt2) +
  scale_x_continuous(n.breaks =  7)

ggsave(filename = "fig3_sim_total_ibd_and_segments_by_relationship.pdf",
       width = 8)
