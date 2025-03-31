rm(list=ls())

source("relationships.R")

require(ibdsegments)

relationship_names <- paste0(1:5, "C")
number_of_simulations <- 1e4

sim_dfs <- list()
for (i_true_relationship in seq_along(relationship_names)){
  set.seed(i_true_relationship)

  true_relationship <- relationship_names[i_true_relationship]
  sim <- r_cibd(n = number_of_simulations,
                pedigree = peds[[true_relationship]],
                ids = persons[[true_relationship]], chromosome_length = L)

  sim_df <- data.frame(true_relationship = rep(true_relationship,
                                               number_of_simulations))


  for (relationship in relationship_names){
    cat(true_relationship, " ll: ", relationship, "\n")
    ll <- d_cibd(sim, pedigree = peds[[relationship]],
                 ids = persons[[relationship]], log10 = TRUE)

    sim_df[[relationship]] <- ll
  }

  sim_dfs[[true_relationship]] <- sim_df
}

sims <- do.call(rbind, sim_dfs)

# readr::write_csv(sims, "cousin_logliks_1e4.csv")
# rm(list=ls())
# sims <- readr::read_csv("cousin_logliks_1e4.csv")

df_log10lrs <- do.call(rbind,lapply(seq_len(length(relationship_names) - 1), function(i_relationship){
  h1 <- relationship_names[i_relationship]
  h2 <- relationship_names[i_relationship + 1]
  hyps <- paste0("H1: ", h1, ", H2: ", h2)

  rbind(data.frame(Category="H1 true",
                   Stat = "Continuous IBD",
                   Hypotheses = hyps,
                   log10lr = sims[[h1]][sims$true_relationship==h1] -
                     sims[[h2]][sims$true_relationship==h1]),
        data.frame(Category="H2 true",
                   Stat = "Continuous IBD",
                   Hypotheses = hyps,
                   log10lr = sims[[h1]][sims$true_relationship==h2] -
                     sims[[h2]][sims$true_relationship==h2]))
}))


require(ggplot2)

ggplot(df_log10lrs) + aes(x = log10lr, y = Hypotheses, col= Category)
  geom_boxplot() +
  theme_bw() +
  xlab(expression(log[10](LR)))

ggsave("fig6_cousins.pdf", width = 8, height = 3.5)

# numerical summary
stats = aggregate(log10lr ~ Category + Hypotheses, data = df_log10lrs, FUN = function(x) mean(x>0))[c(2,1,3)]
names(stats)[3] <- "Pr(LR>1)"
stats[["Pr(LR<1)"]] <- aggregate(log10lr ~ Category + Hypotheses, data = df_log10lrs, FUN = function(x) mean(x<0))$log10lr
stats[["median LR"]] <- aggregate(log10lr ~ Category + Hypotheses, data = df_log10lrs, FUN = median)$log10lr

kableExtra::kbl(data.frame(H1 = rep(paste0(1:4,"C"), each = 2),
                           H2 = rep(paste0(2:5,"C"), each = 2),
                           stats[-1], check.names = FALSE)
                , booktabs = T, digits = 4, format = "latex")


h <- unique(df_log10lrs$Hypotheses)[1]
accuracy <-
  sapply(unique(df_log10lrs$Hypotheses), function(h){
    accuracy_h1 <- mean(df_log10lrs$log10lr[df_log10lrs$Hypotheses==h & df_log10lrs$Category == "H1 true"] > 0)
    accuracy_h2 <- mean(df_log10lrs$log10lr[df_log10lrs$Hypotheses==h & df_log10lrs$Category == "H2 true"] < 0)

    mean(c(accuracy_h1, accuracy_h2))
  })
dput(accuracy)
