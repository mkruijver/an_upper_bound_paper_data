rm(list=ls())

source("relationships.R")

relationship_names_by_inv_kappa <-
  list("2" = c("GP", "HS", "N"),
       "4" = c("GGP", "HN", "GN", "1C"),
       "8" = c("G^2GP", "H1C", "G^2N", "1C1R"))

require(ibdsegments)

inv_kappas <- names(relationship_names_by_inv_kappa)

for (inv_kappa in inv_kappas){

  number_of_simulations <- 1e4 # 1e6 in the paper, runs overnight

  cat("Starting", number_of_simulations, "simulations with kappa_inv", inv_kappa, "\n")

  tic <- Sys.time()
  sim_dfs <- list()

  relationship_names <- relationship_names_by_inv_kappa[[inv_kappa]]

  for (i_true_relationship in seq_along(relationship_names)){
    set.seed(i_true_relationship)

    true_relationship <- relationship_names[i_true_relationship]

    sim <- r_cibd(n = number_of_simulations,
                  pedigree = peds[[true_relationship]],
                  ids = persons[[true_relationship]], chromosome_length = L)

    sim_df <- data.frame(true_relationship = rep(true_relationship,
                                                 number_of_simulations),
                         segment_count = sim$stats$segments,
                         total_ibd = sim$stats$total_length)

    for (relationship in relationship_names){
      cat(true_relationship, " ll: ", relationship, "\n")
      ll <- d_cibd(sim, pedigree = peds[[relationship]],
                   ids = persons[[relationship]], log10 = TRUE)

      sim_df[[relationship]] <- ll
    }

    sim_dfs[[true_relationship]] <- sim_df
  }
  toc <- Sys.time()

  print(toc-tic)
  cat("Finished", number_of_simulations, "\n")


  sims <- do.call(rbind, sim_dfs)

  fn <- gsub("\\+0", "", paste0("kappa_1_", inv_kappa, "_rels_",
                                format(number_of_simulations, scientific = TRUE), ".csv"))

  readr::write_csv(sims, fn)
}

compute_lrs_from_likelihoods <- function(sims, relationship_names){
  number_of_relationships <- length(relationship_names)

  sims <- sims[c("true_relationship", relationship_names)]

  row_max <- apply(sims[-1], 1, max)

  do.call(rbind, lapply(relationship_names, function(relationship_name){
    ll_h1 <- (sims[[relationship_name]] - row_max)
    ll_h2 <- log10(rowMeans(10 ^ (as.matrix(
      sims[-c(1,match(relationship_name, names(sims)))]) - row_max)))

    log10lr <- ll_h1 - ll_h2

    data.frame(Relationship = relationship_name,
               True = ifelse(sims$true_relationship == relationship_name,
                             yes = "H1", no = "H2"),
               Value = log10lr)
  }))
}

sims_by_inv_kappa <- list("2" = readr::read_csv("kappa_1_2_rels_1e4.csv"),
                          "4" = readr::read_csv("kappa_1_4_rels_1e4.csv"),
                          "8" = readr::read_csv("kappa_1_8_rels_1e4.csv"))

relationship_names_by_inv_kappa <-
  list("2" = c("GP", "HS", "N"),
       "4" = c("GGP", "HN", "GN", "1C"),
       "8" = c("G^2GP", "H1C", "G^2N", "1C1R"))


inv_kappas <- c("2", "4", "8")
df_lrs_by_inv_kappa <- lapply(inv_kappas, function(inv_kappa){
  nm <- relationship_names_by_inv_kappa[[inv_kappa]]
  compute_lrs_from_likelihoods(sims_by_inv_kappa[[inv_kappa]][c("true_relationship", nm)], nm)
})

df_lrs <- rbind(data.frame(Kappa="1/2", df_lrs_by_inv_kappa[[1]]),
                data.frame(Kappa="1/4", df_lrs_by_inv_kappa[[2]]),
                data.frame(Kappa="1/8", df_lrs_by_inv_kappa[[3]]))

require(ggplot2)
require(gridExtra)

xlims <- list("1/2" = c(-20, 20),
              "1/4" = c(-10, 10),
              "1/8" = c(-5,5))

plots <- lapply(split(df_lrs, df_lrs$Kappa), function(df_lrs){
  ggplot(df_lrs) + aes(x = Value, col = Relationship, lty = True) +
    stat_density(aes(x=Value), geom= "line", position = "identity", lwd=1) +
    theme_bw() +
    ylab("Density") +
    xlab(if(df_lrs$Kappa[1] == "1/8") expression(log[10](LR)) else NULL ) +
    xlim(xlims[[df_lrs$Kappa[1]]]) +
    facet_wrap(~Kappa, ncol = 1, labeller = label_bquote(kappa[1]==.(Kappa))) + ylab("Density")
})

pdf("fig7_same_kappa_log10_lrs_kde.pdf", width = 8, height = 8)
grid.arrange(plots[[1]], plots[[2]], plots[[3]])
dev.off()

## numerical summary

stats_df <-
  do.call(rbind,lapply(split(df_lrs, df_lrs$Kappa), function(df_lrs){
    stats = aggregate(Value ~ True + Relationship, data = df_lrs, FUN = function(x) mean(x>0))[c(2,1,3)]
    names(stats)[3] <- "Pr(LR>1)"

    stats[["Pr(LR<1)"]] <- aggregate(Value ~ True + Relationship, data = df_lrs, FUN = function(x) mean(x<0))$Value
    stats[["median LR"]] <- aggregate(Value ~ True + Relationship, data = df_lrs, FUN = median)$Value

    accuracy <-
      sapply(unique(df_lrs$Relationship), function(r){
        accuracy_h1 <- mean(df_lrs$Value[df_lrs$Relationship==r & df_lrs$True == "H1"] > 0)
        accuracy_h2 <- mean(df_lrs$Value[df_lrs$Relationship==r & df_lrs$True == "H2"] < 0)

        mean(c(accuracy_h1, accuracy_h2))
      })
    stats$Accuracy <- accuracy[stats$Relationship]

    data.frame(Kappa = df_lrs$Kappa[1],
               H1 = stats$Relationship,
               H2 = paste0("Not ",stats$Relationship),
               stats[-1], check.names = FALSE)
  }))


kableExtra::kbl(stats_df, row.names = FALSE,
                booktabs = T, digits = 4, format = "latex")
