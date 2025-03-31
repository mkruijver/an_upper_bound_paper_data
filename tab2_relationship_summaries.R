## For all relationships, list out kappa_1 and the moments of total IBD

source("relationships.R")

df_sharing <- data.frame(Relationship = names(peds))

df_sharing$kappa1 <- sapply(names(peds), function(relationship_name){

  paste0("1/",
         as.character(1/ibdsegments::d_ibd(ibd = 1,
                                           pedigree = peds[[relationship_name]],
                                           ids = persons[[relationship_name]],
                                           states = "kappa")))
})

moments_by_rel <- sapply(names(peds), function(relationship_name){
  ibdsegments::total_ibd_dist_moments(pedigree = peds[[relationship_name]],
                                      ids = persons[[relationship_name]],
                                      chromosome_length = L)
})

df_sharing$total_shared_mean <- unlist(moments_by_rel["mean", ])
df_sharing$total_shared_sd <- unlist(moments_by_rel["sd", ])

compute_p0 <- function(relationship_name){
  prod(sapply(L, function(l){
    ibdsegments::d_cibd(x = l, ibd = 0, pedigree = peds[[relationship_name]],
           ids = persons[[relationship_name]])
  }))
}

df_sharing$p0 <- sapply(names(peds), compute_p0)

df_sharing$p0 <- ifelse(df_sharing$p0 < 1e-4,
                        format(df_sharing$p0, scientific = TRUE, digits = 4),
                        format(round(df_sharing$p0, digits = 4)))

pmf_num0 <- function(relationship_name){

  pmf_num0_by_chromosome <- lapply(L, function(l) {

    p0_l <- ibdsegments::d_cibd(x = l, ibd = 0,
                                pedigree = peds[[relationship_name]],
                                ids = persons[[relationship_name]])

    list(x = c(0, 1),
         px = c(1-p0_l, p0_l))
  })

  Reduce(add_dists, pmf_num0_by_chromosome)
}

add_dists <- function(a,b){

  ibdsegments:::pmf_of_sum(x1 = a$x,
                           p1 = a$px,
                           x2 = b$x,
                           p2 = b$px, eps = 0)
}

n0_pmfs_by_relationship_name <- sapply(names(peds), pmf_num0, simplify = FALSE)

pmf_exp <- function(pmf) sum(pmf$x * pmf$px)
pmf_exp2 <- function(pmf) sum(pmf$x^2 * pmf$px)

df_sharing$n0_expected <- sapply(n0_pmfs_by_relationship_name, pmf_exp)
df_sharing$n0_sd <- sqrt(sapply(n0_pmfs_by_relationship_name, pmf_exp2) -
                           sapply(n0_pmfs_by_relationship_name, pmf_exp)^2)

o <- df_sharing$Relationship[order(df_sharing$total_shared_mean, df_sharing$total_shared_sd, decreasing = TRUE)]

df_sharing_sorted <- df_sharing[match(o, df_sharing$Relationship),]

setNames(df_sharing_sorted,
         c("Relationship", "$kappa_1$", "Total IBD (mean)", "Total IBD (sd)", "Pr(Total IBD = 0)", "N0 (mean)", "N0 (sd)")) |>
  xtable::xtable() |>
  print(type = "latex", include.rownames = FALSE, booktabs = TRUE)

