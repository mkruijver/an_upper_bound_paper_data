
set.seed(1)

x <- ibdsegments::r_cibd(n = 1,
                    pedigree = pedtools::cousinPed(),
                    chromosome_length = 100)

x$samples[-c(1,2)] |>
xtable::xtable() |>
  print(type = "latex", include.rownames = FALSE, booktabs = TRUE)
