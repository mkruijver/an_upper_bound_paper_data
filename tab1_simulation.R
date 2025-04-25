ped <- pedtools::cousinPed()

i_v <- inheritance_space(pedigree = ped, ids = pedtools::leaves(ped),
                         states = "v", exploit_symmetries = FALSE)
i_ibd <- inheritance_space(pedigree = ped, ids = pedtools::leaves(ped),
                           exploit_symmetries = FALSE)

chromosome_length <- 100.0

set.seed(2)

r_v <- ibdsegments:::random_ibd(n = 1L, chromosome_length = chromosome_length,
           ibd_state_by_v = i_v$ibd_state_by_v,
           number_of_transmissions = i_v$number_of_relevant_transmissions,
           fixed_transmission_masks = i_v$relevant_masks,
           state_stats = 1L)$samples
dec2bin <- Vectorize(function(x) paste(as.integer(rev(intToBits(x)[1:8])), collapse = " "))
r_v$state <- paste0(r_v$state, ": [", dec2bin(r_v$state), "]")

set.seed(2)
r_ibd <- ibdsegments:::random_ibd(n = 1L, chromosome_length = chromosome_length,
                               ibd_state_by_v = i_ibd$ibd_state_by_v,
                               number_of_transmissions = i_ibd$number_of_relevant_transmissions,
                               fixed_transmission_masks = i_ibd$relevant_masks,
                               state_stats = 1L)$samples

r_v
r_ibd

r_v[-c(1,2)] |>
xtable::xtable() |>
  print(type = "latex", include.rownames = FALSE, booktabs = TRUE)



r_ibd[-c(1,2)] |>
  xtable::xtable() |>
  print(type = "latex", include.rownames = FALSE, booktabs = TRUE)


