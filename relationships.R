# this script defines all relationships used in the study


# sex-averaged chromosome lengths
L <- c(267.77274, 251.72703, 218.30851, 202.88815, 197.07529, 186.02051,
       178.39839, 161.53634, 157.34513, 169.28027, 154.50387, 165.48962,
       127.22724, 116.00254, 117.31752, 126.58805, 129.52754, 116.51781,
       106.35181, 107.75898, 62.87959, 70.8374)

require(pedtools)
peds <- list()
persons <- list()

# linear relationships
for(d in 2:11){
  if (d==2) nm <- "GP"
  if (d==3) nm <- "GGP"
  if (d>3) nm <- paste0("G^", d-2, "GP")

  peds[[nm]] <- linearPed(n= d)
  persons[[nm]] <- c("1", pedtools::leaves(peds[[nm]]))
}

# nephew (uncle type) relationships
peds[["N"]] <- avuncularPed()
peds[["GN"]] <- avuncularPed(removal = 2)

for (n in 3:10){
  nm <- paste0("G^", n-1, "N")

  peds[[nm]] <- avuncularPed(removal = n)
}

# half sibs and down to (great) half nieces
peds[["HS"]] <- halfSibPed()
peds[["HN"]] <- avuncularPed(half = TRUE)

# cousins
for(cousin_deg in 1:5){
  nm <- paste0(cousin_deg, "C")
  peds[[nm]] <- cousinPed(degree = cousin_deg)

  max_removed <- if(cousin_deg<5) 1 else 0

  if (max_removed > 0){
    for (num_removed in 1:max_removed){
      nm <- paste0(cousin_deg, "C", num_removed, "R")

      peds[[nm]] <- cousinPed(degree = cousin_deg, removal = num_removed)
    }
  }
}

# half cousins
for(cousin_deg in 1:4){
  nm <- paste0("H", cousin_deg, "C")
  peds[[nm]] <- halfCousinPed(degree = cousin_deg)

  max_removed <- 1

  if (max_removed > 0){

    for (num_removed in 1:max_removed){
      nm <- paste0("H", cousin_deg, "C", num_removed, "R")
      peds[[nm]] <- halfCousinPed(degree = cousin_deg, removal = num_removed)
    }
  }
}

for (nm in names(peds)){
  if (is.null(persons[[nm]])){
    persons[[nm]] <- leaves(peds[[nm]])
  }
}
