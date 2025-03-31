require(pedtools)

draw_link <- function(id1, id2, p, ped, label, delta_control = p$boxh,
                      anchor1 = "r", anchor2 = "l"){
  par(mar = p$scaling$mar, xpd=TRUE)

  idx1 <- match(id1, ped$ID)
  idx2 <- match(id2, ped$ID)

  w <- p$boxw
  h <- p$boxh

  ax <- setNames(c(-1, 1, 0, 0, -1, 1), nm = c("l", "r", "b", "t", "bl", "br"))
  ay <- setNames(c(1, 1, 2, -2, 2, 2), nm = c("l", "r", "b", "t", "bl", "br"))

  x1 <- p$x[idx1] + ax[anchor1] * w/2
  y1 <- p$y[idx1] + ay[anchor1] * h/2

  x2 <- p$x[idx2] + ax[anchor2] * w/2
  y2 <- p$y[idx2] + ay[anchor2] * h/2

  bezier_curve <- function(t, p0, p1, delta_control){
    p_control <- (p0 + p1) / 2 + delta_control
    (1-t)^2 * p0 + 2*(1-t)*t * p_control + t^2 * p1
  }

  t_plot <- seq(from = 0, to = 1, length = 101)

  x_plot <- bezier_curve(t = t_plot, p0 = x1, p1 = x2, delta_control = 0)
  y_plot <- bezier_curve(t = t_plot, p0 = y1, p1 = y2, delta_control = delta_control)

  lines(x_plot, y_plot, lty=2)

  x_mid <- x_plot[50]
  y_mid <- y_plot[50]

  rect(x_mid - w, y_mid - h/2,
       x_mid + w, y_mid + h/2, col = "lightblue", border = "blue")
  text(x_mid, y_mid, label, col = "blue")
}


close_ped <- halfSibPed(nch2 = 2, sex1 = 2) |> addChildren(father = 6) |> addChildren(father = 8) |>
  addChildren(mother = 4, sex = 2) |> addChildren(mother = 12, sex = 2)


close_ped <- relabel(close_ped, new = c(1, 2, 3, 5, 6, 7, 8, 11, 12, 14, 4, 10 ,9, 13))

pdf("fig1_peds.pdf", width = 8.5, height = 7)

  par(mfrow=c(1,3))

  p <- plot(close_ped, draw=TRUE)#, title="Siblings and nephews/nieces")

  draw_link("6", "7", p, close_ped, "FS", delta_control = -0.225)
  draw_link("5", "6", p, close_ped, "HS", delta_control = -0.225)
  draw_link("10", "6", p, close_ped, "HN", delta_control = -0.5)
  draw_link("11", "6", p, close_ped, "N", delta_control = -0.5, anchor1 = "l", anchor2 = "r")

  draw_link("13", "6", p, close_ped, "HGN", anchor1 = "r", anchor2 = "bl", delta_control = +0.4)
  draw_link("14", "6", p, close_ped, "GN", anchor1="l", anchor2 = "br", delta_control = +0.4)

  # illustrate the cousin peds
  ped_cousins <- cousinPed(degree = 5, symmetric = TRUE)
  ped <- ped_cousins
  p <- plot(ped_cousins, draw = TRUE)#, title = "Cousins")
  draw_link("8", "9", p, ped, "1C")
  draw_link("12", "13", p, ped, "2C")
  draw_link("16", "17", p, ped, "3C")
  draw_link("20", "21", p, ped, "4C")
  draw_link("23", "24", p, ped, "5C")

  # illustrate removal
  ped_1c_removal<- cousinPed(degree = 1, removal = 4, symmetric = TRUE)
  p <- plot(ped_1c_removal, draw = TRUE)#, title = "Cousins removed")
  draw_link("7", "16", p, ped_1c_removal, "1C4R", delta_control = p$boxh * 10)
  draw_link("7", "14", p, ped_1c_removal, "1C3R", delta_control = p$boxh * 6)
  draw_link("7", "12", p, ped_1c_removal, "1C2R", delta_control = p$boxh * 3)
  draw_link("7", "10", p, ped_1c_removal, "1C1R", delta_control = p$boxh * 1)
  draw_link("7", "8", p, ped_1c_removal, "1C", delta_control = p$boxh * 0)
dev.off()
