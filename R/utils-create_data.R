denormalize <- function(unit_interval_points, range) {
  return(unit_interval_points * (range[2]-range[1]) + range[1])
}

# fictitious data for K beacons
fix_beacons <- function(num_beacons,
                        general_position = FALSE,
                        X = c(0,   10),
                        Y = c(0,    7),
                        Z = c(0, -0.5)){
  stopifnot(num_beacons>=3)
  p <- if (general_position) 3 else 2
  initial_design <- MaxPro::MaxProLHD(n = num_beacons, p = p)$Design
  design <- MaxPro::MaxPro(InitialDesign = initial_design, s = 2)$Design
  if (p == 2) design <- cbind(design, rep(0, num_beacons))
  colnames(design) <- c("x", "y", "z")
  rownames(design) <- paste0("b", 1:num_beacons)
  design[, "x"] <- denormalize(unit_interval_points = design[, "x"], range = X)
  design[, "y"] <- denormalize(unit_interval_points = design[, "y"], range = Y)
  design[, "z"] <- denormalize(unit_interval_points = design[, "z"], range = Z)
  return(design)
}

distance <- function(P1, P2){
  if(is.vector(P1)) {
    P1 <- matrix(P1, nrow = 1)
  }
  if(is.vector(P2)) {
    P2 <- matrix(P2, nrow = 1)
  }
  D <- as.matrix(dist(rbind(P1, P2), diag = TRUE))
  D <- matrix(D[1:nrow(P1), (nrow(P1)+1):nrow(D)], nrow = nrow(P1))
  colnames(D) <- rownames(P2)
  rownames(D) <- rownames(P1)
  return(D)
}

change_origin <- function(points , origin){
  points - matrix(rep(origin, nrow(points)), nrow=nrow(points), byrow = TRUE)
}
