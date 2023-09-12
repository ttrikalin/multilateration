
p <- c(0, 0, -3)
B <- fix_beacons(4)
r <- distance(P1 = p, P2 = B)

tri_params <- trilateration_params(distances = r, beacons = B, base_point_index = 1)

x<- c(0, 0, 0)

norm((tri_params$A %*% x - tri_params$c_values))
