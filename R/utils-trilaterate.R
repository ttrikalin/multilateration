trilateration_params <- function(distances, beacons, base_point_index = 1){
  A <- change_origin(points = B, origin = B[base_point_index,])
  c_values <- lapply(
    X = 1:nrow(B),
    FUN = function(j) {
      0.5*(r[base_point_index]^2 -
             r[j]^2 +
             norm(as.matrix(B[j,])) -
             norm(as.matrix(B[base_point_index,])))
    })
  c_values <- unlist(c_values)
  list(
    base_point_index = base_point_index,
    A = A,
    c_values = c_values
  )
}



library(gurobi)
library(Matrix)

