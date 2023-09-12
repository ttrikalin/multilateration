set.seed(12345)

# example data for P(t)
transmitter_ss <- function(t, period = 10, mean = -10, amplitude = 12.5) {
  amplitude*sin(2*pi*t/period)/2 + mean
}
delta_esp32   <- -20
delta         <- -10
eta_p         <- 1.9
sigma_p       <- .5
tmax          <- 120
distance_monitors <- c(.25, 1.5, 3 )

receiver_ss <- function(distance, transmitter_logpower, delta, eta = eta_p, sigma = sigma_p) {
  transmitter_logpower + delta - eta*log(distance) + rnorm(n = length(distance), mean = 0, sd = sigma)
}
phone_position1 <- function(t, xmin = 0.1, xmax = 2.5, tmax= 120) {
  xmin + (xmax-xmin) *t / tmax
}
phone_position2 <- function(t, xmin = 0.1, xmax = 2.5, tmax = 120) {
  x <- phone_position1(t=t , xmin=xmin, xmax = xmax, tmax=tmax)
  2*(x * (t<=(tmax/2)) + (xmax - x + xmin) * (t>(tmax/2)))
}

phone_position <- phone_position2


# generate data for three monitors
# measurements times M1 to M3
t1 <- seq.int(from = 0, to = tmax, length.out = tmax)
t2 <- seq.int(from = runif(n=1), by = 1, length.out = tmax)
t3 <- seq.int(from = runif(n=1), by = 1, length.out = tmax)
m1 <- receiver_ss(distance = distance_monitors[1],
                  transmitter_logpower = transmitter_ss(t=t1),
                  delta = delta_esp32)
m2 <- receiver_ss(distance = distance_monitors[2],
                  transmitter_logpower = transmitter_ss(t=t2),
                  delta = delta_esp32)
m3 <- receiver_ss(distance = distance_monitors[3],
                  transmitter_logpower = transmitter_ss(t=t3),
                  delta = delta_esp32)

# generate phone data
t0 <- seq.int(from = 0, to = tmax , by =2)
phone <- receiver_ss(distance = phone_position(t=t0),
                     transmitter_logpower = transmitter_ss(t=t0),
                     delta = delta)




#### Analysis
library(laGP)


### estimate the eta_p(t)
### here I need smoothers of the m1, m2, m3
### but i use approxfuns
m1.est <- approxfun(x= t1, y= m1, rule = 2)
m2.est <- approxfun(x= t2, y= m2, rule = 2)
m3.est <- approxfun(x= t3, y= m3, rule = 2)

m.diff21 <- m2.est(t.obs) - m1.est(t.obs)
m.diff31 <- m3.est(t.obs) - m1.est(t.obs)
eta.est31 <- - m.diff31 / (log(distance_monitors[3] / distance_monitors[1]))
eta.est21 <- - m.diff21 / (log(distance_monitors[2] / distance_monitors[1]))
eta.est <- mean(c(eta.est21, eta.est31))

#### Estimate the power function over time
P.obs1 <- m1  -delta_esp32 + eta.est * log(distance_monitors[1])
P.obs2 <- m2  -delta_esp32 + eta.est * log(distance_monitors[2])
P.obs3 <- m3  -delta_esp32 + eta.est * log(distance_monitors[3])

P.obs <- c(P.obs1, P.obs2, P.obs3)
t.obs <- c(t1, t2, t3)[order(c(t1, t2, t3))]
P.obs <- P.obs[order(c(t1, t2, t3))]

tt <- as.matrix(seq(0, tmax, length.out = 1000))

P.alc <- laGP::aGP(X = as.matrix(t.obs),
             Z = P.obs,
             XX = tt,
             omp.threads = 1,
             verb = 0,
             d= 2)

plot(t.obs, P.obs)
lines(y=P.alc$mean, x=tt, col = "red")


### Maybe not needed, but can package the P.alc predictions in
### a function
P.est <- approxfun(x= tt, y = P.alc$mean, rule = 2)






### Now that we have the P.est model we estimate the position function.
log.x <- ( P.est(t0) - phone + delta)/eta.est


### fit a GP (or a lowess or another smoother, or a Kalman filter) to the
### position function -- here GP.
### Rationale: Position does not change too fast.

lx.alc <- laGP::aGP(X = as.matrix(t0),
                   Z = log.x,
                   XX = tt,
                   verb = 0)
plot(y = log.x, x = t0)
lines(y= lx.alc$mean, x= tt)
lx.est <- approxfun(x= tt, y = lx.alc$mean, rule = 2)
lx.se.est <- approxfun(x= tt, y = sqrt(lx.alc$var), rule = 2)


# This is the position with 50% confidence region
plot(t0, phone_position(t0))
lines(tt, exp(lx.est(tt)))
lines(tt, exp(lx.est(tt) + qnorm(0.25) *lx.se.est(tt)))
lines(tt, exp(lx.est(tt) - qnorm(0.25) *lx.se.est(tt)))

