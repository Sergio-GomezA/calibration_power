cdf1 <- function(x) 0.4 + 0.2 * (1 - exp(-sqrt(6 * x)))

inv_cdf1 <- function(u) (-log(1 - 5 * (u - 0.4)))^2 / 6

inv_cdf1(cdf1(3.5))

inv_cdf1_t <- function(u) (-log(5 * (u - 0.4)))^2 / 6

inv_cdf1_t(cdf1(3.5))


cdf2 <- function(x) 0.8 + 0.2 * (1 - exp(-sqrt(6 * x)))

inv_cdf2 <- function(u) (-log(1 - 5 * (u - 0.8)))^2 / 6

inv_cdf2(cdf2(5))

inv_cdf2_t <- function(u) (-log(5 * (1 - u)))^2 / 6

inv_cdf2_t(cdf2(5))


inv_cdf2_t <- function(u) (-log(5 * (1 - u)))^2 / 6


inv_cdf2_t(.9280)
inv_cdf2_t(5 * .9280 - 4)


U <- runif(1000, 0, 1)

V <- (U[U > 0.8] - 0.8) / 0.2

plot(density(V))

x <- seq(0, 1, length.out = 100)

plot(sapply(x, \(x) sum(V <= x)))

cor(U[U > 0.8], V)

(-log(5 - 5 * 0.83754864))^2 / 6

f_q2 <- function(u, beta = 2, gamma = 3) gamma + beta * log(u / (1 - u))

f_q2(0.8807)

f_q2(0.1760700)


f_q3 <- function(u) (-log(1 - u))^2 / 6
f_q3(0.8930) * 60
f_q3(0.1148) * 60
f_q2(0.747, gamma = 5)

f_q4 <- function(x) (sin(2 * x) + cos(2 * x) + 2) / 2 / pi

f_q4(2.55483605)


f_q4(1.9161)
k <- (2 + sqrt(2)) / 2 / pi

0.375 * k


(-log(1 - (0.92802 - 0.8) / 0.2))^2 / 6


t <- sqrt(2) - 1

# Create grid
x <- seq(-1.5, 1.5, length.out = 1000)
y <- seq(-1.5, 1.5, length.out = 1000)

grid <- expand.grid(x = x, y = y)

# Check inequalities
grid$inside <- (abs(grid$x) + t * abs(grid$y) <= 1) &
    (t * abs(grid$x) + abs(grid$y) <= 1)

grid$inside <- pmax(
    abs(grid$x) + t * abs(grid$y),
    t * abs(grid$x) + abs(grid$y)
) <=
    1
grid$inside <- (abs(grid$x) + abs(grid$y) <= 1 + t) &
    (abs(grid$x) <= 1) &
    (abs(grid$y) <= 1)

grid$inside <- (abs(grid$x) + abs(grid$y) <= 1 + (4 - sqrt(2)) / 7) &
    (abs(grid$x) <= 1) &
    (abs(grid$y) <= 1)

grid$inside <- (abs(grid$x) + abs(grid$y) <= 1 + (sqrt(2) - 1) / 2) &
    (abs(grid$x) <= 1) &
    (abs(grid$y) <= 1)

grid$inside <- (abs(grid$x) + abs(grid$y) <= 1 + 1 / sqrt(2)) &
    (abs(grid$x) <= 1) &
    (abs(grid$y) <= 1)
# Plot
plot(
    grid$x,
    grid$y,
    col = ifelse(grid$inside, "steelblue", "white"),
    pch = 15,
    cex = 0.3,
    asp = 1,
    xlab = "x",
    ylab = "y",
    main = "Region defined by the inequalities"
)

box()


plot(cos((1 + 2 * (0:7)) * pi / 8), sin((1 + 2 * (0:7)) * pi / 8))

sqrt(2) - 1
