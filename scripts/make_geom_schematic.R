# Figure 2 (Section 2.3): the model curve, velocity and acceleration, the
# caustic and the Maxwell set, for xi = 1.2 and phi = 70 degrees.
this_dir <- dirname(gsub("~+~", " ", sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)), fixed = TRUE))
root <- dirname(this_dir)
args <- commandArgs(TRUE)
out <- if (length(args)) args[1] else file.path(root, "figures", "geom_schematic.pdf")
xi <- 1.2; phi <- 70*pi/180; cs <- cos(phi); sn <- sin(phi)
u <- function(t) t + 0.5*t^2*xi*cs; v <- function(t) 0.5*t^2*xi*sn
tt <- seq(-4, 4, length = 801)
dx <- 1 + tt*xi*cs; dy <- tt*xi*sn; sp <- sqrt(dx^2 + dy^2)
kap <- xi*sn/sp^3
ex <- u(tt) - dy/sp/kap; ey <- v(tt) + dx/sp/kap
ts <- -cs/xi
i0 <- which.min(abs(tt - ts)); cusp <- c(ex[i0], ey[i0])
z <- cusp + 1.3*c(cs, sn)
r <- sort(Re(polyroot(c(-z[1], 1 - xi*cs*z[1] - xi*sn*z[2], 1.5*xi*cs, 0.5*xi^2))))[c(1, 3)]
# three-root region from the discriminant (11)
gx <- seq(-3.4, 3.4, length = 400); gy <- seq(-1.8, 4.4, length = 400)
A <- xi^2/2; B <- 1.5*xi*cs
disc <- outer(gx, gy, function(a, b) {
  C <- 1 - xi*cs*a - xi*sn*b; K <- -a
  18*A*B*C*K - 4*B^3*K + B^2*C^2 - 4*A*C^3 - 27*A^2*K^2 })
# same sans-serif font, palette and axis style as the ggplot figures
COL_CURVE <- "#0072B2"; COL_CAUSTIC <- "#E69F00"; COL_MAXWELL <- "#A50F15"; COL_N <- "#009E73"
pdf(out, width = 6.2, height = 5, family = "Helvetica")
par(mar = c(3.4, 3.4, 0.5, 0.5), mgp = c(2.1, 0.6, 0), las = 1, tcl = -0.3)
plot(NA, xlim = c(-3.2, 3.2), ylim = c(-1.6, 4.2), asp = 1,
     xlab = expression(z[u] ~ "(along the tangent at" ~ c[0] * ")"), ylab = expression(z[v] ~ "(along the normal at" ~ c[0] * ")"))
# opaque tint (the caustic colour at 25% on white): transparency renders unevenly across PDF viewers
.filled <- function() image(gx, gy, (disc > 0) + 0, col = c(NA, "#F9E7BF"), add = TRUE)
.filled()
contour(gx, gy, disc, levels = 0, drawlabels = FALSE, add = TRUE, col = COL_CAUSTIC, lwd = 1.5)
th <- seq(0, 2*pi, length = 200)
for (R in 1:2) lines(R*cos(th), R*sin(th), col = "grey45", lty = 2, lwd = 0.8)
text(c(.72, 1.43), c(-.72, -1.43), c("r = 1", "r = 2"), cex = .65, col = "grey45", pos = 4)
abline(h = 0, lty = 2, col = "grey50"); abline(v = 0, lty = 2, col = "grey50")
lines(u(tt), v(tt), col = COL_CURVE, lwd = 2.5)
text(u(-2.6) - .1, v(-2.6), expression(gamma[L](c)), pos = 2, col = COL_CURVE)
segments(cusp[1], cusp[2], cusp[1] + 3.4*cs, cusp[2] + 3.4*sn, col = COL_MAXWELL, lwd = 2)
segments(z[1], z[2], u(r), v(r), col = "grey30")
points(u(r), v(r), pch = 21, bg = "white"); points(z[1], z[2], pch = 19, cex = .8)
arrows(0, 0, 1, 0, length = .08, lwd = 2)
arrows(0, 0, xi*cs, xi*sn, length = .08, lwd = 2, col = COL_N)
a <- seq(0, phi, length = 30); lines(.45*cos(a), .45*sin(a), col = COL_N)
text(.6*cos(phi/2), .6*sin(phi/2), expression(varphi), col = COL_N)
text(1.05, -.22, "T", font = 2); text(xi*cs + .2, xi*sn + .05, "N", col = COL_N, font = 2)
points(0, 0, pch = 19); text(0, 0, expression(gamma[L](c[0])), pos = 1, offset = .9, cex = .8)
text(cusp[1] + 3.1*cs, cusp[2] + 3.1*sn, "Maxwell set", pos = 4, cex = .75, col = COL_MAXWELL)
text(2.6, 3.0, "caustic", cex = .75, col = COL_CAUSTIC)
text(z[1] + .05, z[2], "z", pos = 4)
dev.off()
cat(sprintf("wrote %s\n", out))
