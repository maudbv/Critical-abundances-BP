## circular stats:

th  <- rnorm(100, 1, 4) %% (2*pi)
err <- rnorm(100, mean = 0, sd = 0.8)
icp <- 10

bc <- 2
bs <- 3

plot(th, cos(th))
plot(th, sin(th))
plot(abs(th-3.11), 34*cos(th) + 0.01*sin(th) )


y   <- icp + bc * cos(th) + bs * sin(th) + err

plot(y, cos(th))

plot(th, y)

x = db$ASPECT
z = db$SRali
plot(cos(x), z)
plot(abs(x - 180), z)

plot(x, abs(x - 180))


plot(x, cos(x))
plot(abs(x-180), cos(x))
plot(x, sin(x))
plot(abs(x-180), sin(x))


plot(abs(x - 90), x)

## Converting to radians in order to calculate a shift in aspect with 0 = East
radians = degrees × (π / 180°)
xx = x * (pi/180)
# shift by -π/2 to transform E = 0, S = π/2, W = π and N = 2π/4
yy = xx - pi/2

# convert back to degrees:
degrees = radians / ( π / 180°)  = radians * 180°/ π
y = yy*180/pi

plot(x, xx)
plot(x, xx - pi/4)

plot(cos(xx),cos(xx + pi/4))

zo <- cos(x) + sin(x)

plot(zo, z)

envplot@data$Northern <- cos(envplot@data$ASPECT* (pi/180))
envplot@data$Eastern <- sin(envplot@data$ASPECT* (pi/180))
spplot(envplot, zcol = "Northern", pretty = T, cuts= 6, colorkey = TRUE, col.regions = terrain.colors(6))
spplot(envplot, zcol = "Eastern", pretty = T, cuts= 6, colorkey = TRUE, col.regions = terrain.colors(6))

spplot(envplot, zcol = c("Eastern", "Northern"), pretty = T, cuts= 6, colorkey = TRUE, col.regions = terrain.colors(6))


hist(envplot@data$Northern)
hist( abs(envplot@data$ASPECT-180))
plot(envplot@data$ASPECT, envplot@data$Northern)
plot(abs(envplot@data$ASPECT-180),envplot@data$Northern)
plot(abs(envplot@data$ASPECT-180), sin(envplot@data$ASPECT))
plot(envplot@data$ASPECT, cos(envplot@data$ASPECT - pi/2))

theta.i <- envplot@data$ASPECT*pi/180
atan2(sum(sin(theta.i)),sum(cos(theta.i))) * 180/pi

x = 0:359
plot(x, abs(x-180))
plot(x, cos(x*pi/180))
plot(x,  sin(x*pi/180))
plot(x, cos(x*pi/180) + sin(x*pi/180))

plot(abs(x-180), cos(x*pi/180))
