lower.x <- 2.92

upper.x <- 100

step <- (upper.x - lower.x) / 100

sigma <- 1

df <- 2

bounds <- c(0, 0+3*sigma)

cord.x <- c(lower.x,seq(lower.x,upper.x,step),upper.x)

cord.y <- c(0,dt(seq(lower.x,upper.x,step),df),0)

curve(dt(x,df),xlim=bounds, xlab="value", ylab="density") 

polygon(cord.x,cord.y,col='skyblue')

abline(v=0)
