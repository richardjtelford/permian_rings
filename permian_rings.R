#import packages
library("readxl")
library("ggplot2")
library("ggfortify")
library("dplyr")
library("tidyr")
library("zoo")

#import data
permian_rings <- read_excel("data/L&R2017.data.xlsx")

names(permian_rings)[1:2] <- c("index", "mean_curve")
table(names(permian_rings))#k6050 repeated

plot(permian_rings[,  names(permian_rings) == "K6050"])

#plot data
gather(permian_rings, key = tree, value = width, -index, -mean_curve) %>%
  ggplot(aes(x = index, y = width, colour = tree)) + 
    geom_path() +
    geom_path(data = permian_rings, aes(y = mean_curve), size = 2, colour = "black") +
    geom_path(data = cbind(permian_rings, mean = rowMeans(permian_rings[, -(1:2)], na.rm = TRUE)), aes(y = mean), size = 2, colour = "red")

## simple spectrum
spec <- spectrum(permian_rings$mean_curve, plot = FALSE)

g <- autoplot(spec) +
  geom_vline(xintercept = 1/11, linetype = "dashed", colour = "grey60") +
  geom_vline(xintercept = 1/3, linetype = "dashed", colour = "red")
g  

#acf
autoplot(acf(permian_rings$mean_curve))
ar(permian_rings$mean_curve)


# white noise smoothed
tst <- rnorm(nrow(permian_rings))# * 100)
tst_3smooth <- rollmean(x = tst, k = 3)#rolling mean
spec3 <- spectrum(tst_3smooth, plot = FALSE)
g %+% fortify(spec3)

##comparison with sunspots
sun_spec <- spectrum(sunspot.year, plot = FALSE)
g %+% fortify(sun_spec)
