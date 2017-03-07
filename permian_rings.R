#import packages
library("readxl")
library("ggplot2")
library("ggfortify")
library("dplyr")
library("tidyr")
library("zoo")
library("assertthat")
library("GGally")

#import data
permian_rings <- read_excel("data/L&R2017.data.xlsx")

names(permian_rings)[1:2] <- c("index", "mean_curve")

assert_that(!anyDuplicated(names(permian_rings)))#Check no duplicate column names

#check standardised
tolerance = .Machine$double.eps^0.5 # assumes Matlab has same precision

assert_that(all(colMeans(permian_rings[, -(1:2)], na.rm = TRUE) < tolerance))

permian_rings[, -1] <- scale(permian_rings[, -1])
assert_that(all(sapply(permian_rings[, -1], sd, na.rm = TRUE) - 1 < tolerance))



inMean <- data_frame(
  tree = c("K6049.1", "K6052", "K3733", "K4960.2", "K349", "HOG-01.2", "K6044.2", "K621b", "KH0025", "KH0067_2", "K6004a", "K6050", "K2542", "K4842", "K6046.1", "K1177", "K4613", "K6047", "K3257"), 
  used = c(rep(TRUE, 11), rep(FALSE, 8))
)

stdMean <- rowMeans(permian_rings[, inMean$tree[inMean$used]], na.rm = TRUE)

#correlation with original mean curve
cor(stdMean, permian_rings$mean_curve)
cor(stdMean, permian_rings$mean_curve)^2


#plot data
gather(permian_rings, key = tree, value = width, -index, -mean_curve) %>%
  left_join(inMean, by = c("tree" = "tree")) %>%
  ggplot(aes(x = index, y = width, colour = tree)) + 
    geom_path(show.legend = FALSE) +
    geom_path(data = permian_rings, aes(y = mean_curve), size = 2, colour = "black") +
    geom_path(data = data_frame(index = permian_rings$index, stdMean), aes(y = stdMean), size = 2, colour = "red") + 
  facet_wrap(~used)

gather(permian_rings, key = tree, value = width, -index, -mean_curve) %>%
  left_join(inMean, by = c("tree" = "tree")) %>%
  filter(!is.na(width)) %>% 
  ggplot(aes(x = index, y = tree, fill = width)) + 
  geom_raster() +
  scale_fill_gradient2(mid = "grey95") +
  facet_wrap(~used, ncol = 1, scales = "free_y", strip.position = "right") +
  theme_bw()

## sampling depth
samp_depth <- data_frame(
  index = permian_rings$index, 
  depth = rowSums(!is.na(permian_rings[, inMean$tree[inMean$used]]))
)

ggplot(samp_depth, aes(x = index, y = depth)) + 
  geom_bar(stat = "identity", width = 1)

##agreement

samp_depth$agreement <- apply(permian_rings[, inMean$tree[inMean$used]], 1, sd, na.rm = TRUE)

ggplot(samp_depth, aes(x = index, y = agreement, colour = depth)) + 
  geom_line()


## crossdating
geom_ccf <- function(data, mapping, ..., lim = NA){  
  x  <- eval(mapping$x, data)
  y <-  eval(mapping$y, data)
  CCF <- ccf(x, y, na.action = na.omit, plot = FALSE, ...)
  autoplot(CCF) + ylim(-lim, lim)
}
pr <- permian_rings[, inMean$tree[inMean$used]]

max_ccf <- function(dat, ...){
  res <- sapply(2:ncol(dat), function(i){
    sapply(1:(ncol(dat) - 1), function(j){
      if(i > j){
        CCF <- ccf(dat[,i], dat[,j], na.action = na.omit, plot = FALSE, ...)
        max(abs(CCF$acf), na.rm = TRUE)
       } else{
        NA
      }
    })
  })
  max(res, na.rm = TRUE)
}

gpairs_lower <- function(g){
  g$plots <- g$plots[-(1:g$nrow)]
  g$yAxisLabels <- g$yAxisLabels[-1]
  g$nrow <- g$nrow -1
  
  g$plots <- g$plots[-(seq(g$ncol, length(g$plots), by = g$ncol))]
  g$xAxisLabels <- g$xAxisLabels[-g$ncol]
  g$ncol <- g$ncol - 1
  
  g
}

g <- ggpairs(setNames(pr, make.names(names(pr))), 
        lower = list(continuous = wrap(geom_ccf, lag.max = 20, lim = max_ccf(pr, max.lag = 20))),
        upper = list(continuous = "blank"),
        diag = list(continuous = 'blankDiag')
        )

gpairs_lower(g)

## simple spectrum
spec <- spectrum(permian_rings$mean_curve, plot = FALSE)

g <- autoplot(spec) +
  geom_vline(xintercept = 1/11, linetype = "dashed", colour = "grey60") +
  geom_vline(xintercept = 1/3, linetype = "dashed", colour = "red")
g  


spec2 <- spectrum(stdMean, plot = FALSE)

g %+% fortify(spec2)  

bothSpec <- rbind(
  cbind(curve = "original", fortify(spec)),
  cbind(curve = "recalc", fortify(spec2))
)  

g %+% bothSpec + aes(colour = curve)

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


##correlations with the mean record
cors <- sapply(select(permian_rings, -mean_curve, -index), cor, permian_rings$mean_curve, use="pair")
ggplot(data.frame(cors, used = inMean$used), aes(x = cors, fill = used)) + geom_histogram(bins = 10)

split(cors, inMean$used)


