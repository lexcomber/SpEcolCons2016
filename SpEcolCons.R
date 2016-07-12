## R Analysis
## Lex Comber, a.comber@leeds.ac.uk
## July 2016
## Created specially / spatially for the 
# Spatial Ecology & Conservation 4
# University of Bristol (UK)
# 12 to 14 July 2016

## Keynote talk: 'R in spatial analysis for conservation: space is special'

## Note
# as a reminder a general introduction to R can be found here: https://cran.r-project.org/doc/contrib/Owen-TheRGuide.pdf - work up to page 28 and don't worry if you cannot answer the questions on matrix algebra
# if you really want to go for it then try this excellent book: https://uk.sagepub.com/en-gb/eur/an-introduction-to-r-for-spatial-analysis-and-mapping/book241031 
 
# You will need to the libraries the first time you use them
# However, if you have run the R scripts previously 
# then can jump to the 'library(GISTools)' line
install.packages('GISTools', dep = T) # comment this line out when this has run with a '#'
install.packages('rgdal', dep = T) # comment this line out when this has run with a '#'
install.packages('GWmodel', dep = T) # comment this line out when this has run with a '#'
install.packages('repmis', dep = T) # comment this line out when this has run with a '#'

# you will only need to run the lines below after 
# you have installed the packages on your computer
library(GISTools)
library(rgdal)
library(GWmodel)
library(repmis)

# Part 1. load and prepare the data

# 1.1 load prepared data files
source_data("https://github.com/lexcomber/SpEcolCons2016/blob/master/SpEcolCons.RData?raw=True")

# and have a look at what has been loaded

ls()

# Acknowldgements Rainfall
# https://catalogue.ceh.ac.uk/documents/f2856ee8-da6e-4b67-bedb-590520c77b3c
# You must cite: Tanguy, M.; Dixon, H.; Prosdocimi, I.; Morris, D. G.; Keller, V. D. J. (2015). Gridded estimates of daily and monthly areal rainfall for the United Kingdom (1890-2014) [CEH-GEAR]. NERC Environmental Information Data Centre. http://doi.org/10.5285/f2856ee8-da6e-4b67-bedb-590520c77b3c

# Acknowldgements PET
#http://csi.cgiar.org/Aridity/
#Trabucco, A., and Zomer, R.J. 2009. Global Aridity Index (Global-Aridity) and Global Potential Evapo-Transpiration (Global-PET) Geospatial Database. CGIAR Consortium for Spatial Information. Published online, available from the CGIAR-CSI GeoPortal at: http://www.csi.cgiar.org.
# clip PET to OS

# Acknowledgement Wildness Quality Index
# Kuiters, A. T., van Eupen, M., Carver, S., Fisher, M., Kun, Z., & Vancura, V. (2011). Wilderness register and indicator for Europe. Final Report - http://ec.europa.eu/environment/nature/natura2000/wilderness/pdf/Wilderness_register_indicator.pdf   

# 1.2 interpolate the over OS 10k grid
# wqi at 1k 
# pet at 30 arc second (1k at equator) 
# rain at 1k 
# everything projected to OS

# Aggregate the single variable layers
wqi.agg <- aggregate(wqi[,1], os10, mean)
pet.agg <- aggregate(pet[,1], os10, mean)
rain.agg <- aggregate(rain[,1], os10, mean)

# convert to points
wqi.agg <- SpatialPointsDataFrame(coordinates(wqi.agg), 
	data = wqi.agg@data, proj4string = os.proj)
pet.agg <- SpatialPointsDataFrame(coordinates(pet.agg), 
	data = pet.agg@data, proj4string = os.proj)
rain.agg <- SpatialPointsDataFrame(coordinates(rain.agg), 
	data = rain.agg@data, proj4string = os.proj)

# Aggregate the rainfall data from 'data'
# define a summary function
year.func <- function(x){
	sum(x, na.rm=T)/mean(x, na.rm = T)
}
years <- sort(unique(data.sp$year))
res.vec <- matrix(data = 0, ncol = length(years), nrow = nrow(os10))
for (i in 1:length(years)){
	year.i <- years[i]
	index.i <- which(data.sp$year == year.i)
	data.i <- data.sp[index.i,]
	data.i <- spTransform(data.i, os.proj)
	agg.i <- aggregate(data.i["year"], os10, year.func)
	res.vec[, i] <- agg.i$year
}
colnames(res.vec) <- years
index <- is.na(res.vec)
res.vec[index] <- 0
# create the aggregated variable
tree.agg <- SpatialPointsDataFrame(os10, 
	data.frame(tree = rowSums(res.vec)), 
	proj4string = os.proj)

# check dimensions
dim(wqi.agg)
dim(rain.agg)
dim(pet.agg)
dim(tree.agg)

# 1.3 create a single dataset for use in the analyses
df <- data.frame(tree = tree.agg@data, pet = pet.agg@data, 
	rain = rain.agg@data, wqi = wqi.agg@data)

d.a <- SpatialPointsDataFrame(coordinates(os10), data = data.frame(df))
summary(d.a)
# get rind of any NAs
index<- !is.na(rowSums(d.a@data))
d.a <- d.a[index,]

## Part 2. Analysis 

# 2.1 Initial look at data

# A Initial figures / maps of data

if (.Platform$GUI == "AQUA") {
	quartz(w=5) } else  {x11(w=5) }
par(mar = c(1,1,1,1))
plot(os10m, lwd = 0.5)
plot(data.sp[os10m,], pch = 1, cex = 0.2, col = "#2525254C", add = T)
title("Data points")

if (.Platform$GUI == "AQUA") {
	quartz(w=5) } else  {x11(w=5) }
par(mar = c(1,1,1,1))
plot(os10m, lwd = 0.5)
sh <-auto.shading(d.a$tree[d.a$tree>0], n = 7, cols=brewer.pal(7,'Reds'))
choropleth(d.a, v = d.a$tree, sh, pch = 19, cex = 0.4, add = T)
choro.legend(px = "topright", sh = sh)
title("Trees per 10km sq")

if (.Platform$GUI == "AQUA") {
	quartz(w=5) } else  {x11(w=5) }
par(mar = c(1,1,1,1))
plot(os10m, lwd = 0.5)
sh <-auto.shading(d.a$rain[d.a$rain > 0 ], n = 7,
	cols=brewer.pal(7,'Blues'))
choropleth(d.a, v = d.a$rain, sh, pch = 19, cex = 0.4, add = T)
choro.legend(px = "topright", sh = sh, title = "mm", )
title("Mean monthly Autumn Rain")

if (.Platform$GUI == "AQUA") {
	quartz(w=5) } else  {x11(w=5) }
par(mar = c(1,1,1,1))
plot(os10m, lwd = 0.5)
sh <-auto.shading(d.a$pet[!is.na(d.a$pet)], n = 7,
	cols=brewer.pal(7,'YlOrRd'))
choropleth(d.a, v = d.a$pet, sh, pch = 19, cex = 0.4, add = T)
choro.legend(px = "topright", sh = sh)
title("mean annual PET")

if (.Platform$GUI == "AQUA") {
	quartz(w=5) } else  {x11(w=5) }
par(mar = c(1,1,1,1))
plot(os10m, lwd = 0.5)
sh <-auto.shading(d.a$wqi[!is.na(d.a$wqi)], n = 7,
	cols=brewer.pal(7,'GnBu'))
choropleth(d.a, v = d.a$wqi, sh, pch = 19, cex = 0.4, add = T)
choro.legend(px = "topright", sh = sh, title = "mm")
title("Mean Wildness Quality Index")

# Some descriptive stats
d.a@data[sample(1:nrow(d.a), 6),]
round(cor(d.a@data, 
	use = "pairwise.complete.obs"), 3)

# 2.2 OLS Regression

m <- lm(tree~pet+rain+wqi, data = data.frame(d.a))
summary(m)
coef(m)

# 2.3 GW model with the 'GWmodel' package

# use bw.gwr to optimise the bandwidth
bw2 <- bw.gwr(d.a$tree~d.a$pet+d.a$rain+d.a$wqi, data = d.a, kernel = "bisquare", adaptive = TRUE, approach = "CV") 
bw2

# create GW regression points
hg <- spsample(os10m,3500,'hexagonal',offset=c(0.5,0.5)) # 7500 later?

# create gwr.mod 
m.g2 <- gwr.basic(d.a$tree~d.a$pet+d.a$rain+d.a$wqi, data = d.a, regression.points = hg, bw = bw2, kernel = "bisquare", adaptive = TRUE)

# have a look at the model summary
m.g2

# convert the gw.model SDF to a Spatial Points Data Frame
m.gr2 = SpatialPointsDataFrame(m.g2$SDF,data.frame(m.g2$SDF), proj4string = os.proj)

# 2.4 map the regression coefficient estimates
# PET and GW coefs for PET 
if (.Platform$GUI == "AQUA") {
	quartz(w=8, h = 5.5) } else  {x11(w=8, h=5.5) }
par(mar = c(1,1,1,1))
par(mfrow = c(1,3)) 
cx <- 0.4
sh <-auto.shading(d.a$pet[!is.na(d.a$pet)], n = 6,
	cols=brewer.pal(7,'YlOrRd')[2:7])
plot(os10m, lwd = 0.5)
choropleth(d.a, v = d.a$pet, sh, pch = 19, cex = cx, add = T)
choro.legend(px = "topleft", sh = sh, box.col = NA, cex = 0.7)
txt=expression(paste("Mean annual PET"))
title(txt)

plot(os10m, lwd = 0.5)
plot(data.sp[os10m,], add = T, cex = 0.4, pch = 1, col = "#2525254C")
title("Data Points")

sh = auto.shading(m.gr2$d.a.pet, n = 7, 
	cols = brewer.pal(7, "Spectral"))
plot(os10m, lwd = 0.5, add = F)
choropleth(m.gr2, m.gr2$d.a.pet, sh, pch = 19, cex = cx, add = T)
txt=expression(paste("GW Coeff: ",B[PET]))
choro.legend(px = "topright", sh = sh, box.col = NA, fmt = "%2.2f", cex = 0.7, title = txt)
tit <- sprintf("PET Global Coeff: %s",round(coef(m)[2], 3))
title(tit)

# WQI and GW coefs for WQI
if (.Platform$GUI == "AQUA") {
	quartz(w=8, h = 5.5) } else  {x11(w=8, h=5.5) }
par(mar = c(1,1,1,1))
par(mfrow = c(1,3)) 
cx <- 0.4
sh <-auto.shading(d.a$wqi[!is.na(d.a$wqi)], n = 6,
	cols=brewer.pal(7,'YlOrRd')[2:7])
plot(os10m, lwd = 0.5)
choropleth(d.a, v = d.a$wqi, sh, pch = 19, cex = cx, add = T)
choro.legend(px = "topleft", sh = sh, box.col = NA, cex = 0.7)
txt=expression(paste("Wildness Quality Index"))
title(txt)

plot(os10m, lwd = 0.5)
plot(data.sp[os10m,], add = T, cex = 0.4, pch = 1, col = "#2525254C")
title("Data Points")

sh = auto.shading(m.gr2$d.a.wqi, n = 7, 
	cols = brewer.pal(7, "Spectral"))
plot(os10m, lwd = 0.5, add = F)
choropleth(m.gr2, m.gr2$d.a.wqi, sh, pch = 19, cex = cx, add = T)
txt=expression(paste("GW WQI: ", B[PET]))
choro.legend(px = "topright", sh = sh, box.col = NA, fmt = "%2.2f", cex = 0.7, title = txt)
tit <- sprintf("WQI Global Coeff: %s",round(coef(m)[4], 3))
title(tit)
plot(os10m, lwd = 0.5, add = T)

# 2.5 Examine LOCAL Collinearity
gw.col <- gwr.collin.diagno(d.a$tree~d.a$pet+d.a$rain+d.a$wqi, data = d.a, bw = bw2, kernel = "bisquare", adaptive = TRUE)
### VIFs greater than 10 -> suggestive of collinearity
### CNs greater than 30 -> suggestive of collinearity
summary(gw.col$SDF)
# VIFs & CNs for the variables considered
summary(gw.col$SDF[,1:4])

if (.Platform$GUI == "AQUA") {
	quartz(w=8, h = 5.5) } else  {x11(w=8, h=5.5) }
par(mar = c(1,1,1,1))
par(mfrow = c(1,3)) 
cx <- 0.5

sh = auto.shading(gw.col$SDF@data[,4], n = 7, 
	cols = brewer.pal(7, "Reds"))
plot(os10m, lwd = 0.5, add = F)
choropleth(gw.col$SDF, v = gw.col$SDF@data[,4], sh, pch = 19, cex = cx, add = T)
title("Local CN: all > 30!")
choro.legend(px = "topright", sh = sh, box.col = NA, fmt = "%2.2f", cex = 0.7, title = "CN")
plot(os10m, lwd = 0.5, add = T)

sh = auto.shading(gw.col$SDF@data[,1], n = 5, 
	cols = brewer.pal(5, "Blues"))
sh$breaks <- c(2,5,7,10)
plot(os10m, lwd = 0.5, add = F)
choropleth(gw.col$SDF, v = gw.col$SDF@data[,1], sh, pch = 19, cex = cx, add = T)
title("Local VIF (PET): some > 10")
plot(os10m, lwd = 0.5, add = T)
choro.legend(px = "topright", sh = sh, box.col = NA, fmt = "%2.2f", cex = 0.7, title = "CN")

sh = auto.shading(gw.col$SDF@data[,3], n = 5, 
	cols = brewer.pal(5, "Greens"))
sh$breaks <- c(2,5,7,10)
plot(os10m, lwd = 0.5, add = F)
choropleth(gw.col$SDF, v = gw.col$SDF@data[,3], sh, pch = 19, cex = cx, add = T)
title("Local VIF (WQI): some > 10")
plot(os10m, lwd = 0.5, add = T)
choro.legend(px = "topright", sh = sh, box.col = NA, fmt = "%2.2f", cex = 0.7, title = "CN")

# 2.6 Compensate for local collinearity
# using GW locally comnsated ridge regression function

# First determine the new bandwidth
# This tries to find the bandwidth to fit a ridge term in every location where the CN was greater than 30 
lcr.bw <- bw.gwr.lcr(d.a$tree~d.a$pet+d.a$rain+d.a$wqi, data = d.a, kernel = "bisquare", adaptive = TRUE, lambda.adjust = TRUE, cn.thresh = 30)
lcr.bw

# Then run the locally compensated GWR 
lcr.m <- gwr.lcr(d.a$tree~d.a$pet+d.a$rain+d.a$wqi, data = d.a, regression.points = hg, bw = lcr.bw, kernel = "bisquare", adaptive = TRUE, lambda.adjust = TRUE, cn.thresh = 30)
summary(lcr.m$SDF$Local_CN)
summary(lcr.m$SDF)
# extract the SPDF for analysis
m.gr3 = SpatialPointsDataFrame(lcr.m$SDF,data.frame(lcr.m$SDF), proj4string = os.proj)

# Map Original GWR coefficient estimates and LCR GWR ones
# PET and GW coefs for PET 
if (.Platform$GUI == "AQUA") {
	quartz(w=8, h = 5.5) } else  {x11(w=8, h=5.5) }
par(mar = c(1,1,1,1))
par(mfrow = c(1,3)) 
cx <- 0.4
sh = auto.shading(m.gr2$d.a.pet, n = 7, 
	cols = brewer.pal(7, "Blues"))
plot(os10m, lwd = 0.5, add = F)
choropleth(m.gr2, m.gr2$d.a.pet, sh, pch = 19, cex = cx, add = T)
txt=expression(paste("GW Coeff: ",B[PET]))
choro.legend(px = "topright", sh = sh, box.col = NA, fmt = "%2.2f", cex = 0.7, title = txt)
tit <- sprintf("GWR PET Coeffs")
title(tit)

plot(os10m, lwd = 0.5, add = F)
choropleth(m.gr3, m.gr3$d.a.pet, sh, pch = 19, cex = cx, add = T)
txt=expression(paste("LCR GW Coeff: ",B[PET]))
choro.legend(px = "topright", sh = sh, box.col = NA, fmt = "%2.2f", cex = 0.7, title = txt)
tit <- sprintf("LCR GWR PET Coeffs")
title(tit)

val <- (m.gr2$d.a.pet - m.gr3$d.a.pet)
sh = auto.shading(val, n = 7, 
	cols = brewer.pal(7, "Spectral"))
plot(os10m, lwd = 0.5, add = F)
choropleth(m.gr2, val, sh, pch = 19, cex = cx, add = T)
txt= "Diff"
choro.legend(px = "topright", sh = sh, box.col = NA, fmt = "%2.2f", cex = 0.7, title = txt)
tit <- sprintf("LCR GWR vs GWR PET Coeffs")
title(tit)

# WQI and GW coefs for WQI
if (.Platform$GUI == "AQUA") {
	quartz(w=8, h = 5.5) } else  {x11(w=8, h=5.5) }
par(mar = c(1,1,1,1))
par(mfrow = c(1,3)) 
cx <- 0.4
sh = auto.shading(m.gr2$d.a.wqi, n = 7, 
	cols = brewer.pal(7, "Greens"))
plot(os10m, lwd = 0.5, add = F)
choropleth(m.gr2, m.gr2$d.a.wqi, sh, pch = 19, cex = cx, add = T)
txt=expression(paste("GW Coeff: ",B[PET]))
choro.legend(px = "topright", sh = sh, box.col = NA, fmt = "%2.2f", cex = 0.7, title = txt)
tit <- sprintf("GWR WQI Coeffs")
title(tit)

plot(os10m, lwd = 0.5, add = F)
choropleth(m.gr3, m.gr3$d.a.wqi, sh, pch = 19, cex = cx, add = T)
txt=expression(paste("LCR GW Coeff: ",B[PET]))
choro.legend(px = "topright", sh = sh, box.col = NA, fmt = "%2.2f", cex = 0.7, title = txt)
tit <- sprintf("LCR GWR WQI Coeffs")
title(tit)

val <- (m.gr2$d.a.wqi - m.gr3$d.a.wqi)
sh = auto.shading(val, n = 7, 
	cols = brewer.pal(7, "Spectral"))
plot(os10m, lwd = 0.5, add = F)
choropleth(m.gr2, val, sh, pch = 19, cex = cx, add = T)
txt= "Diff"
choro.legend(px = "topright", sh = sh, box.col = NA, fmt = "%2.2f", cex = 0.7, title = txt)
tit <- sprintf("LCR GWR vs GWR WQI Coeffs")
title(tit)

#### END