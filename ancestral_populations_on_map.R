# calculating the number of ancestral populations

library(LEA)

vcf2geno("olfusa.filtered.vcf", "olfusa.geno")
pc = pca("olfusa.geno", scale = TRUE)

tw = tracy.widom(pc)
tw$pvalues[1:5]
plot(tw$percentage)


project = NULL
project = snmf("olfusa.geno",
               K = 6:11,
               entropy = TRUE,
               repetitions = 2500,
               seed = 1236,
               CPU = 4,
               alpha = 100,
               project = "new")

par(family="Times")
plot(project, col = "#d73027", pch = 19, cex = 1.2)

library(dplyr)

# selecting the best run for K = 8
best = which.min(cross.entropy(project, K = 8))
my.colors <- c("#ffffbf", "#abd9e9", "#313695", "#f46d43", 
              "#d73027", "#a50026", "#fee090", "#74add1")
barchart(project, K = 8, run = best,
         border = NA, space = 0, sort.by.Q = TRUE,
         col = my.colors,
         xlab = "Individuals",
         ylab = "Ancestry proportions",
         main = "Ancestry matrix (K=8)") 

# exporting matrix
matrix8 <- Q(project, K=8, run=best)
write.csv(matrix8, "matrix8")

#######################################################################################
# plotting ancestral populations on map
par(mfrow=c(1,1))

# spatial packages
library(raster)
library(rgeos)
library(rgdal)
library(ggmap)
library(sp)
library(terra)

# colors
library(colorspace)

# Read in the shapefile (obtained from the Cartographic Service of Iceland at https://www.lmi.is/)
IS1<-readOGR(dsn="./Maps/", layer="IS1")
IS2<-readOGR(dsn="./Maps/", layer="IS2")
IS3<-readOGR(dsn="./Maps/", layer="IS3")

zones_clipped_1 <- raster::crop(IS1, extent(1569500,1650000,170000,230000))
zones_clipped_2 <- raster::crop(IS2, extent(1569500,1650000,170000,230000))
zones_clipped_3 <- raster::crop(IS3, extent(1569500,1650000,170000,230000))

qmatrix = Q(project, K=8, run=best)
popmap <- read.table("popmap.olfusa.clean.QGIS.tsv", header=TRUE)
k<- data.frame(qmatrix)
k$ID <- popmap$POP

k %>%
        group_by(k$ID) %>%
        summarize(mean_V1 = mean(V1, na.rm=TRUE),
                  mean_V2 = mean(V2, na.rm=TRUE),
                  mean_V3 = mean(V3, na.rm=TRUE),
                  mean_V4 = mean(V4, na.rm=TRUE),
                  mean_V5 = mean(V5, na.rm=TRUE),
                  mean_V6 = mean(V6, na.rm=TRUE),
                  mean_V7 = mean(V7, na.rm=TRUE),
                  mean_V8 = mean(V8, na.rm=TRUE))
                  
k_summary <- k %>%
        group_by(k$ID) %>%
        summarize(mean_V1 = mean(V1, na.rm=TRUE),
                  mean_V2 = mean(V2, na.rm=TRUE),
                  mean_V3 = mean(V3, na.rm=TRUE),
                  mean_V4 = mean(V4, na.rm=TRUE),
                  mean_V5 = mean(V5, na.rm=TRUE),
                  mean_V6 = mean(V6, na.rm=TRUE),
                  mean_V7 = mean(V7, na.rm=TRUE),
                  mean_V8 = mean(V8, na.rm=TRUE))
                  
coord <- data.frame(popmap$POP, popmap$LAT, popmap$LONG)

coord_summary <- coord %>%
        group_by(popmap.POP) %>%
        summarize(mean_LAT = mean(popmap.LAT, na.rm=TRUE),
                  mean_LONG = mean(popmap.LONG, na.rm=TRUE))

k_number <- k %>%
        group_by(k$ID) %>%
        tally()


k_summary$number <- k_number$n

data <- data.frame(k_summary)

data$LAT <- coord_summary$mean_LAT
data$LONG <- coord_summary$mean_LONG

y<- data$LAT
x<- data$LONG

data$lon=x
data$lat=y
coordinates(data) <- c("lon", "lat")
proj4string(data) <- CRS("+init=epsg:4326") # WGS 84
CRS.new <- CRS("+proj=lcc +lat_1=64.25 +lat_2=65.75 +lat_0=65 +lon_0=-19 +x_0=1700000 +y_0=300000 +ellps=GRS80 +units=m +no_defs")
d <- spTransform(data, CRS.new)
d1 <- data.frame(d)

par(mai=c(0.1,0.1,0.1,0.1))
plot(zones_clipped_1, yaxs = 'i', xaxs = 'i', lwd=2)
plot(zones_clipped_2, add=T, col="light gray",border="dark gray", yaxs = 'i', xaxs = 'i')
plot(zones_clipped_3, add=T, col="light gray", yaxs = 'i', xaxs = 'i')

plot(d, add=TRUE, pch=1, cex=0.01, alpha=0)

library(maps)
library(plotrix)
library(scales)
library(seqinr)

points(d1$lon, d1$lat, cex = d1$number/1800, col='white', pch=19)

my.colors <- c("#ffffbf", "#abd9e9", "#313695", "#f46d43", 
               "#d73027", "#a50026", "#fee090", "#74add1")

for (i in 1:(2)) {my.colors[i]<-col2alpha(color=my.colors[i], alpha=1)}

#for same size of pies
for (x in 1:nrow(d)){floating.pie(d1$lon[x],d1$lat[x], c(d1$mean_V1[x],d1$mean_V2[x],d1$mean_V3[x],d1$mean_V4[x],
                                                         d1$mean_V5[x], d1$mean_V6[x],d1$mean_V7[x],d1$mean_V8[x]),
                                  radius=1200, col =my.colors)}
