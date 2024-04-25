#### Code to Run on HPC
# Bring in packages
library(pacman)
pacman::p_load(Hmsc, tidyverse, corrplot, raster, spThin,
               terra, mapdata, sf, ggplot2, sp, parallel,
               lattice, gridExtra, gstat)


# Load model data
load("/home/nhaydt/NRT_Project/bact2_lit_out_19Ap24.Rdata")
load("/home/nhaydt/NRT_Project/fung2_lit_out_19Ap24.Rdata")

# Load predictor raster stack
preds <- rast("/home/nhaydt/NRT_Project/preds_Feb24.tif")
preds <- stack(preds)

## Bring in grid points
grid <- read.csv("/home/nhaydt/NRT_Project/data/US_grid_250thous_random.csv", stringsAsFactors=FALSE,
                 header = TRUE) ## Across US

# colnames(grid) <- c("x", "y", "precip", "max_temp", "min_temp",
#                     "et", "veg_class", "soil_pH",
#                     "soil_water", "human_mod")

xy.grid <- as.matrix(cbind(grid$x,grid$y))
XData.grid <- grid[ ,c(3:ncol(grid))]

# Convert the environmental and spatial predictors into a format that can be used as input for the predict function
Gradient_blit <- prepareGradient(full_blit_mod, 
                                 XDataNew = XData.grid, 
                                 sDataNew = list(Route=xy.grid))
Gradient_flit <- prepareGradient(full_flit_mod,
                                 XDataNew = XData.grid,
                                 sDataNew = list(Route=xy.grid))

# Compute the posterior predictive distribution
# BACTERIA MODEL
start_time <- Sys.time()
predY_blit <-  predict(full_blit_mod, Gradient=Gradient_blit,
                       expected = TRUE,
                       nParallel = 1,
                       predictEtaMean = TRUE)
end_time <- Sys.time()
end_time - start_time
# FUNGI MODEL
start_time <- Sys.time()
predY_flit <-  predict(full_flit_mod, Gradient=Gradient_flit,
                       expected = TRUE,
                       nParallel = 1,
                       predictEtaMean = TRUE)
end_time <- Sys.time()
end_time - start_time

# Get matrix of posterior mean occurrence probabilities
# gc()
EpredY <- Reduce("+",predY_blit)/length(predY_blit)
pred.df.blit <- as.data.frame(EpredY)
write.csv(pred.df.blit, file = "/home/nhaydt/NRT_Project/mean_bact_preds_250thous_US.csv")

EpredY <- Reduce("+",predY_flit )/length(predY_flit )
pred.df.flit <- as.data.frame(EpredY)
write.csv(pred.df.flit, file = "/home/nhaydt/NRT_Project/mean_fung_preds_250thous_US.csv")