# Chapter 10 scripts, By Danlin Yu
# # Set up your workspace to where you store your data:

setwd("D:/OneDrive - montclair.edu/MyDocs/Collaboration/SpatialAnalysisBookProposal")
setwd("C:/Users/cjkno/Documents/My Documents/Classes - '23 Spring/Research/Paper (Side) #1/R Project- Side Paper #1/Chapter10")

# setwd("/mnt/d/OneDrive - montclair.edu/MyDocs/Collaboration/SpatialAnalysisBookProposal")

# set the font for any output figures:
library(extrafont)
loadfonts(device = "win")

# set the font for any output figures:

windowsFonts(FontStyle=windowsFont("Times New Roman"))

# load necessary libraries:

library(spdep)
library(sp)
library(ggplot2)
library(spatstat)
library(classInt)
library(dplyr)
library(plm)
library(splm)

spatial_data <- st_read(".", "gblntrans")
# construct the neighbors list

# Arrange spatial_data by 'COUNTYNAME' alphabetically, this is critical for panel data analysis because the panel data transformation
# Will re-arrange the spatial units alphabetically:

spatial_data <- spatial_data %>% arrange(COUNTYNAME)

spatial_coord <- st_coordinates(st_centroid(spatial_data))


# Import necessary packages
library(reshape2)
library(sf)
library(tidyverse)

# Convert sf object to data.frame
spatial_data_df <- as.data.frame(spatial_data)

# Extract only the needed data for panel regression:
spatial_data_df <- spatial_data_df[, c("COUNTYNAME","GDPPC95","GDPPC96","GDPPC97","GDPPC98","GDPPC99","GDPPC00","GDPPC01",
									   "FININCPC95","FININCPC96","FININCPC97","FININCPC98","FININCPC99","FININCPC00","FININCPC01",
									   "FDIPC95","FDIPC96","FDIPC97","FDIPC98","FDIPC99","FDIPC00","FDIPC01",
									   "FIXINVPC95","FIXINVPC96","FIXINVPC97","FIXINVPC98","FIXINVPC99","FIXINVPC00","FIXINVPC01",
									   "URB95","URB96","URB97","URB98","URB99","URB00","URB01")]


# Reshape the data to long format
spatial_data_long <- spatial_data_df %>%
  pivot_longer(-COUNTYNAME,
               names_to = c(".value", "Year"),
               names_pattern = "([A-Za-z]+)([0-9]+)",
               values_drop_na = TRUE)


# Correct the year values
spatial_data_long$Year <- ifelse(as.numeric(spatial_data_long$Year) < 50,
                                  paste0("20", spatial_data_long$Year),
                                  paste0("19", spatial_data_long$Year))

# Creating a pdata.frame
spatial_data_p <- pdata.frame(spatial_data_long, index = c("COUNTYNAME", "Year"))

# Set the formula

gb.fm.p <- GDPPC ~ FININCPC + FDIPC + FIXINVPC + URB

# Fit Fixed Effects, Random Effects, and Two-Ways Models
fe_model <- plm(gb.fm.p, data = spatial_data_p, model = "within")
re_model <- plm(gb.fm.p, data = spatial_data_p, model = "random")
tw_model <- plm(gb.fm.p, data = spatial_data_p, effect = "twoways")

# Conduct Hausman Test between Fixed Effects and Random Effects Models
hausman_test <- phtest(fe_model, re_model)

# Print models summaries
summary(fe_model)
summary(re_model)
summary(tw_model)

# Print Hausman Test result
print(hausman_test)

# Testing for individual and/or time effects:

# Estimate a pooled OLS model

pool.mod <- plm(gb.fm.p, data = spatial_data_p, model = "pooling")

# Test for individual effects
individual_effects_test <- plmtest(pool.mod, effect = "individual")

# Test for time effects
time_effects_test <- plmtest(pool.mod, effect = "time")

# Test for Two-Ways:
twoways_effects_test <- plmtest(pool.mod, effect = "twoways")

# print test results
print(individual_effects_test)
print(time_effects_test)
print(twoways_effects_test)

# Start the spatial panel regression estimation:

# 1. Spatial individual fixed effect model (lag model):

# Load the splm library if it is not yet loaded:

library(splm)

# Create the weight matrix list:

gb.listw <- nb2listw(poly2nb(spatial_data))

# The lag switch must be TRUE and spatial.error must be set to be "none" so the spml() will return a spatial lag panel model.
# Otherwise it will return a spatial autocorrelation panel model.

gb.splm.individual.fixed <- spml (gb.fm.p, data = spatial_data_p, listw = gb.listw, model = "within", effect = "individual", lag = TRUE, spatial.error = "none")

# Summarize the result, pay attention to the spatial lag parameter:

summary(gb.splm.individual.fixed)

# 2. Spatial individual fixed effect model (error model):

# Load the splm library if it is not yet loaded:

library(splm)

# Create the weight matrix list:

gb.listw <- nb2listw(poly2nb(spatial_data))

# The lag switch must be TRUE and spatial.error must be set to be "none" so the spml() will return a spatial lag panel model.
# Otherwise it will return a spatial autocorrelation panel model.

gb.spem.individual.fixed <- spml (gb.fm.p, data = spatial_data_p, listw = gb.listw, model = "within", effect = "individual", lag = FALSE, spatial.error = "kkp")

# Summarize the result, pay attention to the spatial lag parameter:

summary(gb.spem.individual.fixed)


# 3. Spatial individual fixed effect model (autocorrelation model):

# Load the splm library if it is not yet loaded:

library(splm)

# Create the weight matrix list:

gb.listw <- nb2listw(poly2nb(spatial_data))

# The lag switch must be TRUE and spatial.error must be set to be "kpp" so the spml() will return a spatial autocorrelation panel model.

gb.sapm.individual.fixed <- spml (gb.fm.p, data = spatial_data_p, listw = gb.listw, model = "within", effect = "individual", lag = TRUE, spatial.error = "kkp")

# Summarize the result, pay attention to the spatial lag parameter:

summary(gb.sapm.individual.fixed)

# LM test for spatial lag or spatial error models:

# Fixed effect models

lag.test.fixed <- slmtest (gb.fm.p, data = spatial_data_p, listw = gb.listw, model = "within", test = "lml")
error.test.fixed <- slmtest (gb.fm.p, data = spatial_data_p, listw = gb.listw, model = "within", test = "lme")
rlag.test.fixed <- slmtest (gb.fm.p, data = spatial_data_p, listw = gb.listw, model = "within", test = "rlml")
rerror.test.fixed <- slmtest (gb.fm.p, data = spatial_data_p, listw = gb.listw, model = "within", test = "rlme")

# Show the results:

lag.test.fixed
error.test.fixed
rlag.test.fixed
rerror.test.fixed

# Random effect models:

lag.test.random <- slmtest (gb.fm.p, data = spatial_data_p, listw = gb.listw, model = "random", test = "lml")
error.test.random <- slmtest (gb.fm.p, data = spatial_data_p, listw = gb.listw, model = "random", test = "lme")
rlag.test.random <- slmtest (gb.fm.p, data = spatial_data_p, listw = gb.listw, model = "random", test = "rlml")
rerror.test.random <- slmtest (gb.fm.p, data = spatial_data_p, listw = gb.listw, model = "random", test = "rlme")

# Show the results:

lag.test.random
error.test.random
rlag.test.random
rerror.test.random

# Spatial panel Durbin model:

# Create the spatial lag terms for the independent variables:

FININCPC.lag <- slag(spatial_data_p$FININCPC, gb.listw)
FDIPC.lag <- slag(spatial_data_p$FDIPC, gb.listw)
FIXINVPC.lag <- slag(spatial_data_p$FIXINVPC, gb.listw)
URB.lag <- slag(spatial_data_p$URB, gb.listw)

spatial_data_p_dm <- data.frame (spatial_data_p, FININCPC.lag, FDIPC.lag, FIXINVPC.lag, URB.lag)

# Creating a pdata.frame
spatial_data_p_dm <- pdata.frame(spatial_data_p_dm, index = c("COUNTYNAME", "Year"))


gb.fm.spdm <- GDPPC ~ FININCPC + FDIPC + FIXINVPC + URB + FININCPC.lag + FDIPC.lag + FIXINVPC.lag + URB.lag

# Spatial Durbin model (non-spatial) fixed/random:

gb.spdm.fixed <- plm (gb.fm.spdm, data = spatial_data_p_dm, model = "within", effect = "individual")
gb.spdm.random <- plm (gb.fm.spdm, data = spatial_data_p_dm, model = "random", effect = "individual")

summary(gb.spdm.fixed)
summary(gb.spdm.random)

# Spatial lag/error Durbin model fixed/random:

gb.spdm.lag.fixed <- spml (gb.fm.spdm, data = spatial_data_p_dm, listw = gb.listw, model = "within", effect = "individual", lag = TRUE, spatial.error = "none")
gb.spdm.error.fixed <- update(gb.spdm.lag.fixed, lag = FALSE, spatial.error = "kkp")
gb.spdm.lag.random <- update (gb.spdm.lag.fixed, model = "random")
gb.spdm.error.random <- update (gb.spdm.error.fixed, model = "random")

# Summarize the results:

summary(gb.spdm.lag.fixed)
summary(gb.spdm.error.fixed)
summary(gb.spdm.lag.random)
summary(gb.spdm.error.random)


# GWPR analysis using the Greater Beijing Area data:

# Load the GWPR scripts:
source("GWPR_unbalanced.r")

# Obtain the annually changing nearest neighbors using the adaptive strategy for each year, the search criterion is AIC:
gb.bd.p<-gwpr.bdwt(gb.fm.p, data = spatial_data_p, coords = spatial_coord, method = "aic", gweight = gwr.adaptive)

# You can check the varying nearest neighbors by:
gb.bd.p

# Now calibrate the GWPR model, with individual fixed effects:
gb.gwpr<- gwpr.unbalanced(gb.fm.p, data = spatial_data_p, coords = spatial_coord, nearestnb = gb.bd.p, effect = "individual", model = "within", inst.method = "baltagi")

# Map one of the spatially varying coefficients - FININCPC
# Convert row names to a column

gb.gwpr$gwpr.b$COUNTYNAME <- rownames(gb.gwpr$gwpr.b)
gb.gwpr$localpv$COUNTYNAME <- rownames(gb.gwpr$localpv)

# Join the dataframes
gb.df <- spatial_data %>%
  dplyr::left_join(gb.gwpr$gwpr.b, by = "COUNTYNAME") %>% 
  dplyr::left_join(gb.gwpr$localpv, by = "COUNTYNAME")

# Filter based on the p-value

gb.df$FININCPC.x[gb.df$FININCPC.y >=0.05] <-NA

# Now, make the map
p10.1<-ggplot() +
  geom_sf(data = gb.df, aes(fill = FININCPC.x)) +
  theme_minimal() +
  scale_fill_gradient(low = "gray90", high = "black", name = "Financial Income\nPer capita", na.value = "white") +
  labs(title = "Varying coefficients of Financial income per capita on GDP per capita") + 
  theme(text = element_text(family = "Times New Roman"))
  

p10.1
ggsave(file ="figure10.1.jpg", p10.1, width = 7.25, height = 7, dpi = 600)



# Starting with a New R session:
# Do a cross sectional GWR analysis and compare:

setwd("/yourworkspace/")


# set the font for any output figures:
library(extrafont)
loadfonts(device = "win")

# set the font for any output figures:

windowsFonts(FontStyle=windowsFont("Times New Roman"))

# load necessary libraries:

library(spdep)
library(sp)
library(ggplot2)
library(spatstat)
library(spgwr)


spatial_data <- st_read(".", "gblntrans")
# construct the neighbors list

# Arrange spatial_data by 'COUNTYNAME' alphabetically, this is critical for panel data analysis because the panel data transformation
# Will re-arrange the spatial units alphabetically:

spatial_data <- spatial_data %>% arrange(COUNTYNAME)

spatial_coord <- st_coordinates(st_centroid(spatial_data))

gb.fm <- GDPPC01 ~ FININCPC01 + FDIPC01 + FIXINVPC01 + URB01

gb.bd <- gwr.sel(gb.fm, data = spatial_data, coords = spatial_coord, adapt = T, gweight = gwr.bisquare, method = "aic")

gb.gwr.crs <- gwr(gb.fm, data = spatial_data, coords = spatial_coord, gweight = gwr.bisquare, adapt = gb.bd)

gb.out <- st_as_sf(gb.gwr.crs$SDF)
gb.out <- st_set_crs(gb.out, st_crs(spatial_data)) # Make sure the coordinate systems are the same

gb.out.poly <- st_join(spatial_data, gb.out, join = st_contains) # Spatial Join

# Map the Financial Income per capita:

p10.2<-ggplot() +
  geom_sf(data = gb.out.poly, aes(fill = FININCPC01.y)) +
  scale_fill_gradient(low = "gray90", high = "black") +
  theme_minimal() +
  labs(title = "Varying coefficients of Financial Income per capita on \nGDP per capita 2001", fill = "Financial Income\nPer capita") +
  theme(text = element_text(family = "Times New Roman"))
  
p10.2

ggsave (file ="figure10.2.jpg", p10.2, width = 7.25, height = 7, dpi = 600)



# The next block of codes attempt to evaluate the varying coefficients through RESF_VC routine:

library(vegan)
library(fields)
source("./resf/resf_panel.R")
source("./resf/resf_vc_panel.R")
source("./resf/resf_p.R")
source("./resf/resf_vc_p.R")
source("./resf/meigen2.R")
source("./resf/meigen_f2.R")

# Extract the spatial unit ID and time ID, this is required for resf panel vc model.

s_id <- spatial_data_p$COUNTYNAME
t_id <- spatial_data_p$Year

# Extract the dependent and independent variables:

lmod <- lm(gb.fm.p, data = spatial_data_p)
dep <- lmod$model[, 1]
indp <- lmod$model[, -1]

# The ids for each spatial units used to construct the Eigenvector:

cid<-seq(1,dim(spatial_coord)[1])

# All coefficients are assumed to be spatially varying. If there is strong theoretical guidance that some of the independent variables are not spatially varying, then they need to be specified using "xconst =" argument

gb.meig  <-meigen2(coords=spatial_coord, s_id=cid)      ## Moran's eigenvectors

################ fixed effect spatially varying coefficient model (the result could be unstable because the individual effects might very well absorb the spatially varying effect.
gb.resf.svc.fe<- resf_vc_panel(y=dep, x=indp, s_id=s_id, t_id=t_id, meig=gb.meig , pmodel="within", effect="individual")


################ random effects spatially varying model (I personally recommend these models, which seem stable)
gb.resf.svc.rd<- resf_vc_panel(y=dep, x=indp,s_id=s_id, t_id=t_id, meig=gb.meig, pmodel="random",effect="individual")

# First check to see which variables are actually varying across space, and which are not:

gb.resf.svc.fe$vc
gb.resf.svc.rd$vc

# 1 represents the coefficients of the variable is spatially varying, and 0 is not.

# Now extract the spatially varying coefficients and corresponding p-values:
# Because of how the RESF SVC model's output structure, the output coefficients and p-values are identical set stacked over the temporal periods.
# That means only the first time periods's values need to be extracted:

gb.out.fe.b <- gb.resf.svc.fe$b_vc [1: length(cid),]
gb.out.fe.p <- gb.resf.svc.fe$p_vc [1: length(cid),]

gb.out.rd.b <- gb.resf.svc.rd$b_vc [1: length(cid),]
gb.out.rd.p <- gb.resf.svc.rd$p_vc [1: length(cid),]

# Map the varying coefficients of FININCPC:

# Rename the columns

names(gb.out.fe.p)[names(gb.out.fe.p) == "FININCPC"] <- "FININCPC.p"
names(gb.out.fe.b)[names(gb.out.fe.b) == "FININCPC"] <- "FININCPC.b"

names(gb.out.rd.p)[names(gb.out.rd.p) == "FININCPC"] <- "FININCPC.p"
names(gb.out.rd.b)[names(gb.out.rd.b) == "FININCPC"] <- "FININCPC.b"

# Combine data by columns
spatial_data.fe <- st_bind_cols(spatial_data, gb.out.fe.p["FININCPC.p"], gb.out.fe.b["FININCPC.b"])

spatial_data.rd <- st_bind_cols(spatial_data, gb.out.rd.p["FININCPC.p"], gb.out.rd.b["FININCPC.b"])


spatial_data.fe$FININCPC.b[spatial_data.fe$FININCPC.p >=0.05] <-NA

spatial_data.rd$FININCPC.b[spatial_data.rd$FININCPC.p >=0.05] <-NA

# Now, make the maps

p10.3<-ggplot() +
  geom_sf(data = spatial_data.fe, aes(fill = FININCPC.b)) +
  theme_minimal() +
  scale_fill_gradient(low = "gray90", high = "black", name = "Financial Income\nPer capita", na.value = "white") +
  labs(title = "Varying coefficients of Financial income per capita \non GDP per capita based on RESF fixed effect model") + 
  theme(text = element_text(family = "Times New Roman"))
  

p10.3
ggsave(file ="figure10.3.jpg", p10.3, width = 7.25, height = 7, dpi = 600)


p10.4<-ggplot() +
  geom_sf(data = spatial_data.rd, aes(fill = FININCPC.b)) +
  theme_minimal() +
  scale_fill_gradient(low = "gray90", high = "black", name = "Financial Income\nPer capita", na.value = "white") +
  labs(title = "Varying coefficients of Financial income per capita \non GDP per capita based on RESF random effect model") + 
  theme(text = element_text(family = "Times New Roman"))
  

p10.4
ggsave(file ="figure10.4.jpg", p10.4, width = 7.25, height = 7, dpi = 600)

