# Spatial Panel Regression 

library(spdep)
library(sp)
library(ggplot2)
library(spatstat)
library(classInt)
library(dplyr)
library(plm)
library(splm)
library(reshape2)
library(sf)
library(tidyverse)

spatial_data <- st_read("C:/Users/cjkno/Documents/My Documents/Classes - '23 Spring/Research/Paper (Side) #1/Shapefile with All Data 10July2024", "DMAs_CirculatoryER_CensusData_GTrends_Simp")

# Sort Names Alphabetically 
spatial_data <- spatial_data %>% arrange(NAME)

# Select a model

# Extract only the needed data for panel regression:

#McD
gb.fm.p <- CircWtd ~ BIPOC + Hsptls + McD # Individual, Time, ~pooled, Lag
spatial_data <- spatial_data[, c("NAME",
                                     "Hsptls16",	"Hsptls17",	"Hsptls18",	"Hsptls19",	"Hsptls20",	"Hsptls21",	"Hsptls22",
                                     "CircWtd16",	"CircWtd17",	"CircWtd18", "CircWtd19", "CircWtd20", "CircWtd21", "CircWtd22",
                                     "P_BIPOC16",	"P_BIPOC17",	"P_BIPOC18",	"P_BIPOC19", "P_BIPOC20", "P_BIPOC21", "P_BIPOC22",
                                     "McD16",	"McD17",	"McD18",	"McD19",	"McD20",	"McD21",	"McD22")]

# # TB
# gb.fm.p <- CircWtd ~ BIPOC + Hsptls + TB # Individual, Time, ~pooled
# spatial_data <- spatial_data[, c("NAME",
#                                      "Hsptls14",	"Hsptls16",	"Hsptls17",	"Hsptls18",	"Hsptls19",	"Hsptls20",	"Hsptls21",	"Hsptls22",
#                                      "CircWtd14",	"CircWtd16",	"CircWtd17",	"CircWtd18", "CircWtd19", "CircWtd20", "CircWtd21", "CircWtd22",
#                                      "P_BIPOC14",	"P_BIPOC16",	"P_BIPOC17",	"P_BIPOC18",	"P_BIPOC19", "P_BIPOC20", "P_BIPOC21", "P_BIPOC22",
#                                      "TB14",	"TB16",	"TB17",	"TB18",	"TB19",	"TB20",	"TB21",	"TB22")]
# 
# # BK
# gb.fm.p <- CircWtd ~ BIPOC + Hsptls + BK # Individual, Time, ~pooled
# spatial_data <- spatial_data[, c("NAME",
#                                      "Hsptls17",	"Hsptls18",	"Hsptls19",	"Hsptls20",	"Hsptls21",	"Hsptls22",
#                                      "CircWtd17",	"CircWtd18", "CircWtd19", "CircWtd20", "CircWtd21", "CircWtd22",
#                                      "P_BIPOC17",	"P_BIPOC18",	"P_BIPOC19", "P_BIPOC20", "P_BIPOC21", "P_BIPOC22",
#                                      "BK17",	"BK18",	"BK19",	"BK20",	"BK21",	"BK22")]
# 
# # JITB
# gb.fm.p <- CircWtd ~ BIPOC + Hsptls + JITB # Random, None
# spatial_data <- spatial_data[, c("NAME",
#                                      "Hsptls19",	"Hsptls20",	"Hsptls21",	"Hsptls22",
#                                      "CircWtd19", "CircWtd20", "CircWtd21", "CircWtd22",
#                                      "P_BIPOC19", "P_BIPOC20", "P_BIPOC21", "P_BIPOC22",
#                                      "JITB19",	"JITB20",	"JITB21",	"JITB22")]
# 
# # KFC
# gb.fm.p <- CircWtd ~ BIPOC + Hsptls + KFC # Individual, Time
# spatial_data <- spatial_data[, c("NAME",
#                                      "Hsptls18",	"Hsptls19",	"Hsptls20",	"Hsptls21",	"Hsptls22",
#                                      "CircWtd18", "CircWtd19", "CircWtd20", "CircWtd21", "CircWtd22",
#                                      "P_BIPOC18",	"P_BIPOC19", "P_BIPOC20", "P_BIPOC21", "P_BIPOC22",
#                                      "KFC18",	"KFC19",	"KFC20",	"KFC21",	"KFC22")]


# Convert sf object to data.frame
spatial_coord <- st_coordinates(st_centroid(spatial_data))

spatial_data_df <- as.data.frame(spatial_data)


# Reshape the data to long format
spatial_data_long <- spatial_data_df %>%
  pivot_longer(-NAME,
               names_to = c(".value", "Year"),
               names_pattern = "([A-Za-z]+)([0-9]+)",
               values_drop_na = TRUE)


# Correct the year values
spatial_data_long$Year <- paste0("20", spatial_data_long$Year)

# Drop year with no data for circulatory   
# spatial_data_long <- subset(spatial_data_long, Year != 2015)

# Creating a pdata.frame
spatial_data_p <- pdata.frame(spatial_data_long, index = c("NAME", "Year"))

# Set the formula

# Fit Fixed Effects, Random Effects, and Two-Ways Models
fe_model <- plm(gb.fm.p, data = spatial_data_p, model = "within") # fe_model
re_model <- plm(gb.fm.p, data = spatial_data_p, model = "random") # re_model
tw_model <- plm(gb.fm.p, data = spatial_data_p, effect = "twoways") # tw_model

# Conduct Hausman Test between Fixed Effects and Random Effects Models
hausman_test <- phtest(fe_model, re_model)
print(hausman_test) #Null hypothesis is that a random effect model should be used, alternative is a fixed effects model 

# Print models summaries
summary(fe_model)
summary(re_model)
summary(tw_model)


### Testing for individual and/or time effects:

# Select a model

# Estimate a pooled OLS model
pool.mod <- plm(gb.fm.p, data = spatial_data_p, model = "pooling")

# Test for individual effects
plmtest(pool.mod, effect = "individual") # individual_effects_test

# Test for time effects
plmtest(pool.mod, effect = "time") # time_effects_test

# Test for Two-Ways:
plmtest(pool.mod, effect = "twoways") # twoways_effects_test



### Testing spatial models 

# LM test for spatial lag or spatial error models:

# Fixed effect models:
slmtest (gb.fm.p, data = spatial_data_p, listw = gb.listw, model = "within", test = "lml") # Lag Fixed
slmtest (gb.fm.p, data = spatial_data_p, listw = gb.listw, model = "within", test = "lme") # Error Fixed
slmtest (gb.fm.p, data = spatial_data_p, listw = gb.listw, model = "within", test = "rlml") # R Lag Fixed 
slmtest (gb.fm.p, data = spatial_data_p, listw = gb.listw, model = "within", test = "rlme") # R Error Fixed

# Random effect models:
slmtest (gb.fm.p, data = spatial_data_p, listw = gb.listw, model = "random", test = "lml") # Lag Random
slmtest (gb.fm.p, data = spatial_data_p, listw = gb.listw, model = "random", test = "lme") # Error Random
slmtest (gb.fm.p, data = spatial_data_p, listw = gb.listw, model = "random", test = "rlml") # R Lag Random
slmtest (gb.fm.p, data = spatial_data_p, listw = gb.listw, model = "random", test = "rlme") # R Error Random


# 1. Spatial individual fixed effect model (lag model):
# The lag switch must be TRUE and spatial.error must be set to be "none" so the spml() will return a spatial lag panel model.
# Otherwise it will return a spatial autocorrelation panel model.

gb.listw <- nb2listw(poly2nb(spatial_data), style = "W")

gb.splm.individual.fixed <- spml(formula = gb.fm.p, 
                                 data = spatial_data_p, 
                                 listw = gb.listw, 
                                 model = "within", 
                                 effect = "individual", 
                                 lag = TRUE, 
                                 spatial.error = "none")

summary(gb.splm.individual.fixed)

#Results_McD <- gb.splm.individual.fixed


# 2. Spatial individual fixed effect model (error model):
# The lag switch must be TRUE and spatial.error must be set to be "none" so the spml() will return a spatial lag panel model.
# Otherwise it will return a spatial autocorrelation panel model.

gb.spem.individual.fixed <- spml (gb.fm.p, data = spatial_data_p, 
                                  listw = gb.listw, 
                                  model = "within", 
                                  effect = "individual", 
                                  lag = FALSE, 
                                  spatial.error = "kkp")

summary(gb.spem.individual.fixed)


# 3. Spatial individual fixed effect model (autocorrelation model):
# The lag switch must be TRUE and spatial.error must be set to be "kpp" so the spml() will return a spatial autocorrelation panel model.

gb.sapm.individual.fixed <- spml (gb.fm.p, data = spatial_data_p, 
                                  listw = gb.listw, 
                                  model = "within", 
                                  effect = "individual", 
                                  lag = TRUE, spatial.error = "kkp")

summary(gb.sapm.individual.fixed)






######################################################################################
# # Extract only the needed data for panel regression:
# spatial_data_df <- spatial_data_df[, c("NAME",	"TtlHsptls", "Circ10",	"Circ11",	"Circ12",	"Circ13",	"Circ14",	"Circ15",	"Circ16",	"Circ17",	"Circ18",	"Circ19",	"Circ20",	"Circ21",	"Circ22",
#                                        "CircWtd19",	"CircWtd10",	"CircWtd11",	"CircWtd12",	"CircWtd13",	"CircWtd14",	"CircWtd15",	"CircWtd16",	"CircWtd17",	"CircWtd18",
#                                        "TTL10",	"TTL11",	"TTL12",	"TTL13",	"TTL14",	"TTL15",	"TTL16",	"TTL17",	"TTL18",	"TTL19",
#                                        "P_WHTE10",	"P_WHTE11",	"P_WHTE12",	"P_WHTE13",	"P_WHTE14",	"P_WHTE15",	"P_WHTE16",	"P_WHTE17",	"P_WHTE18",	"P_WHTE19",
#                                        "P_BIPOC10",	"P_BIPOC11",	"P_BIPOC12",	"P_BIPOC13",	"P_BIPOC14",	"P_BIPOC15",	"P_BIPOC16",	"P_BIPOC17",	"P_BIPOC18",	"P_BIPOC19",
#                                        "AGE10",	"AGE11",	"AGE12",	"AGE13",	"AGE14",	"AGE15",	"AGE16",	"AGE17",	"AGE18",	"AGE19",
#                                        "P_AGE19",	"P_AGE10",	"P_AGE11",	"P_AGE12",	"P_AGE13",	"P_AGE14",	"P_AGE15",	"P_AGE16",	"P_AGE17",	"P_AGE18",
#                                        "McD10",	"McD11",	"McD12",	"McD13",	"McD14",	"McD15",	"McD16",	"McD17",	"McD18",	"McD19",	"McD20",	"McD21",	"McD22",	
#                                        "TB10",	"TB11",	"TB12",	"TB13",	"TB14",	"TB15",	"TB16",	"TB17",	"TB18",	"TB19",	"TB20",	"TB21",	"TB22",	
#                                        "BK10",	"BK11",	"BK12",	"BK13",	"BK14",	"BK15",	"BK16",	"BK17",	"BK18",	"BK19",	"BK20",	"BK21",	"BK22",	
#                                        "CJ10",	"CJ11",	"CJ12",	"CJ13",	"CJ14",	"CJ15",	"CJ16",	"CJ17",	"CJ18",	"CJ19",	"CJ20",	"CJ21",	"CJ22",	
#                                        "JITB10",	"JITB11",	"JITB12",	"JITB13",	"JITB14",	"JITB15",	"JITB16",	"JITB17",	"JITB18",	"JITB19",	"JITB20",	"JITB21",	"JITB22",	
#                                        "KFC10",	"KFC11",	"KFC12",	"KFC13",	"KFC14",	"KFC15",	"KFC16",	"KFC17",	"KFC18",	"KFC19",	"KFC20",	"KFC21",	"KFC22")]
# 
# Extract only the needed data for panel regression:
# spatial_data_df <- spatial_data_df[, c("NAME",
#                                        "Hsptls10",	"Hsptls11",	"Hsptls12",	"Hsptls13",	"Hsptls14",	"Hsptls16",	"Hsptls17",	"Hsptls18",	"Hsptls19",	"Hsptls20",	"Hsptls21",	"Hsptls22",
#                                        "CircWtd10",	"CircWtd11",	"CircWtd12",	"CircWtd13",	"CircWtd14",	"CircWtd16",	"CircWtd17",	"CircWtd18", "CircWtd19", "CircWtd20", "CircWtd21", "CircWtd22",
#                                        "P_BIPOC10",	"P_BIPOC11",	"P_BIPOC12",	"P_BIPOC13",	"P_BIPOC14",	"P_BIPOC16",	"P_BIPOC17",	"P_BIPOC18",	"P_BIPOC19", "P_BIPOC20", "P_BIPOC21", "P_BIPOC22",
#                                        "McD16",	"McD17",	"McD18",	"McD19",	"McD20",	"McD21",	"McD22",	
#                                        "TB14",	"TB16",	"TB17",	"TB18",	"TB19",	"TB20",	"TB21",	"TB22",	
#                                        "BK17",	"BK18",	"BK19",	"BK20",	"BK21",	"BK22",	
#                                        "JITB19",	"JITB20",	"JITB21",	"JITB22",	
#                                        "KFC18",	"KFC19",	"KFC20",	"KFC21",	"KFC22")]



# # Convert to numeric fields 
# spatial_data_df$P_BIPOC18 <- as.numeric(spatial_data_df$P_BIPOC18)
# 
# indexval <- which(colnames(spatial_data_df) == "McD15")
# 
# for (i in indexval:ncol(spatial_data_df)){
#   spatial_data_df[[i]] <- as.numeric(spatial_data_df[[i]]) # Remember, using single brackets[] returns a data frame, double brackets [[]] returns a vector 
# }