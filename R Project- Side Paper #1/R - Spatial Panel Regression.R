# Spatial Panel Regression 

library(spdep)
library(sp)
library(ggplot2)
library(spatstat)
library(classInt)
library(dplyr)
library(plm)
library(splm)


spatial_data <- st_read("C:/Users/cjkno/Documents/My Documents/Classes - '23 Spring/Research/Paper (Side) #1/Shapefile with All Data 08July2024", "DMAs_CirculatoryER_CensusData_GTrends")
