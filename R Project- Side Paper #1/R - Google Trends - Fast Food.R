### GOOGLE TRENDS DATA

#devtools::install_github("PMassicotte/gtrendsR")
#devtools::install_github("mamut86/GTT")

#Load packages 
library(readr) 
library(gtrendsR) 
library(purrr) 
library(GTT)
library(dplyr)

setwd("C:/Users/cjkno/Documents/My Documents/Classes - '23 Spring/Research/Paper (Side) #1/DATA/Google Trends")

################################ Pulling Codes ################################################
### Function to get DMA code

DMA_code <- function(geo){
  codes <- unique(countries$sub_code[substr(countries$sub_code, 1,5) == geo])
  if(length(codes) > 1){
    countries[countries$sub_code %in% codes[2:length(codes)], 2:3]
  } else{
    message('No DMA code for this geo')
  }
}

View(countries) # List of all codes

code <- DMA_code('US-CA') # Narrow list to one state

#Newark is 450 US-NJ-450
#Nielsen DMA for NY is 501


### GTT for getting topic 

kg <- kgraph(keyword = "McDonald's") # Enter keyword, then check output to find descriptions
kg$entities[[1]]$id


################################ Testing Formatting ###############################################

### GTRENDSR - Testing Keywords and Formatting 

trends <- gtrends(keyword = c('McDonalds', 'IT Job'), 
                  geo = c('US-CA'), 
                  time = "2015-01-01 2015-12-31",
                  gprop = "web",
                  low_search_volume = TRUE,
                  cookie_url = "http://trends.google.com/Cookies/NID")

DMA <- trends$interest_by_DMA

plot(trends[["interest_over_time"]][["date"]], 
     trends[["interest_over_time"]][["hits"]],
     type = "b", 
     col = "red")

results.time <- data.frame(trends$interest_over_time)
results.DMA <- data.frame(trends$interest_by_DMA)

# results <- results %>% dplyr::select(hits)
# colnames(results)[1] <- 'lead'

################################################################################

### GTRENDSR  

# Do simple search to create the empty data frame

trends <- gtrends(keyword = c('McDonalds', 'IT Job'), 
                  geo = c('US-CA'), 
                  time = "2010-01-01 2022-12-31",
                  gprop = "web",
                  low_search_volume = TRUE,
                  cookie_url = "http://trends.google.com/Cookies/NID")

#TIME: Clean 
results.time <- data.frame(trends$interest_over_time)
results.time <- results.time %>% dplyr::select(-geo, -time, -gprop, -category)

# TIME: Split query of interest and standardizing query into different columns
results.time.IT <- results.time %>%  subset(keyword == "IT Job")
results.time.IT <- data.frame(results.time.IT, row.names = TRUE)
results.time <- results.time %>%  subset(keyword != "IT Job")
results.time <- data.frame(results.time, row.names = TRUE)
results.time <- cbind(results.time, results.time.IT)

# TIME: Divide query of interest by standardizing query 
results.time$KeywordStandardized <- results.time[,1]/results.time[,3]
colnames(results.time)[5] <- paste0(results.time[1, 2])
results.time <- results.time %>% dplyr::select(-keyword, -hits)


# DMA: Clean
results.DMA <- data.frame(trends$interest_by_dma)
results.DMA <- results.DMA %>% dplyr::select(-geo, -gprop)

# DMA: Split query of interest and standardizing query into different columns
results.DMA.IT <- results.DMA %>%  subset(keyword == "IT Job")
results.DMA <- results.DMA %>%  subset(keyword != "IT Job")
results.DMA <- full_join(results.DMA, results.DMA.IT, by = "location")

# DMA: Divide query of interest by standardizing query 
results.DMA$KeywordStandardized <- results.DMA[,2]/results.DMA[,4]
colnames(results.DMA)[6] <-paste0(results.DMA[1, 3])
results.DMA <- results.DMA %>% dplyr::select(-keyword.x, -keyword.y, -hits.y, -hits.x)

### For Loop 

# Import keywords list 
kwlist <- read.csv("Fast Food Keywords.csv")

# Loop gtrends functions and save output to results data frame 
for (i in 1:nrow(kwlist)){
    trends <- gtrends(keyword = c(kwlist[i, 1], "IT Job"), 
                      geo = c('US-CA'), 
                      time = "2010-01-01 2022-12-31",
                      gprop = "web",
                      low_search_volume = TRUE,
                      cookie_url = "http://trends.google.com/Cookies/NID")
    
    #Clean
    trendsbytime <- data.frame(trends$interest_over_time)
    trendsbytime <- trendsbytime %>% dplyr::select(-geo, -time, -gprop, -category) 
    
    # TIME: Split query of interest and standardizing query into different columns
    trendsbytime.IT <- trendsbytime %>%  subset(keyword == "IT Job")
    trendsbytime.IT <- data.frame(trendsbytime.IT, row.names = TRUE)
    trendsbytime <- trendsbytime %>%  subset(keyword != "IT Job")
    trendsbytime <- data.frame(trendsbytime, row.names = TRUE)
    trendsbytime <- cbind(trendsbytime, trendsbytime.IT)
    
    # Divide query of interest by standardizing query 
    trendsbytime$KeywordStandardized <- trendsbytime[,1]/trendsbytime[,3]
    colnames(trendsbytime)[5] <- paste0(kwlist[i, 1])
    trendsbytime <- trendsbytime %>% dplyr::select(-keyword, -hits)
    
    # Clean
    trendsbyDMA <- data.frame(trends$interest_by_DMA)
    trendsbyDMA <- trendsbyDMA %>% dplyr::select(-geo, -gprop)
    
    # DMA: Split query of interest and standardizing query into different columns
    trendsbyDMA.IT <- trendsbyDMA %>%  subset(keyword == "IT Job")
    trendsbyDMA <- trendsbyDMA %>%  subset(keyword != "IT Job")
    trendsbyDMA <- full_join(trendsbyDMA, trendsbyDMA.IT, by = "location")
    
    # Divide query of interest by standardizing query 
    trendsbyDMA$KeywordStandardized <- trendsbyDMA[,2]/trendsbyDMA[,4]
    colnames(trendsbyDMA)[6] <- paste0(kwlist[i, 1])
    trendsbyDMA <- trendsbyDMA %>% dplyr::select(-keyword.x, -keyword.y, -hits.y, -hits.x)
    

    results.time.final <- cbind(results.time, trendsbytime)
    results.DMA.final <- full_join(results.DMA, trendsbyDMA, by = "location")
    

}

# # Loop gtrends functions and save output to results data frame 
# for (i in 1:nrow(kwlist)){
#   trends <- gtrends(keyword = c(kwlist[i, 1], "IT Job"), 
#                     geo = c('US-CA'), 
#                     time = "2010-01-01 2022-12-31",
#                     gprop = "web",
#                     low_search_volume = TRUE,
#                     cookie_url = "http://trends.google.com/Cookies/NID")
#   
#   trendsbytime <- data.frame(trends$interest_over_time)
#   trendsbytime <- trendsbytime %>% dplyr::select(-geo, -time, -gprop, -category)
#   trendsbytime <- trendsbytime %>%  subset(keyword != "IT Job")
#   colnames(trendsbytime)[2] <- paste0(kwlist[i, 1])
#   trendsbytime <- data.frame(trendsbytime, row.names = TRUE)
#   trendsbytime <- trendsbytime %>% dplyr::select(-keyword)
#   
#   
#   trendsbyDMA <- data.frame(trends$interest_by_DMA)
#   trendsbyDMA <- trendsbyDMA %>% dplyr::select(-geo, -gprop)
#   trendsbyDMA <- trendsbyDMA %>%  subset(keyword != "IT Job")
#   colnames(trendsbyDMA)[2] <- paste0(kwlist[i, 1])
#   trendsbyDMA <- trendsbyDMA %>% dplyr::select(-keyword)
#   
#   
#   results.time <- cbind(results.time, trendsbytime)
#   results.DMA <- full_join(results.DMA, trendsbyDMA, by = "location")
#   
#   #results.DMA <- cbind(results.DMA, trendsbyDMA)
#   
#   
# }


write.csv(results.time.final, file = "GTrends Mined Data - Time.csv", row.names = T)
write.csv(results.DMA.final, file = "GTrends Mined Data - DMA.csv", row.names = T)




