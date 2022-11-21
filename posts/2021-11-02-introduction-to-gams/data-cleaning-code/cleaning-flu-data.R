library(dplyr)
library(lubridate)

unzip("data/raw-flu-data/california-flu-data.zip")
flu_data_raw <- read.csv("clinical-sentinel-laboratory-influenza-and-other-respiratory-virus-surveillance-data-by-region-and-influenza-season.csv")

flu_data <- flu_data_raw %>%
  filter(season %in% c( "2009-2010","2010-2011","2011-2012", "2012-2013"),
         region =="California",
         Respiratory_Virus =="Total_Influenza",
         !Specimens_Tested ==0)%>%
  mutate(
    weekending_date = lubridate::as_date(weekending, format = "%m/%d/%Y"),
    week = lubridate::week(weekending_date),
    year = lubridate::year(weekending_date))%>%
  group_by(season)%>%
  mutate(year_min = min(year),
         week_centered =ifelse(year==year_min, week-52,week))%>%
  transmute(season= factor(season), 
            year=year, 
            week =week, 
            
            week_centered = week_centered, 
            tests_pos = Number_Positive,
            tests_total = Specimens_Tested)

write.csv(flu_data,"data/california-flu-data.csv",row.names = FALSE)  
         
         
