## Extract Difference Between EWS Appearance and Derivative Change ##

extract.ews.difference <- function(ews.data, deriv.data,data_source,sensitivity){
  
  require(tidyverse)
  
  out.dat <- ews.data %>%
    select(-cases)%>% # drop columns that will be duplicated below
    left_join(deriv.data, by = "true.date")%>% # merge EWS data and model derivative data
    #mutate(deriv.direction =as.numeric(levels(deriv.direction))[deriv.direction]) %>% # convert factor to numeric, keeping value levels
    group_by(metric.code,period)%>%
    summarise(first.deriv = true.date[rep(rle(deriv.direction)$lengths > 7 & rle(deriv.direction)$value == 1 & lag(rle(deriv.direction)$value) == 0,times = rle(deriv.direction)$lengths)][1],
              first.ews = true.date[rep(rle(threshold.crossed)$lengths >= sensitivity & rle(threshold.crossed)$value == 1,times = rle(threshold.crossed)$lengths)][1])%>%
    #^ extract the date index of the first period of positive first derivative (i.e. when 'deriv.direction' = 1) that persists for longer than a week ('lengths > 7')
    # and extract the date index of the first detection of an EWS (i.e. when 'threshold.crossed' = 1)  that persists for length 'sensitivity' 
    mutate(day.prior = base::difftime(first.deriv, first.ews, units = "days"), # calculate difference between date of derivative change and first EWS detection
           data.source= data_source)  %>%
    ungroup()%>%
    mutate(prediction = ifelse(is.na(day.prior) & is.na(first.deriv) & !is.na(first.ews),"unknown",
                               ifelse(is.na(day.prior) & !is.na(first.deriv) & is.na(first.ews), "miss",   
                                      ifelse((is.na(day.prior) & is.na(first.deriv) & is.na(first.ews) )|day.prior == 0, "match", 
                                             ifelse(day.prior > 0, "prior", "post"))))) %>% # categorise predictions
    mutate(prediction = factor(prediction, levels = c("match","miss", "post", "prior","unknown"))) # ensure all possible levels coded in to factor structure
  
  return(out.dat)
}
