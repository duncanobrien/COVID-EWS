################################################################################
          ### Expanding/Iterative GAM Function ###
################################################################################

expanding.gam <- function(mod.dat,sensitivity,train.length){
  #@mod.dat = data frame consisting of number of COVID-19 cases, date (numeric) and day-of-week (integer)
  #@sensitivity = length of time a stationary period persists for after a decrease of the same length for the decrease to be acknowledged
  #@train.length = minimum number of data points allocated to model fitting
  
  library(mgcv) #GAM fitting
  library(gratia) #derivative estimation
  set.seed(123)
  
  #curve.direction fn for identifying direction of change in GAM derivatives from
  #Burthe et al. 2016 Do early warning indicators consistently predict 
  #nonlinear change in long-term ecological data? J. Appl. Ecol. 53, 
  #666–676. (doi:https://doi.org/10.1111/1365-2664.12519)
  curve.direction <- function(x=0, L.CI, U.CI){
    pos <- ifelse(x<U.CI, 1, -1) 
    negs<-ifelse(x>L.CI, -1, 1)
    return(pos+negs)}
  
  success <-FALSE #initiate success tag
  stat.date <- data.frame(end.date = NA, start.date = NA,end.index = NA,start.index = NA) #initiate out frame
  end.pt <- 1 #initiate window end point
  
  while(end.pt < (dim(mod.dat)[1]-sensitivity)){ #check model won't exceed length of data but continue loop if smaller
    
    start <- end.pt #re-initiate window start point based previous window
    end.pt <- start + train.length -1 #re-initiate window end point based previous window
    
    while(isFALSE(success) & end.pt < (dim(mod.dat)[1])){ #if a stationary period is not detected and model doesn't exceed length of data, fit GAM up to current end.pt
      sub.dat <- na.omit(mod.dat[start:end.pt,]) #subset dataframe
      
      if(sum(sub.dat$cases) == 0){ #if all cases 0 (i.e. at start of monitoring), classify as fail and add next data point
        success <- FALSE
        end.pt = end.pt +1
      }else{
        tmp.gam <-  suppressWarnings(mgcv::gam(abs(cases)~ s(Weekday,bs = "cs",k=7)+s(Date, bs="tp", k=dim(sub.dat)[1]/4) , data = sub.dat,
                                               family = gaussian(), method = "REML")) 
        #fit generic GAM to subsetted dataframe. Warnings suppressed as can break function loop
        tmp.deriv <- data.frame(gratia::derivatives(tmp.gam, term = "s(Date)", interval = "confidence",n= dim(sub.dat)[1]),
                                "deriv.direction" = curve.direction(0, gratia::derivatives(tmp.gam, term = "s(Date)", interval = "confidence",n= dim(sub.dat)[1])$lower, gratia::derivatives(tmp.gam, term = "s(Date)", interval = "confidence",n= dim(sub.dat)[1])$upper)/2)
        #estimate derivatives and direction of change
        r <-tail(tmp.deriv$data[rep(rle(tmp.deriv$deriv.direction)$lengths >= sensitivity & rle(tmp.deriv$deriv.direction)$value == 0 & lag(rle(tmp.deriv$deriv.direction)$value) == -1,times = rle(tmp.deriv$deriv.direction)$lengths >= sensitivity)])
        #calculate whether a period of 'sensitivity' length of 0s (no change) follows a period of -1
        if(length(r) !=0){
          if(!is.na(r)[1]){
            if(((rle(tmp.deriv$deriv.direction)$lengths >= sensitivity & rle(tmp.deriv$deriv.direction)$value == -1 ) | (rle(tmp.deriv$deriv.direction)$lengths < sensitivity & rle(tmp.deriv$deriv.direction)$value == 0))[1]){ 
              #if((rle(tmp.deriv$deriv.direction)$lengths >= sensitivity & rle(tmp.deriv$deriv.direction)$value == -1)[1]){ 
              stat.date$end.date[dim(stat.date)[1]] <- max(r,na.rm=T)
              stat.date$end.index[dim(stat.date)[1]] <- end.pt
              success <- TRUE #if the first 'sensitivity' elements are '-1', or '0' for a period shorter than 'sensitivity', append end.pt to previous, and toggle success tag to 'TRUE' to break 'while' loop
            }else{
              stat.date <- rbind(stat.date,cbind(end.date=max(r,na.rm=T),start.date=mod.dat$Date[start],end.index = end.pt,start.index = start))
              success <- TRUE #if rle is met, record dates of window and toggle success tag to 'TRUE' to break 'while' loop
            }}
        }else{
          success <- FALSE #if rle not met, ensure success tag remains 'FALSE'
        }
        end.pt = end.pt +1 #update window end point
      }
    }
    success <- FALSE #reset success tag to 'FALSE' when starting new window and assessment (second 'while' loop)
  }
  return(stat.date[-1,]) # drop init row and return dates of interest 
}