###########################################################################################
### Figure 1 - Sequential UK COVID-19 EWS ### 
###########################################################################################

### Preamble ###
source("Code/composite_EWS_wrapper_fn.R")
source("Code/curve_direction_fn.R")
source("Code/expanding_GAM_fn.R")
require(tidyverse) #dplyr, ggplot etc
require(ggthemes) #theme_clean
require(pbmcapply) #paralleled lapply for EWS assessment
require(mgcv) #gam fitting
require(gratia) #derivative estimation
require(lubridate) #working with dates functions
require(data.table) #rbindlist
require(scales) #plot scales functions
require(egg) #ggarrange

### Load data ###
cov.uk.dat <- read.csv("Data/data_2021-Jun-09.csv") %>%
  mutate(true.date = as.Date(date,format ="%Y-%m-%d")) %>% #ensure data col in Date format
  arrange(true.date)%>% #sort by increasing date
  mutate(Date = lubridate::decimal_date(true.date)) %>% #convert to decimal for ease of EWS assessment
  mutate(cases = as.numeric(newCasesBySpecimenDate))%>% #ensure case col is numeric
  mutate(Country = "UK")%>% #set county category as "UK"
  mutate(Weekday = factor(weekdays(as.Date(true.date)),
                          levels = c("Monday", "Tuesday", "Wednesday","Thursday", "Friday", "Saturday", "Sunday"),
                          ordered = TRUE), 
         Weekday = as.integer(Weekday)) %>% #create weekday factor from true.date col
  select(Country,Date,cases,true.date,Weekday)%>% #select on pertinent cols
  slice(24:n()) #only keep timeseries after first case detected (confounds EWS assessment if timeseries begins with 0s) 
#https://coronavirus.data.gov.uk/details/cases

metrics <- c("skew","SD","acf") # metrics of interest (Skewness, variance and autocorrelation)

#### Expanding GAM ###
exp.uk <- expanding.gam(cov.uk.dat,sensitivity = 7,train.length = 10) #estimate wave onset using expanding_gam function.
# @sensitivity refers to number of decreasing time points required for a wave to 'subside'
# @train.length refers to number of time points used to fit first GAM (minimum required for convergence)

ukgam1 <-mgcv::gam(abs(cases)~ s(Weekday,bs = "cs",k=7)+s(Date, bs="tp", k=20), 
                   data = cov.uk.dat[1:exp.uk$end.index[1],],
                   family = gaussian(), method = "REML") #fit GAM for first 'wave'
plot.gam(ukgam1,pages=1) #check smooths
uk.deriv1 <- data.frame(gratia::derivatives(ukgam1, term = "s(Date)", interval = "confidence",n= dim(cov.uk.dat[1:exp.uk$end.index[1],])[1]),
                        "deriv.direction" = curve.direction(0, gratia::derivatives(ukgam1, term = "s(Date)", interval = "confidence",n= dim(cov.uk.dat[1:exp.uk$end.index[1],])[1])$lower, gratia::derivatives(ukgam1, term = "s(Date)", interval = "confidence",n= dim(cov.uk.dat[1:exp.uk$end.index[1],])[1])$upper)/2)
#estimate GAM derivatives

ukgam2 <-mgcv::gam(abs(cases)~ s(Weekday,bs = "cs",k=7) + s(Date, bs="tp", k=20) , 
                   data =  cov.uk.dat[exp.uk$start.index[2]:exp.uk$end.index[3],],
                   family = gaussian(), method = "REML")
uk.deriv2 <- data.frame(gratia::derivatives(ukgam2, term = "s(Date)", interval = "confidence",n= dim(cov.uk.dat[exp.uk$start.index[2]:exp.uk$end.index[3],])[1]),
                        "deriv.direction" = curve.direction(0, gratia::derivatives(ukgam2, term = "s(Date)", interval = "confidence",n= dim(cov.uk.dat[exp.uk$start.index[2]:exp.uk$end.index[3],])[1])$lower, gratia::derivatives(ukgam2, term = "s(Date)", interval = "confidence",n= dim(cov.uk.dat[exp.uk$start.index[2]:exp.uk$end.index[3],])[1])$upper)/2)

ukgam3 <-mgcv::gam(abs(cases)~ s(Weekday,bs = "cs",k=7)+s(Date, bs="tp", k=5) , 
                   data = cov.uk.dat[exp.uk$start.index[4]:dim(cov.uk.dat)[1],],
                   family = gaussian(), method = "REML")
uk.deriv3 <- data.frame(gratia::derivatives(ukgam3, term = "s(Date)", interval = "confidence",n= dim(cov.uk.dat[exp.uk$start.index[4]:dim(cov.uk.dat)[1],])[1]),
                        "deriv.direction" = curve.direction(0, gratia::derivatives(ukgam3, term = "s(Date)", interval = "confidence",n= dim(cov.uk.dat[exp.uk$start.index[4]:dim(cov.uk.dat)[1],])[1])$lower, gratia::derivatives(ukgam3, term = "s(Date)", interval = "confidence",n= dim(cov.uk.dat[exp.uk$start.index[4]:dim(cov.uk.dat)[1],])[1])$upper)/2)

### Perform EWS assessment ###
cutoff.uk.cov.ews.multigam <- pbmclapply(cov.uk.dat$Date[1:(length(cov.uk.dat$Date)-13)],FUN = function(x){
  
  data.cut <- cov.uk.dat[as.numeric(cov.uk.dat$Date) >=  as.numeric(x),]
  # subset dataset to start date x
  
  tmp <- comp_EWS_wrapper(data.frame(timedat = as.numeric(data.cut$Date),
                                     biomass = data.cut$cases),
                          metrics = metrics,  threshold = 2, burn_in =7,
                          plotIt = F, ggplotIt = F, tail.direction = "one.tailed",
                          interpolate = F, method = "w_comp")
  # @data.frame consists of timedat and case data (biomass byproduct of primary usage in other disciplines)
  # @metrics refers EWS indicators of interest
  # @threshold = threshold*sigma value for EWS indicator to trangress to constitute a 'signal'
  # @burn_in = number of time points to train expanding indicator mean and standard error
  # @plotIt and ggplotIt whether to plot EWS strength trends
  # @tail.direction = one.tailed/two.tailed (if one.tailed, only positive threshold surpassing considered a 'signal', two.tailed is both direction)
  # @interpolate whether missing values should be interpolated
  # @method = "w_comp"/"dakos". "w_comp" is expanding whereas "dakos" is rolling window
  
  tmp$cutoff <- paste(x) # paste cutoff date for reference
  return(tmp)
  
}, mc.cores =3 )
#perform EWS assessment repeated for each time point (the timeseries is initiated from each time point)
names(cutoff.uk.cov.ews.multigam) <- cov.uk.dat$true.date[1:(length(cov.uk.dat$true.date)-13)] #name list elements to start date of assessment
cutoff.uk.cov.ews.multigam <- data.table::rbindlist(cutoff.uk.cov.ews.multigam, idcol = "start.date")%>%
  mutate(period = ifelse(time <= exp.uk$end.date[1], "first", #define periods of constancy prior to waves
                         ifelse(time >= exp.uk$start.date[2] & time <= exp.uk$end.date[3], "second",
                                ifelse(time > exp.uk$end.date[3],"third","fourth")))) # classify EWS results dependent on which wave it falls in
save(cutoff.uk.cov.ews.multigam,file = "Results/UK/uk.multi.gam.ews.RData")
load("Results/UK/uk.multi.gam.ews.RData")

uk.multifirst.plot.dat <- cutoff.uk.cov.ews.multigam[cutoff.uk.cov.ews.multigam$start.date == "2020-02-22" & cutoff.uk.cov.ews.multigam$period == "first",]%>%
  left_join(data.frame("true.date" = cov.uk.dat$true.date,"time" = cov.uk.dat$Date), by = c("time"))
# subset EWS assessment data to first wave only and rbind 'true.date' for plotting efficency later

uk.multisecond.plot.dat <- cutoff.uk.cov.ews.multigam[cutoff.uk.cov.ews.multigam$start.date == cov.uk.dat$true.date[exp.uk$start.index[2]] & cutoff.uk.cov.ews.multigam$period == "second",]%>%
  left_join(data.frame("true.date" = cov.uk.dat$true.date,"time" = cov.uk.dat$Date), by = c("time"))
# subset EWS assessment data to second wave only

uk.multithird.plot.dat <- cutoff.uk.cov.ews.multigam[cutoff.uk.cov.ews.multigam$start.date == cov.uk.dat$true.date[exp.uk$start.index[4]] & cutoff.uk.cov.ews.multigam$period == "third",]%>%
  left_join(data.frame("true.date" = cov.uk.dat$true.date,"time" = cov.uk.dat$Date), by = c("time"))
# subset EWS assessment data to third wave only

total.multiukplot.dat <- rbind(uk.multifirst.plot.dat,uk.multisecond.plot.dat,uk.multithird.plot.dat) %>%
  rename(cases = count.used) # merge EWS subsets in to singular dataframe

multi.deriv.cov.uk <- cov.uk.dat[,2:4]%>%
  cbind(rbind(uk.deriv1,uk.deriv2,uk.deriv3)) # merge GAM subsets in to singular dataframe

multi.uk.p2 <- ggplot(data = total.multiukplot.dat, aes(x=true.date,y=str)) +
  geom_hline(yintercept = 2, linetype="solid", color = "grey", size=1)+
  geom_line(aes(col= metric.code))+
  geom_point(data =total.multiukplot.dat[which(total.multiukplot.dat$threshold.crossed==1)], aes(x=true.date, y = str, col = metric.code,alpha= "Detected")) +
  scale_alpha_manual(values = c(1),
                     breaks = c("Detected"),
                     guide = guide_legend(override.aes = list(linetype = c(0),shape = c(16)))) +
  scale_colour_manual(values = scales::hue_pal()(7),guide = guide_legend(override.aes = list(linetype = rep(1,7),shape=NA))) +
  ggthemes::theme_clean() + xlab("Date") + ylab("Strength of EWS") +
  geom_vline(xintercept = total.multiukplot.dat$true.date[total.multiukplot.dat$period != lag(total.multiukplot.dat$period)], linetype="dashed", 
             color = "black", size=1)+
  scale_x_date(date_breaks = "3 months", date_labels = "%b-%Y")+
  labs(color='EWS Indicator\nStrength',alpha ="EWS Detection") +
  theme(plot.margin = margin(c(10, 8, 0, 10)),
        legend.key.height = unit(0.3,"cm"),
        legend.key.width = unit(0.5,"cm"),
        legend.title = element_text(size = 10),
        legend.text = element_text(size=10))

multi.uk.p1 <- ggplot(multi.deriv.cov.uk, aes(x=true.date, y=cases)) +
  aes(group=NA)+
  geom_path(aes(col = as.factor(deriv.direction)))+
  geom_point(data= total.multiukplot.dat[which(total.multiukplot.dat$threshold.crossed==1 & total.multiukplot.dat$metric.code == "acf + SD + skew")], 
             aes(y =-5000, alpha = ""),size = 3,pch= "|",col = "#53B400")+
  ylab("COVID-19 Cases") + 
  xlab("Date")+
  geom_vline(xintercept = total.multiukplot.dat$true.date[total.multiukplot.dat$period != lag(total.multiukplot.dat$period)], 
             linetype="dashed",  color = "black", size=1)+
  scale_colour_manual(values=c("#22B4F5","black","#F07589"),name = "Predicted Nonlinear\nCase Trend", labels = c("Negative","No trend","Positive"))+
  scale_y_continuous(labels = scales::number_format(accuracy = 1000)) + 
  scale_x_date(date_breaks = "3 months", date_labels = "%b-%Y") +
  ggthemes::theme_clean() + ggtitle("UK Daily COVID Cases: Three Waves")+  
  annotate("label", x = as.Date("2020-06-01"), y =max(total.multiukplot.dat$cases)*0.8 , label = "EWS indicator: acf + SD + skewness")+
  scale_alpha_manual(values = 0.8, labels = "EWS detected",
                     guide = guide_legend(override.aes = list(size = 5))) +
  labs(alpha = "EWS")  +
  theme(plot.margin = margin(c(10, 8, 0, 10)),
        legend.key.height = unit(0.3,"cm" ),
        legend.key.width = unit(0.5,"cm"),
        legend.title = element_text(size = 10),
        legend.text = element_text(size=10))

pdf(file="Results/UK/multigam_UK_fig.pdf",
    width=9, height = 6,onefile=F)
egg::ggarrange(multi.uk.p1,multi.uk.p2,nrow = 2,heights = c(2, 2), labels = c("a","b"))
dev.off() # plot GAM predicted trends and EWS trends
