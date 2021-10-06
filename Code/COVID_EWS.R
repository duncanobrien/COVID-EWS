require(tidyverse) #dplyr, ggplot etc
require(ggthemes) #theme_clean
require(pbmcapply) #paralleled lapply for EWS assessment
require(mgcv) #gam fitting
require(gratia) #derivative estimation

source("Code/composite_EWS_wrapper_fn.R")
source("Code/curve_direction_fn.R")
source("Code/expanding_GAM_fn.R")
source("Code/ews_difference_fn.R")

###########################################################################################
### Read in Data ###
###########################################################################################
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

cov.jpn.dat <- read.csv("Data/pcr_positive_daily_JPN.csv") %>%
  mutate(true.date = as.Date(Date,format ="%Y/%m/%d")) %>% #slightly different date formatting in raw data
  arrange(true.date)%>%
  mutate(Date = lubridate::decimal_date(true.date)) %>%
  mutate(Country = "Japan")%>%
  mutate(Weekday = factor(weekdays(as.Date(true.date)),
                          levels = c("Monday", "Tuesday", "Wednesday","Thursday", "Friday", "Saturday", "Sunday"),
                          ordered = TRUE), 
         Weekday = as.integer(Weekday)) %>%
  select(Country,Date,cases,true.date,Weekday)%>%
  slice(28:n())

cov.usa.dat <- read.csv("Data/data_table_for_daily_case_trends_united_states.csv") %>%
  mutate(true.date = as.Date(Date,format ="%B %d %Y")) %>% #slightly different date formatting in raw data
  arrange(true.date)%>%
  mutate(Date = lubridate::decimal_date(true.date)) %>%
  mutate(New.Cases = gsub(",","",New.Cases))%>% #remove commas from case col
  mutate(cases = as.numeric(as.character(New.Cases)))%>%
  mutate(Country = "USA")%>%
  mutate(Weekday = factor(weekdays(as.Date(true.date)),
                          levels = c("Monday", "Tuesday", "Wednesday","Thursday", "Friday", "Saturday", "Sunday"),
                          ordered = TRUE), 
         Weekday = as.integer(Weekday)) %>%
  select(Country,Date,cases,true.date,Weekday)

cov.chl.dat <- as.data.frame(t(read.csv("Data/daily_cov_CHL.csv",header=FALSE))) %>%
  janitor::row_to_names(row_number = 1)%>%
  mutate(true.date = as.Date(Fecha,format ="%Y-%m-%d")) %>% #slightly different date formatting in raw data
  arrange(true.date)%>%
  mutate(Date = lubridate::decimal_date(true.date)) %>%
  mutate(cases = as.numeric(as.character(`Casos nuevos totales`)))%>%
  mutate(Country = "Chile")%>%
  mutate(Weekday = factor(weekdays(as.Date(true.date)),
                          levels = c("Monday", "Tuesday", "Wednesday","Thursday", "Friday", "Saturday", "Sunday"),
                          ordered = TRUE), 
         Weekday = as.integer(Weekday)) %>%
  select(Country,Date,cases,true.date,Weekday)%>%
  slice(3:n())

cov.WHO.dat <- read.csv("Data/WHO-COVID-19-global-data.csv") %>%
  mutate(true.date = as.Date(Date_reported,format ="%Y-%m-%d")) %>%
  group_by(Country) %>%
  arrange(Country,true.date) %>%
  mutate(Date = lubridate::decimal_date(true.date))%>%
  mutate(cases = abs(as.numeric(as.character(New_cases))))%>% #convert to numeric and assume negative values are typos
  mutate(Weekday = factor(weekdays(as.Date(true.date)),
                          levels = c("Monday", "Tuesday", "Wednesday","Thursday", "Friday", "Saturday", "Sunday"),
                          ordered = TRUE), 
         Weekday = as.integer(Weekday)) %>%
  select(Country,Date,cases,true.date,Weekday)

metrics <- c("skew","SD","acf") # metrics of interest (Skewness, variance and autocorrelation)

###########################################################################################
### UK ### Exemplar country, the workflow of which is repeated for all other countries
###########################################################################################
#https://coronavirus.data.gov.uk/details/cases
metrics <- c("skew","SD","acf") # metrics of interest (Skewness, variance and autocorrelation)

#### Expanding GAM
exp.uk <- expanding.gam(cov.uk.dat,sensitivity = 7,train.length = 10) #estimate wave onset using expanding_gam function.
          # @sensitivity refers to number of decreasing time points required for a wave to 'subside'
          # @train.length refers to number of time points used to fit first GAM (minimum required for convergence)

plot(cov.uk.dat$cases)
abline(v=exp.uk$end.index) #visualise 'wave' groups

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
plot.gam(ukgam2,pages=1)
uk.deriv2 <- data.frame(gratia::derivatives(ukgam2, term = "s(Date)", interval = "confidence",n= dim(cov.uk.dat[exp.uk$start.index[2]:exp.uk$end.index[3],])[1]),
                        "deriv.direction" = curve.direction(0, gratia::derivatives(ukgam2, term = "s(Date)", interval = "confidence",n= dim(cov.uk.dat[exp.uk$start.index[2]:exp.uk$end.index[3],])[1])$lower, gratia::derivatives(ukgam2, term = "s(Date)", interval = "confidence",n= dim(cov.uk.dat[exp.uk$start.index[2]:exp.uk$end.index[3],])[1])$upper)/2)

ukgam3 <-mgcv::gam(abs(cases)~ s(Weekday,bs = "cs",k=7)+s(Date, bs="tp", k=5) , 
                   data = cov.uk.dat[exp.uk$start.index[4]:dim(cov.uk.dat)[1],],
                   family = gaussian(), method = "REML")
plot.gam(ukgam3,pages =1)
uk.deriv3 <- data.frame(gratia::derivatives(ukgam3, term = "s(Date)", interval = "confidence",n= dim(cov.uk.dat[exp.uk$start.index[4]:dim(cov.uk.dat)[1],])[1]),
                        "deriv.direction" = curve.direction(0, gratia::derivatives(ukgam3, term = "s(Date)", interval = "confidence",n= dim(cov.uk.dat[exp.uk$start.index[4]:dim(cov.uk.dat)[1],])[1])$lower, gratia::derivatives(ukgam3, term = "s(Date)", interval = "confidence",n= dim(cov.uk.dat[exp.uk$start.index[4]:dim(cov.uk.dat)[1],])[1])$upper)/2)

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

png(file="Results/UK/multigam_UK_fig.png",
    width=3000, height = 2500 ,res=300)
egg::ggarrange(multi.uk.p1,multi.uk.p2,nrow = 2,heights = c(2, 2), labels = c("a","b"))
dev.off() # plot GAM predicted trends and EWS trends

uk.multiews.diff2 <- extract.ews.difference(total.multiukplot.dat,multi.deriv.cov.uk,"UK",2)
# calculate difference between first observed EWS warning and derivative increase.
# A warning is defined as two consecutive EWSs

#### Continuous EWS (repeat EWS assessment as above but with no restarting for each wave)
no.restart.uk.cov.ews <- comp_EWS_wrapper(data.frame(timedat = as.numeric(cov.uk.dat$Date),
                                                     biomass = cov.uk.dat$cases),
                                          metrics = metrics,  threshold = 2, burn_in =7,
                                          plotIt = F, ggplotIt = F, tail.direction = "one.tailed",
                                          interpolate = F, method = "w_comp")
save(no.restart.uk.cov.ews,file = "Results/UK/continuous_uk.ews.RData")
load("Results/UK/continuous_uk.ews.RData")

no.restart.ukplot.dat <- no.restart.uk.cov.ews%>%
  left_join(data.frame("true.date" = cov.uk.dat$true.date,"time" = cov.uk.dat$Date), by = c("time"))%>%
  rename(cases = count.used)

uk.no.restart.p2 <- ggplot(data = no.restart.ukplot.dat, aes(x=true.date,y=str)) +
  geom_hline(yintercept = 2, linetype="solid", color = "grey", size=1)+
  geom_line(aes(col= metric.code))+
  geom_point(data =no.restart.ukplot.dat[which(no.restart.ukplot.dat$threshold.crossed==1)], aes(x=true.date, y = str, col = metric.code,alpha= "Detected")) +
  scale_alpha_manual(values = c(1),
                     breaks = c("Detected"),
                     guide = guide_legend(override.aes = list(linetype = c(0),shape = c(16)))) +
  scale_colour_manual(values = scales::hue_pal()(7),guide = guide_legend(override.aes = list(linetype = rep(1,7),shape=NA))) +
  ggthemes::theme_clean() + xlab("Date") + ylab("Strength of EWS") +
  scale_x_date(date_breaks = "3 months", date_labels = "%b-%Y")+
  labs(color='EWS Indicator\nStrength',alpha ="EWS Detection") +
  theme(plot.margin = margin(c(10, 8, 0, 10)),
        legend.key.height = unit(0.3,"cm"),
        legend.key.width = unit(0.5,"cm"),
        legend.title = element_text(size = 10),
        legend.text = element_text(size=10))

uk.no.restart.p1 <- ggplot(multi.deriv.cov.uk, aes(x=true.date, y=cases)) +
  aes(group=NA)+
  geom_path(aes(col = as.factor(deriv.direction)))+
  geom_point(data= no.restart.ukplot.dat[which(no.restart.ukplot.dat$threshold.crossed==1 & no.restart.ukplot.dat$metric.code == "acf + SD + skew")], 
             aes(y =-5000, alpha = ""),size = 3,pch= "|",col = "#53B400")+
  ylab("COVID-19 Cases") + 
  xlab("Date")+
  scale_colour_manual(values=c("#22B4F5","black","#F07589"),name = "Predicted Nonlinear\nCase Trend", labels = c("Negative","No trend","Positive"))+
  scale_y_continuous(labels = scales::number_format(accuracy = 1000)) + 
  scale_x_date(date_breaks = "3 months", date_labels = "%b-%Y") +
  ggthemes::theme_clean() + ggtitle("UK Daily COVID Cases: Continuous Assessment")+  
  annotate("label", x = as.Date("2020-05-01"), y =max(no.restart.ukplot.dat$cases)*0.8 , label = "EWS indicator: acf + SD + skewness")+
  scale_alpha_manual(values = 0.8, labels = "EWS detected") +
  labs(alpha = "EWS")  +
  theme(plot.margin = margin(c(10, 8, 0, 10)),
        legend.key.height = unit(0.3,"cm" ),
        legend.key.width = unit(0.1,"cm"),
        legend.title = element_text(size = 10),
        legend.text = element_text(size=10))

png(file="Results/UK/UK_continuous_fig.png",
    width=3000, height = 2500 ,res=300)
egg::ggarrange(uk.no.restart.p1,uk.no.restart.p2,nrow = 2,heights = c(2, 2), labels = c("a","b"))
dev.off()

### Log GAM (repeat EWS assessment as in expanding GAM methodolgy but using log tranformed case data)
trans.log.cutoff.uk.cov.ews <- pbmclapply(cov.uk.dat$Date[1:(length(cov.uk.dat$Date)-13)],FUN = function(x){
  
  data.cut <- cov.uk.dat[as.numeric(cov.uk.dat$Date) >=  as.numeric(x),]
  
  tmp <- comp_EWS_wrapper(data.frame(timedat = as.numeric(data.cut$Date),
                                     biomass = log1p(data.cut$cases)),
                          metrics = metrics,  threshold = 2, burn_in =7,
                          plotIt = F, ggplotIt = F, tail.direction = "one.tailed",
                          interpolate = F, method = "w_comp")
  tmp$cutoff <- paste(x)
  return(tmp)
  
}, mc.cores =3 )
names(trans.log.cutoff.uk.cov.ews) <- cov.uk.dat$true.date[1:(length(cov.uk.dat$true.date)-13)]
temp <- trans.log.cutoff.uk.cov.ews
trans.log.cutoff.uk.cov.ews <- data.table::rbindlist(trans.log.cutoff.uk.cov.ews, idcol = "start.date")%>%
  #rename(start.date = date) %>%
  mutate(period = ifelse(time < 2020.505, "first", #define periods of constancy prior to waves
                         ifelse(time >= 2020.505 & time < 2021.321, "second",
                                ifelse(time >= 2021.321,"third","fourth"))))
save(trans.log.cutoff.uk.cov.ews,file = "Results/UK/trans.uklog.ews.RData")
load("Results/UK/trans.uklog.ews.RData")

uk.first.logplot.dat <- trans.log.cutoff.uk.cov.ews[trans.log.cutoff.uk.cov.ews$start.date == "2020-02-22" & trans.log.cutoff.uk.cov.ews$period == "first",]%>%
  left_join(data.frame("true.date" = cov.uk.dat$true.date,"time" = cov.uk.dat$Date), by = c("time"))

uk.second.logplot.dat <- trans.log.cutoff.uk.cov.ews[trans.log.cutoff.uk.cov.ews$start.date == "2020-07-04" & trans.log.cutoff.uk.cov.ews$period == "second",]%>%
  left_join(data.frame("true.date" = cov.uk.dat$true.date,"time" = cov.uk.dat$Date), by = c("time"))

uk.third.logplot.dat <- trans.log.cutoff.uk.cov.ews[trans.log.cutoff.uk.cov.ews$start.date == "2021-04-28" & trans.log.cutoff.uk.cov.ews$period == "third",]%>%
  left_join(data.frame("true.date" = cov.uk.dat$true.date,"time" = cov.uk.dat$Date), by = c("time"))

total.logukplot.dat <- rbind(uk.first.logplot.dat,uk.second.logplot.dat,uk.third.logplot.dat) %>%
  rename(cases = count.used)

uk.logtotal.p2<-ggplot(data = total.logukplot.dat, aes(x=true.date,y=str)) +
  geom_hline(yintercept = 2, linetype="solid", color = "grey", size=1)+
  geom_line(aes(col= metric.code,alpha= "Indicator strength"))+
  geom_point(data =total.logukplot.dat[which(total.logukplot.dat$threshold.crossed==1)], aes(x=true.date, y = str, col = metric.code,alpha= "EWS detected")) +
  scale_alpha_manual(values = c(1, 1),
                     breaks = c("Indicator strength", "EWS detected"),
                     guide = guide_legend(override.aes = list(linetype = c(1, 0),shape = c(NA, 16)))) +
  #scale_colour_manual(values = scales::hue_pal()(7)) +
  ggthemes::theme_clean() + xlab("Date") + ylab("Strength of EWS") +
  geom_vline(xintercept = total.logukplot.dat$true.date[total.logukplot.dat$period != lag(total.logukplot.dat$period)], linetype="dashed", 
             color = "black", size=1)+
  scale_x_date(date_breaks = "3 months", date_labels = "%b-%Y")+
  labs(color='Indicator Code',alpha ="EWS") +
  theme(plot.margin = margin(c(10, 8, 0, 10)),
        legend.key.height = unit(0.3,"cm"),
        legend.key.width = unit(0.1,"cm"),
        legend.title = element_text(size = 10),
        legend.text = element_text(size=10))

uk.logtotal.p1 <-ggplot(deriv.logcov.uk, aes(x=true.date, y=log(cases))) +
  aes(group=NA)+
  geom_path(aes(col = as.factor(deriv.direction)))+
  geom_point(data= total.logukplot.dat[which(total.logukplot.dat$threshold.crossed==1 & total.logukplot.dat$metric.code == "acf + SD + skew")], 
             aes(y =-1, alpha = ""),size = 3,pch= "|",col = "#53B400")+
  ylab("Log COVID-19 Cases") + 
  xlab("Date")+
  geom_vline(xintercept = total.logukplot.dat$true.date[total.logukplot.dat$period != lag(total.logukplot.dat$period)], 
             linetype="dashed",  color = "black", size=1)+
  scale_colour_manual(values=c("#22B4F5","black","#F07589"),name = "Predicted Nonlinear\nCase Trend", labels = c("Negative","No trend","Positive"))+
  scale_y_continuous(labels = scales::number_format(accuracy = 1)) + 
  scale_x_date(date_breaks = "3 months", date_labels = "%b-%Y") +
  ggthemes::theme_clean() + ggtitle("Log UK Daily COVID Cases: Three Waves")+  
  annotate("label", x = as.Date("2020-06-01"), y =max(log(total.logukplot.dat$cases))*0.8 , label = "EWS indicator: acf + SD + skewness")+
  scale_alpha_manual(values = 0.8, labels = "EWS detected",
                     guide = guide_legend(override.aes = list(size = 5))) +
  labs(alpha = "EWS")  +
  theme(plot.margin = margin(c(10, 8, 0, 10)),
        legend.key.height = unit(0.3,"cm" ),
        legend.key.width = unit(0.5,"cm"),
        legend.title = element_text(size = 10),
        legend.text = element_text(size=10))

png(file="Results/UK/UK_logGAM.png",
    width=3000, height = 2500 ,res=300)
egg::ggarrange(uk.logtotal.p1,uk.logtotal.p2,nrow = 2,heights = c(2, 2), labels = c("a","b"))
dev.off()

###########################################################################################
### World - Argentina ###
###########################################################################################
#https://covid19.who.int/info/
cov.arg.dat <- cov.WHO.dat %>%
  filter(Country == "Argentina")%>%
  slice(67:n())

ggplot(cov.arg.dat,aes(x=Date, y = cases)) + geom_path() + 
  scale_x_continuous(breaks=seq(2020,2022,0.2)) + 
  scale_y_continuous( labels = scales::number_format(accuracy = 100))+
  ggtitle("Daily Argentina COVID Cases")+ 
  ylab("Cases") + theme_bw()

#### Expanding GAM
exp.arg <- expanding.gam(cov.arg.dat,sensitivity = 7,train.length = 10)
plot(cov.arg.dat$cases)
abline(v=exp.arg$end.index)
arggam1 <-mgcv::gam(abs(cases)~ s(Weekday,bs = "cs",k=7)+s(Date, bs="tp", k=20), 
                    data = cov.arg.dat[exp.arg$start.index[1]:exp.arg$end.index[3],],
                    family = gaussian(), method = "REML")
plot.gam(arggam1,pages=1)
arg.deriv1 <- data.frame(gratia::derivatives(arggam1, term = "s(Date)", interval = "confidence",n= dim(cov.arg.dat[exp.arg$start.index[1]:exp.arg$end.index[3],])[1]),
                         "deriv.direction" = curve.direction(0, gratia::derivatives(arggam1, term = "s(Date)", interval = "confidence",n= dim(cov.arg.dat[exp.arg$start.index[1]:exp.arg$end.index[3],])[1])$lower, gratia::derivatives(arggam1, term = "s(Date)", interval = "confidence",n= dim(cov.arg.dat[exp.arg$start.index[1]:exp.arg$end.index[3],])[1])$upper)/2)
plot(arg.deriv1$deriv.direction)

arggam2 <-mgcv::gam(abs(cases)~ s(Weekday,bs = "cs",k=7) + s(Date, bs="tp", k=15) , 
                    data =  cov.arg.dat[exp.arg$start.index[4]:dim(cov.arg.dat)[1],],
                    family = gaussian(), method = "REML")
plot.gam(arggam2,pages=1)
arg.deriv2 <- data.frame(gratia::derivatives(arggam2, term = "s(Date)", interval = "confidence",n= dim(cov.arg.dat[exp.arg$start.index[4]:dim(cov.arg.dat)[1],])[1]),
                         "deriv.direction" = curve.direction(0, gratia::derivatives(arggam2, term = "s(Date)", interval = "confidence",n= dim(cov.arg.dat[exp.arg$start.index[4]:dim(cov.arg.dat)[1],])[1])$lower, gratia::derivatives(arggam2, term = "s(Date)", interval = "confidence",n= dim(cov.arg.dat[exp.arg$start.index[4]:dim(cov.arg.dat)[1],])[1])$upper)/2)
plot(arg.deriv2$deriv.direction)

cutoff.arg.cov.ews.multigam <- pbmclapply(cov.arg.dat$Date[1:(length(cov.arg.dat$Date)-13)],FUN = function(x){
  
  data.cut <- cov.arg.dat[as.numeric(cov.arg.dat$Date) >=  as.numeric(x),]
  
  tmp <- comp_EWS_wrapper(data.frame(timedat = as.numeric(data.cut$Date),
                                     biomass = data.cut$cases),
                          metrics = metrics,  threshold = 2, burn_in =7,
                          plotIt = F, ggplotIt = F, tail.direction = "one.tailed",
                          interpolate = F, method = "w_comp")
  tmp$cutoff <- paste(x)
  return(tmp)
  
}, mc.cores =3 )
names(cutoff.arg.cov.ews.multigam) <- cov.arg.dat$true.date[1:(length(cov.arg.dat$true.date)-13)]
temp <- cutoff.arg.cov.ews.multigam
cutoff.arg.cov.ews.multigam <- temp 
cutoff.arg.cov.ews.multigam <- data.table::rbindlist(cutoff.arg.cov.ews.multigam, idcol = "start.date")%>%
  mutate(period = ifelse(time <= exp.arg$end.date[3], "first", #define periods of constancy prior to waves
                         ifelse(time >= exp.arg$start.date[4] , "second","third")))
save(cutoff.arg.cov.ews.multigam,file = "Results/Argentina/arg.multi.gam.ews.RData")
load("Results/Argentina/arg.multi.gam.ews.RData")

arg.multifirst.plot.dat <- cutoff.arg.cov.ews.multigam[cutoff.arg.cov.ews.multigam$start.date == "2020-03-09" & cutoff.arg.cov.ews.multigam$period == "first",]%>%
  left_join(data.frame("true.date" = cov.arg.dat$true.date,"time" = cov.arg.dat$Date), by = c("time"))

arg.multisecond.plot.dat <- cutoff.arg.cov.ews.multigam[cutoff.arg.cov.ews.multigam$start.date == cov.arg.dat$true.date[exp.arg$start.index[4]] & cutoff.arg.cov.ews.multigam$period == "second",]%>%
  left_join(data.frame("true.date" = cov.arg.dat$true.date,"time" = cov.arg.dat$Date), by = c("time"))

total.multiargplot.dat <- rbind(arg.multifirst.plot.dat,arg.multisecond.plot.dat) %>%
  rename(cases = count.used)

multi.deriv.cov.arg <- cov.arg.dat[,2:4]%>%
  cbind(rbind(arg.deriv1,arg.deriv2))


multi.arg.p2 <- ggplot(data = total.multiargplot.dat, aes(x=true.date,y=str)) +
  geom_hline(yintercept = 2, linetype="solid", color = "grey", size=1)+
  geom_line(aes(col= metric.code,alpha= "Indicator strength"))+
  geom_point(data =total.multiargplot.dat[which(total.multiargplot.dat$threshold.crossed==1)], aes(x=true.date, y = str, col = metric.code,alpha= "EWS detected")) +
  scale_alpha_manual(values = c(1, 1),
                     breaks = c("Indicator strength", "EWS detected"),
                     guide = guide_legend(override.aes = list(linetype = c(1, 0),shape = c(NA, 16)))) +
  #scale_colour_manual(values = scales::hue_pal()(7)) +
  ggthemes::theme_clean() + xlab("Date") + ylab("Strength of EWS") +
  geom_vline(xintercept = total.multiargplot.dat$true.date[total.multiargplot.dat$period != lag(total.multiargplot.dat$period)], linetype="dashed", 
             color = "black", size=1)+
  scale_x_date(date_breaks = "3 months", date_labels = "%b-%Y")+
  labs(color='EWS Indicator',alpha ="EWS") +
  theme(plot.margin = margin(c(10, 8, 0, 10)),
        legend.key.height = unit(0.3,"cm"),
        legend.key.width = unit(0.1,"cm"),
        legend.title = element_text(size = 10),
        legend.text = element_text(size=10))


multi.arg.p1 <- ggplot(multi.deriv.cov.arg, aes(x=true.date, y=cases)) +
  aes(group=NA)+
  geom_path(aes(col = as.factor(deriv.direction)))+
  geom_point(data= total.multiargplot.dat[which(total.multiargplot.dat$threshold.crossed==1 & total.multiargplot.dat$metric.code == "acf + SD + skew")], 
             aes(y =-5000, alpha = ""),size = 3,pch= "|",col = "#53B400")+
  ylab("COVID-19 Cases") + 
  xlab("Date")+
  geom_vline(xintercept = total.multiargplot.dat$true.date[total.multiargplot.dat$period != lag(total.multiargplot.dat$period)], 
             linetype="dashed",  color = "black", size=1)+
  scale_colour_manual(values=c("#22B4F5","black","#F07589"),name = "Predicted Nonlinear\nCase Trend", labels = c("Negative","No trend","Positive"))+
  scale_y_continuous(labels = scales::number_format(accuracy = 1000)) + 
  scale_x_date(date_breaks = "3 months", date_labels = "%b-%Y") +
  ggthemes::theme_clean() + ggtitle("Argentina Daily COVID Cases: Three Waves")+  
  annotate("label", x = as.Date("2020-06-01"), y =max(total.multiargplot.dat$cases)*0.8 , label = "EWS indicator: acf + SD + skewness")+
  scale_alpha_manual(values = 0.8, labels = "EWS detected") +
  labs(alpha = "EWS")  +
  theme(plot.margin = margin(c(10, 8, 0, 10)),
        legend.key.height = unit(0.3,"cm" ),
        legend.key.width = unit(0.1,"cm"),
        legend.title = element_text(size = 10),
        legend.text = element_text(size=10))

png(file="Results/Argentina/multigam_arg_fig.png",
    width=3000, height = 2500 ,res=300)
egg::ggarrange(multi.arg.p1,multi.arg.p2,nrow = 2,heights = c(2, 2), labels = c("a","b"))
dev.off()

arg.multiews.diff2 <- extract.ews.difference(total.multiargplot.dat,multi.deriv.cov.arg,"Argentina",2)

###########################################################################################
### World - Belgium ###
###########################################################################################
#https://covid19.who.int/info/
cov.bel.dat <- cov.WHO.dat %>%
  filter(Country == "Belgium")%>%
  slice(34:n())


ggplot(cov.bel.dat,aes(x=Date, y = cases)) + geom_path() + 
  scale_x_continuous(breaks=seq(2020,2022,0.2)) + 
  scale_y_continuous( labels = scales::number_format(accuracy = 100))+
  ggtitle("Daily Belgium COVID Cases")+ 
  ylab("Cases") + theme_bw()

#### Expanding GAM
exp.bel <- expanding.gam(cov.bel.dat,sensitivity = 7,train.length = 10)
plot(cov.bel.dat$cases)
abline(v=exp.bel$end.index)

belgam1 <-mgcv::gam(abs(cases)~ s(Weekday,bs = "cs",k=7)+s(Date, bs="tp", k=30), 
                    data = cov.bel.dat[exp.bel$start.index[1]:exp.bel$end.index[1],],
                    family = gaussian(), method = "REML")
plot.gam(belgam1,pages=1)
bel.deriv1 <- data.frame(gratia::derivatives(belgam1, term = "s(Date)", interval = "confidence",n= dim(cov.bel.dat[exp.bel$start.index[1]:exp.bel$end.index[1],])[1]),
                         "deriv.direction" = curve.direction(0, gratia::derivatives(belgam1, term = "s(Date)", interval = "confidence",n= dim(cov.bel.dat[exp.bel$start.index[1]:exp.bel$end.index[1],])[1])$lower, gratia::derivatives(belgam1, term = "s(Date)", interval = "confidence",n= dim(cov.bel.dat[exp.bel$start.index[1]:exp.bel$end.index[1],])[1])$upper)/2)
plot(bel.deriv1$deriv.direction)

belgam2 <-mgcv::gam(abs(cases)~ s(Weekday,bs = "cs",k=7) + s(Date, bs="tp", k=15) , 
                    data =  cov.bel.dat[exp.bel$start.index[2]:exp.bel$end.index[5],],
                    family = gaussian(), method = "REML")
plot.gam(belgam2,pages=1)
bel.deriv2 <- data.frame(gratia::derivatives(belgam2, term = "s(Date)", interval = "confidence",n= dim(cov.bel.dat[exp.bel$start.index[2]:exp.bel$end.index[5],])[1]),
                         "deriv.direction" = curve.direction(0, gratia::derivatives(belgam2, term = "s(Date)", interval = "confidence",n= dim(cov.bel.dat[exp.bel$start.index[2]:exp.bel$end.index[5],])[1])$lower, gratia::derivatives(belgam2, term = "s(Date)", interval = "confidence",n= dim(cov.bel.dat[exp.bel$start.index[2]:exp.bel$end.index[5],])[1])$upper)/2)
plot(bel.deriv2$deriv.direction)

belgam3 <-mgcv::gam(abs(cases)~ s(Weekday,bs = "cs",k=7) + s(Date, bs="tp", k=20) , 
                    data =  cov.bel.dat[exp.bel$start.index[6]:dim(cov.bel.dat)[1],],
                    family = gaussian(), method = "REML")
plot.gam(belgam3,pages=1)
bel.deriv3 <- data.frame(gratia::derivatives(belgam3, term = "s(Date)", interval = "confidence",n= dim(cov.bel.dat[exp.bel$start.index[6]:dim(cov.bel.dat)[1],])[1]),
                         "deriv.direction" = curve.direction(0, gratia::derivatives(belgam3, term = "s(Date)", interval = "confidence",n= dim(cov.bel.dat[exp.bel$start.index[6]:dim(cov.bel.dat)[1],])[1])$lower, gratia::derivatives(belgam3, term = "s(Date)", interval = "confidence",n= dim(cov.bel.dat[exp.bel$start.index[6]:dim(cov.bel.dat)[1],])[1])$upper)/2)
plot(bel.deriv3$deriv.direction)

cutoff.bel.cov.ews.multigam <- pbmclapply(cov.bel.dat$Date[1:(length(cov.bel.dat$Date)-13)],FUN = function(x){
  
  data.cut <- cov.bel.dat[as.numeric(cov.bel.dat$Date) >=  as.numeric(x),]
  
  tmp <- comp_EWS_wrapper(data.frame(timedat = as.numeric(data.cut$Date),
                                     biomass = data.cut$cases),
                          metrics = metrics,  threshold = 2, burn_in =7,
                          plotIt = F, ggplotIt = F, tail.direction = "one.tailed",
                          interpolate = F, method = "w_comp")
  tmp$cutoff <- paste(x)
  return(tmp)
  
}, mc.cores =3 )
names(cutoff.bel.cov.ews.multigam) <- cov.bel.dat$true.date[1:(length(cov.bel.dat$true.date)-13)]
temp <- cutoff.bel.cov.ews.multigam
cutoff.bel.cov.ews.multigam <- data.table::rbindlist(cutoff.bel.cov.ews.multigam, idcol = "start.date")%>%
  mutate(period = ifelse(time <= exp.bel$end.date[1], "first", #define periods of constancy prior to waves
                         ifelse(time >= exp.bel$start.date[2] & time <= exp.bel$end.date[5], "second",
                                ifelse(time >= exp.bel$start.date[6], "third","fourth"))))
save(cutoff.bel.cov.ews.multigam,file = "Results/Belgium/bel.multi.gam.ews.RData")
load("Results/Belgium/bel.multi.gam.ews.RData")

bel.multifirst.plot.dat <- cutoff.bel.cov.ews.multigam[cutoff.bel.cov.ews.multigam$start.date == "2020-02-05" & cutoff.bel.cov.ews.multigam$period == "first",]%>%
  left_join(data.frame("true.date" = cov.bel.dat$true.date,"time" = cov.bel.dat$Date), by = c("time"))

bel.multisecond.plot.dat <- cutoff.bel.cov.ews.multigam[cutoff.bel.cov.ews.multigam$start.date == cov.bel.dat$true.date[exp.bel$start.index[2]] & cutoff.bel.cov.ews.multigam$period == "second",]%>%
  left_join(data.frame("true.date" = cov.bel.dat$true.date,"time" = cov.bel.dat$Date), by = c("time"))

bel.multithird.plot.dat <- cutoff.bel.cov.ews.multigam[cutoff.bel.cov.ews.multigam$start.date == cov.bel.dat$true.date[exp.bel$start.index[6]] & cutoff.bel.cov.ews.multigam$period == "third",]%>%
  left_join(data.frame("true.date" = cov.bel.dat$true.date,"time" = cov.bel.dat$Date), by = c("time"))

total.multibelplot.dat <- rbind(bel.multifirst.plot.dat,bel.multisecond.plot.dat,bel.multithird.plot.dat) %>%
  rename(cases = count.used)

multi.deriv.cov.bel <- cov.bel.dat[,2:4]%>%
  cbind(rbind(bel.deriv1,bel.deriv2,bel.deriv3))


multi.bel.p2 <- ggplot(data = total.multibelplot.dat, aes(x=true.date,y=str)) +
  geom_hline(yintercept = 2, linetype="solid", color = "grey", size=1)+
  geom_line(aes(col= metric.code,alpha= "Indicator strength"))+
  geom_point(data =total.multibelplot.dat[which(total.multibelplot.dat$threshold.crossed==1)], aes(x=true.date, y = str, col = metric.code,alpha= "EWS detected")) +
  scale_alpha_manual(values = c(1, 1),
                     breaks = c("Indicator strength", "EWS detected"),
                     guide = guide_legend(override.aes = list(linetype = c(1, 0),shape = c(NA, 16)))) +
  #scale_colour_manual(values = scales::hue_pal()(7)) +
  ggthemes::theme_clean() + xlab("Date") + ylab("Strength of EWS") +
  geom_vline(xintercept = total.multibelplot.dat$true.date[total.multibelplot.dat$period != lag(total.multibelplot.dat$period)], linetype="dashed", 
             color = "black", size=1)+
  scale_x_date(date_breaks = "3 months", date_labels = "%b-%Y")+
  labs(color='EWS Indicator',alpha ="EWS") +
  theme(plot.margin = margin(c(10, 8, 0, 10)),
        legend.key.height = unit(0.3,"cm"),
        legend.key.width = unit(0.1,"cm"),
        legend.title = element_text(size = 10),
        legend.text = element_text(size=10))


multi.bel.p1 <- ggplot(multi.deriv.cov.bel, aes(x=true.date, y=cases)) +
  aes(group=NA)+
  geom_path(aes(col = as.factor(deriv.direction)))+
  geom_point(data= total.multibelplot.dat[which(total.multibelplot.dat$threshold.crossed==1 & total.multibelplot.dat$metric.code == "acf + SD + skew")], 
             aes(y =-1000, alpha = ""),size = 3,pch= "|",col = "#53B400")+
  ylab("COVID-19 Cases") + 
  xlab("Date")+
  geom_vline(xintercept = total.multibelplot.dat$true.date[total.multibelplot.dat$period != lag(total.multibelplot.dat$period)], 
             linetype="dashed",  color = "black", size=1)+
  scale_colour_manual(values=c("#22B4F5","black","#F07589"),name = "Predicted Nonlinear\nCase Trend", labels = c("Negative","No trend","Positive"))+
  scale_y_continuous(labels = scales::number_format(accuracy = 1000)) + 
  scale_x_date(date_breaks = "3 months", date_labels = "%b-%Y") +
  ggthemes::theme_clean() + ggtitle("Belgium Daily COVID Cases: Three Waves")+  
  annotate("label", x = as.Date("2020-05-01"), y =max(total.multibelplot.dat$cases)*0.8 , label = "EWS indicator: acf + SD + skewness")+
  scale_alpha_manual(values = 0.8, labels = "EWS detected") +
  labs(alpha = "EWS")  +
  theme(plot.margin = margin(c(10, 8, 0, 10)),
        legend.key.height = unit(0.3,"cm" ),
        legend.key.width = unit(0.1,"cm"),
        legend.title = element_text(size = 10),
        legend.text = element_text(size=10))

png(file="Results/Belgium/multigam_bel_fig.png",
    width=3000, height = 2500 ,res=300)
egg::ggarrange(multi.bel.p1,multi.bel.p2,nrow = 2,heights = c(2, 2), labels = c("a","b"))
dev.off()

bel.multiews.diff2 <- extract.ews.difference(total.multibelplot.dat,multi.deriv.cov.bel,"Belgium",2)

###########################################################################################
### World - Brazil ###
###########################################################################################
#https://covid19.who.int/info/
cov.bzl.dat <- cov.WHO.dat %>%
  filter(Country == "Brazil")%>%
  slice(61:n())

ggplot(cov.bzl.dat,aes(x=Date, y = cases)) + geom_path() + 
  scale_x_continuous(breaks=seq(2020,2022,0.2)) + 
  scale_y_continuous( labels = scales::number_format(accuracy = 1000))+
  ggtitle("Daily Brazil COVID Cases")+ 
  ylab("Cases") + theme_bw()

#### Expanding GAM
exp.bzl <- expanding.gam(cov.bzl.dat,sensitivity = 7,train.length = 10)
plot(cov.bzl.dat$Date,cov.bzl.dat$cases)
abline(v=exp.bzl$end.date)
bzlgam1 <-mgcv::gam(abs(cases)~ s(Weekday,bs = "cs",k=7)+s(Date, bs="tp", k=15), 
                    data = cov.bzl.dat[exp.bzl$start.index[1]:exp.bzl$end.index[2],],
                    family = gaussian(), method = "REML")
plot.gam(bzlgam1,pages=1)
bzl.deriv1 <- data.frame(gratia::derivatives(bzlgam1, term = "s(Date)", interval = "confidence",n= dim(cov.bzl.dat[exp.bzl$start.index[1]:exp.bzl$end.index[2],])[1]),
                         "deriv.direction" = curve.direction(0, gratia::derivatives(bzlgam1, term = "s(Date)", interval = "confidence",n= dim(cov.bzl.dat[exp.bzl$start.index[1]:exp.bzl$end.index[2],])[1])$lower, gratia::derivatives(bzlgam1, term = "s(Date)", interval = "confidence",n= dim(cov.bzl.dat[exp.bzl$start.index[1]:exp.bzl$end.index[2],])[1])$upper)/2)
plot(bzl.deriv1$deriv.direction)
bzlgam2 <-mgcv::gam(abs(cases)~ s(Weekday,bs = "cs",k=7) + s(Date, bs="tp", k=20) , 
                    data =  cov.bzl.dat[exp.bzl$start.index[3]:dim(cov.bzl.dat)[1],],
                    family = gaussian(), method = "REML")
plot.gam(bzlgam2,pages=1)
bzl.deriv2 <- data.frame(gratia::derivatives(bzlgam2, term = "s(Date)", interval = "confidence",n= dim(cov.bzl.dat[exp.bzl$start.index[3]:dim(cov.bzl.dat)[1],])[1]),
                         "deriv.direction" = curve.direction(0, gratia::derivatives(bzlgam2, term = "s(Date)", interval = "confidence",n= dim(cov.bzl.dat[exp.bzl$start.index[3]:dim(cov.bzl.dat)[1],])[1])$lower, gratia::derivatives(bzlgam2, term = "s(Date)", interval = "confidence",n= dim(cov.bzl.dat[exp.bzl$start.index[3]:dim(cov.bzl.dat)[1],])[1])$upper)/2)
plot(bzl.deriv2$deriv.direction)

cutoff.bzl.cov.ews.multigam <- pbmclapply(cov.bzl.dat$Date[1:(length(cov.bzl.dat$Date)-13)],FUN = function(x){
  
  data.cut <- cov.bzl.dat[cov.bzl.dat$Date >=  as.numeric(x),]
  
  tmp <- comp_EWS_wrapper(data.frame(timedat = as.numeric(data.cut$Date),
                                     biomass = data.cut$cases),
                          metrics = metrics,  threshold = 2, burn_in =7,
                          plotIt = F, ggplotIt = F, tail.direction = "one.tailed",
                          interpolate = F, method = "w_comp")
  tmp$cutoff <- paste(x)
  return(tmp)
  
}, mc.cores =3 )
names(cutoff.bzl.cov.ews.multigam) <- cov.bzl.dat$true.date[1:(length(cov.bzl.dat$true.date)-13)]
temp <- cutoff.bzl.cov.ews.multigam
cutoff.bzl.cov.ews.multigam <- data.table::rbindlist(cutoff.bzl.cov.ews.multigam, idcol = "start.date")%>%
  mutate(period = ifelse(time <= exp.bzl$end.date[2], "first", #define periods of constancy prior to waves
                         ifelse(time >= exp.bzl$start.date[3], "second","third")))
save(cutoff.bzl.cov.ews.multigam,file = "Results/Brazil/bzl.multi.gam.ews.RData")
load("Results/Brazil/bzl.multi.gam.ews.RData")

bzl.multifirst.plot.dat <- cutoff.bzl.cov.ews.multigam[cutoff.bzl.cov.ews.multigam$start.date == "2020-03-03" & cutoff.bzl.cov.ews.multigam$period == "first",]%>%
  left_join(data.frame("true.date" = cov.bzl.dat$true.date,"time" = cov.bzl.dat$Date), by = c("time"))

bzl.multisecond.plot.dat <- cutoff.bzl.cov.ews.multigam[cutoff.bzl.cov.ews.multigam$start.date == cov.bzl.dat$true.date[exp.bzl$start.index[3]] & cutoff.bzl.cov.ews.multigam$period == "second",]%>%
  left_join(data.frame("true.date" = cov.bzl.dat$true.date,"time" = cov.bzl.dat$Date), by = c("time"))

total.multibzlplot.dat <- rbind(bzl.multifirst.plot.dat,bzl.multisecond.plot.dat) %>%
  rename(cases = count.used)

multi.deriv.cov.bzl <- cov.bzl.dat[,2:4]%>%
  cbind(rbind(bzl.deriv1,bzl.deriv2))

multi.bzl.p2 <- ggplot(data = total.multibzlplot.dat, aes(x=true.date,y=str)) +
  geom_hline(yintercept = 2, linetype="solid", color = "grey", size=1)+
  geom_line(aes(col= metric.code,alpha= "Indicator strength"))+
  geom_point(data =total.multibzlplot.dat[which(total.multibzlplot.dat$threshold.crossed==1)], aes(x=true.date, y = str, col = metric.code,alpha= "EWS detected")) +
  scale_alpha_manual(values = c(1, 1),
                     breaks = c("Indicator strength", "EWS detected"),
                     guide = guide_legend(override.aes = list(linetype = c(1, 0),shape = c(NA, 16)))) +
  #scale_colour_manual(values = scales::hue_pal()(7)) +
  ggthemes::theme_clean() + xlab("Date") + ylab("Strength of EWS") +
  geom_vline(xintercept = total.multibzlplot.dat$true.date[total.multibzlplot.dat$period != lag(total.multibzlplot.dat$period)], linetype="dashed", 
             color = "black", size=1)+
  scale_x_date(date_breaks = "3 months", date_labels = "%b-%Y")+
  labs(color='EWS Indicator',alpha ="EWS") +
  theme(plot.margin = margin(c(10, 8, 0, 10)),
        legend.key.height = unit(0.3,"cm"),
        legend.key.width = unit(0.1,"cm"),
        legend.title = element_text(size = 10),
        legend.text = element_text(size=10))


multi.bzl.p1 <- ggplot(multi.deriv.cov.bzl, aes(x=true.date, y=cases)) +
  aes(group=NA)+
  geom_path(aes(col = as.factor(deriv.direction)))+
  geom_point(data= total.multibzlplot.dat[which(total.multibzlplot.dat$threshold.crossed==1 & total.multibzlplot.dat$metric.code == "acf + SD + skew")], 
             aes(y =-5000, alpha = ""),size = 3,pch= "|",col = "#53B400")+
  ylab("COVID-19 Cases") + 
  xlab("Date")+
  geom_vline(xintercept = total.multibzlplot.dat$true.date[total.multibzlplot.dat$period != lag(total.multibzlplot.dat$period)], 
             linetype="dashed",  color = "black", size=1)+
  scale_colour_manual(values=c("#22B4F5","black","#F07589"),name = "Predicted Nonlinear\nCase Trend", labels = c("Negative","No trend","Positive"))+
  scale_y_continuous(labels = scales::number_format(accuracy = 1000)) + 
  scale_x_date(date_breaks = "3 months", date_labels = "%b-%Y") +
  ggthemes::theme_clean() + ggtitle("Brazil Daily COVID Cases: Two Waves")+  
  annotate("label", x = as.Date("2020-06-01"), y =max(total.multibzlplot.dat$cases)*0.8 , label = "EWS indicator: acf + SD + skewness")+
  scale_alpha_manual(values = 0.8, labels = "EWS detected") +
  labs(alpha = "EWS")  +
  theme(plot.margin = margin(c(10, 8, 0, 10)),
        legend.key.height = unit(0.3,"cm" ),
        legend.key.width = unit(0.1,"cm"),
        legend.title = element_text(size = 10),
        legend.text = element_text(size=10))

png(file="Results/Brazil/multigam_bzl_fig.png",
    width=3000, height = 2500 ,res=300)
egg::ggarrange(multi.bzl.p1,multi.bzl.p2,nrow = 2,heights = c(2, 2), labels = c("a","b"))
dev.off()

bzl.multiews.diff2 <- extract.ews.difference(total.multibzlplot.dat,multi.deriv.cov.bzl,"Brazil",2)

###########################################################################################
### World - Canada ###
###########################################################################################
#https://covid19.who.int/info/
cov.can.dat <- cov.WHO.dat %>%
  filter(Country == "Canada")%>%
  slice(56:n())


ggplot(cov.can.dat,aes(x=Date, y = cases)) + geom_path() + 
  scale_x_continuous(breaks=seq(2020,2022,0.2)) + 
  scale_y_continuous( labels = scales::number_format(accuracy = 100))+
  ggtitle("Daily Canada COVID Cases")+ 
  ylab("Cases") + theme_bw()

#### Expanding GAM
exp.can <- expanding.gam(cov.can.dat,sensitivity = 7,train.length = 10)
plot(cov.can.dat$cases)
abline(v=exp.can$end.index)

cangam1 <-mgcv::gam(abs(cases)~ s(Weekday,bs = "cs",k=7)+s(Date, bs="tp", k=40), 
                    data = cov.can.dat[exp.can$start.index[1]:exp.can$end.index[1],],
                    family = gaussian(), method = "REML")
plot.gam(cangam1,pages=1)
can.deriv1 <- data.frame(gratia::derivatives(cangam1, term = "s(Date)", interval = "confidence",n= dim(cov.can.dat[exp.can$start.index[1]:exp.can$end.index[1],])[1]),
                         "deriv.direction" = curve.direction(0, gratia::derivatives(cangam1, term = "s(Date)", interval = "confidence",n= dim(cov.can.dat[exp.can$start.index[1]:exp.can$end.index[1],])[1])$lower, gratia::derivatives(cangam1, term = "s(Date)", interval = "confidence",n= dim(cov.can.dat[exp.can$start.index[1]:exp.can$end.index[1],])[1])$upper)/2)
plot(can.deriv1$deriv.direction)

cangam2 <-mgcv::gam(abs(cases)~ s(Weekday,bs = "cs",k=7) + s(Date, bs="tp", k=20) , 
                    data =  cov.can.dat[exp.can$start.index[2]:exp.can$end.index[2],],
                    family = gaussian(), method = "REML")
plot.gam(cangam2,pages=1)
can.deriv2 <- data.frame(gratia::derivatives(cangam2, term = "s(Date)", interval = "confidence",n= dim(cov.can.dat[exp.can$start.index[2]:exp.can$end.index[2],])[1]),
                         "deriv.direction" = curve.direction(0, gratia::derivatives(cangam2, term = "s(Date)", interval = "confidence",n= dim(cov.can.dat[exp.can$start.index[2]:exp.can$end.index[2],])[1])$lower, gratia::derivatives(cangam2, term = "s(Date)", interval = "confidence",n= dim(cov.can.dat[exp.can$start.index[2]:exp.can$end.index[2],])[1])$upper)/2)
plot(can.deriv2$deriv.direction)

cangam3 <-mgcv::gam(abs(cases)~ s(Weekday,bs = "cs",k=7) + s(Date, bs="tp", k=20) , 
                    data =  cov.can.dat[exp.can$start.index[3]:dim(cov.can.dat)[1],],
                    family = gaussian(), method = "REML")
plot.gam(cangam3,pages=1)
can.deriv3 <- data.frame(gratia::derivatives(cangam3, term = "s(Date)", interval = "confidence",n= dim(cov.can.dat[exp.can$start.index[3]:dim(cov.can.dat)[1],])[1]),
                         "deriv.direction" = curve.direction(0, gratia::derivatives(cangam3, term = "s(Date)", interval = "confidence",n= dim(cov.can.dat[exp.can$start.index[3]:dim(cov.can.dat)[1],])[1])$lower, gratia::derivatives(cangam3, term = "s(Date)", interval = "confidence",n= dim(cov.can.dat[exp.can$start.index[3]:dim(cov.can.dat)[1],])[1])$upper)/2)
plot(can.deriv3$deriv.direction)

cutoff.can.cov.ews.multigam <- pbmclapply(cov.can.dat$Date[1:(length(cov.can.dat$Date)-13)],FUN = function(x){
  
  data.cut <- cov.can.dat[as.numeric(cov.can.dat$Date) >=  as.numeric(x),]
  
  tmp <- comp_EWS_wrapper(data.frame(timedat = as.numeric(data.cut$Date),
                                     biomass = data.cut$cases),
                          metrics = metrics,  threshold = 2, burn_in =7,
                          plotIt = F, ggplotIt = F, tail.direction = "one.tailed",
                          interpolate = F, method = "w_comp")
  tmp$cutoff <- paste(x)
  return(tmp)
  
}, mc.cores =3 )
names(cutoff.can.cov.ews.multigam) <- cov.can.dat$true.date[1:(length(cov.can.dat$true.date)-13)]
temp <- cutoff.can.cov.ews.multigam
cutoff.can.cov.ews.multigam <- temp
cutoff.can.cov.ews.multigam <- data.table::rbindlist(cutoff.can.cov.ews.multigam, idcol = "start.date")%>%
  mutate(period = ifelse(time <= exp.can$end.date[1], "first", #define periods of constancy prior to waves
                         ifelse(time >= exp.can$start.date[2] & time <= exp.can$end.date[2], "second",
                                ifelse(time >= exp.can$start.date[3], "third","fourth"))))
save(cutoff.can.cov.ews.multigam,file = "Results/Canada/can.multi.gam.ews.RData")
load("Results/Canada/can.multi.gam.ews.RData")

can.multifirst.plot.dat <- cutoff.can.cov.ews.multigam[cutoff.can.cov.ews.multigam$start.date == "2020-03-04" & cutoff.can.cov.ews.multigam$period == "first",]%>%
  left_join(data.frame("true.date" = cov.can.dat$true.date,"time" = cov.can.dat$Date), by = c("time"))

can.multisecond.plot.dat <- cutoff.can.cov.ews.multigam[cutoff.can.cov.ews.multigam$start.date == cov.can.dat$true.date[exp.can$start.index[2]] & cutoff.can.cov.ews.multigam$period == "second",]%>%
  left_join(data.frame("true.date" = cov.can.dat$true.date,"time" = cov.can.dat$Date), by = c("time"))

can.multithird.plot.dat <- cutoff.can.cov.ews.multigam[cutoff.can.cov.ews.multigam$start.date == cov.can.dat$true.date[exp.can$start.index[3]] & cutoff.can.cov.ews.multigam$period == "third",]%>%
  left_join(data.frame("true.date" = cov.can.dat$true.date,"time" = cov.can.dat$Date), by = c("time"))

total.multicanplot.dat <- rbind(can.multifirst.plot.dat,can.multisecond.plot.dat,can.multithird.plot.dat) %>%
  rename(cases = count.used)

multi.deriv.cov.can <- cov.can.dat[,2:4]%>%
  cbind(rbind(can.deriv1,can.deriv2,can.deriv3))


multi.can.p2 <- ggplot(data = total.multicanplot.dat, aes(x=true.date,y=str)) +
  geom_hline(yintercept = 2, linetype="solid", color = "grey", size=1)+
  geom_line(aes(col= metric.code,alpha= "Indicator strength"))+
  geom_point(data =total.multicanplot.dat[which(total.multicanplot.dat$threshold.crossed==1)], aes(x=true.date, y = str, col = metric.code,alpha= "EWS detected")) +
  scale_alpha_manual(values = c(1, 1),
                     breaks = c("Indicator strength", "EWS detected"),
                     guide = guide_legend(override.aes = list(linetype = c(1, 0),shape = c(NA, 16)))) +
  #scale_colour_manual(values = scales::hue_pal()(7)) +
  ggthemes::theme_clean() + xlab("Date") + ylab("Strength of EWS") +
  geom_vline(xintercept = total.multicanplot.dat$true.date[total.multicanplot.dat$period != lag(total.multicanplot.dat$period)], linetype="dashed", 
             color = "black", size=1)+
  scale_x_date(date_breaks = "3 months", date_labels = "%b-%Y")+
  labs(color='EWS Indicator',alpha ="EWS") +
  theme(plot.margin = margin(c(10, 8, 0, 10)),
        legend.key.height = unit(0.3,"cm"),
        legend.key.width = unit(0.1,"cm"),
        legend.title = element_text(size = 10),
        legend.text = element_text(size=10))


multi.can.p1 <- ggplot(multi.deriv.cov.can, aes(x=true.date, y=cases)) +
  aes(group=NA)+
  geom_path(aes(col = as.factor(deriv.direction)))+
  geom_point(data= total.multicanplot.dat[which(total.multicanplot.dat$threshold.crossed==1 & total.multicanplot.dat$metric.code == "acf + SD + skew")], 
             aes(y =-1000, alpha = ""),size = 3,pch= "|",col = "#53B400")+
  ylab("COVID-19 Cases") + 
  xlab("Date")+
  geom_vline(xintercept = total.multicanplot.dat$true.date[total.multicanplot.dat$period != lag(total.multicanplot.dat$period)], 
             linetype="dashed",  color = "black", size=1)+
  scale_colour_manual(values=c("#22B4F5","black","#F07589"),name = "Predicted Nonlinear\nCase Trend", labels = c("Negative","No trend","Positive"))+
  scale_y_continuous(labels = scales::number_format(accuracy = 1000)) + 
  scale_x_date(date_breaks = "3 months", date_labels = "%b-%Y") +
  ggthemes::theme_clean() + ggtitle("Canada Daily COVID Cases: Three Waves")+  
  annotate("label", x = as.Date("2020-06-01"), y =max(total.multicanplot.dat$cases)*0.8 , label = "EWS indicator: acf + SD + skewness")+
  scale_alpha_manual(values = 0.8, labels = "EWS detected") +
  labs(alpha = "EWS")  +
  theme(plot.margin = margin(c(10, 8, 0, 10)),
        legend.key.height = unit(0.3,"cm" ),
        legend.key.width = unit(0.1,"cm"),
        legend.title = element_text(size = 10),
        legend.text = element_text(size=10))

png(file="Results/Canada/multigam_can_fig.png",
    width=3000, height = 2500 ,res=300)
egg::ggarrange(multi.can.p1,multi.can.p2,nrow = 2,heights = c(2, 2), labels = c("a","b"))
dev.off()

can.multiews.diff2 <- extract.ews.difference(total.multicanplot.dat,multi.deriv.cov.can,"Canada",2)

###########################################################################################
### Chile ###
###########################################################################################
#https://www.minciencia.gob.cl/covid19/
ggplot(cov.chl.dat,aes(x=Date, y = cases)) + geom_path() + 
  scale_x_continuous(breaks=seq(2020,2022,0.2)) + 
  scale_y_continuous( labels = scales::number_format(accuracy = 1000))+
  ggtitle("Daily Chile COVID Cases")+ 
  ylab("Cases") + theme_bw()


#### Expanding GAM
exp.chl <- expanding.gam(cov.chl.dat,sensitivity = 7,train.length = 10)
plot(cov.chl.dat$cases)
abline(v=exp.chl$end.index)
chlgam1 <-mgcv::gam(abs(cases)~ s(Weekday,bs = "cs",k=7)+s(Date, bs="tp", k=10), 
                    data = cov.chl.dat[exp.chl$start.index[1]:exp.chl$end.index[1],],
                    family = gaussian(), method = "REML")
plot.gam(chlgam1,pages=1)
chl.deriv1 <- data.frame(gratia::derivatives(chlgam1, term = "s(Date)", interval = "confidence",n= dim(cov.chl.dat[exp.chl$start.index[1]:exp.chl$end.index[1],])[1]),
                         "deriv.direction" = curve.direction(0, gratia::derivatives(chlgam1, term = "s(Date)", interval = "confidence",n= dim(cov.chl.dat[exp.chl$start.index[1]:exp.chl$end.index[1],])[1])$lower, gratia::derivatives(chlgam1, term = "s(Date)", interval = "confidence",n= dim(cov.chl.dat[exp.chl$start.index[1]:exp.chl$end.index[1],])[1])$upper)/2)
plot(chl.deriv1$deriv.direction)

chlgam2 <-mgcv::gam(abs(cases)~ s(Weekday,bs = "cs",k=7) + s(Date, bs="tp", k=20) , 
                    data =  cov.chl.dat[exp.chl$start.index[2]:dim(cov.chl.dat)[1],],
                    family = gaussian(), method = "REML")
plot.gam(chlgam2,pages=1)
chl.deriv2 <- data.frame(gratia::derivatives(chlgam2, term = "s(Date)", interval = "confidence",n= dim(cov.chl.dat[exp.chl$start.index[2]:dim(cov.chl.dat)[1],])[1]),
                         "deriv.direction" = curve.direction(0, gratia::derivatives(chlgam2, term = "s(Date)", interval = "confidence",n= dim(cov.chl.dat[exp.chl$start.index[2]:dim(cov.chl.dat)[1],])[1])$lower, gratia::derivatives(chlgam2, term = "s(Date)", interval = "confidence",n= dim(cov.chl.dat[exp.chl$start.index[2]:dim(cov.chl.dat)[1],])[1])$upper)/2)
plot(chl.deriv2$deriv.direction)

cutoff.chl.cov.ews.multigam <- pbmclapply(cov.chl.dat$Date[1:(length(cov.chl.dat$Date)-13)],FUN = function(x){
  
  data.cut <- cov.chl.dat[as.numeric(cov.chl.dat$Date) >=  as.numeric(x),]
  
  tmp <- comp_EWS_wrapper(data.frame(timedat = as.numeric(data.cut$Date),
                                     biomass = data.cut$cases),
                          metrics = metrics,  threshold = 2, burn_in =7,
                          plotIt = F, ggplotIt = F, tail.direction = "one.tailed",
                          interpolate = F, method = "w_comp")
  tmp$cutoff <- paste(x)
  return(tmp)
  
}, mc.cores =3 )
names(cutoff.chl.cov.ews.multigam) <- cov.chl.dat$true.date[1:(length(cov.chl.dat$true.date)-13)]
temp <- cutoff.chl.cov.ews.multigam
cutoff.chl.cov.ews.multigam <- temp
cutoff.chl.cov.ews.multigam <- data.table::rbindlist(cutoff.chl.cov.ews.multigam, idcol = "start.date")%>%
  mutate(period = ifelse(time <= exp.chl$end.date[1], "first", #define periods of constancy prior to waves
                         ifelse(time > exp.chl$end.date[1] , "second","third")))
save(cutoff.chl.cov.ews.multigam,file = "Results/Chile/chl.multi.gam.ews.RData")
load("Results/Chile/chl.multi.gam.ews.RData")

chl.multifirst.plot.dat <- cutoff.chl.cov.ews.multigam[cutoff.chl.cov.ews.multigam$start.date == "2020-03-04" & cutoff.chl.cov.ews.multigam$period == "first",]%>%
  left_join(data.frame("true.date" = cov.chl.dat$true.date,"time" = cov.chl.dat$Date), by = c("time"))

chl.multisecond.plot.dat <- cutoff.chl.cov.ews.multigam[cutoff.chl.cov.ews.multigam$start.date == cov.chl.dat$true.date[exp.chl$start.index[2]] & cutoff.chl.cov.ews.multigam$period == "second",]%>%
  left_join(data.frame("true.date" = cov.chl.dat$true.date,"time" = cov.chl.dat$Date), by = c("time"))

total.multichlplot.dat <- rbind(chl.multifirst.plot.dat,chl.multisecond.plot.dat) %>%
  rename(cases = count.used)

multi.deriv.cov.chl <- cov.chl.dat[,2:4]%>%
  cbind(rbind(chl.deriv1,chl.deriv2))

multi.chl.p2 <- ggplot(data = total.multichlplot.dat, aes(x=true.date,y=str)) +
  geom_hline(yintercept = 2, linetype="solid", color = "grey", size=1)+
  geom_line(aes(col= metric.code,alpha= "Indicator strength"))+
  geom_point(data =total.multichlplot.dat[which(total.multichlplot.dat$threshold.crossed==1)], aes(x=true.date, y = str, col = metric.code,alpha= "EWS detected")) +
  scale_alpha_manual(values = c(1, 1),
                     breaks = c("Indicator strength", "EWS detected"),
                     guide = guide_legend(override.aes = list(linetype = c(1, 0),shape = c(NA, 16)))) +
  #scale_colour_manual(values = scales::hue_pal()(7)) +
  ggthemes::theme_clean() + xlab("Date") + ylab("Strength of EWS") +
  geom_vline(xintercept = total.multichlplot.dat$true.date[total.multichlplot.dat$period != lag(total.multichlplot.dat$period)], linetype="dashed", 
             color = "black", size=1)+
  scale_x_date(date_breaks = "3 months", date_labels = "%b-%Y")+
  labs(color='EWS Indicator',alpha ="EWS") +
  theme(plot.margin = margin(c(10, 8, 0, 10)),
        legend.key.height = unit(0.3,"cm"),
        legend.key.width = unit(0.1,"cm"),
        legend.title = element_text(size = 10),
        legend.text = element_text(size=10))


multi.chl.p1 <- ggplot(multi.deriv.cov.chl, aes(x=true.date, y=cases)) +
  aes(group=NA)+
  geom_path(aes(col = as.factor(deriv.direction)))+
  geom_point(data= total.multichlplot.dat[which(total.multichlplot.dat$threshold.crossed==1 & total.multichlplot.dat$metric.code == "acf + SD + skew")], 
                          aes(y =-1000, alpha = ""),size = 3,pch= "|",col = "#53B400")+
  ylab("COVID-19 Cases") + 
  xlab("Date")+
  geom_vline(xintercept = total.multichlplot.dat$true.date[total.multichlplot.dat$period != lag(total.multichlplot.dat$period)], 
             linetype="dashed",  color = "black", size=1)+
  scale_colour_manual(values=c("#22B4F5","black","#F07589"),name = "Predicted Nonlinear\nCase Trend", labels = c("Negative","No trend","Positive"))+
  scale_y_continuous(labels = scales::number_format(accuracy = 1000)) + 
  scale_x_date(date_breaks = "3 months", date_labels = "%b-%Y") +
  ggthemes::theme_clean() + ggtitle("Chile Daily COVID Cases: Two Waves")+  
  annotate("label", x = as.Date("2020-06-01"), y =max(total.multichlplot.dat$cases)*0.8 , label = "EWS indicator: acf + SD + skewness")+
  scale_alpha_manual(values = 0.8, labels = "EWS detected") +
  labs(alpha = "EWS")  +
  theme(plot.margin = margin(c(10, 8, 0, 10)),
        legend.key.height = unit(0.3,"cm" ),
        legend.key.width = unit(0.1,"cm"),
        legend.title = element_text(size = 10),
        legend.text = element_text(size=10))

png(file="Results/Chile/multigam_chl_fig.png",
    width=3000, height = 2500 ,res=300)
egg::ggarrange(multi.chl.p1,multi.chl.p2,nrow = 2,heights = c(2, 2), labels = c("a","b"))
dev.off()

chl.multiews.diff2 <- extract.ews.difference(total.multichlplot.dat,multi.deriv.cov.chl,"Chile",2)

###########################################################################################
### World - Colombia ###
###########################################################################################
#https://covid19.who.int/info/
cov.col.dat <- cov.WHO.dat %>%
  filter(Country == "Colombia")%>%
  slice(85:n())


ggplot(cov.col.dat,aes(x=Date, y = cases)) + geom_path() + 
  scale_x_continuous(breaks=seq(2020,2022,0.2)) + 
  scale_y_continuous( lacols = scales::number_format(accuracy = 100))+
  ggtitle("Daily Colombia COVID Cases")+ 
  ylab("Cases") + theme_bw()

#### Expanding GAM
exp.col <- expanding.gam(cov.col.dat,sensitivity = 7,train.length = 10)
plot(cov.col.dat$cases)
abline(v=exp.col$end.index)

colgam1 <-mgcv::gam(abs(cases)~ s(Weekday,bs = "cs",k=7)+s(Date, bs="tp", k=30), 
                    data = cov.col.dat[exp.col$start.index[1]:exp.col$end.index[4],],
                    family = gaussian(), method = "REML")
plot.gam(colgam1,pages=1)
col.deriv1 <- data.frame(gratia::derivatives(colgam1, term = "s(Date)", interval = "confidence",n= dim(cov.col.dat[exp.col$start.index[1]:exp.col$end.index[4],])[1]),
                         "deriv.direction" = curve.direction(0, gratia::derivatives(colgam1, term = "s(Date)", interval = "confidence",n= dim(cov.col.dat[exp.col$start.index[1]:exp.col$end.index[4],])[1])$lower, gratia::derivatives(colgam1, term = "s(Date)", interval = "confidence",n= dim(cov.col.dat[exp.col$start.index[1]:exp.col$end.index[4],])[1])$upper)/2)
plot(col.deriv1$deriv.direction)

colgam2 <-mgcv::gam(abs(cases)~ s(Weekday,bs = "cs",k=7) + s(Date, bs="tp", k=20) , 
                    data =  cov.col.dat[exp.col$start.index[5]:dim(cov.col.dat)[1],],
                    family = gaussian(), method = "REML")
plot.gam(colgam2,pages=1)
col.deriv2 <- data.frame(gratia::derivatives(colgam2, term = "s(Date)", interval = "confidence",n= dim(cov.col.dat[exp.col$start.index[5]:dim(cov.col.dat)[1],])[1]),
                         "deriv.direction" = curve.direction(0, gratia::derivatives(colgam2, term = "s(Date)", interval = "confidence",n= dim(cov.col.dat[exp.col$start.index[5]:dim(cov.col.dat)[1],])[1])$lower, gratia::derivatives(colgam2, term = "s(Date)", interval = "confidence",n= dim(cov.col.dat[exp.col$start.index[5]:dim(cov.col.dat)[1],])[1])$upper)/2)
plot(col.deriv2$deriv.direction)

cutoff.col.cov.ews.multigam <- pbmclapply(cov.col.dat$Date[1:(length(cov.col.dat$Date)-13)],FUN = function(x){
  
  data.cut <- cov.col.dat[as.numeric(cov.col.dat$Date) >=  as.numeric(x),]
  
  tmp <- comp_EWS_wrapper(data.frame(timedat = as.numeric(data.cut$Date),
                                     biomass = data.cut$cases),
                          metrics = metrics,  threshold = 2, burn_in =7,
                          plotIt = F, ggplotIt = F, tail.direction = "one.tailed",
                          interpolate = F, method = "w_comp")
  tmp$cutoff <- paste(x)
  return(tmp)
  
}, mc.cores =3 )
names(cutoff.col.cov.ews.multigam) <- cov.col.dat$true.date[1:(length(cov.col.dat$true.date)-13)]
temp <- cutoff.col.cov.ews.multigam
cutoff.col.cov.ews.multigam <- data.table::rbindlist(cutoff.col.cov.ews.multigam, idcol = "start.date")%>%
  mutate(period = ifelse(time <= exp.col$end.date[4], "first", #define periods of constancy prior to waves
                         ifelse(time >= exp.col$start.date[5] , "second", "third")))
save(cutoff.col.cov.ews.multigam,file = "/Results/Colombia/col.multi.gam.ews.RData")
load("Results/Colombia/col.multi.gam.ews.RData")

col.multifirst.plot.dat <- cutoff.col.cov.ews.multigam[cutoff.col.cov.ews.multigam$start.date == "2020-03-27" & cutoff.col.cov.ews.multigam$period == "first",]%>%
  left_join(data.frame("true.date" = cov.col.dat$true.date,"time" = cov.col.dat$Date), by = c("time"))

col.multisecond.plot.dat <- cutoff.col.cov.ews.multigam[cutoff.col.cov.ews.multigam$start.date == cov.col.dat$true.date[exp.col$start.index[5]] & cutoff.col.cov.ews.multigam$period == "second",]%>%
  left_join(data.frame("true.date" = cov.col.dat$true.date,"time" = cov.col.dat$Date), by = c("time"))

total.multicolplot.dat <- rbind(col.multifirst.plot.dat,col.multisecond.plot.dat) %>%
  rename(cases = count.used)

multi.deriv.cov.col <- cov.col.dat[,2:4]%>%
  cbind(rbind(col.deriv1,col.deriv2))


multi.col.p2 <- ggplot(data = total.multicolplot.dat, aes(x=true.date,y=str)) +
  geom_hline(yintercept = 2, linetype="solid", color = "grey", size=1)+
  geom_line(aes(col= metric.code,alpha= "Indicator strength"))+
  geom_point(data =total.multicolplot.dat[which(total.multicolplot.dat$threshold.crossed==1)], aes(x=true.date, y = str, col = metric.code,alpha= "EWS detected")) +
  scale_alpha_manual(values = c(1, 1),
                     breaks = c("Indicator strength", "EWS detected"),
                     guide = guide_legend(override.aes = list(linetype = c(1, 0),shape = c(NA, 16)))) +
  #scale_colour_manual(values = scales::hue_pal()(7)) +
  ggthemes::theme_clean() + xlab("Date") + ylab("Strength of EWS") +
  geom_vline(xintercept = total.multicolplot.dat$true.date[total.multicolplot.dat$period != lag(total.multicolplot.dat$period)], linetype="dashed", 
             color = "black", size=1)+
  scale_x_date(date_breaks = "3 months", date_labels = "%b-%Y")+
  labs(color='EWS Indicator',alpha ="EWS") +
  theme(plot.margin = margin(c(10, 8, 0, 10)),
        legend.key.height = unit(0.3,"cm"),
        legend.key.width = unit(0.1,"cm"),
        legend.title = element_text(size = 10),
        legend.text = element_text(size=10))


multi.col.p1 <- ggplot(multi.deriv.cov.col, aes(x=true.date, y=cases)) +
  aes(group=NA)+
  geom_path(aes(col = as.factor(deriv.direction)))+
  geom_point(data= total.multicolplot.dat[which(total.multicolplot.dat$threshold.crossed==1 & total.multicolplot.dat$metric.code == "acf + SD + skew")], 
             aes(y =-1000, alpha = ""),size = 3,pch= "|",col = "#53B400")+
  ylab("COVID-19 Cases") + 
  xlab("Date")+
  geom_vline(xintercept = total.multicolplot.dat$true.date[total.multicolplot.dat$period != lag(total.multicolplot.dat$period)], 
             linetype="dashed",  color = "black", size=1)+
  scale_colour_manual(values=c("#22B4F5","black","#F07589"),name = "Predicted Nonlinear\nCase Trend", labels = c("Negative","No trend","Positive"))+
  scale_y_continuous(labels = scales::number_format(accuracy = 1000)) + 
  scale_x_date(date_breaks = "3 months", date_labels = "%b-%Y") +
  ggthemes::theme_clean() + ggtitle("Colombia Daily COVID Cases: Two Waves")+  
  annotate("label", x = as.Date("2020-06-10"), y =max(total.multicolplot.dat$cases)*0.8 , label = "EWS indicator: acf + SD + skewness")+
  scale_alpha_manual(values = 0.8, labels = "EWS detected") +
  labs(alpha = "EWS")  +
  theme(plot.margin = margin(c(10, 8, 0, 10)),
        legend.key.height = unit(0.3,"cm" ),
        legend.key.width = unit(0.1,"cm"),
        legend.title = element_text(size = 10),
        legend.text = element_text(size=10))

png(file="Results/Colombia/multigam_col_fig.png",
    width=3000, height = 2500 ,res=300)
egg::ggarrange(multi.col.p1,multi.col.p2,nrow = 2,heights = c(2, 2), labels = c("a","b"))
dev.off()

col.multiews.diff2 <- extract.ews.difference(total.multicolplot.dat,multi.deriv.cov.col,"Colombia",2)

###########################################################################################
### World - France ###
###########################################################################################
#https://covid19.who.int/info/
cov.fr.dat <- cov.WHO.dat %>%
  filter(Country == "France") %>%
  slice(47:n())

ggplot(cov.fr.dat,aes(x=Date, y = cases)) + geom_path() + 
  scale_x_continuous(breaks=seq(2020,2022,0.2)) + 
  scale_y_continuous( labels = scales::number_format(accuracy = 100))+
  ggtitle("Daily France COVID Cases")+ 
  ylab("Cases") + theme_bw()


#### Expanding GAM
exp.fr <- expanding.gam(cov.fr.dat,sensitivity = 7,train.length = 10)
plot(cov.fr.dat$cases)
abline(v=exp.fr$end.index)
frgam1 <-mgcv::gam(abs(cases)~ s(Weekday,bs = "cs",k=7)+s(Date, bs="tp", k=30), 
                    data = cov.fr.dat[exp.fr$start.index[1]:exp.fr$end.index[1],],
                    family = gaussian(), method = "REML")
plot.gam(frgam1,pages=1)
fr.deriv1 <- data.frame(gratia::derivatives(frgam1, term = "s(Date)", interval = "confidence",n= dim(cov.fr.dat[exp.fr$start.index[1]:exp.fr$end.index[1],])[1]),
                         "deriv.direction" = curve.direction(0, gratia::derivatives(frgam1, term = "s(Date)", interval = "confidence",n= dim(cov.fr.dat[exp.fr$start.index[1]:exp.fr$end.index[1],])[1])$lower, gratia::derivatives(frgam1, term = "s(Date)", interval = "confidence",n= dim(cov.fr.dat[exp.fr$start.index[1]:exp.fr$end.index[1],])[1])$upper)/2)

frgam2 <-mgcv::gam(abs(cases)~ s(Weekday,bs = "cs",k=7) + s(Date, bs="tp", k=15) , 
                    data =  cov.fr.dat[exp.fr$start.index[2]:dim(cov.fr.dat)[1],],
                    family = gaussian(), method = "REML")
plot.gam(frgam2,pages=1)
fr.deriv2 <- data.frame(gratia::derivatives(frgam2, term = "s(Date)", interval = "confidence",n= dim(cov.fr.dat[exp.fr$start.index[2]:dim(cov.fr.dat)[1],])[1]),
                         "deriv.direction" = curve.direction(0, gratia::derivatives(frgam2, term = "s(Date)", interval = "confidence",n= dim(cov.fr.dat[exp.fr$start.index[2]:dim(cov.fr.dat)[1],])[1])$lower, gratia::derivatives(frgam2, term = "s(Date)", interval = "confidence",n= dim(cov.fr.dat[exp.fr$start.index[2]:dim(cov.fr.dat)[1],])[1])$upper)/2)
plot(fr.deriv2$deriv.direction)

cutoff.fr.cov.ews.multigam <- pbmclapply(cov.fr.dat$Date[1:(length(cov.fr.dat$Date)-13)],FUN = function(x){
  
  data.cut <- cov.fr.dat[as.numeric(cov.fr.dat$Date) >=  as.numeric(x),]
  
  tmp <- comp_EWS_wrapper(data.frame(timedat = as.numeric(data.cut$Date),
                                     biomass = data.cut$cases),
                          metrics = metrics,  threshold = 2, burn_in =7,
                          plotIt = F, ggplotIt = F, tail.direction = "one.tailed",
                          interpolate = F, method = "w_comp")
  tmp$cutoff <- paste(x)
  return(tmp)
  
}, mc.cores =3 )
names(cutoff.fr.cov.ews.multigam) <- cov.fr.dat$true.date[1:(length(cov.fr.dat$true.date)-13)]
temp <- cutoff.fr.cov.ews.multigam
cutoff.fr.cov.ews.multigam <- data.table::rbindlist(cutoff.fr.cov.ews.multigam, idcol = "start.date")%>%
  mutate(period = ifelse(time <= exp.fr$end.date[1], "first", #define periods of constancy prior to waves
                         ifelse(time >= exp.fr$start.date[2], "second","third")))
save(cutoff.fr.cov.ews.multigam,file = "Results/France/fr.multi.gam.ews.RData")
load("/Results/France/fr.multi.gam.ews.RData")

fr.multifirst.plot.dat <- cutoff.fr.cov.ews.multigam[cutoff.fr.cov.ews.multigam$start.date == "2020-02-18" & cutoff.fr.cov.ews.multigam$period == "first",]%>%
  left_join(data.frame("true.date" = cov.fr.dat$true.date,"time" = cov.fr.dat$Date), by = c("time"))

fr.multisecond.plot.dat <- cutoff.fr.cov.ews.multigam[cutoff.fr.cov.ews.multigam$start.date == cov.fr.dat$true.date[exp.fr$start.index[2]] & cutoff.fr.cov.ews.multigam$period == "second",]%>%
  left_join(data.frame("true.date" = cov.fr.dat$true.date,"time" = cov.fr.dat$Date), by = c("time"))

total.multifrplot.dat <- rbind(fr.multifirst.plot.dat,fr.multisecond.plot.dat) %>%
  rename(cases = count.used)

multi.deriv.cov.fr <- cov.fr.dat[,2:4]%>%
  cbind(rbind(fr.deriv1,fr.deriv2))


multi.fr.p2 <- ggplot(data = total.multifrplot.dat, aes(x=true.date,y=str)) +
  geom_hline(yintercept = 2, linetype="solid", color = "grey", size=1)+
  geom_line(aes(col= metric.code,alpha= "Indicator strength"))+
  geom_point(data =total.multifrplot.dat[which(total.multifrplot.dat$threshold.crossed==1)], aes(x=true.date, y = str, col = metric.code,alpha= "EWS detected")) +
  scale_alpha_manual(values = c(1, 1),
                     breaks = c("Indicator strength", "EWS detected"),
                     guide = guide_legend(override.aes = list(linetype = c(1, 0),shape = c(NA, 16)))) +
  #scale_colour_manual(values = scales::hue_pal()(7)) +
  ggthemes::theme_clean() + xlab("Date") + ylab("Strength of EWS") +
  geom_vline(xintercept = total.multifrplot.dat$true.date[total.multifrplot.dat$period != lag(total.multifrplot.dat$period)], linetype="dashed", 
             color = "black", size=1)+
  scale_x_date(date_breaks = "3 months", date_labels = "%b-%Y")+
  labs(color='EWS Indicator',alpha ="EWS") +
  theme(plot.margin = margin(c(10, 8, 0, 10)),
        legend.key.height = unit(0.3,"cm"),
        legend.key.width = unit(0.1,"cm"),
        legend.title = element_text(size = 10),
        legend.text = element_text(size=10))


multi.fr.p1 <- ggplot(multi.deriv.cov.fr, aes(x=true.date, y=cases)) +
  aes(group=NA)+
  geom_path(aes(col = as.factor(deriv.direction)))+
  geom_point(data= total.multifrplot.dat[which(total.multifrplot.dat$threshold.crossed==1 & total.multifrplot.dat$metric.code == "acf + SD + skew")], 
             aes(y =-5000, alpha = ""),size = 3,pch= "|",col = "#53B400")+
  ylab("COVID-19 Cases") + 
  xlab("Date")+
  geom_vline(xintercept = total.multifrplot.dat$true.date[total.multifrplot.dat$period != lag(total.multifrplot.dat$period)], 
             linetype="dashed",  color = "black", size=1)+
  scale_colour_manual(values=c("#22B4F5","black","#F07589"),name = "Predicted Nonlinear\nCase Trend", labels = c("Negative","No trend","Positive"))+
  scale_y_continuous(labels = scales::number_format(accuracy = 1000)) + 
  scale_x_date(date_breaks = "3 months", date_labels = "%b-%Y") +
  ggthemes::theme_clean() + ggtitle("France Daily COVID Cases: Two Waves")+  
  annotate("label", x = as.Date("2020-05-01"), y =max(total.multifrplot.dat$cases)*0.8 , label = "EWS indicator: acf + SD + skewness")+
  scale_alpha_manual(values = 0.8, labels = "EWS detected") +
  labs(alpha = "EWS")  +
  theme(plot.margin = margin(c(10, 8, 0, 10)),
        legend.key.height = unit(0.3,"cm" ),
        legend.key.width = unit(0.1,"cm"),
        legend.title = element_text(size = 10),
        legend.text = element_text(size=10))

png(file="Results/France/multigam_fr_fig.png",
    width=3000, height = 2500 ,res=300)
egg::ggarrange(multi.fr.p1,multi.fr.p2,nrow = 2,heights = c(2, 2), labels = c("a","b"))
dev.off()

fr.multiews.diff2 <- extract.ews.difference(total.multifrplot.dat,multi.deriv.cov.fr,"France",2)

###########################################################################################
### World - Germany ###
###########################################################################################
#https://covid19.who.int/info/
cov.ger.dat <- cov.WHO.dat %>%
  filter(Country == "Germany")%>%
  slice(55:n())

ggplot(cov.ger.dat,aes(x=Date, y = cases)) + geom_path() + 
  scale_x_continuous(breaks=seq(2020,2022,0.2)) + 
  scale_y_continuous( labels = scales::number_format(accuracy = 100))+
  ggtitle("Daily Germany COVID Cases")+ 
  ylab("Cases") + theme_bw()

#### Expanding GAM
exp.ger <- expanding.gam(cov.ger.dat,sensitivity = 7,train.length = 10)
plot(cov.ger.dat$cases)
abline(v=exp.ger$end.index)
gergam1 <-mgcv::gam(abs(cases)~ s(Weekday,bs = "cs",k=7)+s(Date, bs="tp", k=20), 
                    data = cov.ger.dat[exp.ger$start.index[1]:exp.ger$end.index[1],],
                    family = gaussian(), method = "REML")
plot.gam(gergam1,pages=1)
ger.deriv1 <- data.frame(gratia::derivatives(gergam1, term = "s(Date)", interval = "confidence",n= dim(cov.ger.dat[exp.ger$start.index[1]:exp.ger$end.index[1],])[1]),
                         "deriv.direction" = curve.direction(0, gratia::derivatives(gergam1, term = "s(Date)", interval = "confidence",n= dim(cov.ger.dat[exp.ger$start.index[1]:exp.ger$end.index[1],])[1])$lower, gratia::derivatives(gergam1, term = "s(Date)", interval = "confidence",n= dim(cov.ger.dat[exp.ger$start.index[1]:exp.ger$end.index[1],])[1])$upper)/2)
plot(ger.deriv1$deriv.direction)

gergam2 <-mgcv::gam(abs(cases)~ s(Weekday,bs = "cs",k=7) + s(Date, bs="tp", k=15) , 
                    data =  cov.ger.dat[exp.ger$start.index[2]:dim(cov.ger.dat)[1],],
                    family = gaussian(), method = "REML")
plot.gam(gergam2,pages=1)
ger.deriv2 <- data.frame(gratia::derivatives(gergam2, term = "s(Date)", interval = "confidence",n= dim(cov.ger.dat[exp.ger$start.index[2]:dim(cov.ger.dat)[1],])[1]),
                         "deriv.direction" = curve.direction(0, gratia::derivatives(gergam2, term = "s(Date)", interval = "confidence",n= dim(cov.ger.dat[exp.ger$start.index[2]:dim(cov.ger.dat)[1],])[1])$lower, gratia::derivatives(gergam2, term = "s(Date)", interval = "confidence",n= dim(cov.ger.dat[exp.ger$start.index[2]:dim(cov.ger.dat)[1],])[1])$upper)/2)
plot(ger.deriv2$deriv.direction)

cutoff.ger.cov.ews.multigam <- pbmclapply(cov.ger.dat$Date[1:(length(cov.ger.dat$Date)-13)],FUN = function(x){
  
  data.cut <- cov.ger.dat[as.numeric(cov.ger.dat$Date) >=  as.numeric(x),]
  
  tmp <- comp_EWS_wrapper(data.frame(timedat = as.numeric(data.cut$Date),
                                     biomass = data.cut$cases),
                          metrics = metrics,  threshold = 2, burn_in =7,
                          plotIt = F, ggplotIt = F, tail.direction = "one.tailed",
                          interpolate = F, method = "w_comp")
  tmp$cutoff <- paste(x)
  return(tmp)
  
}, mc.cores =3 )
names(cutoff.ger.cov.ews.multigam) <- cov.ger.dat$true.date[1:(length(cov.ger.dat$true.date)-13)]
temp <- cutoff.ger.cov.ews.multigam
cutoff.ger.cov.ews.multigam <- data.table::rbindlist(cutoff.ger.cov.ews.multigam, idcol = "start.date")%>%
  mutate(period = ifelse(time <= exp.ger$end.date[1], "first", #define periods of constancy prior to waves
                         ifelse(time >= exp.ger$start.date[2], "second","third")))
save(cutoff.ger.cov.ews.multigam,file = "Results/Germany/ger.multi.gam.ews.RData")
load("Results/Germany/ger.multi.gam.ews.RData")

ger.multifirst.plot.dat <- cutoff.ger.cov.ews.multigam[cutoff.ger.cov.ews.multigam$start.date == "2020-02-26" & cutoff.ger.cov.ews.multigam$period == "first",]%>%
  left_join(data.frame("true.date" = cov.ger.dat$true.date,"time" = cov.ger.dat$Date), by = c("time"))

ger.multisecond.plot.dat <- cutoff.ger.cov.ews.multigam[cutoff.ger.cov.ews.multigam$start.date == cov.ger.dat$true.date[exp.ger$start.index[2]] & cutoff.ger.cov.ews.multigam$period == "second",]%>%
  left_join(data.frame("true.date" = cov.ger.dat$true.date,"time" = cov.ger.dat$Date), by = c("time"))

total.multigerplot.dat <- rbind(ger.multifirst.plot.dat,ger.multisecond.plot.dat) %>%
  rename(cases = count.used)

multi.deriv.cov.ger <- cov.ger.dat[,2:4]%>%
  cbind(rbind(ger.deriv1,ger.deriv2))


multi.ger.p2 <- ggplot(data = total.multigerplot.dat, aes(x=true.date,y=str)) +
  geom_hline(yintercept = 2, linetype="solid", color = "grey", size=1)+
  geom_line(aes(col= metric.code,alpha= "Indicator strength"))+
  geom_point(data =total.multigerplot.dat[which(total.multigerplot.dat$threshold.crossed==1)], aes(x=true.date, y = str, col = metric.code,alpha= "EWS detected")) +
  scale_alpha_manual(values = c(1, 1),
                     breaks = c("Indicator strength", "EWS detected"),
                     guide = guide_legend(override.aes = list(linetype = c(1, 0),shape = c(NA, 16)))) +
  #scale_colour_manual(values = scales::hue_pal()(7)) +
  ggthemes::theme_clean() + xlab("Date") + ylab("Strength of EWS") +
  geom_vline(xintercept = total.multigerplot.dat$true.date[total.multigerplot.dat$period != lag(total.multigerplot.dat$period)], linetype="dashed", 
             color = "black", size=1)+
  scale_x_date(date_breaks = "3 months", date_labels = "%b-%Y")+
  labs(color='EWS Indicator',alpha ="EWS") +
  theme(plot.margin = margin(c(10, 8, 0, 10)),
        legend.key.height = unit(0.3,"cm"),
        legend.key.width = unit(0.1,"cm"),
        legend.title = element_text(size = 10),
        legend.text = element_text(size=10))


multi.ger.p1 <- ggplot(multi.deriv.cov.ger, aes(x=true.date, y=cases)) +
  aes(group=NA)+
  geom_path(aes(col = as.factor(deriv.direction)))+
  geom_point(data= total.multigerplot.dat[which(total.multigerplot.dat$threshold.crossed==1 & total.multigerplot.dat$metric.code == "acf + SD + skew")], 
             aes(y =-5000, alpha = ""),size = 3,pch= "|",col = "#53B400")+
  ylab("COVID-19 Cases") + 
  xlab("Date")+
  geom_vline(xintercept = total.multigerplot.dat$true.date[total.multigerplot.dat$period != lag(total.multigerplot.dat$period)], 
             linetype="dashed",  color = "black", size=1)+
  scale_colour_manual(values=c("#22B4F5","black","#F07589"),name = "Predicted Nonlinear\nCase Trend", labels = c("Negative","No trend","Positive"))+
  scale_y_continuous(labels = scales::number_format(accuracy = 1000)) + 
  scale_x_date(date_breaks = "3 months", date_labels = "%b-%Y") +
  ggthemes::theme_clean() + ggtitle("Germany Daily COVID Cases: Two Waves")+  
  annotate("label", x = as.Date("2020-05-01"), y =max(total.multigerplot.dat$cases)*0.8 , label = "EWS indicator: acf + SD + skewness")+
  scale_alpha_manual(values = 0.8, labels = "EWS detected") +
  labs(alpha = "EWS")  +
  theme(plot.margin = margin(c(10, 8, 0, 10)),
        legend.key.height = unit(0.3,"cm" ),
        legend.key.width = unit(0.1,"cm"),
        legend.title = element_text(size = 10),
        legend.text = element_text(size=10))

png(file="Results/Germany/multigam_ger_fig.png",
    width=3000, height = 2500 ,res=300)
egg::ggarrange(multi.ger.p1,multi.ger.p2,nrow = 2,heights = c(2, 2), labels = c("a","b"))
dev.off()

ger.multiews.diff2 <- extract.ews.difference(total.multigerplot.dat,multi.deriv.cov.ger,"Germany",2)

###########################################################################################
### World - India ###
###########################################################################################
#https://prsindia.org/covid-19/cases
cov.ind.WHO.dat <- cov.WHO.dat %>%
  filter(Country == "India")%>%
  slice(61:n())

ggplot(cov.ind.WHO.dat,aes(x=Date, y = cases)) + geom_path() + 
  scale_x_continuous(breaks=seq(2020,2022,0.2)) + 
  scale_y_continuous( labels = scales::number_format(accuracy = 100000))+
  ggtitle("Daily India COVID Cases")+ 
  ylab("Cases") + theme_bw()

#### Expanding GAM
exp.ind <- expanding.gam(cov.ind.WHO.dat,sensitivity = 7,train.length = 10)
plot(cov.ind.WHO.dat$cases)
abline(v=exp.ind$end.index)
indgam1 <-mgcv::gam(abs(cases)~ s(Weekday,bs = "cs",k=7)+s(Date, bs="tp", k=20), 
                    data = cov.ind.WHO.dat[exp.ind$start.index[1]:exp.ind$end.index[1],],
                    family = gaussian(), method = "REML")
plot.gam(indgam1,pages=1)
ind.deriv1 <- data.frame(gratia::derivatives(indgam1, term = "s(Date)", interval = "confidence",n= dim(cov.ind.WHO.dat[exp.ind$start.index[1]:exp.ind$end.index[1],])[1]),
                         "deriv.direction" = curve.direction(0, gratia::derivatives(indgam1, term = "s(Date)", interval = "confidence",n= dim(cov.ind.WHO.dat[exp.ind$start.index[1]:exp.ind$end.index[1],])[1])$lower, gratia::derivatives(indgam1, term = "s(Date)", interval = "confidence",n= dim(cov.ind.WHO.dat[exp.ind$start.index[1]:exp.ind$end.index[1],])[1])$upper)/2)
plot(ind.deriv1$deriv.direction)
indgam2 <-mgcv::gam(abs(cases)~ s(Weekday,bs = "cs",k=7) + s(Date, bs="tp", k=20) , 
                    data =  cov.ind.WHO.dat[(exp.ind$end.index[1]+1):(dim(cov.ind.WHO.dat)[1]),],
                    family = gaussian(), method = "REML")
plot.gam(indgam2,pages=1)
ind.deriv2 <- data.frame(gratia::derivatives(indgam2, term = "s(Date)", interval = "confidence",n= dim(cov.ind.WHO.dat[(exp.ind$end.index[1]+1):(dim(cov.ind.WHO.dat)[1]),])[1]),
                         "deriv.direction" = curve.direction(0, gratia::derivatives(indgam2, term = "s(Date)", interval = "confidence",n= dim(cov.ind.WHO.dat[(exp.ind$end.index[1]+1):(dim(cov.ind.WHO.dat)[1]),])[1])$lower, gratia::derivatives(indgam2, term = "s(Date)", interval = "confidence",n= dim(cov.ind.WHO.dat[(exp.ind$end.index[1]+1):(dim(cov.ind.WHO.dat)[1]),])[1])$upper)/2)
plot(ind.deriv2$deriv.direction)

cutoff.ind.cov.ews.multigam <- pbmclapply(cov.ind.WHO.dat$Date[1:(length(cov.ind.WHO.dat$Date)-13)],FUN = function(x){
  
  data.cut <- cov.ind.WHO.dat[as.numeric(cov.ind.WHO.dat$Date) >=  as.numeric(x),]
  
  tmp <- comp_EWS_wrapper(data.frame(timedat = as.numeric(data.cut$Date),
                                     biomass = data.cut$cases),
                          metrics = metrics,  threshold = 2, burn_in =7,
                          plotIt = F, ggplotIt = F, tail.direction = "one.tailed",
                          interpolate = F, method = "w_comp")
  tmp$cutoff <- paste(x)
  return(tmp)
  
}, mc.cores =3 )
names(cutoff.ind.cov.ews.multigam) <- cov.ind.WHO.dat$true.date[1:(length(cov.ind.WHO.dat$true.date)-13)]
temp <- cutoff.ind.cov.ews.multigam
cutoff.ind.cov.ews.multigam <- temp
cutoff.ind.cov.ews.multigam <- data.table::rbindlist(cutoff.ind.cov.ews.multigam, idcol = "start.date")%>%
  mutate(period = ifelse(time <= exp.ind$end.date[1], "first", #define periods of constancy prior to waves
                         ifelse(time >exp.ind$end.date[1] , "second","third")))
save(cutoff.ind.cov.ews.multigam,file = "Results/India/ind.multi.gam.ews.RData")
load("Results/India/ind.multi.gam.ews.RData")

ind.multifirst.plot.dat <- cutoff.ind.cov.ews.multigam[cutoff.ind.cov.ews.multigam$start.date == "2020-03-03" & cutoff.ind.cov.ews.multigam$period == "first",]%>%
  left_join(data.frame("true.date" = cov.ind.WHO.dat$true.date,"time" = cov.ind.WHO.dat$Date), by = c("time"))

ind.multisecond.plot.dat <- cutoff.ind.cov.ews.multigam[cutoff.ind.cov.ews.multigam$start.date == cov.ind.WHO.dat$true.date[exp.ind$end.index[1]+1] & cutoff.ind.cov.ews.multigam$period == "second",]%>%
  left_join(data.frame("true.date" = cov.ind.WHO.dat$true.date,"time" = cov.ind.WHO.dat$Date), by = c("time"))

total.multiindplot.dat <- rbind(ind.multifirst.plot.dat,ind.multisecond.plot.dat) %>%
  rename(cases = count.used)

multi.deriv.cov.ind <- cov.ind.WHO.dat[,2:4]%>%
  cbind(rbind(ind.deriv1,ind.deriv2))


multi.ind.p2 <- ggplot(data = total.multiindplot.dat, aes(x=true.date,y=str)) +
  geom_hline(yintercept = 2, linetype="solid", color = "grey", size=1)+
  geom_line(aes(col= metric.code,alpha= "Indicator strength"))+
  geom_point(data =total.multiindplot.dat[which(total.multiindplot.dat$threshold.crossed==1)], aes(x=true.date, y = str, col = metric.code,alpha= "EWS detected")) +
  scale_alpha_manual(values = c(1, 1),
                     breaks = c("Indicator strength", "EWS detected"),
                     guide = guide_legend(override.aes = list(linetype = c(1, 0),shape = c(NA, 16)))) +
  #scale_colour_manual(values = scales::hue_pal()(7)) +
  ggthemes::theme_clean() + xlab("Date") + ylab("Strength of EWS") +
  geom_vline(xintercept = total.multiindplot.dat$true.date[total.multiindplot.dat$period != lag(total.multiindplot.dat$period)], linetype="dashed", 
             color = "black", size=1)+
  scale_x_date(date_breaks = "3 months", date_labels = "%b-%Y")+
  labs(color='EWS Indicator',alpha ="EWS") +
  theme(plot.margin = margin(c(10, 8, 0, 10)),
        legend.key.height = unit(0.3,"cm"),
        legend.key.width = unit(0.1,"cm"),
        legend.title = element_text(size = 10),
        legend.text = element_text(size=10))


multi.ind.p1 <- ggplot(multi.deriv.cov.ind, aes(x=true.date, y=cases)) +
  aes(group=NA)+
  geom_path(aes(col = as.factor(deriv.direction)))+
  geom_point(data= total.multiindplot.dat[which(total.multiindplot.dat$threshold.crossed==1 & total.multiindplot.dat$metric.code == "acf + SD + skew")], 
             aes(y =-20000, alpha = ""),size = 3,pch= "|",col = "#53B400")+
  ylab("COVID-19 Cases") + 
  xlab("Date")+
  geom_vline(xintercept = total.multiindplot.dat$true.date[total.multiindplot.dat$period != lag(total.multiindplot.dat$period)], 
             linetype="dashed",  color = "black", size=1)+
  scale_colour_manual(values=c("#22B4F5","black","#F07589"),name = "Predicted Nonlinear\nCase Trend", labels = c("Negative","No trend","Positive"))+
  scale_y_continuous(labels = scales::number_format(accuracy = 1000)) + 
  scale_x_date(date_breaks = "3 months", date_labels = "%b-%Y") +
  ggthemes::theme_clean() + ggtitle("India Daily COVID Cases: Two Waves")+  
  annotate("label", x = as.Date("2020-06-01"), y =max(total.multiindplot.dat$cases)*0.8 , label = "EWS indicator: acf + SD + skewness")+
  scale_alpha_manual(values = 0.8, labels = "EWS detected") +
  labs(alpha = "EWS")  +
  theme(plot.margin = margin(c(10, 8, 0, 10)),
        legend.key.height = unit(0.3,"cm" ),
        legend.key.width = unit(0.1,"cm"),
        legend.title = element_text(size = 10),
        legend.text = element_text(size=10))

png(file="Results/India/multigam_ind_fig.png",
    width=3000, height = 2500 ,res=300)
egg::ggarrange(multi.ind.p1,multi.ind.p2,nrow = 2,heights = c(2, 2), labels = c("a","b"))
dev.off()

ind.multiews.diff2 <- extract.ews.difference(total.multiindplot.dat,multi.deriv.cov.ind,"India",2)

###########################################################################################
### World - Israel ###
###########################################################################################
#https://covid19.who.int/info/
cov.isr.dat <- cov.WHO.dat %>%
  filter(Country == "Israel")%>%
  slice(64:n())


ggplot(cov.isr.dat,aes(x=Date, y = cases)) + geom_path() + 
  scale_x_continuous(breaks=seq(2020,2022,0.2)) + 
  scale_y_continuous( labels = scales::number_format(accuracy = 100))+
  ggtitle("Daily Israel COVID Cases")+ 
  ylab("Cases") + theme_bw()

#### Expanding GAM
exp.isr <- expanding.gam(cov.isr.dat,sensitivity = 7,train.length = 10)
plot(cov.isr.dat$cases)
abline(v=exp.isr$end.index)

isrgam1 <-mgcv::gam(abs(cases)~ s(Weekday,bs = "cs",k=7)+s(Date, bs="tp", k=15), 
                    data = cov.isr.dat[exp.isr$start.index[1]:exp.isr$end.index[1],],
                    family = gaussian(), method = "REML")
plot.gam(isrgam1,pages=1)
isr.deriv1 <- data.frame(gratia::derivatives(isrgam1, term = "s(Date)", interval = "confidence",n= dim(cov.isr.dat[exp.isr$start.index[1]:exp.isr$end.index[1],])[1]),
                         "deriv.direction" = curve.direction(0, gratia::derivatives(isrgam1, term = "s(Date)", interval = "confidence",n= dim(cov.isr.dat[exp.isr$start.index[1]:exp.isr$end.index[1],])[1])$lower, gratia::derivatives(isrgam1, term = "s(Date)", interval = "confidence",n= dim(cov.isr.dat[exp.isr$start.index[1]:exp.isr$end.index[1],])[1])$upper)/2)

isrgam2 <-mgcv::gam(abs(cases)~ s(Weekday,bs = "cs",k=7) + s(Date, bs="tp", k=15) , 
                    data =  cov.isr.dat[exp.isr$start.index[2]:exp.isr$end.index[3],],
                    family = gaussian(), method = "REML")
plot.gam(isrgam2,pages=1)
isr.deriv2 <- data.frame(gratia::derivatives(isrgam2, term = "s(Date)", interval = "confidence",n= dim(cov.isr.dat[exp.isr$start.index[2]:exp.isr$end.index[3],])[1]),
                         "deriv.direction" = curve.direction(0, gratia::derivatives(isrgam2, term = "s(Date)", interval = "confidence",n= dim(cov.isr.dat[exp.isr$start.index[2]:exp.isr$end.index[3],])[1])$lower, gratia::derivatives(isrgam2, term = "s(Date)", interval = "confidence",n= dim(cov.isr.dat[exp.isr$start.index[2]:exp.isr$end.index[3],])[1])$upper)/2)
plot(isr.deriv2$deriv.direction)

isrgam3 <-mgcv::gam(abs(cases)~ s(Weekday,bs = "cs",k=7) + s(Date, bs="tp", k=20) , 
                    data =  cov.isr.dat[exp.isr$start.index[4]:dim(cov.isr.dat)[1],],
                    family = gaussian(), method = "REML")
plot.gam(isrgam3,pages=1)
isr.deriv3 <- data.frame(gratia::derivatives(isrgam3, term = "s(Date)", interval = "confidence",n= dim(cov.isr.dat[exp.isr$start.index[4]:dim(cov.isr.dat)[1],])[1]),
                         "deriv.direction" = curve.direction(0, gratia::derivatives(isrgam3, term = "s(Date)", interval = "confidence",n= dim(cov.isr.dat[exp.isr$start.index[4]:dim(cov.isr.dat)[1],])[1])$lower, gratia::derivatives(isrgam3, term = "s(Date)", interval = "confidence",n= dim(cov.isr.dat[exp.isr$start.index[4]:dim(cov.isr.dat)[1],])[1])$upper)/2)
plot(isr.deriv3$deriv.direction)

cutoff.isr.cov.ews.multigam <- pbmclapply(cov.isr.dat$Date[1:(length(cov.isr.dat$Date)-13)],FUN = function(x){
  
  data.cut <- cov.isr.dat[as.numeric(cov.isr.dat$Date) >=  as.numeric(x),]
  
  tmp <- comp_EWS_wrapper(data.frame(timedat = as.numeric(data.cut$Date),
                                     biomass = data.cut$cases),
                          metrics = metrics,  threshold = 2, burn_in =7,
                          plotIt = F, ggplotIt = F, tail.direction = "one.tailed",
                          interpolate = F, method = "w_comp")
  tmp$cutoff <- paste(x)
  return(tmp)
  
}, mc.cores =3 )
names(cutoff.isr.cov.ews.multigam) <- cov.isr.dat$true.date[1:(length(cov.isr.dat$true.date)-13)]
temp <- cutoff.isr.cov.ews.multigam
cutoff.isr.cov.ews.multigam <- data.table::rbindlist(cutoff.isr.cov.ews.multigam, idcol = "start.date")%>%
  mutate(period = ifelse(time <= exp.isr$end.date[1], "first", #define periods of constancy prior to waves
                         ifelse(time >= exp.isr$start.date[2] & time <= exp.isr$end.date[3], "second",
                                ifelse(time >= exp.isr$start.date[4],"third","fourth"))))
save(cutoff.isr.cov.ews.multigam,file = "Results/Israel/isr.multi.gam.ews.RData")
load("Results/Israel/isr.multi.gam.ews.RData")

isr.multifirst.plot.dat <- cutoff.isr.cov.ews.multigam[cutoff.isr.cov.ews.multigam$start.date == "2020-03-06" & cutoff.isr.cov.ews.multigam$period == "first",]%>%
  left_join(data.frame("true.date" = cov.isr.dat$true.date,"time" = cov.isr.dat$Date), by = c("time"))

isr.multisecond.plot.dat <- cutoff.isr.cov.ews.multigam[cutoff.isr.cov.ews.multigam$start.date == cov.isr.dat$true.date[exp.isr$start.index[2]] & cutoff.isr.cov.ews.multigam$period == "second",]%>%
  left_join(data.frame("true.date" = cov.isr.dat$true.date,"time" = cov.isr.dat$Date), by = c("time"))

isr.multithird.plot.dat <- cutoff.isr.cov.ews.multigam[cutoff.isr.cov.ews.multigam$start.date == cov.isr.dat$true.date[exp.isr$start.index[4]] & cutoff.isr.cov.ews.multigam$period == "third",]%>%
  left_join(data.frame("true.date" = cov.isr.dat$true.date,"time" = cov.isr.dat$Date), by = c("time"))

total.multiisrplot.dat <- rbind(isr.multifirst.plot.dat,isr.multisecond.plot.dat,isr.multithird.plot.dat) %>%
  rename(cases = count.used)

multi.deriv.cov.isr <- cov.isr.dat[,2:4]%>%
  cbind(rbind(isr.deriv1,isr.deriv2,isr.deriv3))


multi.isr.p2 <- ggplot(data = total.multiisrplot.dat, aes(x=true.date,y=str)) +
  geom_hline(yintercept = 2, linetype="solid", color = "grey", size=1)+
  geom_line(aes(col= metric.code,alpha= "Indicator strength"))+
  geom_point(data =total.multiisrplot.dat[which(total.multiisrplot.dat$threshold.crossed==1)], aes(x=true.date, y = str, col = metric.code,alpha= "EWS detected")) +
  scale_alpha_manual(values = c(1, 1),
                     breaks = c("Indicator strength", "EWS detected"),
                     guide = guide_legend(override.aes = list(linetype = c(1, 0),shape = c(NA, 16)))) +
  #scale_colour_manual(values = scales::hue_pal()(7)) +
  ggthemes::theme_clean() + xlab("Date") + ylab("Strength of EWS") +
  geom_vline(xintercept = total.multiisrplot.dat$true.date[total.multiisrplot.dat$period != lag(total.multiisrplot.dat$period)], linetype="dashed", 
             color = "black", size=1)+
  scale_x_date(date_breaks = "3 months", date_labels = "%b-%Y")+
  labs(color='EWS Indicator',alpha ="EWS") +
  theme(plot.margin = margin(c(10, 8, 0, 10)),
        legend.key.height = unit(0.3,"cm"),
        legend.key.width = unit(0.1,"cm"),
        legend.title = element_text(size = 10),
        legend.text = element_text(size=10))


multi.isr.p1 <- ggplot(multi.deriv.cov.isr, aes(x=true.date, y=cases)) +
  aes(group=NA)+
  geom_path(aes(col = as.factor(deriv.direction)))+
  geom_point(data= total.multiisrplot.dat[which(total.multiisrplot.dat$threshold.crossed==1 & total.multiisrplot.dat$metric.code == "acf + SD + skew")], 
                          aes(y =-1000, alpha = ""),size = 3,pch= "|",col = "#53B400")+
  ylab("COVID-19 Cases") + 
  xlab("Date")+
  geom_vline(xintercept = total.multiisrplot.dat$true.date[total.multiisrplot.dat$period != lag(total.multiisrplot.dat$period)], 
             linetype="dashed",  color = "black", size=1)+
  scale_colour_manual(values=c("#22B4F5","black","#F07589"),name = "Predicted Nonlinear\nCase Trend", labels = c("Negative","No trend","Positive"))+
  scale_y_continuous(labels = scales::number_format(accuracy = 1000)) + 
  scale_x_date(date_breaks = "3 months", date_labels = "%b-%Y") +
  ggthemes::theme_clean() + ggtitle("Israel Daily COVID Cases: Three Waves")+  
  annotate("label", x = as.Date("2020-06-01"), y =max(total.multiisrplot.dat$cases)*0.8 , label = "EWS indicator: acf + SD + skewness")+
  scale_alpha_manual(values = 0.8, labels = "EWS detected") +
  labs(alpha = "EWS")  +
  theme(plot.margin = margin(c(10, 8, 0, 10)),
        legend.key.height = unit(0.3,"cm" ),
        legend.key.width = unit(0.1,"cm"),
        legend.title = element_text(size = 10),
        legend.text = element_text(size=10))

png(file="Results/Israel/multigam_isr_fig.png",
    width=3000, height = 2500 ,res=300)
egg::ggarrange(multi.isr.p1,multi.isr.p2,nrow = 2,heights = c(2, 2), labels = c("a","b"))
dev.off()

isr.multiews.diff2 <- extract.ews.difference(total.multiisrplot.dat,multi.deriv.cov.isr,"Israel",2)

###########################################################################################
### World - Italy ###
###########################################################################################
#https://covid19.who.int/info/
cov.itl.dat <- cov.WHO.dat %>%
  filter(Country == "Italy") %>%
  slice(51:n())

ggplot(cov.itl.dat,aes(x=Date, y = cases)) + geom_path() + 
  scale_x_continuous(breaks=seq(2020,2022,0.2)) + 
  scale_y_continuous( labels = scales::number_format(accuracy = 100))+
  ggtitle("Daily Italy COVID Cases")+ 
  ylab("Cases") + theme_bw()

#### Expanding GAM
exp.itl <- expanding.gam(cov.itl.dat,sensitivity = 7,train.length = 10)
plot(cov.itl.dat$cases)
abline(v=exp.itl$end.index)
itlgam1 <-mgcv::gam(abs(cases)~ s(Weekday,bs = "cs",k=7)+s(Date, bs="tp", k=15), 
                   data = cov.itl.dat[exp.itl$start.index[1]:exp.itl$end.index[1],],
                   family = gaussian(), method = "REML")
plot.gam(itlgam1,pages=1)
itl.deriv1 <- data.frame(gratia::derivatives(itlgam1, term = "s(Date)", interval = "confidence",n= dim(cov.itl.dat[exp.itl$start.index[1]:exp.itl$end.index[1],])[1]),
                        "deriv.direction" = curve.direction(0, gratia::derivatives(itlgam1, term = "s(Date)", interval = "confidence",n= dim(cov.itl.dat[exp.itl$start.index[1]:exp.itl$end.index[1],])[1])$lower, gratia::derivatives(itlgam1, term = "s(Date)", interval = "confidence",n= dim(cov.itl.dat[exp.itl$start.index[1]:exp.itl$end.index[1],])[1])$upper)/2)

itlgam2 <-mgcv::gam(abs(cases)~ s(Weekday,bs = "cs",k=7) + s(Date, bs="tp", k=20) , 
                   data =  cov.itl.dat[exp.itl$start.index[2]:dim(cov.itl.dat)[1],],
                   family = gaussian(), method = "REML")
plot.gam(itlgam2,pages=1)
itl.deriv2 <- data.frame(gratia::derivatives(itlgam2, term = "s(Date)", interval = "confidence",n= dim(cov.itl.dat[exp.itl$start.index[2]:dim(cov.itl.dat)[1],])[1]),
                        "deriv.direction" = curve.direction(0, gratia::derivatives(itlgam2, term = "s(Date)", interval = "confidence",n= dim(cov.itl.dat[exp.itl$start.index[2]:dim(cov.itl.dat)[1],])[1])$lower, gratia::derivatives(itlgam2, term = "s(Date)", interval = "confidence",n= dim(cov.itl.dat[exp.itl$start.index[2]:dim(cov.itl.dat)[1],])[1])$upper)/2)
plot(itl.deriv2$deriv.direction)

cutoff.itl.cov.ews.multigam <- pbmclapply(cov.itl.dat$Date[1:(length(cov.itl.dat$Date)-13)],FUN = function(x){
  
  data.cut <- cov.itl.dat[as.numeric(cov.itl.dat$Date) >=  as.numeric(x),]
  
  tmp <- comp_EWS_wrapper(data.frame(timedat = as.numeric(data.cut$Date),
                                     biomass = data.cut$cases),
                          metrics = metrics,  threshold = 2, burn_in =7,
                          plotIt = F, ggplotIt = F, tail.direction = "one.tailed",
                          interpolate = F, method = "w_comp")
  tmp$cutoff <- paste(x)
  return(tmp)
  
}, mc.cores =3 )
names(cutoff.itl.cov.ews.multigam) <- cov.itl.dat$true.date[1:(length(cov.itl.dat$true.date)-13)]
temp <- cutoff.itl.cov.ews.multigam
cutoff.itl.cov.ews.multigam <- data.table::rbindlist(cutoff.itl.cov.ews.multigam, idcol = "start.date")%>%
  mutate(period = ifelse(time <= exp.itl$end.date[1], "first", #define periods of constancy prior to waves
                         ifelse(time >= exp.itl$start.date[2] , "second","third")))
save(cutoff.itl.cov.ews.multigam,file = "Results/Italy/itl.multi.gam.ews.RData")
load("Results/Italy/itl.multi.gam.ews.RData")

itl.multifirst.plot.dat <- cutoff.itl.cov.ews.multigam[cutoff.itl.cov.ews.multigam$start.date == "2020-02-22" & cutoff.itl.cov.ews.multigam$period == "first",]%>%
  left_join(data.frame("true.date" = cov.itl.dat$true.date,"time" = cov.itl.dat$Date), by = c("time"))

itl.multisecond.plot.dat <- cutoff.itl.cov.ews.multigam[cutoff.itl.cov.ews.multigam$start.date == cov.itl.dat$true.date[exp.itl$start.index[2]] & cutoff.itl.cov.ews.multigam$period == "second",]%>%
  left_join(data.frame("true.date" = cov.itl.dat$true.date,"time" = cov.itl.dat$Date), by = c("time"))

total.multiitlplot.dat <- rbind(itl.multifirst.plot.dat,itl.multisecond.plot.dat) %>%
  rename(cases = count.used)

multi.deriv.cov.itl <- cov.itl.dat[,2:4]%>%
  cbind(rbind(itl.deriv1,itl.deriv2))


multi.itl.p2 <- ggplot(data = total.multiitlplot.dat, aes(x=true.date,y=str)) +
  geom_hline(yintercept = 2, linetype="solid", color = "grey", size=1)+
  geom_line(aes(col= metric.code,alpha= "Indicator strength"))+
  geom_point(data =total.multiitlplot.dat[which(total.multiitlplot.dat$threshold.crossed==1)], aes(x=true.date, y = str, col = metric.code,alpha= "EWS detected")) +
  scale_alpha_manual(values = c(1, 1),
                     breaks = c("Indicator strength", "EWS detected"),
                     guide = guide_legend(override.aes = list(linetype = c(1, 0),shape = c(NA, 16)))) +
  #scale_colour_manual(values = scales::hue_pal()(7)) +
  ggthemes::theme_clean() + xlab("Date") + ylab("Strength of EWS") +
  geom_vline(xintercept = total.multiitlplot.dat$true.date[total.multiitlplot.dat$period != lag(total.multiitlplot.dat$period)], linetype="dashed", 
             color = "black", size=1)+
  scale_x_date(date_breaks = "3 months", date_labels = "%b-%Y")+
  labs(color='EWS Indicator',alpha ="EWS") +
  theme(plot.margin = margin(c(10, 8, 0, 10)),
        legend.key.height = unit(0.3,"cm"),
        legend.key.width = unit(0.1,"cm"),
        legend.title = element_text(size = 10),
        legend.text = element_text(size=10))


multi.itl.p1 <- ggplot(multi.deriv.cov.itl, aes(x=true.date, y=cases)) +
  aes(group=NA)+
  geom_path(aes(col = as.factor(deriv.direction)))+
  geom_point(data= total.multiitlplot.dat[which(total.multiitlplot.dat$threshold.crossed==1 & total.multiitlplot.dat$metric.code == "acf + SD + skew")], 
             aes(y =-5000, alpha = ""),size = 3,pch= "|",col = "#53B400")+
  ylab("COVID-19 Cases") + 
  xlab("Date")+
  geom_vline(xintercept = total.multiitlplot.dat$true.date[total.multiitlplot.dat$period != lag(total.multiitlplot.dat$period)], 
             linetype="dashed",  color = "black", size=1)+
  scale_colour_manual(values=c("#22B4F5","black","#F07589"),name = "Predicted Nonlinear\nCase Trend", labels = c("Negative","No trend","Positive"))+
  scale_y_continuous(labels = scales::number_format(accuracy = 1000)) + 
  scale_x_date(date_breaks = "3 months", date_labels = "%b-%Y") +
  ggthemes::theme_clean() + ggtitle("Italy Daily COVID Cases: Two Waves")+  
  annotate("label", x = as.Date("2020-05-01"), y =max(total.multiitlplot.dat$cases)*0.8 , label = "EWS indicator: acf + SD + skewness")+
  scale_alpha_manual(values = 0.8, labels = "EWS detected") +
  labs(alpha = "EWS")  +
  theme(plot.margin = margin(c(10, 8, 0, 10)),
        legend.key.height = unit(0.3,"cm" ),
        legend.key.width = unit(0.1,"cm"),
        legend.title = element_text(size = 10),
        legend.text = element_text(size=10))

png(file="Results/Italy/multigam_itl_fig.png",
    width=3000, height = 2500 ,res=300)
egg::ggarrange(multi.itl.p1,multi.itl.p2,nrow = 2,heights = c(2, 2), labels = c("a","b"))
dev.off()

itl.multiews.diff2 <- extract.ews.difference(total.multiitlplot.dat,multi.deriv.cov.itl,"Italy",2)

###########################################################################################
### Japan ###
###########################################################################################
#https://www.mhlw.go.jp/stf/covid-19/open-data.html
ggplot(cov.jpn.dat,aes(x=Date, y = cases)) + geom_path() + 
  scale_x_continuous(breaks=seq(2020,2022,0.2)) + 
  scale_y_continuous( labels = scales::number_format(accuracy = 1000))+
  ggtitle("Daily Japan COVID Cases")+ 
  ylab("Cases") + theme_bw()


#### Expanding GAM
exp.jpn <- expanding.gam(cov.jpn.dat,sensitivity = 7,train.length = 10)
plot(cov.jpn.dat$cases)
abline(v=exp.jpn$end.index)
jpngam1 <-mgcv::gam(abs(cases)~ s(Weekday,bs = "cs",k=7)+s(Date, bs="tp", k=20), 
                    data = cov.jpn.dat[exp.jpn$start.index[1]:exp.jpn$end.index[1],],
                    family = gaussian(), method = "REML")
plot.gam(jpngam1,pages=1)
jpn.deriv1 <- data.frame(gratia::derivatives(jpngam1, term = "s(Date)", interval = "confidence",n= dim(cov.jpn.dat[exp.jpn$start.index[1]:exp.jpn$end.index[1],])[1]),
                          "deriv.direction" = curve.direction(0, gratia::derivatives(jpngam1, term = "s(Date)", interval = "confidence",n= dim(cov.jpn.dat[exp.jpn$start.index[1]:exp.jpn$end.index[1],])[1])$lower, gratia::derivatives(jpngam1, term = "s(Date)", interval = "confidence",n= dim(cov.jpn.dat[exp.jpn$start.index[1]:exp.jpn$end.index[1],])[1])$upper)/2)
plot(jpn.deriv1$deriv.direction)

jpngam2 <-mgcv::gam(abs(cases)~ s(Weekday,bs = "cs",k=7) + s(Date, bs="tp", k=10) , 
                    data =  cov.jpn.dat[exp.jpn$start.index[2]:exp.jpn$end.index[3],],
                    family = gaussian(), method = "REML")
plot.gam(jpngam2,pages=1)
jpn.deriv2 <- data.frame(gratia::derivatives(jpngam2, term = "s(Date)", interval = "confidence",n= dim(cov.jpn.dat[exp.jpn$start.index[2]:exp.jpn$end.index[3],])[1]),
                          "deriv.direction" = curve.direction(0, gratia::derivatives(jpngam2, term = "s(Date)", interval = "confidence",n= dim(cov.jpn.dat[exp.jpn$start.index[2]:exp.jpn$end.index[3],])[1])$lower, gratia::derivatives(jpngam2, term = "s(Date)", interval = "confidence",n= dim(cov.jpn.dat[exp.jpn$start.index[2]:exp.jpn$end.index[3],])[1])$upper)/2)
plot(jpn.deriv2$deriv.direction)
jpngam3 <-mgcv::gam(abs(cases)~ s(Weekday,bs = "cs",k=7)+s(Date, bs="tp", k=10) , 
                    data = cov.jpn.dat[exp.jpn$start.index[4]:exp.jpn$end.index[4],],
                    family = gaussian(), method = "REML")
plot.gam(jpngam3,pages =1)
jpn.deriv3 <- data.frame(gratia::derivatives(jpngam3, term = "s(Date)", interval = "confidence",n= dim(cov.jpn.dat[exp.jpn$start.index[4]:exp.jpn$end.index[4],])[1]),
                          "deriv.direction" = curve.direction(0, gratia::derivatives(jpngam3, term = "s(Date)", interval = "confidence",n= dim(cov.jpn.dat[exp.jpn$start.index[4]:exp.jpn$end.index[4],])[1])$lower, gratia::derivatives(jpngam3, term = "s(Date)", interval = "confidence",n= dim(cov.jpn.dat[exp.jpn$start.index[4]:exp.jpn$end.index[4],])[1])$upper)/2)

jpngam4 <-mgcv::gam(abs(cases)~ s(Weekday,bs = "cs",k=7)+s(Date, bs="tp", k=10) , 
                    data = cov.jpn.dat[exp.jpn$start.index[5]:dim(cov.jpn.dat)[1],],
                    family = gaussian(), method = "REML")

plot.gam(jpngam4,pages =1)
jpn.deriv4 <- data.frame(gratia::derivatives(jpngam4, term = "s(Date)", interval = "confidence",n= dim(cov.jpn.dat[exp.jpn$start.index[5]:dim(cov.jpn.dat)[1],])[1]),
                         "deriv.direction" = curve.direction(0, gratia::derivatives(jpngam4, term = "s(Date)", interval = "confidence",n= dim(cov.jpn.dat[exp.jpn$start.index[5]:dim(cov.jpn.dat)[1],])[1])$lower, gratia::derivatives(jpngam4, term = "s(Date)", interval = "confidence",n= dim(cov.jpn.dat[exp.jpn$start.index[5]:dim(cov.jpn.dat)[1],])[1])$upper)/2)



cutoff.jpn.cov.ews.multigam <- pbmclapply(cov.jpn.dat$Date[1:(length(cov.jpn.dat$Date)-13)],FUN = function(x){
  
  data.cut <- cov.jpn.dat[as.numeric(cov.jpn.dat$Date) >=  as.numeric(x),]
  
  tmp <- comp_EWS_wrapper(data.frame(timedat = as.numeric(data.cut$Date),
                                      biomass = data.cut$cases),
                          metrics = metrics,  threshold = 2, burn_in =7,
                          plotIt = F, ggplotIt = F, tail.direction = "one.tailed",
                          interpolate = F, method = "w_comp")
  tmp$cutoff <- paste(x)
  return(tmp)
  
}, mc.cores =3 )
names(cutoff.jpn.cov.ews.multigam) <- cov.jpn.dat$true.date[1:(length(cov.jpn.dat$true.date)-13)]
temp <- cutoff.jpn.cov.ews.multigam
cutoff.jpn.cov.ews.multigam <- temp
cutoff.jpn.cov.ews.multigam <- data.table::rbindlist(cutoff.jpn.cov.ews.multigam, idcol = "start.date")%>%
  mutate(period = ifelse(time <= exp.jpn$end.date[1], "first", #define periods of constancy prior to waves
                         ifelse(time >= exp.jpn$start.date[2] & time <= exp.jpn$end.date[3], "second",
                                ifelse(time >= exp.jpn$start.date[4] & time <= exp.jpn$end.date[4],"third",
                                       ifelse(time >= exp.jpn$start.date[5], "fourth","fifth")))))
save(cutoff.jpn.cov.ews.multigam,file = "Results/Japan/jpn.multi.gam.ews.RData")
load("Results/Japan/jpn.multi.gam.ews.RData")

jpn.multifirst.plot.dat <- cutoff.jpn.cov.ews.multigam[cutoff.jpn.cov.ews.multigam$start.date == "2020-02-12" & cutoff.jpn.cov.ews.multigam$period == "first",]%>%
  left_join(data.frame("true.date" = cov.jpn.dat$true.date,"time" = cov.jpn.dat$Date), by = c("time"))

jpn.multisecond.plot.dat <- cutoff.jpn.cov.ews.multigam[cutoff.jpn.cov.ews.multigam$start.date == cov.jpn.dat$true.date[exp.jpn$start.index[2]] & cutoff.jpn.cov.ews.multigam$period == "second",]%>%
  left_join(data.frame("true.date" = cov.jpn.dat$true.date,"time" = cov.jpn.dat$Date), by = c("time"))

jpn.multithird.plot.dat <- cutoff.jpn.cov.ews.multigam[cutoff.jpn.cov.ews.multigam$start.date == cov.jpn.dat$true.date[exp.jpn$start.index[4]] & cutoff.jpn.cov.ews.multigam$period == "third",]%>%
  left_join(data.frame("true.date" = cov.jpn.dat$true.date,"time" = cov.jpn.dat$Date), by = c("time"))

jpn.multifourth.plot.dat <- cutoff.jpn.cov.ews.multigam[cutoff.jpn.cov.ews.multigam$start.date == cov.jpn.dat$true.date[exp.jpn$start.index[5]] & cutoff.jpn.cov.ews.multigam$period == "fourth",]%>%
  left_join(data.frame("true.date" = cov.jpn.dat$true.date,"time" = cov.jpn.dat$Date), by = c("time"))

total.multijpnplot.dat <- rbind(jpn.multifirst.plot.dat,jpn.multisecond.plot.dat,jpn.multithird.plot.dat,jpn.multifourth.plot.dat) %>%
  rename(cases = count.used)

multi.deriv.cov.jpn <- cov.jpn.dat[,2:4]%>%
  cbind(rbind(jpn.deriv1,jpn.deriv2,jpn.deriv3,jpn.deriv4))


multi.jpn.p2 <- ggplot(data = total.multijpnplot.dat, aes(x=true.date,y=str)) +
  geom_hline(yintercept = 2, linetype="solid", color = "grey", size=1)+
  geom_line(aes(col= metric.code,alpha= "Indicator strength"))+
  geom_point(data =total.multijpnplot.dat[which(total.multijpnplot.dat$threshold.crossed==1)], aes(x=true.date, y = str, col = metric.code,alpha= "EWS detected")) +
  scale_alpha_manual(values = c(1, 1),
                     breaks = c("Indicator strength", "EWS detected"),
                     guide = guide_legend(override.aes = list(linetype = c(1, 0),shape = c(NA, 16)))) +
  #scale_colour_manual(values = scales::hue_pal()(7)) +
  ggthemes::theme_clean() + xlab("Date") + ylab("Strength of EWS") +
  geom_vline(xintercept = total.multijpnplot.dat$true.date[total.multijpnplot.dat$period != lag(total.multijpnplot.dat$period)], linetype="dashed", 
             color = "black", size=1)+
  scale_x_date(date_breaks = "3 months", date_labels = "%b-%Y")+
  labs(color='EWS Indicator',alpha ="EWS") +
  theme(plot.margin = margin(c(10, 8, 0, 10)),
        legend.key.height = unit(0.3,"cm"),
        legend.key.width = unit(0.1,"cm"),
        legend.title = element_text(size = 10),
        legend.text = element_text(size=10))


multi.jpn.p1 <- ggplot(multi.deriv.cov.jpn, aes(x=true.date, y=cases)) +
  aes(group=NA)+
  geom_path(aes(col = as.factor(deriv.direction)))+
  geom_point(data= total.multijpnplot.dat[which(total.multijpnplot.dat$threshold.crossed==1 & total.multijpnplot.dat$metric.code == "acf + SD + skew")], 
             aes(y =-1000, alpha = ""),size = 3,pch= "|",col = "#53B400")+
  ylab("COVID-19 Cases") + 
  xlab("Date")+
  geom_vline(xintercept = total.multijpnplot.dat$true.date[total.multijpnplot.dat$period != lag(total.multijpnplot.dat$period)], 
             linetype="dashed",  color = "black", size=1)+
  scale_colour_manual(values=c("#22B4F5","black","#F07589"),name = "Predicted Nonlinear\nCase Trend", labels = c("Negative","No trend","Positive"))+
  scale_y_continuous(labels = scales::number_format(accuracy = 1000)) + 
  scale_x_date(date_breaks = "3 months", date_labels = "%b-%Y") +
  ggthemes::theme_clean() + ggtitle("Japan Daily COVID Cases: Four Waves")+  
  annotate("label", x = as.Date("2020-04-01"), y =max(total.multijpnplot.dat$cases)*0.8 , label = "EWS indicator: acf + SD + skewness")+
  scale_alpha_manual(values = 0.8, labels = "EWS detected") +
  labs(alpha = "EWS")  +
  theme(plot.margin = margin(c(10, 8, 0, 10)),
        legend.key.height = unit(0.3,"cm" ),
        legend.key.width = unit(0.1,"cm"),
        legend.title = element_text(size = 10),
        legend.text = element_text(size=10))

png(file="Results/Japan/multigam_jpn_fig.png",
    width=3000, height = 2500 ,res=300)
egg::ggarrange(multi.jpn.p1,multi.jpn.p2,nrow = 2,heights = c(2, 2), labels = c("a","b"))
dev.off()

jpn.multiews.diff2 <- extract.ews.difference(total.multijpnplot.dat,multi.deriv.cov.jpn,"Japan",2)

###########################################################################################
### World - Mexico ###
###########################################################################################
#https://covid19.who.int/info/
cov.mex.dat <- cov.WHO.dat %>%
  filter(Country == "Mexico")%>%
  slice(67:n())%>%
  mutate(cases = ifelse(true.date == "2020-10-10",3468,cases))


ggplot(cov.mex.dat,aes(x=Date, y = cases)) + geom_path() + 
  scale_x_continuous(breaks=seq(2020,2022,0.2)) + 
  scale_y_continuous( labels = scales::number_format(accuracy = 100))+
  ggtitle("Daily Mexico COVID Cases")+ 
  ylab("Cases") + theme_bw()

#### Expanding GAM
exp.mex <- expanding.gam(cov.mex.dat,sensitivity = 7,train.length = 10)
plot(cov.mex.dat$cases)
abline(v=exp.mex$end.index)

mexgam1 <-mgcv::gam(abs(cases)~ s(Weekday,bs = "cs",k=7)+s(Date, bs="tp", k=20), 
                    data = cov.mex.dat[exp.mex$start.index[1]:exp.mex$end.index[1],],
                    family = gaussian(), method = "REML")
plot.gam(mexgam1,pages=1)
mex.deriv1 <- data.frame(gratia::derivatives(mexgam1, term = "s(Date)", interval = "confidence",n= dim(cov.mex.dat[exp.mex$start.index[1]:exp.mex$end.index[1],])[1]),
                         "deriv.direction" = curve.direction(0, gratia::derivatives(mexgam1, term = "s(Date)", interval = "confidence",n= dim(cov.mex.dat[exp.mex$start.index[1]:exp.mex$end.index[1],])[1])$lower, gratia::derivatives(mexgam1, term = "s(Date)", interval = "confidence",n= dim(cov.mex.dat[exp.mex$start.index[1]:exp.mex$end.index[1],])[1])$upper)/2)
plot(mex.deriv1$deriv.direction)
mexgam2 <-mgcv::gam(abs(cases)~ s(Weekday,bs = "cs",k=7) + s(Date, bs="tp", k=20) , 
                    data =  cov.mex.dat[exp.mex$start.index[2]:dim(cov.mex.dat)[1],],
                    family = gaussian(), method = "REML")
plot.gam(mexgam2,pages=1)
mex.deriv2 <- data.frame(gratia::derivatives(mexgam2, term = "s(Date)", interval = "confidence",n= dim(cov.mex.dat[exp.mex$start.index[2]:dim(cov.mex.dat)[1],])[1]),
                         "deriv.direction" = curve.direction(0, gratia::derivatives(mexgam2, term = "s(Date)", interval = "confidence",n= dim(cov.mex.dat[exp.mex$start.index[2]:dim(cov.mex.dat)[1],])[1])$lower, gratia::derivatives(mexgam2, term = "s(Date)", interval = "confidence",n= dim(cov.mex.dat[exp.mex$start.index[2]:dim(cov.mex.dat)[1],])[1])$upper)/2)
plot(mex.deriv2$deriv.direction)

cutoff.mex.cov.ews.multigam <- pbmclapply(cov.mex.dat$Date[1:(length(cov.mex.dat$Date)-13)],FUN = function(x){
  
  data.cut <- cov.mex.dat[as.numeric(cov.mex.dat$Date) >=  as.numeric(x),]
  
  tmp <- comp_EWS_wrapper(data.frame(timedat = as.numeric(data.cut$Date),
                                     biomass = data.cut$cases),
                          metrics = metrics,  threshold = 2, burn_in =7,
                          plotIt = F, ggplotIt = F, tail.direction = "one.tailed",
                          interpolate = F, method = "w_comp")
  tmp$cutoff <- paste(x)
  return(tmp)
  
}, mc.cores =3 )
names(cutoff.mex.cov.ews.multigam) <- cov.mex.dat$true.date[1:(length(cov.mex.dat$true.date)-13)]
temp <- cutoff.mex.cov.ews.multigam
cutoff.mex.cov.ews.multigam <- temp
cutoff.mex.cov.ews.multigam <- data.table::rbindlist(cutoff.mex.cov.ews.multigam, idcol = "start.date")%>%
  mutate(period = ifelse(time <= exp.mex$end.date[1], "first", #define periods of constancy prior to waves
                         ifelse(time >= exp.mex$start.date[2] ,"second","third")))
save(cutoff.mex.cov.ews.multigam,file = "Results/Mexico/mex.multi.gam.ews.RData")
load("/Results/Mexico/mex.multi.gam.ews.RData")

mex.multifirst.plot.dat <- cutoff.mex.cov.ews.multigam[cutoff.mex.cov.ews.multigam$start.date == "2020-03-09" & cutoff.mex.cov.ews.multigam$period == "first",]%>%
  left_join(data.frame("true.date" = cov.mex.dat$true.date,"time" = cov.mex.dat$Date), by = c("time"))

mex.multisecond.plot.dat <- cutoff.mex.cov.ews.multigam[cutoff.mex.cov.ews.multigam$start.date == cov.mex.dat$true.date[exp.mex$start.index[2]] & cutoff.mex.cov.ews.multigam$period == "second",]%>%
  left_join(data.frame("true.date" = cov.mex.dat$true.date,"time" = cov.mex.dat$Date), by = c("time"))

total.multimexplot.dat <- rbind(mex.multifirst.plot.dat,mex.multisecond.plot.dat) %>%
  rename(cases = count.used)

multi.deriv.cov.mex <- cov.mex.dat[,2:4]%>%
  cbind(rbind(mex.deriv1,mex.deriv2))


multi.mex.p2 <- ggplot(data = total.multimexplot.dat, aes(x=true.date,y=str)) +
  geom_hline(yintercept = 2, linetype="solid", color = "grey", size=1)+
  geom_line(aes(col= metric.code,alpha= "Indicator strength"))+
  geom_point(data =total.multimexplot.dat[which(total.multimexplot.dat$threshold.crossed==1)], aes(x=true.date, y = str, col = metric.code,alpha= "EWS detected")) +
  scale_alpha_manual(values = c(1, 1),
                     breaks = c("Indicator strength", "EWS detected"),
                     guide = guide_legend(override.aes = list(linetype = c(1, 0),shape = c(NA, 16)))) +
  #scale_colour_manual(values = scales::hue_pal()(7)) +
  ggthemes::theme_clean() + xlab("Date") + ylab("Strength of EWS") +
  geom_vline(xintercept = total.multimexplot.dat$true.date[total.multimexplot.dat$period != lag(total.multimexplot.dat$period)], linetype="dashed", 
             color = "black", size=1)+
  scale_x_date(date_breaks = "3 months", date_labels = "%b-%Y")+
  labs(color='EWS Indicator',alpha ="EWS") +
  theme(plot.margin = margin(c(10, 8, 0, 10)),
        legend.key.height = unit(0.3,"cm"),
        legend.key.width = unit(0.1,"cm"),
        legend.title = element_text(size = 10),
        legend.text = element_text(size=10))


multi.mex.p1 <- ggplot(multi.deriv.cov.mex, aes(x=true.date, y=cases)) +
  aes(group=NA)+
  geom_path(aes(col = as.factor(deriv.direction)))+
  geom_point(data= total.multimexplot.dat[which(total.multimexplot.dat$threshold.crossed==1 & total.multimexplot.dat$metric.code == "acf + SD + skew")], 
             aes(y =-5000, alpha = ""),size = 3,pch= "|",col = "#53B400")+
  ylab("COVID-19 Cases") + 
  xlab("Date")+
  geom_vline(xintercept = total.multimexplot.dat$true.date[total.multimexplot.dat$period != lag(total.multimexplot.dat$period)], 
             linetype="dashed",  color = "black", size=1)+
  scale_colour_manual(values=c("#22B4F5","black","#F07589"),name = "Predicted Nonlinear\nCase Trend", labels = c("Negative","No trend","Positive"))+
  scale_y_continuous(labels = scales::number_format(accuracy = 100)) + 
  scale_x_date(date_breaks = "3 months", date_labels = "%b-%Y") +
  ggthemes::theme_clean() + ggtitle("Mexico Daily COVID Cases: Two Waves")+  
  annotate("label", x = as.Date("2020-06-01"), y =max(total.multimexplot.dat$cases)*0.8 , label = "EWS indicator: acf + SD + skewness")+
  scale_alpha_manual(values = 0.8, labels = "EWS detected") +
  labs(alpha = "EWS")  +
  theme(plot.margin = margin(c(10, 8, 0, 10)),
        legend.key.height = unit(0.3,"cm" ),
        legend.key.width = unit(0.1,"cm"),
        legend.title = element_text(size = 10),
        legend.text = element_text(size=10))

png(file="Results/Mexico/multigam_mex_fig.png",
    width=3000, height = 2500 ,res=300)
egg::ggarrange(multi.mex.p1,multi.mex.p2,nrow = 2,heights = c(2, 2), labels = c("a","b"))
dev.off()

mex.multiews.diff2 <- extract.ews.difference(total.multimexplot.dat,multi.deriv.cov.mex,"Mexico",2)

###########################################################################################
### World - Netherlands ###
###########################################################################################
#https://covid19.who.int/info/
cov.ned.dat <- cov.WHO.dat %>%
  filter(Country == "Netherlands")%>%
  slice(51:n())

ggplot(cov.ned.dat,aes(x=Date, y = cases)) + geom_path() + 
  scale_x_continuous(breaks=seq(2020,2022,0.2)) + 
  scale_y_continuous( labels = scales::number_format(accuracy = 100))+
  ggtitle("Daily Netherlands COVID Cases")+ 
  ylab("Cases") + theme_bw()

#### Expanding GAM
exp.ned <- expanding.gam(cov.ned.dat,sensitivity = 7,train.length = 10)
plot(cov.ned.dat$cases)
abline(v=exp.ned$end.index)

nedgam1 <-mgcv::gam(abs(cases)~ s(Weekday,bs = "cs",k=7)+s(Date, bs="tp", k=30), 
                    data = cov.ned.dat[exp.ned$start.index[1]:exp.ned$end.index[1],],
                    family = gaussian(), method = "REML")
plot.gam(nedgam1,pages=1)
ned.deriv1 <- data.frame(gratia::derivatives(nedgam1, term = "s(Date)", interval = "confidence",n= dim(cov.ned.dat[exp.ned$start.index[1]:exp.ned$end.index[1],])[1]),
                         "deriv.direction" = curve.direction(0, gratia::derivatives(nedgam1, term = "s(Date)", interval = "confidence",n= dim(cov.ned.dat[exp.ned$start.index[1]:exp.ned$end.index[1],])[1])$lower, gratia::derivatives(nedgam1, term = "s(Date)", interval = "confidence",n= dim(cov.ned.dat[exp.ned$start.index[1]:exp.ned$end.index[1],])[1])$upper)/2)
plot(ned.deriv1$deriv.direction)

nedgam2 <-mgcv::gam(abs(cases)~ s(Weekday,bs = "cs",k=7) + s(Date, bs="tp", k=20) , 
                    data =  cov.ned.dat[exp.ned$start.index[2]:dim(cov.ned.dat)[1],],
                    family = gaussian(), method = "REML")
plot.gam(nedgam2,pages=1)
ned.deriv2 <- data.frame(gratia::derivatives(nedgam2, term = "s(Date)", interval = "confidence",n= dim(cov.ned.dat[exp.ned$start.index[2]:dim(cov.ned.dat)[1],])[1]),
                         "deriv.direction" = curve.direction(0, gratia::derivatives(nedgam2, term = "s(Date)", interval = "confidence",n= dim(cov.ned.dat[exp.ned$start.index[2]:dim(cov.ned.dat)[1],])[1])$lower, gratia::derivatives(nedgam2, term = "s(Date)", interval = "confidence",n= dim(cov.ned.dat[exp.ned$start.index[2]:dim(cov.ned.dat)[1],])[1])$upper)/2)
plot(ned.deriv2$deriv.direction)

cutoff.ned.cov.ews.multigam <- pbmclapply(cov.ned.dat$Date[1:(length(cov.ned.dat$Date)-13)],FUN = function(x){
  
  data.cut <- cov.ned.dat[as.numeric(cov.ned.dat$Date) >=  as.numeric(x),]
  
  tmp <- comp_EWS_wrapper(data.frame(timedat = as.numeric(data.cut$Date),
                                     biomass = data.cut$cases),
                          metrics = metrics,  threshold = 2, burn_in =7,
                          plotIt = F, ggplotIt = F, tail.direction = "one.tailed",
                          interpolate = F, method = "w_comp")
  tmp$cutoff <- paste(x)
  return(tmp)
  
}, mc.cores =3 )
names(cutoff.ned.cov.ews.multigam) <- cov.ned.dat$true.date[1:(length(cov.ned.dat$true.date)-13)]
temp <- cutoff.ned.cov.ews.multigam
cutoff.ned.cov.ews.multigam <- data.table::rbindlist(cutoff.ned.cov.ews.multigam, idcol = "start.date")%>%
  mutate(period = ifelse(time <= exp.ned$end.date[1], "first", #define periods of constancy prior to waves
                         ifelse(time >= exp.ned$start.date[2] , "second","third")))
save(cutoff.ned.cov.ews.multigam,file = "Results/Netherlands/ned.multi.gam.ews.RData")
load("Results/Netherlands/ned.multi.gam.ews.RData")

ned.multifirst.plot.dat <- cutoff.ned.cov.ews.multigam[cutoff.ned.cov.ews.multigam$start.date == "2020-02-22" & cutoff.ned.cov.ews.multigam$period == "first",]%>%
  left_join(data.frame("true.date" = cov.ned.dat$true.date,"time" = cov.ned.dat$Date), by = c("time"))

ned.multisecond.plot.dat <- cutoff.ned.cov.ews.multigam[cutoff.ned.cov.ews.multigam$start.date == cov.ned.dat$true.date[exp.ned$start.index[2]] & cutoff.ned.cov.ews.multigam$period == "second",]%>%
  left_join(data.frame("true.date" = cov.ned.dat$true.date,"time" = cov.ned.dat$Date), by = c("time"))

total.multinedplot.dat <- rbind(ned.multifirst.plot.dat,ned.multisecond.plot.dat) %>%
  rename(cases = count.used)

multi.deriv.cov.ned <- cov.ned.dat[,2:4]%>%
  cbind(rbind(ned.deriv1,ned.deriv2))


multi.ned.p2 <- ggplot(data = total.multinedplot.dat, aes(x=true.date,y=str)) +
  geom_hline(yintercept = 2, linetype="solid", color = "grey", size=1)+
  geom_line(aes(col= metric.code,alpha= "Indicator strength"))+
  geom_point(data =total.multinedplot.dat[which(total.multinedplot.dat$threshold.crossed==1)], aes(x=true.date, y = str, col = metric.code,alpha= "EWS detected")) +
  scale_alpha_manual(values = c(1, 1),
                     breaks = c("Indicator strength", "EWS detected"),
                     guide = guide_legend(override.aes = list(linetype = c(1, 0),shape = c(NA, 16)))) +
  #scale_colour_manual(values = scales::hue_pal()(7)) +
  ggthemes::theme_clean() + xlab("Date") + ylab("Strength of EWS") +
  geom_vline(xintercept = total.multinedplot.dat$true.date[total.multinedplot.dat$period != lag(total.multinedplot.dat$period)], linetype="dashed", 
             color = "black", size=1)+
  scale_x_date(date_breaks = "3 months", date_labels = "%b-%Y")+
  labs(color='EWS Indicator',alpha ="EWS") +
  theme(plot.margin = margin(c(10, 8, 0, 10)),
        legend.key.height = unit(0.3,"cm"),
        legend.key.width = unit(0.1,"cm"),
        legend.title = element_text(size = 10),
        legend.text = element_text(size=10))


multi.ned.p1 <- ggplot(multi.deriv.cov.ned, aes(x=true.date, y=cases)) +
  aes(group=NA)+
  geom_path(aes(col = as.factor(deriv.direction)))+
  geom_point(data= total.multinedplot.dat[which(total.multinedplot.dat$threshold.crossed==1 & total.multinedplot.dat$metric.code == "acf + SD + skew")], 
             aes(y =-1000, alpha = ""),size = 3,pch= "|",col = "#53B400")+
  ylab("COVID-19 Cases") + 
  xlab("Date")+
  geom_vline(xintercept = total.multinedplot.dat$true.date[total.multinedplot.dat$period != lag(total.multinedplot.dat$period)], 
             linetype="dashed",  color = "black", size=1)+
  scale_colour_manual(values=c("#22B4F5","black","#F07589"),name = "Predicted Nonlinear\nCase Trend", labels = c("Negative","No trend","Positive"))+
  scale_y_continuous(labels = scales::number_format(accuracy = 1000)) + 
  scale_x_date(date_breaks = "3 months", date_labels = "%b-%Y") +
  ggthemes::theme_clean() + ggtitle("Netherlands Daily COVID Cases: Two Waves")+  
  annotate("label", x = as.Date("2020-06-01"), y =max(total.multinedplot.dat$cases)*0.8 , label = "EWS indicator: acf + SD + skewness")+
  scale_alpha_manual(values = 0.8, labels = "EWS detected") +
  labs(alpha = "EWS")  +
  theme(plot.margin = margin(c(10, 8, 0, 10)),
        legend.key.height = unit(0.3,"cm" ),
        legend.key.width = unit(0.1,"cm"),
        legend.title = element_text(size = 10),
        legend.text = element_text(size=10))

png(file="Results/Netherlands/multigam_ned_fig.png",
    width=3000, height = 2500 ,res=300)
egg::ggarrange(multi.ned.p1,multi.ned.p2,nrow = 2,heights = c(2, 2), labels = c("a","b"))
dev.off()

ned.multiews.diff2 <- extract.ews.difference(total.multinedplot.dat,multi.deriv.cov.ned,"Netherlands",2)

###########################################################################################
### World - Pakistan ###
###########################################################################################
#https://covid19.who.int/info/
cov.pk.dat <- cov.WHO.dat %>%
  filter(Country == "Pakistan")%>%
  slice(72:n())

ggplot(cov.pk.dat,aes(x=Date, y = cases)) + geom_path() + 
  scale_x_continuous(breaks=seq(2020,2022,0.2)) + 
  scale_y_continuous( lapks = scales::number_format(accuracy = 100))+
  ggtitle("Daily Pakistan COVID Cases")+ 
  ylab("Cases") + theme_bw()

#### Expanding GAM
exp.pk <- expanding.gam(cov.pk.dat,sensitivity = 7,train.length = 10)
plot(cov.pk.dat$cases)
abline(v=exp.pk$end.index)

pkgam1 <-mgcv::gam(abs(cases)~ s(Weekday,bs = "cs",k=7)+s(Date, bs="tp", k=20), 
                   data = cov.pk.dat[exp.pk$start.index[1]:exp.pk$end.index[2],],
                   family = gaussian(), method = "REML")
plot.gam(pkgam1,pages=1)
pk.deriv1 <- data.frame(gratia::derivatives(pkgam1, term = "s(Date)", interval = "confidence",n= dim(cov.pk.dat[exp.pk$start.index[1]:exp.pk$end.index[2],])[1]),
                        "deriv.direction" = curve.direction(0, gratia::derivatives(pkgam1, term = "s(Date)", interval = "confidence",n= dim(cov.pk.dat[exp.pk$start.index[1]:exp.pk$end.index[2],])[1])$lower, gratia::derivatives(pkgam1, term = "s(Date)", interval = "confidence",n= dim(cov.pk.dat[exp.pk$start.index[1]:exp.pk$end.index[2],])[1])$upper)/2)
plot(pk.deriv1$deriv.direction)

pkgam2 <-mgcv::gam(abs(cases)~ s(Weekday,bs = "cs",k=7) + s(Date, bs="tp", k=20) , 
                   data =  cov.pk.dat[exp.pk$start.index[3]:exp.pk$end.index[4],],
                   family = gaussian(), method = "REML")
plot.gam(pkgam2,pages=1)
pk.deriv2 <- data.frame(gratia::derivatives(pkgam2, term = "s(Date)", interval = "confidence",n= dim(cov.pk.dat[exp.pk$start.index[3]:exp.pk$end.index[4],])[1]),
                        "deriv.direction" = curve.direction(0, gratia::derivatives(pkgam2, term = "s(Date)", interval = "confidence",n= dim(cov.pk.dat[exp.pk$start.index[3]:exp.pk$end.index[4],])[1])$lower, gratia::derivatives(pkgam2, term = "s(Date)", interval = "confidence",n= dim(cov.pk.dat[exp.pk$start.index[3]:exp.pk$end.index[4],])[1])$upper)/2)
plot(pk.deriv2$deriv.direction)

pkgam3 <-mgcv::gam(abs(cases)~ s(Weekday,bs = "cs",k=7) + s(Date, bs="tp", k=15) , 
                   data =  cov.pk.dat[(exp.pk$end.index[4]+1):dim(cov.pk.dat)[1],],
                   family = gaussian(), method = "REML")
plot.gam(pkgam3,pages=1)
pk.deriv3 <- data.frame(gratia::derivatives(pkgam3, term = "s(Date)", interval = "confidence",n= dim(cov.pk.dat[(exp.pk$end.index[4]+1):dim(cov.pk.dat)[1],])[1]),
                        "deriv.direction" = curve.direction(0, gratia::derivatives(pkgam3, term = "s(Date)", interval = "confidence",n= dim(cov.pk.dat[(exp.pk$end.index[4]+1):dim(cov.pk.dat)[1],])[1])$lower, gratia::derivatives(pkgam3, term = "s(Date)", interval = "confidence",n= dim(cov.pk.dat[(exp.pk$end.index[4]+1):dim(cov.pk.dat)[1],])[1])$upper)/2)
plot(pk.deriv3$deriv.direction)

cutoff.pk.cov.ews.multigam <- pbmclapply(cov.pk.dat$Date[1:(length(cov.pk.dat$Date)-13)],FUN = function(x){
  
  data.cut <- cov.pk.dat[as.numeric(cov.pk.dat$Date) >=  as.numeric(x),]
  
  tmp <- comp_EWS_wrapper(data.frame(timedat = as.numeric(data.cut$Date),
                                     biomass = data.cut$cases),
                          metrics = metrics,  threshold = 2, burn_in =7,
                          plotIt = F, ggplotIt = F, tail.direction = "one.tailed",
                          interpolate = F, method = "w_comp")
  tmp$cutoff <- paste(x)
  return(tmp)
  
}, mc.cores =3 )
names(cutoff.pk.cov.ews.multigam) <- cov.pk.dat$true.date[1:(length(cov.pk.dat$true.date)-13)]
temp <- cutoff.pk.cov.ews.multigam
cutoff.pk.cov.ews.multigam <- temp
cutoff.pk.cov.ews.multigam <- data.table::rbindlist(cutoff.pk.cov.ews.multigam, idcol = "start.date")%>%
  mutate(period = ifelse(time <= exp.pk$end.date[2], "first", #define periods of constancy prior to waves
                         ifelse(time >= exp.pk$start.date[3] & time <= exp.pk$end.date[4], "second",
                                ifelse(time > exp.pk$end.date[4], "third","fourth"))))
save(cutoff.pk.cov.ews.multigam,file = "Results/Pakistan/pk.multi.gam.ews.RData")
load("Results/Pakistan/pk.multi.gam.ews.RData")

pk.multifirst.plot.dat <- cutoff.pk.cov.ews.multigam[cutoff.pk.cov.ews.multigam$start.date == "2020-03-14" & cutoff.pk.cov.ews.multigam$period == "first",]%>%
  left_join(data.frame("true.date" = cov.pk.dat$true.date,"time" = cov.pk.dat$Date), by = c("time"))

pk.multisecond.plot.dat <- cutoff.pk.cov.ews.multigam[cutoff.pk.cov.ews.multigam$start.date == cov.pk.dat$true.date[exp.pk$start.index[3]] & cutoff.pk.cov.ews.multigam$period == "second",]%>%
  left_join(data.frame("true.date" = cov.pk.dat$true.date,"time" = cov.pk.dat$Date), by = c("time"))

pk.multithird.plot.dat <- cutoff.pk.cov.ews.multigam[cutoff.pk.cov.ews.multigam$start.date == cov.pk.dat$true.date[exp.pk$end.index[4]+1] & cutoff.pk.cov.ews.multigam$period == "third",]%>%
  left_join(data.frame("true.date" = cov.pk.dat$true.date,"time" = cov.pk.dat$Date), by = c("time"))

total.multipkplot.dat <- rbind(pk.multifirst.plot.dat,pk.multisecond.plot.dat,pk.multithird.plot.dat) %>%
  rename(cases = count.used)

multi.deriv.cov.pk <- cov.pk.dat[,2:4]%>%
  cbind(rbind(pk.deriv1,pk.deriv2,pk.deriv3))


multi.pk.p2 <- ggplot(data = total.multipkplot.dat, aes(x=true.date,y=str)) +
  geom_hline(yintercept = 2, linetype="solid", color = "grey", size=1)+
  geom_line(aes(col= metric.code,alpha= "Indicator strength"))+
  geom_point(data =total.multipkplot.dat[which(total.multipkplot.dat$threshold.crossed==1)], aes(x=true.date, y = str, col = metric.code,alpha= "EWS detected")) +
  scale_alpha_manual(values = c(1, 1),
                     breaks = c("Indicator strength", "EWS detected"),
                     guide = guide_legend(override.aes = list(linetype = c(1, 0),shape = c(NA, 16)))) +
  #scale_colour_manual(values = scales::hue_pal()(7)) +
  ggthemes::theme_clean() + xlab("Date") + ylab("Strength of EWS") +
  geom_vline(xintercept = total.multipkplot.dat$true.date[total.multipkplot.dat$period != lag(total.multipkplot.dat$period)], linetype="dashed", 
             color = "black", size=1)+
  scale_x_date(date_breaks = "3 months", date_labels = "%b-%Y")+
  labs(color='EWS Indicator',alpha ="EWS") +
  theme(plot.margin = margin(c(10, 8, 0, 10)),
        legend.key.height = unit(0.3,"cm"),
        legend.key.width = unit(0.1,"cm"),
        legend.title = element_text(size = 10),
        legend.text = element_text(size=10))


multi.pk.p1 <- ggplot(multi.deriv.cov.pk, aes(x=true.date, y=cases)) +
  aes(group=NA)+
  geom_path(aes(col = as.factor(deriv.direction)))+
  geom_point(data= total.multipkplot.dat[which(total.multipkplot.dat$threshold.crossed==1 & total.multipkplot.dat$metric.code == "acf + SD + skew")], 
             aes(y =-1000, alpha = ""),size = 3,pch= "|",col = "#53B400")+
  ylab("COVID-19 Cases") + 
  xlab("Date")+
  geom_vline(xintercept = total.multipkplot.dat$true.date[total.multipkplot.dat$period != lag(total.multipkplot.dat$period)], 
             linetype="dashed",  color = "black", size=1)+
  scale_colour_manual(values=c("#22B4F5","black","#F07589"),name = "Predicted Nonlinear\nCase Trend", labels = c("Negative","No trend","Positive"))+
  scale_y_continuous(labels = scales::number_format(accuracy = 1000)) + 
  scale_x_date(date_breaks = "3 months", date_labels = "%b-%Y") +
  ggthemes::theme_clean() + ggtitle("Pakistan Daily COVID Cases: Three Waves")+  
  annotate("label", x = as.Date("2020-06-01"), y =max(total.multipkplot.dat$cases)*0.8 , label = "EWS indicator: acf + SD + skewness")+
  scale_alpha_manual(values = 0.8, labels = "EWS detected") +
  labs(alpha = "EWS")  +
  theme(plot.margin = margin(c(10, 8, 0, 10)),
        legend.key.height = unit(0.3,"cm" ),
        legend.key.width = unit(0.1,"cm"),
        legend.title = element_text(size = 10),
        legend.text = element_text(size=10))

png(file="Results/Pakistan/multigam_pk_fig.png",
    width=3000, height = 2500 ,res=300)
egg::ggarrange(multi.pk.p1,multi.pk.p2,nrow = 2,heights = c(2, 2), labels = c("a","b"))
dev.off()

pk.multiews.diff2 <- extract.ews.difference(total.multipkplot.dat,multi.deriv.cov.pk,"Pakistan",2)

###########################################################################################
### World - Peru ###
###########################################################################################
#https://covid19.who.int/info/
cov.per.dat <- cov.WHO.dat %>%
  filter(Country == "Peru")%>%
  slice(67:n())


ggplot(cov.per.dat,aes(x=Date, y = cases)) + geom_path() + 
  scale_x_continuous(breaks=seq(2020,2022,0.2)) + 
  scale_y_continuous( lapers = scales::number_format(accuracy = 100))+
  ggtitle("Daily Peru COVID Cases")+ 
  ylab("Cases") + theme_bw()

#### Expanding GAM
exp.per <- expanding.gam(cov.per.dat,sensitivity = 7,train.length = 10)
plot(cov.per.dat$cases)
abline(v=exp.per$end.index)

pergam1 <-mgcv::gam(abs(cases)~ s(Weekday,bs = "cs",k=7)+s(Date, bs="tp", k=30), 
                    data = cov.per.dat[exp.per$start.index[1]:exp.per$end.index[4],],
                    family = gaussian(), method = "REML")
plot.gam(pergam1,pages=1)
per.deriv1 <- data.frame(gratia::derivatives(pergam1, term = "s(Date)", interval = "confidence",n= dim(cov.per.dat[exp.per$start.index[1]:exp.per$end.index[4],])[1]),
                         "deriv.direction" = curve.direction(0, gratia::derivatives(pergam1, term = "s(Date)", interval = "confidence",n= dim(cov.per.dat[exp.per$start.index[1]:exp.per$end.index[4],])[1])$lower, gratia::derivatives(pergam1, term = "s(Date)", interval = "confidence",n= dim(cov.per.dat[exp.per$start.index[1]:exp.per$end.index[4],])[1])$upper)/2)
plot(per.deriv1$deriv.direction)

pergam2 <-mgcv::gam(abs(cases)~ s(Weekday,bs = "cs",k=7) + s(Date, bs="tp", k=15) , 
                    data =  cov.per.dat[exp.per$start.index[5]:dim(cov.per.dat)[1],],
                    family = gaussian(), method = "REML")
plot.gam(pergam2,pages=1)
per.deriv2 <- data.frame(gratia::derivatives(pergam2, term = "s(Date)", interval = "confidence",n= dim(cov.per.dat[exp.per$start.index[5]:dim(cov.per.dat)[1],])[1]),
                         "deriv.direction" = curve.direction(0, gratia::derivatives(pergam2, term = "s(Date)", interval = "confidence",n= dim(cov.per.dat[exp.per$start.index[5]:dim(cov.per.dat)[1],])[1])$lower, gratia::derivatives(pergam2, term = "s(Date)", interval = "confidence",n= dim(cov.per.dat[exp.per$start.index[5]:dim(cov.per.dat)[1],])[1])$upper)/2)
plot(per.deriv2$deriv.direction)

cutoff.per.cov.ews.multigam <- pbmclapply(cov.per.dat$Date[1:(length(cov.per.dat$Date)-13)],FUN = function(x){
  
  data.cut <- cov.per.dat[as.numeric(cov.per.dat$Date) >=  as.numeric(x),]
  
  tmp <- comp_EWS_wrapper(data.frame(timedat = as.numeric(data.cut$Date),
                                     biomass = data.cut$cases),
                          metrics = metrics,  threshold = 2, burn_in =7,
                          plotIt = F, ggplotIt = F, tail.direction = "one.tailed",
                          interpolate = F, method = "w_comp")
  tmp$cutoff <- paste(x)
  return(tmp)
  
}, mc.cores =3 )
names(cutoff.per.cov.ews.multigam) <- cov.per.dat$true.date[1:(length(cov.per.dat$true.date)-13)]
temp <- cutoff.per.cov.ews.multigam
cutoff.per.cov.ews.multigam <- data.table::rbindlist(cutoff.per.cov.ews.multigam, idcol = "start.date")%>%
  mutate(period = ifelse(time <= exp.per$end.date[4], "first", #define periods of constancy prior to waves
                         ifelse(time >= exp.per$start.date[5] , "second", "third")))
save(cutoff.per.cov.ews.multigam,file = "Results/Peru/per.multi.gam.ews.RData")
load("Results/Peru/per.multi.gam.ews.RData")

per.multifirst.plot.dat <- cutoff.per.cov.ews.multigam[cutoff.per.cov.ews.multigam$start.date == "2020-03-09" & cutoff.per.cov.ews.multigam$period == "first",]%>%
  left_join(data.frame("true.date" = cov.per.dat$true.date,"time" = cov.per.dat$Date), by = c("time"))

per.multisecond.plot.dat <- cutoff.per.cov.ews.multigam[cutoff.per.cov.ews.multigam$start.date == cov.per.dat$true.date[exp.per$start.index[5]] & cutoff.per.cov.ews.multigam$period == "second",]%>%
  left_join(data.frame("true.date" = cov.per.dat$true.date,"time" = cov.per.dat$Date), by = c("time"))

total.multiperplot.dat <- rbind(per.multifirst.plot.dat,per.multisecond.plot.dat) %>%
  rename(cases = count.used)

multi.deriv.cov.per <- cov.per.dat[,2:4]%>%
  cbind(rbind(per.deriv1,per.deriv2))


multi.per.p2 <- ggplot(data = total.multiperplot.dat, aes(x=true.date,y=str)) +
  geom_hline(yintercept = 2, linetype="solid", color = "grey", size=1)+
  geom_line(aes(col= metric.code,alpha= "Indicator strength"))+
  geom_point(data =total.multiperplot.dat[which(total.multiperplot.dat$threshold.crossed==1)], aes(x=true.date, y = str, col = metric.code,alpha= "EWS detected")) +
  scale_alpha_manual(values = c(1, 1),
                     breaks = c("Indicator strength", "EWS detected"),
                     guide = guide_legend(override.aes = list(linetype = c(1, 0),shape = c(NA, 16)))) +
  #scale_colour_manual(values = scales::hue_pal()(7)) +
  ggthemes::theme_clean() + xlab("Date") + ylab("Strength of EWS") +
  geom_vline(xintercept = total.multiperplot.dat$true.date[total.multiperplot.dat$period != lag(total.multiperplot.dat$period)], linetype="dashed", 
             color = "black", size=1)+
  scale_x_date(date_breaks = "3 months", date_labels = "%b-%Y")+
  labs(color='EWS Indicator',alpha ="EWS") +
  theme(plot.margin = margin(c(10, 8, 0, 10)),
        legend.key.height = unit(0.3,"cm"),
        legend.key.width = unit(0.1,"cm"),
        legend.title = element_text(size = 10),
        legend.text = element_text(size=10))


multi.per.p1 <- ggplot(multi.deriv.cov.per, aes(x=true.date, y=cases)) +
  aes(group=NA)+
  geom_path(aes(col = as.factor(deriv.direction)))+
  geom_point(data= total.multiperplot.dat[which(total.multiperplot.dat$threshold.crossed==1 & total.multiperplot.dat$metric.code == "acf + SD + skew")], 
             aes(y =-1000, alpha = ""),size = 3,pch= "|",col = "#53B400")+
  ylab("COVID-19 Cases") + 
  xlab("Date")+
  geom_vline(xintercept = total.multiperplot.dat$true.date[total.multiperplot.dat$period != lag(total.multiperplot.dat$period)], 
             linetype="dashed",  color = "black", size=1)+
  scale_colour_manual(values=c("#22B4F5","black","#F07589"),name = "Predicted Nonlinear\nCase Trend", labels = c("Negative","No trend","Positive"))+
  scale_y_continuous(labels = scales::number_format(accuracy = 1000)) + 
  scale_x_date(date_breaks = "3 months", date_labels= "%b-%Y") +
  ggthemes::theme_clean() + ggtitle("Peru Daily COVID Cases: Two Waves")+  
  annotate("label", x = as.Date("2020-06-01"), y =max(total.multiperplot.dat$cases)*0.8 , label = "EWS indicator: acf + SD + skewness")+
  scale_alpha_manual(values = 0.8, labels = "EWS detected") +
  labs(alpha = "EWS")  +
  theme(plot.margin = margin(c(10, 8, 0, 10)),
        legend.key.height = unit(0.3,"cm" ),
        legend.key.width = unit(0.1,"cm"),
        legend.title = element_text(size = 10),
        legend.text = element_text(size=10))

png(file="Results/Peru/multigam_per_fig.png",
    width=3000, height = 2500 ,res=300)
egg::ggarrange(multi.per.p1,multi.per.p2,nrow = 2,heights = c(2, 2), labels = c("a","b"))
dev.off()

per.multiews.diff2 <- extract.ews.difference(total.multiperplot.dat,multi.deriv.cov.per,"Peru",2)

###########################################################################################
### World - Philippines ###
###########################################################################################
#https://covid19.who.int/info/
cov.phl.dat <- cov.WHO.dat %>%
  filter(Country == "Philippines")%>%
  slice(63:n())

ggplot(cov.phl.dat,aes(x=Date, y = cases)) + geom_path() + 
  scale_x_continuous(breaks=seq(2020,2022,0.2)) + 
  scale_y_continuous( labels = scales::number_format(accuracy = 100))+
  ggtitle("Daily Philippines COVID Cases")+ 
  ylab("Cases") + theme_bw()

#### Expanding GAM
exp.phl <- expanding.gam(cov.phl.dat,sensitivity = 7,train.length = 10)
plot(cov.phl.dat$cases)
abline(v=exp.phl$end.index)

phlgam1 <-mgcv::gam(abs(cases)~ s(Weekday,bs = "cs",k=7)+s(Date, bs="tp", k=13), 
                    data = cov.phl.dat[exp.phl$start.index[1]:exp.phl$end.index[2],],
                    family = gaussian(), method = "REML")
plot.gam(phlgam1,pages=1)
phl.deriv1 <- data.frame(gratia::derivatives(phlgam1, term = "s(Date)", interval = "confidence",n= dim(cov.phl.dat[exp.phl$start.index[1]:exp.phl$end.index[2],])[1]),
                         "deriv.direction" = curve.direction(0, gratia::derivatives(phlgam1, term = "s(Date)", interval = "confidence",n= dim(cov.phl.dat[exp.phl$start.index[1]:exp.phl$end.index[2],])[1])$lower, gratia::derivatives(phlgam1, term = "s(Date)", interval = "confidence",n= dim(cov.phl.dat[exp.phl$start.index[1]:exp.phl$end.index[2],])[1])$upper)/2)
plot(phl.deriv1$deriv.direction)

phlgam2 <-mgcv::gam(abs(cases)~ s(Weekday,bs = "cs",k=7) + s(Date, bs="tp", k=20) , 
                    data =  cov.phl.dat[exp.phl$start.index[3]:dim(cov.phl.dat)[1],],
                    family = gaussian(), method = "REML")
plot.gam(phlgam2,pages=1)
phl.deriv2 <- data.frame(gratia::derivatives(phlgam2, term = "s(Date)", interval = "confidence",n= dim(cov.phl.dat[exp.phl$start.index[3]:dim(cov.phl.dat)[1],])[1]),
                         "deriv.direction" = curve.direction(0, gratia::derivatives(phlgam2, term = "s(Date)", interval = "confidence",n= dim(cov.phl.dat[exp.phl$start.index[3]:dim(cov.phl.dat)[1],])[1])$lower, gratia::derivatives(phlgam2, term = "s(Date)", interval = "confidence",n= dim(cov.phl.dat[exp.phl$start.index[3]:dim(cov.phl.dat)[1],])[1])$upper)/2)
plot(phl.deriv2$deriv.direction)

cutoff.phl.cov.ews.multigam <- pbmclapply(cov.phl.dat$Date[1:(length(cov.phl.dat$Date)-13)],FUN = function(x){
  
  data.cut <- cov.phl.dat[as.numeric(cov.phl.dat$Date) >=  as.numeric(x),]
  
  tmp <- comp_EWS_wrapper(data.frame(timedat = as.numeric(data.cut$Date),
                                     biomass = data.cut$cases),
                          metrics = metrics,  threshold = 2, burn_in =7,
                          plotIt = F, ggplotIt = F, tail.direction = "one.tailed",
                          interpolate = F, method = "w_comp")
  tmp$cutoff <- paste(x)
  return(tmp)
  
}, mc.cores =3 )
names(cutoff.phl.cov.ews.multigam) <- cov.phl.dat$true.date[1:(length(cov.phl.dat$true.date)-13)]
temp <- cutoff.phl.cov.ews.multigam
cutoff.phl.cov.ews.multigam <- data.table::rbindlist(cutoff.phl.cov.ews.multigam, idcol = "start.date")%>%
  mutate(period = ifelse(time <= exp.phl$end.date[2], "first", #define periods of constancy prior to waves
                         ifelse(time >= exp.phl$start.date[3] , "second","third")))
save(cutoff.phl.cov.ews.multigam,file = "Results/Philippines/phl.multi.gam.ews.RData")
load("Results/Philippines/phl.multi.gam.ews.RData")

phl.multifirst.plot.dat <- cutoff.phl.cov.ews.multigam[cutoff.phl.cov.ews.multigam$start.date == "2020-03-05" & cutoff.phl.cov.ews.multigam$period == "first",]%>%
  left_join(data.frame("true.date" = cov.phl.dat$true.date,"time" = cov.phl.dat$Date), by = c("time"))

phl.multisecond.plot.dat <- cutoff.phl.cov.ews.multigam[cutoff.phl.cov.ews.multigam$start.date == cov.phl.dat$true.date[exp.phl$start.index[3]] & cutoff.phl.cov.ews.multigam$period == "second",]%>%
  left_join(data.frame("true.date" = cov.phl.dat$true.date,"time" = cov.phl.dat$Date), by = c("time"))

total.multiphlplot.dat <- rbind(phl.multifirst.plot.dat,phl.multisecond.plot.dat) %>%
  rename(cases = count.used)

multi.deriv.cov.phl <- cov.phl.dat[,2:4]%>%
  cbind(rbind(phl.deriv1,phl.deriv2))


multi.phl.p2 <- ggplot(data = total.multiphlplot.dat, aes(x=true.date,y=str)) +
  geom_hline(yintercept = 2, linetype="solid", color = "grey", size=1)+
  geom_line(aes(col= metric.code,alpha= "Indicator strength"))+
  geom_point(data =total.multiphlplot.dat[which(total.multiphlplot.dat$threshold.crossed==1)], aes(x=true.date, y = str, col = metric.code,alpha= "EWS detected")) +
  scale_alpha_manual(values = c(1, 1),
                     breaks = c("Indicator strength", "EWS detected"),
                     guide = guide_legend(override.aes = list(linetype = c(1, 0),shape = c(NA, 16)))) +
  #scale_colour_manual(values = scales::hue_pal()(7)) +
  ggthemes::theme_clean() + xlab("Date") + ylab("Strength of EWS") +
  geom_vline(xintercept = total.multiphlplot.dat$true.date[total.multiphlplot.dat$period != lag(total.multiphlplot.dat$period)], linetype="dashed", 
             color = "black", size=1)+
  scale_x_date(date_breaks = "3 months", date_labels = "%b-%Y")+
  labs(color='EWS Indicator',alpha ="EWS") +
  theme(plot.margin = margin(c(10, 8, 0, 10)),
        legend.key.height = unit(0.3,"cm"),
        legend.key.width = unit(0.1,"cm"),
        legend.title = element_text(size = 10),
        legend.text = element_text(size=10))


multi.phl.p1 <- ggplot(multi.deriv.cov.phl, aes(x=true.date, y=cases)) +
  aes(group=NA)+
  geom_path(aes(col = as.factor(deriv.direction)))+
  geom_point(data= total.multiphlplot.dat[which(total.multiphlplot.dat$threshold.crossed==1 & total.multiphlplot.dat$metric.code == "acf + SD + skew")], 
             aes(y =-1000, alpha = ""),size = 3,pch= "|",col = "#53B400")+
  ylab("COVID-19 Cases") + 
  xlab("Date")+
  geom_vline(xintercept = total.multiphlplot.dat$true.date[total.multiphlplot.dat$period != lag(total.multiphlplot.dat$period)], 
             linetype="dashed",  color = "black", size=1)+
  scale_colour_manual(values=c("#22B4F5","black","#F07589"),name = "Predicted Nonlinear\nCase Trend", labels = c("Negative","No trend","Positive"))+
  scale_y_continuous(labels = scales::number_format(accuracy = 1000)) + 
  scale_x_date(date_breaks = "3 months", date_labels = "%b-%Y") +
  ggthemes::theme_clean() + ggtitle("Philippines Daily COVID Cases: Two Waves")+  
  annotate("label", x = as.Date("2020-06-01"), y =max(total.multiphlplot.dat$cases)*0.8 , label = "EWS indicator: acf + SD + skewness")+
  scale_alpha_manual(values = 0.8, labels = "EWS detected") +
  labs(alpha = "EWS")  +
  theme(plot.margin = margin(c(10, 8, 0, 10)),
        legend.key.height = unit(0.3,"cm" ),
        legend.key.width = unit(0.1,"cm"),
        legend.title = element_text(size = 10),
        legend.text = element_text(size=10))

png(file="Results/Philippines/multigam_phl_fig.png",
    width=3000, height = 2500 ,res=300)
egg::ggarrange(multi.phl.p1,multi.phl.p2,nrow = 2,heights = c(2, 2), labels = c("a","b"))
dev.off()

phl.multiews.diff2 <- extract.ews.difference(total.multiphlplot.dat,multi.deriv.cov.phl,"Philippines",2)

###########################################################################################
### World - Poland ###
###########################################################################################
#https://covid19.who.int/info/
cov.pol.dat <- cov.WHO.dat %>%
  filter(Country == "Poland")%>%
  slice(65:n())

ggplot(cov.pol.dat,aes(x=Date, y = cases)) + geom_path() + 
  scale_x_continuous(breaks=seq(2020,2022,0.2)) + 
  scale_y_continuous( labels = scales::number_format(accuracy = 100))+
  ggtitle("Daily Poland COVID Cases")+ 
  ylab("Cases") + theme_bw()

#### Expanding GAM
exp.pol <- expanding.gam(cov.pol.dat,sensitivity = 7,train.length = 10)
plot(cov.pol.dat$cases)
abline(v=exp.pol$end.index)

polgam1 <-mgcv::gam(abs(cases)~ s(Weekday,bs = "cs",k=7)+s(Date, bs="tp", k=30), 
                    data = cov.pol.dat[exp.pol$start.index[1]:exp.pol$end.index[4],],
                    family = gaussian(), method = "REML")
plot.gam(polgam1,pages=1)
pol.deriv1 <- data.frame(gratia::derivatives(polgam1, term = "s(Date)", interval = "confidence",n= dim(cov.pol.dat[exp.pol$start.index[1]:exp.pol$end.index[4],])[1]),
                         "deriv.direction" = curve.direction(0, gratia::derivatives(polgam1, term = "s(Date)", interval = "confidence",n= dim(cov.pol.dat[exp.pol$start.index[1]:exp.pol$end.index[4],])[1])$lower, gratia::derivatives(polgam1, term = "s(Date)", interval = "confidence",n= dim(cov.pol.dat[exp.pol$start.index[1]:exp.pol$end.index[4],])[1])$upper)/2)
plot(pol.deriv1$deriv.direction)

polgam2 <-mgcv::gam(abs(cases)~ s(Weekday,bs = "cs",k=7) + s(Date, bs="tp", k=30) , 
                    data =  cov.pol.dat[exp.pol$start.index[5]:dim(cov.pol.dat)[1],],
                    family = gaussian(), method = "REML")
plot.gam(polgam2,pages=1)
pol.deriv2 <- data.frame(gratia::derivatives(polgam2, term = "s(Date)", interval = "confidence",n= dim(cov.pol.dat[exp.pol$start.index[5]:dim(cov.pol.dat)[1],])[1]),
                         "deriv.direction" = curve.direction(0, gratia::derivatives(polgam2, term = "s(Date)", interval = "confidence",n= dim(cov.pol.dat[exp.pol$start.index[5]:dim(cov.pol.dat)[1],])[1])$lower, gratia::derivatives(polgam2, term = "s(Date)", interval = "confidence",n= dim(cov.pol.dat[exp.pol$start.index[5]:dim(cov.pol.dat)[1],])[1])$upper)/2)
plot(pol.deriv2$deriv.direction)

cutoff.pol.cov.ews.multigam <- pbmclapply(cov.pol.dat$Date[1:(length(cov.pol.dat$Date)-13)],FUN = function(x){
  
  data.cut <- cov.pol.dat[as.numeric(cov.pol.dat$Date) >=  as.numeric(x),]
  
  tmp <- comp_EWS_wrapper(data.frame(timedat = as.numeric(data.cut$Date),
                                     biomass = data.cut$cases),
                          metrics = metrics,  threshold = 2, burn_in =7,
                          plotIt = F, ggplotIt = F, tail.direction = "one.tailed",
                          interpolate = F, method = "w_comp")
  tmp$cutoff <- paste(x)
  return(tmp)
  
}, mc.cores =3 )
names(cutoff.pol.cov.ews.multigam) <- cov.pol.dat$true.date[1:(length(cov.pol.dat$true.date)-13)]
temp <- cutoff.pol.cov.ews.multigam
cutoff.pol.cov.ews.multigam <- data.table::rbindlist(cutoff.pol.cov.ews.multigam, idcol = "start.date")%>%
  mutate(period = ifelse(time <= exp.pol$end.date[4], "first", #define periods of constancy prior to waves
                         ifelse(time >= exp.pol$start.date[5] , "second","third")))
save(cutoff.pol.cov.ews.multigam,file = "Results/Poland/pol.multi.gam.ews.RData")
load("Results/Poland/pol.multi.gam.ews.RData")

pol.multifirst.plot.dat <- cutoff.pol.cov.ews.multigam[cutoff.pol.cov.ews.multigam$start.date == "2020-03-07" & cutoff.pol.cov.ews.multigam$period == "first",]%>%
  left_join(data.frame("true.date" = cov.pol.dat$true.date,"time" = cov.pol.dat$Date), by = c("time"))

pol.multisecond.plot.dat <- cutoff.pol.cov.ews.multigam[cutoff.pol.cov.ews.multigam$start.date == cov.pol.dat$true.date[exp.pol$start.index[5]] & cutoff.pol.cov.ews.multigam$period == "second",]%>%
  left_join(data.frame("true.date" = cov.pol.dat$true.date,"time" = cov.pol.dat$Date), by = c("time"))

total.multipolplot.dat <- rbind(pol.multifirst.plot.dat,pol.multisecond.plot.dat) %>%
  rename(cases = count.used)

multi.deriv.cov.pol <- cov.pol.dat[,2:4]%>%
  cbind(rbind(pol.deriv1,pol.deriv2))


multi.pol.p2 <- ggplot(data = total.multipolplot.dat, aes(x=true.date,y=str)) +
  geom_hline(yintercept = 2, linetype="solid", color = "grey", size=1)+
  geom_line(aes(col= metric.code,alpha= "Indicator strength"))+
  geom_point(data =total.multipolplot.dat[which(total.multipolplot.dat$threshold.crossed==1)], aes(x=true.date, y = str, col = metric.code,alpha= "EWS detected")) +
  scale_alpha_manual(values = c(1, 1),
                     breaks = c("Indicator strength", "EWS detected"),
                     guide = guide_legend(override.aes = list(linetype = c(1, 0),shape = c(NA, 16)))) +
  #scale_colour_manual(values = scales::hue_pal()(7)) +
  ggthemes::theme_clean() + xlab("Date") + ylab("Strength of EWS") +
  geom_vline(xintercept = total.multipolplot.dat$true.date[total.multipolplot.dat$period != lag(total.multipolplot.dat$period)], linetype="dashed", 
             color = "black", size=1)+
  scale_x_date(date_breaks = "3 months", date_labels = "%b-%Y")+
  labs(color='EWS Indicator',alpha ="EWS") +
  theme(plot.margin = margin(c(10, 8, 0, 10)),
        legend.key.height = unit(0.3,"cm"),
        legend.key.width = unit(0.1,"cm"),
        legend.title = element_text(size = 10),
        legend.text = element_text(size=10))


multi.pol.p1 <- ggplot(multi.deriv.cov.pol, aes(x=true.date, y=cases)) +
  aes(group=NA)+
  geom_path(aes(col = as.factor(deriv.direction)))+
  geom_point(data= total.multipolplot.dat[which(total.multipolplot.dat$threshold.crossed==1 & total.multipolplot.dat$metric.code == "acf + SD + skew")], 
             aes(y =-1000, alpha = ""),size = 3,pch= "|",col = "#53B400")+
  ylab("COVID-19 Cases") + 
  xlab("Date")+
  geom_vline(xintercept = total.multipolplot.dat$true.date[total.multipolplot.dat$period != lag(total.multipolplot.dat$period)], 
             linetype="dashed",  color = "black", size=1)+
  scale_colour_manual(values=c("#22B4F5","black","#F07589"),name = "Predicted Nonlinear\nCase Trend", labels = c("Negative","No trend","Positive"))+
  scale_y_continuous(labels = scales::number_format(accuracy = 1000)) + 
  scale_x_date(date_breaks = "3 months", date_labels = "%b-%Y") +
  ggthemes::theme_clean() + ggtitle("Poland Daily COVID Cases: Two Waves")+  
  annotate("label", x = as.Date("2020-06-01"), y =max(total.multipolplot.dat$cases)*0.8 , label = "EWS indicator: acf + SD + skewness")+
  scale_alpha_manual(values = 0.8, labels = "EWS detected") +
  labs(alpha = "EWS")  +
  theme(plot.margin = margin(c(10, 8, 0, 10)),
        legend.key.height = unit(0.3,"cm" ),
        legend.key.width = unit(0.1,"cm"),
        legend.title = element_text(size = 10),
        legend.text = element_text(size=10))

png(file="Results/Poland/multigam_pol_fig.png",
    width=3000, height = 2500 ,res=300)
egg::ggarrange(multi.pol.p1,multi.pol.p2,nrow = 2,heights = c(2, 2), labels = c("a","b"))
dev.off()

pol.multiews.diff2 <- extract.ews.difference(total.multipolplot.dat,multi.deriv.cov.pol,"Poland",2)

###########################################################################################
### World - Portugal ###
###########################################################################################
#https://covid19.who.int/info/
cov.por.dat <- cov.WHO.dat %>%
  filter(Country == "Portugal")%>%
  slice(122:n())

ggplot(cov.por.dat,aes(x=Date, y = cases)) + geom_path() + 
  scale_x_continuous(breaks=seq(2020,2022,0.2)) + 
  scale_y_continuous( labels = scales::number_format(accuracy = 100))+
  ggtitle("Daily Portugal COVID Cases")+ 
  ylab("Cases") + theme_bw()

#### Expanding GAM
exp.por <- expanding.gam(cov.por.dat,sensitivity = 7,train.length = 10)
plot(cov.por.dat$cases)
abline(v=exp.por$end.index)

porgam1 <-mgcv::gam(abs(cases)~ s(Weekday,bs = "cs",k=7)+s(Date, bs="tp", k=25), 
                    data = cov.por.dat[exp.por$start.index[1]:exp.por$end.index[4],],
                    family = gaussian(), method = "REML")
plot.gam(porgam1,pages=1)
por.deriv1 <- data.frame(gratia::derivatives(porgam1, term = "s(Date)", interval = "confidence",n= dim(cov.por.dat[exp.por$start.index[1]:exp.por$end.index[4],])[1]),
                         "deriv.direction" = curve.direction(0, gratia::derivatives(porgam1, term = "s(Date)", interval = "confidence",n= dim(cov.por.dat[exp.por$start.index[1]:exp.por$end.index[4],])[1])$lower, gratia::derivatives(porgam1, term = "s(Date)", interval = "confidence",n= dim(cov.por.dat[exp.por$start.index[1]:exp.por$end.index[4],])[1])$upper)/2)
plot(por.deriv1$deriv.direction)

porgam2 <-mgcv::gam(abs(cases)~ s(Weekday,bs = "cs",k=7) + s(Date, bs="tp", k=30) , 
                    data =  cov.por.dat[exp.por$start.index[5]:dim(cov.por.dat)[1],],
                    family = gaussian(), method = "REML")
plot.gam(porgam2,pages=1)
por.deriv2 <- data.frame(gratia::derivatives(porgam2, term = "s(Date)", interval = "confidence",n= dim(cov.por.dat[exp.por$start.index[5]:dim(cov.por.dat)[1],])[1]),
                         "deriv.direction" = curve.direction(0, gratia::derivatives(porgam2, term = "s(Date)", interval = "confidence",n= dim(cov.por.dat[exp.por$start.index[5]:dim(cov.por.dat)[1],])[1])$lower, gratia::derivatives(porgam2, term = "s(Date)", interval = "confidence",n= dim(cov.por.dat[exp.por$start.index[5]:dim(cov.por.dat)[1],])[1])$upper)/2)
plot(por.deriv2$deriv.direction)

cutoff.por.cov.ews.multigam <- pbmclapply(cov.por.dat$Date[1:(length(cov.por.dat$Date)-13)],FUN = function(x){
  
  data.cut <- cov.por.dat[as.numeric(cov.por.dat$Date) >=  as.numeric(x),]
  
  tmp <- comp_EWS_wrapper(data.frame(timedat = as.numeric(data.cut$Date),
                                     biomass = data.cut$cases),
                          metrics = metrics,  threshold = 2, burn_in =7,
                          plotIt = F, ggplotIt = F, tail.direction = "one.tailed",
                          interpolate = F, method = "w_comp")
  tmp$cutoff <- paste(x)
  return(tmp)
  
}, mc.cores =3 )
names(cutoff.por.cov.ews.multigam) <- cov.por.dat$true.date[1:(length(cov.por.dat$true.date)-13)]
temp <- cutoff.por.cov.ews.multigam
cutoff.por.cov.ews.multigam <- temp
cutoff.por.cov.ews.multigam <- data.table::rbindlist(cutoff.por.cov.ews.multigam, idcol = "start.date")%>%
  mutate(period = ifelse(time > exp.por$start.date[1] & time <= exp.por$end.date[4], "first", 
                         ifelse(time >= exp.por$start.date[5], "second","third")))
save(cutoff.por.cov.ews.multigam,file = "Results/Portugal/por.multi.gam.ews.RData")
load("Results/Portugal/por.multi.gam.ews.RData")

por.multifirst.plot.dat <- cutoff.por.cov.ews.multigam[cutoff.por.cov.ews.multigam$start.date == "2020-05-03" & cutoff.por.cov.ews.multigam$period == "first",]%>%
  left_join(data.frame("true.date" = cov.por.dat$true.date,"time" = cov.por.dat$Date), by = c("time"))

por.multisecond.plot.dat <- cutoff.por.cov.ews.multigam[cutoff.por.cov.ews.multigam$start.date == cov.por.dat$true.date[exp.por$start.index[5]] & cutoff.por.cov.ews.multigam$period == "second",]%>%
  left_join(data.frame("true.date" = cov.por.dat$true.date,"time" = cov.por.dat$Date), by = c("time"))

total.multiporplot.dat <- rbind(por.multifirst.plot.dat,por.multisecond.plot.dat) %>%
  rename(cases = count.used)

multi.deriv.cov.por <- cov.por.dat[,2:4]%>%
  cbind(rbind(por.deriv1,por.deriv2))


multi.por.p2 <- ggplot(data = total.multiporplot.dat, aes(x=true.date,y=str)) +
  geom_hline(yintercept = 2, linetype="solid", color = "grey", size=1)+
  geom_line(aes(col= metric.code,alpha= "Indicator strength"))+
  geom_point(data =total.multiporplot.dat[which(total.multiporplot.dat$threshold.crossed==1)], aes(x=true.date, y = str, col = metric.code,alpha= "EWS detected")) +
  scale_alpha_manual(values = c(1, 1),
                     breaks = c("Indicator strength", "EWS detected"),
                     guide = guide_legend(override.aes = list(linetype = c(1, 0),shape = c(NA, 16)))) +
  #scale_colour_manual(values = scales::hue_pal()(7)) +
  ggthemes::theme_clean() + xlab("Date") + ylab("Strength of EWS") +
  geom_vline(xintercept = total.multiporplot.dat$true.date[total.multiporplot.dat$period != lag(total.multiporplot.dat$period)], linetype="dashed", 
             color = "black", size=1)+
  scale_x_date(date_breaks = "3 months", date_labels = "%b-%Y")+
  labs(color='EWS Indicator',alpha ="EWS") +
  theme(plot.margin = margin(c(10, 8, 0, 10)),
        legend.key.height = unit(0.3,"cm"),
        legend.key.width = unit(0.1,"cm"),
        legend.title = element_text(size = 10),
        legend.text = element_text(size=10))


multi.por.p1 <- ggplot(multi.deriv.cov.por, aes(x=true.date, y=cases)) +
  aes(group=NA)+
  geom_path(aes(col = as.factor(deriv.direction)))+
  geom_point(data= total.multiporplot.dat[which(total.multiporplot.dat$threshold.crossed==1 & total.multiporplot.dat$metric.code == "acf + SD + skew")], 
             aes(y =-1000, alpha = ""),size = 3,pch= "|",col = "#53B400")+
  ylab("COVID-19 Cases") + 
  xlab("Date")+
  geom_vline(xintercept = total.multiporplot.dat$true.date[total.multiporplot.dat$period != lag(total.multiporplot.dat$period)], 
             linetype="dashed",  color = "black", size=1)+
  scale_colour_manual(values=c("#22B4F5","black","#F07589"),name = "Predicted Nonlinear\nCase Trend", labels = c("Negative","No trend","Positive"))+
  scale_y_continuous(labels = scales::number_format(accuracy = 1000)) + 
  scale_x_date(date_breaks = "3 months", date_labels = "%b-%Y") +
  ggthemes::theme_clean() + ggtitle("Portugal Daily COVID Cases: One Wave")+  
  annotate("label", x = as.Date("2020-07-01"), y =max(total.multiporplot.dat$cases)*0.8 , label = "EWS indicator: acf + SD + skewness")+
  scale_alpha_manual(values = 0.8, labels = "EWS detected") +
  labs(alpha = "EWS")  +
  theme(plot.margin = margin(c(10, 8, 0, 10)),
        legend.key.height = unit(0.3,"cm" ),
        legend.key.width = unit(0.1,"cm"),
        legend.title = element_text(size = 10),
        legend.text = element_text(size=10))

png(file="Results/Portugal/multigam_por_fig.png",
    width=3000, height = 2500 ,res=300)
egg::ggarrange(multi.por.p1,multi.por.p2,nrow = 2,heights = c(2, 2), labels = c("a","b"))
dev.off()

por.multiews.diff2 <- extract.ews.difference(total.multiporplot.dat,multi.deriv.cov.por,"Portugal",2)

###########################################################################################
### World - Singapore ###
###########################################################################################
#https://covid19.who.int/info/
cov.sing.dat <- cov.WHO.dat %>%
  filter(Country == "Singapore")%>%
  slice(54:n())

ggplot(cov.sing.dat,aes(x=Date, y = cases)) + geom_path() + 
  scale_x_continuous(breaks=seq(2020,2022,0.2)) + 
  scale_y_continuous( labels = scales::number_format(accuracy = 100))+
  ggtitle("Daily Singapore COVID Cases")+ 
  ylab("Cases") + theme_bw()

#### Expanding GAM
exp.sing <- expanding.gam(cov.sing.dat,sensitivity = 7,train.length = 10)
plot(cov.sing.dat$cases)
abline(v=exp.sing$end.index)

singgam1 <-mgcv::gam(abs(cases)~ s(Weekday,bs = "cs",k=7)+s(Date, bs="tp", k=25), 
                    data = cov.sing.dat[exp.sing$start.index[1]:exp.sing$end.index[3],],
                    family = gaussian(), method = "REML")
plot.gam(singgam1,pages=1)
sing.deriv1 <- data.frame(gratia::derivatives(singgam1, term = "s(Date)", interval = "confidence",n= dim(cov.sing.dat[exp.sing$start.index[1]:exp.sing$end.index[3],])[1]),
                          "deriv.direction" = curve.direction(0, gratia::derivatives(singgam1, term = "s(Date)", interval = "confidence",n= dim(cov.sing.dat[exp.sing$start.index[1]:exp.sing$end.index[3],])[1])$lower, gratia::derivatives(singgam1, term = "s(Date)", interval = "confidence",n= dim(cov.sing.dat[exp.sing$start.index[1]:exp.sing$end.index[3],])[1])$upper)/2)
plot(sing.deriv1$deriv.direction)
singgam2 <-mgcv::gam(abs(cases)~ s(Weekday,bs = "cs",k=7) + s(Date, bs="tp", k=40) , 
                    data =  cov.sing.dat[exp.sing$start.index[4]:dim(cov.sing.dat)[1],],
                    family = gaussian(), method = "REML")
plot.gam(singgam2,pages=1)
sing.deriv2 <- data.frame(gratia::derivatives(singgam2, term = "s(Date)", interval = "confidence",n= dim(cov.sing.dat[exp.sing$start.index[4]:dim(cov.sing.dat)[1],])[1]),
                          "deriv.direction" = curve.direction(0, gratia::derivatives(singgam2, term = "s(Date)", interval = "confidence",n= dim(cov.sing.dat[exp.sing$start.index[4]:dim(cov.sing.dat)[1],])[1])$lower, gratia::derivatives(singgam2, term = "s(Date)", interval = "confidence",n= dim(cov.sing.dat[exp.sing$start.index[4]:dim(cov.sing.dat)[1],])[1])$upper)/2)
plot(sing.deriv2$deriv.direction)

cutoff.sing.cov.ews.multigam <- pbmclapply(cov.sing.dat$Date[1:(length(cov.sing.dat$Date)-13)],FUN = function(x){
  
  data.cut <- cov.sing.dat[as.numeric(cov.sing.dat$Date) >=  as.numeric(x),]
  
  tmp <- comp_EWS_wrapper(data.frame(timedat = as.numeric(data.cut$Date),
                                      biomass = data.cut$cases),
                          metrics = metrics,  threshold = 2, burn_in =7,
                          plotIt = F, ggplotIt = F, tail.direction = "one.tailed",
                          interpolate = F, method = "w_comp")
  tmp$cutoff <- paste(x)
  return(tmp)
  
}, mc.cores =3 )
names(cutoff.sing.cov.ews.multigam) <- cov.sing.dat$true.date[1:(length(cov.sing.dat$true.date)-13)]
temp <- cutoff.sing.cov.ews.multigam
cutoff.sing.cov.ews.multigam <- temp
cutoff.sing.cov.ews.multigam <- data.table::rbindlist(cutoff.sing.cov.ews.multigam, idcol = "start.date")%>%
  mutate(period = ifelse(time <= exp.sing$end.date[3], "first", #define periods of constancy prior to waves
                         ifelse(time >= exp.sing$start.date[4] ,"second","third")))
save(cutoff.sing.cov.ews.multigam,file = "Results/Singapore/sing.multi.gam.ews.RData")
load("Results/Singapore/sing.multi.gam.ews.RData")

sing.multifirst.plot.dat <- cutoff.sing.cov.ews.multigam[cutoff.sing.cov.ews.multigam$start.date == "2020-02-25" & cutoff.sing.cov.ews.multigam$period == "first",]%>%
  left_join(data.frame("true.date" = cov.sing.dat$true.date,"time" = cov.sing.dat$Date), by = c("time"))

sing.multisecond.plot.dat <- cutoff.sing.cov.ews.multigam[cutoff.sing.cov.ews.multigam$start.date == cov.sing.dat$true.date[exp.sing$start.index[4]] & cutoff.sing.cov.ews.multigam$period == "second",]%>%
  left_join(data.frame("true.date" = cov.sing.dat$true.date,"time" = cov.sing.dat$Date), by = c("time"))

total.multisingplot.dat <- rbind(sing.multifirst.plot.dat,sing.multisecond.plot.dat) %>%
  rename(cases = count.used)

multi.deriv.cov.sing <- cov.sing.dat[,2:4]%>%
  cbind(rbind(sing.deriv1,sing.deriv2))


multi.sing.p2 <- ggplot(data = total.multisingplot.dat, aes(x=true.date,y=str)) +
  geom_hline(yintercept = 2, linetype="solid", color = "grey", size=1)+
  geom_line(aes(col= metric.code,alpha= "Indicator strength"))+
  geom_point(data =total.multisingplot.dat[which(total.multisingplot.dat$threshold.crossed==1)], aes(x=true.date, y = str, col = metric.code,alpha= "EWS detected")) +
  scale_alpha_manual(values = c(1, 1),
                     breaks = c("Indicator strength", "EWS detected"),
                     guide = guide_legend(override.aes = list(linetype = c(1, 0),shape = c(NA, 16)))) +
  #scale_colour_manual(values = scales::hue_pal()(7)) +
  ggthemes::theme_clean() + xlab("Date") + ylab("Strength of EWS") +
  geom_vline(xintercept = total.multisingplot.dat$true.date[total.multisingplot.dat$period != lag(total.multisingplot.dat$period)], linetype="dashed", 
             color = "black", size=1)+
  scale_x_date(date_breaks = "3 months", date_labels = "%b-%Y")+
  labs(color='EWS Indicator',alpha ="EWS") +
  theme(plot.margin = margin(c(10, 8, 0, 10)),
        legend.key.height = unit(0.3,"cm"),
        legend.key.width = unit(0.1,"cm"),
        legend.title = element_text(size = 10),
        legend.text = element_text(size=10))


multi.sing.p1 <- ggplot(multi.deriv.cov.sing, aes(x=true.date, y=cases)) +
  aes(group=NA)+
  geom_path(aes(col = as.factor(deriv.direction)))+
  geom_point(data= total.multisingplot.dat[which(total.multisingplot.dat$threshold.crossed==1 & total.multisingplot.dat$metric.code == "acf + SD + skew")], 
             aes(y =-100, alpha = ""),size = 3,pch= "|",col = "#53B400")+
  ylab("COVID-19 Cases") + 
  xlab("Date")+
  geom_vline(xintercept = total.multisingplot.dat$true.date[total.multisingplot.dat$period != lag(total.multisingplot.dat$period)], 
             linetype="dashed",  color = "black", size=1)+
  scale_colour_manual(values=c("#22B4F5","black","#F07589"),name = "Predicted Nonlinear\nCase Trend", labels = c("Negative","No trend","Positive"))+
  scale_y_continuous(labels = scales::number_format(accuracy = 100)) + 
  scale_x_date(date_breaks = "3 months", date_labels = "%b-%Y") +
  ggthemes::theme_clean() + ggtitle("Singapore Daily COVID Cases: One Wave")+  
  annotate("label", x = as.Date("2020-04-01"), y =max(total.multisingplot.dat$cases)*0.8 , label = "EWS indicator: acf + SD + skewness")+
  scale_alpha_manual(values = 0.8, labels = "EWS detected") +
  labs(alpha = "EWS")  +
  theme(plot.margin = margin(c(10, 8, 0, 10)),
        legend.key.height = unit(0.3,"cm" ),
        legend.key.width = unit(0.1,"cm"),
        legend.title = element_text(size = 10),
        legend.text = element_text(size=10))

png(file="Results/Singapore/multigam_sing_fig.png",
    width=3000, height = 2500 ,res=300)
egg::ggarrange(multi.sing.p1,multi.sing.p2,nrow = 2,heights = c(2, 2), labels = c("a","b"))
dev.off()

sing.multiews.diff2 <- extract.ews.difference(total.multisingplot.dat,multi.deriv.cov.sing,"Singapore",2)

###########################################################################################
### World - Spain ###
###########################################################################################
#https://covid19.who.int/info/
cov.esp.dat <- cov.WHO.dat %>%
  filter(Country == "Spain")%>%
  slice(44:n())

ggplot(cov.esp.dat,aes(x=Date, y = cases)) + geom_path() + 
  scale_x_continuous(breaks=seq(2020,2022,0.2)) + 
  scale_y_continuous( labels = scales::number_format(accuracy = 100))+
  ggtitle("Daily Spain COVID Cases")+ 
  ylab("Cases") + theme_bw()

#### Expanding GAM
exp.esp <- expanding.gam(cov.esp.dat,sensitivity = 7,train.length = 10)
plot(cov.esp.dat$cases)
abline(v=exp.esp$end.index)

espgam1 <-mgcv::gam(abs(cases)~ s(Weekday,bs = "cs",k=7)+s(Date, bs="tp", k=15), 
                     data = cov.esp.dat[exp.esp$start.index[1]:exp.esp$end.index[1],],
                     family = gaussian(), method = "REML")
plot.gam(espgam1,pages=1)
esp.deriv1 <- data.frame(gratia::derivatives(espgam1, term = "s(Date)", interval = "confidence",n= dim(cov.esp.dat[exp.esp$start.index[1]:exp.esp$end.index[1],])[1]),
                            "deriv.direction" = curve.direction(0, gratia::derivatives(espgam1, term = "s(Date)", interval = "confidence",n= dim(cov.esp.dat[exp.esp$start.index[1]:exp.esp$end.index[1],])[1])$lower, gratia::derivatives(espgam1, term = "s(Date)", interval = "confidence",n= dim(cov.esp.dat[exp.esp$start.index[1]:exp.esp$end.index[1],])[1])$upper)/2)
espgam2 <-mgcv::gam(abs(cases)~ s(Weekday,bs = "cs",k=7) + s(Date, bs="tp", k=15) , 
                     data =  cov.esp.dat[exp.esp$start.index[2]:dim(cov.esp.dat)[1],],
                     family = gaussian(), method = "REML")
plot.gam(espgam2,pages=1)
esp.deriv2 <- data.frame(gratia::derivatives(espgam2, term = "s(Date)", interval = "confidence",n= dim(cov.esp.dat[exp.esp$start.index[2]:dim(cov.esp.dat)[1],])[1]),
                            "deriv.direction" = curve.direction(0, gratia::derivatives(espgam2, term = "s(Date)", interval = "confidence",n= dim(cov.esp.dat[exp.esp$start.index[2]:dim(cov.esp.dat)[1],])[1])$lower, gratia::derivatives(espgam2, term = "s(Date)", interval = "confidence",n= dim(cov.esp.dat[exp.esp$start.index[2]:dim(cov.esp.dat)[1],])[1])$upper)/2)

plot(esp.deriv2$deriv.direction)

cutoff.esp.cov.ews.multigam <- pbmclapply(cov.esp.dat$Date[1:(length(cov.esp.dat$Date)-13)],FUN = function(x){
  
  data.cut <- cov.esp.dat[as.numeric(cov.esp.dat$Date) >=  as.numeric(x),]
  
  tmp <- comp_EWS_wrapper(data.frame(timedat = as.numeric(data.cut$Date),
                                       biomass = data.cut$cases),
                          metrics = metrics,  threshold = 2, burn_in =7,
                          plotIt = F, ggplotIt = F, tail.direction = "one.tailed",
                          interpolate = F, method = "w_comp")
  tmp$cutoff <- paste(x)
  return(tmp)
  
}, mc.cores =3 )
names(cutoff.esp.cov.ews.multigam) <- cov.esp.dat$true.date[1:(length(cov.esp.dat$true.date)-13)]
temp <- cutoff.esp.cov.ews.multigam
cutoff.esp.cov.ews.multigam <- data.table::rbindlist(cutoff.esp.cov.ews.multigam, idcol = "start.date")%>%
  mutate(period = ifelse(time <= exp.esp$end.date[1], "first", #define periods of constancy prior to waves
                         ifelse(time >= exp.esp$start.date[2] , "second","third")))
save(cutoff.esp.cov.ews.multigam,file = "Results/Spain/esp.multi.gam.ews.RData")
load("Results/Spain/esp.multi.gam.ews.RData")

esp.multifirst.plot.dat <- cutoff.esp.cov.ews.multigam[cutoff.esp.cov.ews.multigam$start.date == "2020-02-15" & cutoff.esp.cov.ews.multigam$period == "first",]%>%
  left_join(data.frame("true.date" = cov.esp.dat$true.date,"time" = cov.esp.dat$Date), by = c("time"))

esp.multisecond.plot.dat <- cutoff.esp.cov.ews.multigam[cutoff.esp.cov.ews.multigam$start.date == cov.esp.dat$true.date[exp.esp$start.index[2]] & cutoff.esp.cov.ews.multigam$period == "second",]%>%
  left_join(data.frame("true.date" = cov.esp.dat$true.date,"time" = cov.esp.dat$Date), by = c("time"))

total.multiespplot.dat <- rbind(esp.multifirst.plot.dat,esp.multisecond.plot.dat) %>%
  rename(cases = count.used)

multi.deriv.cov.esp <- cov.esp.dat[,2:4]%>%
  cbind(rbind(esp.deriv1,esp.deriv2))


multi.esp.p2 <- ggplot(data = total.multiespplot.dat, aes(x=true.date,y=str)) +
  geom_hline(yintercept = 2, linetype="solid", color = "grey", size=1)+
  geom_line(aes(col= metric.code,alpha= "Indicator strength"))+
  geom_point(data =total.multiespplot.dat[which(total.multiespplot.dat$threshold.crossed==1)], aes(x=true.date, y = str, col = metric.code,alpha= "EWS detected")) +
  scale_alpha_manual(values = c(1, 1),
                     breaks = c("Indicator strength", "EWS detected"),
                     guide = guide_legend(override.aes = list(linetype = c(1, 0),shape = c(NA, 16)))) +
  #scale_colour_manual(values = scales::hue_pal()(7)) +
  ggthemes::theme_clean() + xlab("Date") + ylab("Strength of EWS") +
  geom_vline(xintercept = total.multiespplot.dat$true.date[total.multiespplot.dat$period != lag(total.multiespplot.dat$period)], linetype="dashed", 
             color = "black", size=1)+
  scale_x_date(date_breaks = "3 months", date_labels = "%b-%Y")+
  labs(color='EWS Indicator',alpha ="EWS") +
  theme(plot.margin = margin(c(10, 8, 0, 10)),
        legend.key.height = unit(0.3,"cm"),
        legend.key.width = unit(0.1,"cm"),
        legend.title = element_text(size = 10),
        legend.text = element_text(size=10))

multi.esp.p1 <- ggplot(multi.deriv.cov.esp, aes(x=true.date, y=cases)) +
  aes(group=NA)+
  geom_path(aes(col = as.factor(deriv.direction)))+
  geom_point(data= total.multiespplot.dat[which(total.multiespplot.dat$threshold.crossed==1 & total.multiespplot.dat$metric.code == "acf + SD + skew")], 
             aes(y =-5000, alpha = ""),size = 3,pch= "|",col = "#53B400")+
  ylab("COVID-19 Cases") + 
  xlab("Date")+
  geom_vline(xintercept = total.multiespplot.dat$true.date[total.multiespplot.dat$period != lag(total.multiespplot.dat$period)], 
             linetype="dashed",  color = "black", size=1)+
  scale_colour_manual(values=c("#22B4F5","black","#F07589"),name = "Predicted Nonlinear\nCase Trend", labels = c("Negative","No trend","Positive"))+
  scale_y_continuous(labels = scales::number_format(accuracy = 1000)) + 
  scale_x_date(date_breaks = "3 months", date_labels = "%b-%Y") +
  ggthemes::theme_clean() + ggtitle("Spain Daily COVID Cases: Two Waves")+  
  annotate("label", x = as.Date("2020-05-01"), y =max(total.multiespplot.dat$cases)*0.8 , label = "EWS indicator: acf + SD + skewness")+
  scale_alpha_manual(values = 0.8, labels = "EWS detected") +
  labs(alpha = "EWS")  +
  theme(plot.margin = margin(c(10, 8, 0, 10)),
        legend.key.height = unit(0.3,"cm" ),
        legend.key.width = unit(0.1,"cm"),
        legend.title = element_text(size = 10),
        legend.text = element_text(size=10))

png(file="Results/Spain/multigam_esp_fig.png",
    width=3000, height = 2500 ,res=300)
egg::ggarrange(multi.esp.p1,multi.esp.p2,nrow = 2,heights = c(2, 2), labels = c("a","b"))
dev.off()

esp.multiews.diff2 <- extract.ews.difference(total.multiespplot.dat,multi.deriv.cov.esp,"Spain",2)

###########################################################################################
### World - South Africa ###
###########################################################################################
#https://covid19.who.int/info/
cov.sa.dat <- cov.WHO.dat %>%
  filter(Country == "South Africa")%>%
  slice(69:n())

ggplot(cov.sa.dat,aes(x=Date, y = cases)) + geom_path() + 
  scale_x_continuous(breaks=seq(2020,2022,0.2)) + 
  scale_y_continuous( labels = scales::number_format(accuracy = 100))+
  ggtitle("Daily South Africa COVID Cases")+ 
  ylab("Cases") + theme_bw()

#### Expanding GAM
exp.sa <- expanding.gam(cov.sa.dat,sensitivity = 7,train.length = 10)
plot(cov.sa.dat$cases)
abline(v=exp.sa$end.index)

sagam1 <-mgcv::gam(abs(cases)~ s(Weekday,bs = "cs",k=7)+s(Date, bs="tp", k=20), 
                    data = cov.sa.dat[exp.sa$start.index[1]:exp.sa$end.index[3],],
                    family = gaussian(), method = "REML")
plot.gam(sagam1,pages=1)
sa.deriv1 <- data.frame(gratia::derivatives(sagam1, term = "s(Date)", interval = "confidence",n= dim(cov.sa.dat[exp.sa$start.index[1]:exp.sa$end.index[3],])[1]),
                         "deriv.direction" = curve.direction(0, gratia::derivatives(sagam1, term = "s(Date)", interval = "confidence",n= dim(cov.sa.dat[exp.sa$start.index[1]:exp.sa$end.index[3],])[1])$lower, gratia::derivatives(sagam1, term = "s(Date)", interval = "confidence",n= dim(cov.sa.dat[exp.sa$start.index[1]:exp.sa$end.index[3],])[1])$upper)/2)
plot(sa.deriv1$deriv.direction)
sagam2 <-mgcv::gam(abs(cases)~ s(Weekday,bs = "cs",k=7) + s(Date, bs="tp", k=15) , 
                    data =  cov.sa.dat[exp.sa$start.index[4]:exp.sa$end.index[4],],
                    family = gaussian(), method = "REML")
plot.gam(sagam2,pages=1)
sa.deriv2 <- data.frame(gratia::derivatives(sagam2, term = "s(Date)", interval = "confidence",n= dim(cov.sa.dat[exp.sa$start.index[4]:exp.sa$end.index[4],])[1]),
                         "deriv.direction" = curve.direction(0, gratia::derivatives(sagam2, term = "s(Date)", interval = "confidence",n= dim(cov.sa.dat[exp.sa$start.index[4]:exp.sa$end.index[4],])[1])$lower, gratia::derivatives(sagam2, term = "s(Date)", interval = "confidence",n= dim(cov.sa.dat[exp.sa$start.index[4]:exp.sa$end.index[4],])[1])$upper)/2)
plot(sa.deriv2$deriv.direction)
sagam3 <-mgcv::gam(abs(cases)~ s(Weekday,bs = "cs",k=7) + s(Date, bs="tp", k=40) , 
                   data =  cov.sa.dat[exp.sa$start.index[5]:dim(cov.sa.dat)[1],],
                   family = gaussian(), method = "REML")
plot.gam(sagam3,pages=1)
sa.deriv3 <- data.frame(gratia::derivatives(sagam3, term = "s(Date)", interval = "confidence",n= dim(cov.sa.dat[exp.sa$start.index[5]:dim(cov.sa.dat)[1],])[1]),
                        "deriv.direction" = curve.direction(0, gratia::derivatives(sagam3, term = "s(Date)", interval = "confidence",n= dim(cov.sa.dat[exp.sa$start.index[5]:dim(cov.sa.dat)[1],])[1])$lower, gratia::derivatives(sagam3, term = "s(Date)", interval = "confidence",n= dim(cov.sa.dat[exp.sa$start.index[5]:dim(cov.sa.dat)[1],])[1])$upper)/2)
plot(sa.deriv3$deriv.direction)

cutoff.sa.cov.ews.multigam <- pbmclapply(cov.sa.dat$Date[1:(length(cov.sa.dat$Date)-13)],FUN = function(x){
  
  data.cut <- cov.sa.dat[as.numeric(cov.sa.dat$Date) >=  as.numeric(x),]
  
  tmp <- comp_EWS_wrapper(data.frame(timedat = as.numeric(data.cut$Date),
                                     biomass = data.cut$cases),
                          metrics = metrics,  threshold = 2, burn_in =7,
                          plotIt = F, ggplotIt = F, tail.direction = "one.tailed",
                          interpolate = F, method = "w_comp")
  tmp$cutoff <- paste(x)
  return(tmp)
  
}, mc.cores =3 )
names(cutoff.sa.cov.ews.multigam) <- cov.sa.dat$true.date[1:(length(cov.sa.dat$true.date)-13)]
temp <- cutoff.sa.cov.ews.multigam
cutoff.sa.cov.ews.multigam <- temp
cutoff.sa.cov.ews.multigam <- data.table::rbindlist(cutoff.sa.cov.ews.multigam, idcol = "start.date")%>%
  mutate(period = ifelse(time <= exp.sa$end.date[3], "first", #define periods of constancy prior to waves
                         ifelse(time >= exp.sa$start.date[4] & time <= exp.sa$end.date[4] , "second",
                                ifelse(time >= exp.sa$start.date[5],"third","fourth"))))
save(cutoff.sa.cov.ews.multigam,file = "Results/South Africa/sa.multi.gam.ews.RData")
load("Results/South Africa/sa.multi.gam.ews.RData")

sa.multifirst.plot.dat <- cutoff.sa.cov.ews.multigam[cutoff.sa.cov.ews.multigam$start.date == "2020-03-11" & cutoff.sa.cov.ews.multigam$period == "first",]%>%
  left_join(data.frame("true.date" = cov.sa.dat$true.date,"time" = cov.sa.dat$Date), by = c("time"))

sa.multisecond.plot.dat <- cutoff.sa.cov.ews.multigam[cutoff.sa.cov.ews.multigam$start.date == cov.sa.dat$true.date[exp.sa$start.index[4]] & cutoff.sa.cov.ews.multigam$period == "second",]%>%
  left_join(data.frame("true.date" = cov.sa.dat$true.date,"time" = cov.sa.dat$Date), by = c("time"))

sa.multithird.plot.dat <- cutoff.sa.cov.ews.multigam[cutoff.sa.cov.ews.multigam$start.date == cov.sa.dat$true.date[exp.sa$start.index[5]] & cutoff.sa.cov.ews.multigam$period == "third",]%>%
  left_join(data.frame("true.date" = cov.sa.dat$true.date,"time" = cov.sa.dat$Date), by = c("time"))

total.multisaplot.dat <- rbind(sa.multifirst.plot.dat,sa.multisecond.plot.dat,sa.multithird.plot.dat) %>%
  rename(cases = count.used)

multi.deriv.cov.sa <- cov.sa.dat[,2:4]%>%
  cbind(rbind(sa.deriv1,sa.deriv2,sa.deriv3))


multi.sa.p2 <- ggplot(data = total.multisaplot.dat, aes(x=true.date,y=str)) +
  geom_hline(yintercept = 2, linetype="solid", color = "grey", size=1)+
  geom_line(aes(col= metric.code,alpha= "Indicator strength"))+
  geom_point(data =total.multisaplot.dat[which(total.multisaplot.dat$threshold.crossed==1)], aes(x=true.date, y = str, col = metric.code,alpha= "EWS detected")) +
  scale_alpha_manual(values = c(1, 1),
                     breaks = c("Indicator strength", "EWS detected"),
                     guide = guide_legend(override.aes = list(linetype = c(1, 0),shape = c(NA, 16)))) +
  #scale_colour_manual(values = scales::hue_pal()(7)) +
  ggthemes::theme_clean() + xlab("Date") + ylab("Strength of EWS") +
  geom_vline(xintercept = total.multisaplot.dat$true.date[total.multisaplot.dat$period != lag(total.multisaplot.dat$period)], linetype="dashed", 
             color = "black", size=1)+
  scale_x_date(date_breaks = "3 months", date_labels = "%b-%Y")+
  labs(color='EWS Indicator',alpha ="EWS") +
  theme(plot.margin = margin(c(10, 8, 0, 10)),
        legend.key.height = unit(0.3,"cm"),
        legend.key.width = unit(0.1,"cm"),
        legend.title = element_text(size = 10),
        legend.text = element_text(size=10))


multi.sa.p1 <- ggplot(multi.deriv.cov.sa, aes(x=true.date, y=cases)) +
  aes(group=NA)+
  geom_path(aes(col = as.factor(deriv.direction)))+
  geom_point(data= total.multisaplot.dat[which(total.multisaplot.dat$threshold.crossed==1 & total.multisaplot.dat$metric.code == "acf + SD + skew")], 
             aes(y =-5000, alpha = ""),size = 3,pch= "|",col = "#53B400")+
  ylab("COVID-19 Cases") + 
  xlab("Date")+
  geom_vline(xintercept = total.multisaplot.dat$true.date[total.multisaplot.dat$period != lag(total.multisaplot.dat$period)], 
             linetype="dashed",  color = "black", size=1)+
  scale_colour_manual(values=c("#22B4F5","black","#F07589"),name = "Predicted Nonlinear\nCase Trend", labels = c("Negative","No trend","Positive"))+
  scale_y_continuous(labels = scales::number_format(accuracy = 1000)) + 
  scale_x_date(date_breaks = "3 months", date_labels = "%b-%Y") +
  ggthemes::theme_clean() + ggtitle("South Africa Daily COVID Cases: Three Waves")+  
  annotate("label", x = as.Date("2020-06-01"), y =max(total.multisaplot.dat$cases)*0.8 , label = "EWS indicator: acf + SD + skewness")+
  scale_alpha_manual(values = 0.8, labels = "EWS detected") +
  labs(alpha = "EWS")  +
  theme(plot.margin = margin(c(10, 8, 0, 10)),
        legend.key.height = unit(0.3,"cm" ),
        legend.key.width = unit(0.1,"cm"),
        legend.title = element_text(size = 10),
        legend.text = element_text(size=10))

png(file="Results/South Africa/multigam_sa_fig.png",
    width=3000, height = 2500 ,res=300)
egg::ggarrange(multi.sa.p1,multi.sa.p2,nrow = 2,heights = c(2, 2), labels = c("a","b"))
dev.off()

sa.multiews.diff2 <- extract.ews.difference(total.multisaplot.dat,multi.deriv.cov.sa,"South Africa",2)

###########################################################################################
### USA ###
###########################################################################################
#https://covid.cdc.gov/covid-data-tracker/#trends_dailytrendscases

ggplot(cov.usa.dat,aes(x=Date, y = cases)) + geom_path() + 
  scale_x_continuous(breaks=seq(2020,2022,0.2)) + ggtitle("Daily USA COVID Cases")+ 
  scale_y_continuous( labels = scales::number_format(accuracy = 100000))+
  ylab("Cases") + theme_bw()

#### Expanding GAM
exp.usa <- expanding.gam(cov.usa.dat,sensitivity = 7,train.length = 10)
plot(cov.usa.dat$cases)
abline(v=exp.usa$end.index)
usagam1 <-mgcv::gam(abs(cases)~ s(Weekday,bs = "cs",k=7)+s(Date, bs="tp", k=30), 
                    data = cov.usa.dat[exp.usa$start.index[1]:exp.usa$end.index[2],],
                    family = gaussian(), method = "REML")
plot.gam(usagam1,pages=1)
usa.deriv1 <- data.frame(gratia::derivatives(usagam1, term = "s(Date)", interval = "confidence",n= dim(cov.usa.dat[exp.usa$start.index[1]:exp.usa$end.index[2],])[1]),
                          "deriv.direction" = curve.direction(0, gratia::derivatives(usagam1, term = "s(Date)", interval = "confidence",n= dim(cov.usa.dat[exp.usa$start.index[1]:exp.usa$end.index[2],])[1])$lower, gratia::derivatives(usagam1, term = "s(Date)", interval = "confidence",n= dim(cov.usa.dat[exp.usa$start.index[1]:exp.usa$end.index[2],])[1])$upper)/2)
usagam2 <-mgcv::gam(abs(cases)~ s(Weekday,bs = "cs",k=7) + s(Date, bs="tp", k=10) , 
                    data =  cov.usa.dat[exp.usa$start.index[3]:exp.usa$end.index[5],],
                    family = gaussian(), method = "REML")
plot.gam(usagam2,pages=1)
usa.deriv2 <- data.frame(gratia::derivatives(usagam2, term = "s(Date)", interval = "confidence",n= dim(cov.usa.dat[exp.usa$start.index[3]:exp.usa$end.index[5],])[1]),
                          "deriv.direction" = curve.direction(0, gratia::derivatives(usagam2, term = "s(Date)", interval = "confidence",n= dim(cov.usa.dat[exp.usa$start.index[3]:exp.usa$end.index[5],])[1])$lower, gratia::derivatives(usagam2, term = "s(Date)", interval = "confidence",n= dim(cov.usa.dat[exp.usa$start.index[3]:exp.usa$end.index[5],])[1])$upper)/2)
plot(usa.deriv2$deriv.direction)
usagam3 <-mgcv::gam(abs(cases)~ s(Weekday,bs = "cs",k=7)+s(Date, bs="tp", k=15) , 
                    data = cov.usa.dat[(exp.usa$end.index[5]+1):dim(cov.usa.dat)[1],],
                    family = gaussian(), method = "REML")
plot.gam(usagam3,pages =1)
usa.deriv3 <- data.frame(gratia::derivatives(usagam3, term = "s(Date)", interval = "confidence",n= dim(cov.usa.dat[(exp.usa$end.index[5]+1):dim(cov.usa.dat)[1],])[1]),
                          "deriv.direction" = curve.direction(0, gratia::derivatives(usagam3, term = "s(Date)", interval = "confidence",n= dim(cov.usa.dat[(exp.usa$end.index[5]+1):dim(cov.usa.dat)[1],])[1])$lower, gratia::derivatives(usagam3, term = "s(Date)", interval = "confidence",n= dim(cov.usa.dat[(exp.usa$end.index[5]+1):dim(cov.usa.dat)[1],])[1])$upper)/2)


cutoff.usa.cov.ews.multigam <- pbmclapply(cov.usa.dat$Date[1:(length(cov.usa.dat$Date)-13)],FUN = function(x){
  
  data.cut <- cov.usa.dat[as.numeric(cov.usa.dat$Date) >=  as.numeric(x),]
  
  tmp <- comp_EWS_wrapper(data.frame(timedat = as.numeric(data.cut$Date),
                                      biomass = data.cut$cases),
                          metrics = metrics,  threshold = 2, burn_in =7,
                          plotIt = F, ggplotIt = F, tail.direction = "one.tailed",
                          interpolate = F, method = "w_comp")
  tmp$cutoff <- paste(x)
  return(tmp)
  
}, mc.cores =3 )
names(cutoff.usa.cov.ews.multigam) <- cov.usa.dat$true.date[1:(length(cov.usa.dat$true.date)-13)]
temp <- cutoff.usa.cov.ews.multigam
cutoff.usa.cov.ews.multigam <- data.table::rbindlist(cutoff.usa.cov.ews.multigam, idcol = "start.date")%>%
  mutate(period = ifelse(time <= exp.usa$end.date[2], "first", #define periods of constancy prior to waves
                         ifelse(time >= exp.usa$start.date[3] & time <= exp.usa$end.date[5], "second",
                                ifelse(time > exp.usa$end.date[5],"third","fourth"))))
save(cutoff.usa.cov.ews.multigam,file = "Results/USA/usa.multi.gam.ews.RData")
load("Results/USA/usa.multi.gam.ews.RData")

usa.multifirst.plot.dat <- cutoff.usa.cov.ews.multigam[cutoff.usa.cov.ews.multigam$start.date == "2020-01-22" & cutoff.usa.cov.ews.multigam$period == "first",]%>%
  left_join(data.frame("true.date" = cov.usa.dat$true.date,"time" = cov.usa.dat$Date), by = c("time"))

usa.multisecond.plot.dat <- cutoff.usa.cov.ews.multigam[cutoff.usa.cov.ews.multigam$start.date == cov.usa.dat$true.date[exp.usa$start.index[3]] & cutoff.usa.cov.ews.multigam$period == "second",]%>%
  left_join(data.frame("true.date" = cov.usa.dat$true.date,"time" = cov.usa.dat$Date), by = c("time"))

usa.multithird.plot.dat <- cutoff.usa.cov.ews.multigam[cutoff.usa.cov.ews.multigam$start.date == cov.usa.dat$true.date[exp.usa$end.index[5]+1] & cutoff.usa.cov.ews.multigam$period == "third",]%>%
  left_join(data.frame("true.date" = cov.usa.dat$true.date,"time" = cov.usa.dat$Date), by = c("time"))

total.multiusaplot.dat <- rbind(usa.multifirst.plot.dat,usa.multisecond.plot.dat,usa.multithird.plot.dat) %>%
  rename(cases = count.used)

multi.deriv.cov.usa <- cov.usa.dat[,2:4]%>%
  cbind(rbind(usa.deriv1,usa.deriv2,usa.deriv3))


multi.usa.p2 <- ggplot(data = total.multiusaplot.dat, aes(x=true.date,y=str)) +
  geom_hline(yintercept = 2, linetype="solid", color = "grey", size=1)+
  geom_line(aes(col= metric.code,alpha= "Indicator strength"))+
  geom_point(data =total.multiusaplot.dat[which(total.multiusaplot.dat$threshold.crossed==1)], aes(x=true.date, y = str, col = metric.code,alpha= "EWS detected")) +
  scale_alpha_manual(values = c(1, 1),
                     breaks = c("Indicator strength", "EWS detected"),
                     guide = guide_legend(override.aes = list(linetype = c(1, 0),shape = c(NA, 16)))) +
  #scale_colour_manual(values = scales::hue_pal()(7)) +
  ggthemes::theme_clean() + xlab("Date") + ylab("Strength of EWS") +
  geom_vline(xintercept = total.multiusaplot.dat$true.date[total.multiusaplot.dat$period != lag(total.multiusaplot.dat$period)], linetype="dashed", 
             color = "black", size=1)+
  scale_x_date(date_breaks = "3 months", date_labels = "%b-%Y")+
  labs(color='EWS Indicator',alpha ="EWS") +
  theme(plot.margin = margin(c(10, 8, 0, 10)),
        legend.key.height = unit(0.3,"cm"),
        legend.key.width = unit(0.1,"cm"),
        legend.title = element_text(size = 10),
        legend.text = element_text(size=10))

multi.usa.p1 <- ggplot(multi.deriv.cov.usa, aes(x=true.date, y=cases)) +
  aes(group=NA)+
  geom_path(aes(col = as.factor(deriv.direction)))+
  geom_point(data= total.multiusaplot.dat[which(total.multiusaplot.dat$threshold.crossed==1 & total.multiusaplot.dat$metric.code == "acf + SD + skew")], 
             aes(y =-15000, alpha = ""),size = 3,pch= "|",col = "#53B400")+
  ylab("COVID-19 Cases") + 
  xlab("Date")+
  geom_vline(xintercept = total.multiusaplot.dat$true.date[total.multiusaplot.dat$period != lag(total.multiusaplot.dat$period)], 
             linetype="dashed",  color = "black", size=1)+
  scale_colour_manual(values=c("#22B4F5","black","#F07589"),name = "Predicted Nonlinear\nCase Trend", labels = c("Negative","No trend","Positive"))+
  scale_y_continuous(labels = scales::number_format(accuracy = 1000)) + 
  scale_x_date(date_breaks = "3 months", date_labels = "%b-%Y") +
  ggthemes::theme_clean() + ggtitle("USA Daily COVID Cases: Two Waves")+  
  annotate("label", x = as.Date("2020-05-01"), y =max(total.multiusaplot.dat$cases)*0.8 , label = "EWS indicator: acf + SD + skewness")+
  scale_alpha_manual(values = 0.8, labels = "EWS detected") +
  labs(alpha = "EWS")  +
  theme(plot.margin = margin(c(10, 8, 0, 10)),
        legend.key.height = unit(0.3,"cm" ),
        legend.key.width = unit(0.1,"cm"),
        legend.title = element_text(size = 10),
        legend.text = element_text(size=10))

png(file="Results/USA/multigam_usa_fig.png",
    width=3000, height = 2500 ,res=300)
egg::ggarrange(multi.usa.p1,multi.usa.p2,nrow = 2,heights = c(2, 2), labels = c("a","b"))
dev.off()

usa.multiews.diff2 <- extract.ews.difference(total.multiusaplot.dat,multi.deriv.cov.usa,"USA",2)


###########################################################################################
### EWS Difference Comparison  ###
###########################################################################################
continents.eu <- c("UK","Belgium","France","Germany","Italy","Netherlands","Poland","Portugal","Spain")
continents.na <- c("USA","Mexico","Canada")
continents.sa <- c("Argentina","Chile","Brazil","Colombia","Peru")
continents.asia <- c("India","Japan","Philippines","Singapore","Pakistan")

comp.multiews.diff2 <- rbind(uk.multiews.diff2,ind.multiews.diff2,chl.multiews.diff2,
                             bzl.multiews.diff2,usa.multiews.diff2,sing.multiews.diff2,
                             itl.multiews.diff2,jpn.multiews.diff2,esp.multiews.diff2,
                             ger.multiews.diff2,fr.multiews.diff2,sa.multiews.diff2,
                             isr.multiews.diff2,arg.multiews.diff2,mex.multiews.diff2,
                             ned.multiews.diff2,por.multiews.diff2,bel.multiews.diff2,
                             phl.multiews.diff2,pol.multiews.diff2,can.multiews.diff2,
                             col.multiews.diff2,pk.multiews.diff2,per.multiews.diff2)%>%
  filter(period != "fourth") #create EWS-deriv difference data frame and drop the lone fourth wave in Japan
write.csv(comp.multiews.diff2,row.names = FALSE,file="Results/all.multiews.2.csv")
comp.multiews.diff2 <- read.csv(file="Results/all.multiews.2.csv")

region.comp2 <- comp.multiews.diff2 %>% #categorise df into continents 
  mutate(region = ifelse(data.source %in% continents.eu,"Europe",
                         ifelse(data.source %in% continents.na,"N.America",
                                ifelse(data.source %in% continents.sa,"S.America",
                                       ifelse(data.source %in% continents.asia,"Asia","Other"))))) 

median.diff2 <- comp.multiews.diff2 %>% #estimate overall median preemption by EWS 
  mutate(day.prior = as.integer(day.prior,units="days"), #convert from diffdays to numeric
         day.prior = ifelse(abs(day.prior) > 60, NA,day.prior))%>% #exclude estimates > 60 or <-60 days as outliers
  group_by(metric.code,period)%>%
  summarise(mean = mean(day.prior,na.rm = T),
            se = sd(day.prior,na.rm = T)/n(),n=n(),
            n_missed = sum(is.na(day.prior)),
            n_prior = sum(day.prior > 0,na.rm = T))%>%
  mutate(prop_missed = n_missed/n,prop_prior = n_prior/n)
write.table(median.diff2, file = "Results/median_diff2.txt", sep = ",", quote = FALSE, row.names = F)

median.region2 <- region.comp2 %>% #estimate continent-wise median preemption by EWS 
  mutate(day.prior = as.integer(day.prior,units="days"), #convert from diffdays to numeric
         day.prior = ifelse(abs(day.prior) > 60, NA,day.prior))%>% #exclude estimates > 60 or <-60 days as outliers
  group_by(metric.code,period,region)%>%
  summarise(mean = mean(day.prior,na.rm = T),
            se = sd(day.prior,na.rm = T)/n(),n=n(),
            n_missed = sum(is.na(day.prior)),
            n_prior = sum(day.prior > 0,na.rm = T))%>%
  mutate(prop_missed = n_missed/n,prop_prior = n_prior/n)
write.table(median.region2, file = "Results/median_region2.txt", sep = ",", quote = FALSE, row.names = F)


region.scatter2 <- ggplot(data = region.comp2 %>% filter(data.source != "USA"|period !="third") %>% mutate(miss.fac = ifelse(prediction %in% c("miss","unknown"), "Undetected","Detected")), 
                          aes(x = period, y = as.numeric(day.prior))) + 
  geom_point(aes(col=data.source,shape = miss.fac,size = miss.fac),position=position_dodge(width=0.5),alpha=0.5)+
  theme_bw()+
  geom_hline(yintercept = 0,linetype="dashed", color = "black")+
  ggtitle("Two Consecutive EWS") + 
  facet_grid(region~metric.code,scales = "fixed")+
  xlab("Wave") + ylab("Days prior to non-linear case increase") + labs(col='Country',shape="Detection") +
  scale_y_continuous(limits=c(-100,90),breaks = seq(90,-100,by=-60)) +
  scale_x_discrete(labels =c("1", "2","3"))+
  scale_shape_manual(values = c(16, NA),name = "Warning")+
  scale_size_manual(values = c(3, 2))+ 
  guides(size = "none",
         colour = guide_legend(override.aes = list(size = 4),ncol=1),
         shape = guide_legend(override.aes = list(size=c(4,4)))) +
  theme(axis.text = element_text(size = ceiling(12 * 0.7), colour = "black"),
        axis.title = element_text(size = ceiling(12 * 0.8)),
        legend.text = element_text(size = ceiling(12 * 0.9), family = "sans"),
        legend.title = element_text(size = 12,face = "bold",family = "sans"),
        legend.position = "right",
        legend.key = element_rect(fill = "white", colour = NA),
        legend.background = element_rect(colour = "black"),
        plot.title = element_text(size = ceiling(12 * 1.1), face = "bold"),
        panel.grid.major.y = element_line(colour = "gray", linetype = "dotted"),
        panel.grid.major.x = element_line(colour = "gray", linetype = "dotted"),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill="white"))

png(file="Results/region_comparison_2.png",
    width=3000, height = 2500 ,res=300) #visualise scatter of EWS preemption days by continent
region.scatter2
dev.off()
