###########################################################################################
### Supplementary Figure 2 - Continuous UK COVID-19 EWS ### 
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

### Perform EWS assessment ###
metrics <- c("skew","SD","acf") # metrics of interest (Skewness, variance and autocorrelation)

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
