###########################################################################################
### Figure 2 - GLobal EWS Pre-emption ### 
###########################################################################################

### Preamble ###
require(tidyverse) #dplyr, ggplot etc
require(ggthemes) #theme_clean
require(scales) #plot scales functions
require(egg) #ggarrange

comp.multiews.diff2 <- read.csv(file="Results/all.multiews.2.csv")
#load in difference data between EWS detection and non-linear increases estimated from derivatives

continents.eu <- c("UK","Belgium","France","Germany","Italy","Netherlands","Poland","Portugal","Spain")
continents.na <- c("USA","Mexico","Canada")
continents.sa <- c("Argentina","Chile","Brazil","Colombia","Peru")
continents.asia <- c("India","Japan","Philippines","Singapore","Pakistan")
#classify continents

region.comp2 <- comp.multiews.diff2 %>% #categorise df into continents 
  mutate(region = ifelse(data.source %in% continents.eu,"Europe",
                         ifelse(data.source %in% continents.na,"N.America",
                                ifelse(data.source %in% continents.sa,"S.America",
                                       ifelse(data.source %in% continents.asia,"Asia","Other"))))) 

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

pdf(file="Results/region_comparison_2.pdf",
    width=10, height = 8,onefile = F) #visualise scatter of EWS preemption days by continent
region.scatter2
dev.off()
