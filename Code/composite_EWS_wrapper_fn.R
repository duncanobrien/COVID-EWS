### Composite EWS Wrapper fn ###

comp_EWS_wrapper <- function(data,metrics,trait = NULL, threshold = 2,plotIt = F, tail.direction = "one.tailed", burn_in = 5, ggplotIt = T,
                             y_lab = "Generic Indicator Name", trait_lab = "Generic Trait Name",
                             trait_scale = 100000, interpolate = F,data_source = NULL,method = c("w_comp","dakos"), output="raw"){
require(gtools) #"combinations" fn
require(ggplot2)
require(ggthemes) #"theme_clean()" fn
require(egg) #"ggarrange" fn
require(scales) #"pretty_breaks" fn
source("Code/w_composite_ews_fn.R")

  method <- match.arg(method)
  
  if(method == "w_comp"){
   to.test.l<-list()
  for(jj in 1:length(metrics)){
    #jj=1
    to.test.l[[jj]]<-split(combinations(n = length(metrics), r = jj, v = metrics, repeats.allowed = FALSE), seq(nrow(combinations(n = length(metrics), r = jj, v = metrics, repeats.allowed = FALSE))))
  }
  to.test<-unlist(to.test.l, recursive=FALSE)
  
  ##object to store results
  res<-NULL
  
  ##loop through all the metrics
  ##object to store results
  for(i in 1:length(to.test)){
    #i=10
    ##set the weighing to 1 - just required to make the code run, doesnt do anythnig
    W<-data.frame(inds=sort(unlist(to.test[i])), "wei"=1)
    ##run the EWS from clements & ozgul nat comms and save out:
    res[[i]]<-W_composite_ews(dat=data, indicators=sort(unlist(to.test[i])), weights=W, trait = trait, threshold=threshold, burn_in = burn_in, tail.direction = tail.direction, plotIt=plotIt, interpolate = interpolate)
  }
  
  bind.res<-rbindlist(res)
  bind.res$str<-(bind.res$metric.score-bind.res$rolling.mean)/bind.res$rolling.sd
  bind.res$data_source <- data_source
  bind.res<-as.data.table(bind.res)
  
 if(ggplotIt == TRUE){ 
      p<-ggplot(data = bind.res, aes(x=time,y=str)) +geom_line(data =bind.res,aes(col= metric.code))+
    geom_point(data =bind.res[which(bind.res$threshold.crossed==1)], aes(x=time, y = str, col = metric.code)) +
        scale_x_continuous(breaks = scales::pretty_breaks(n = 6)) +
        theme_clean() + xlab("Date") + ylab("Strength of EWS") +
    theme(plot.margin = margin(c(10, 8, 0, 10)),
          legend.key.height = unit(0.3,"cm" ),
          legend.key.width = unit(0.1,"cm"),
          legend.title = element_text(size = 10),
          legend.text = element_text(size=8))
   
  if(is.null(trait)==F){
      plot.dat<-data.frame("timeseries"=bind.res$time[1:(nrow(data)-burn_in)], "value"=bind.res$count.used[1:(nrow(data)-burn_in)],"trait"=trait[burn_in:(nrow(data)-1)])   
      p2 <- ggplot(data = plot.dat, aes(x=timeseries)) +
        geom_line(aes(y=value),linetype=1) + 
        #geom_point(aes(y=(trait*100000)), size =1, shape = 1, col = "blue")+
        geom_line(aes(y=(trait*trait_scale)),linetype=2, size = 0.4, alpha = 0.4,col = "blue") +
        #ylab(y_lab) + 
        scale_y_reverse(y_lab,sec.axis = sec_axis(~./trait_scale, name = trait_lab))+
        theme(axis.title.x=element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank(), 
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              plot.margin = margin(c(0, 0, 10, 10)),
              axis.title.y=element_text(size=10),
              axis.title.y.right = element_text(color = "royalblue1", size=10)) 
      p3 <- ggarrange(p,p2,nrow = 2,heights = c(2, 1))
      print(p3)
  
  }else if(is.null(trait)==T){
    plot.dat<-data.frame("timeseries"=as.numeric(bind.res$time), "value"=bind.res$count.used)   
    p4 <- ggplot(data = plot.dat, aes(x=timeseries, y=value)) +
      geom_line() +
      ylab(y_lab) + 
      theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(), 
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            plot.margin = margin(c(0, 0, 10, 10)),
            axis.title.y=element_text(size=10))+
      scale_y_reverse()
      p5 <- ggarrange(p,p4,nrow = 2,heights = c(2, 1))
      print(p5)
  }
  
}
  }
if(method == "dakos"){
  
  bind.res <- no.plot.ews(timeseries = data, output = output,interpolate = T)
  bind.res$data_source <- data_source
  bind.res<-as.data.table(bind.res)
}  
  
return(bind.res)

}
   