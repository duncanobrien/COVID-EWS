## Curve Direction Fn ##

# taken from Burthe, S.J., Henrys, P.A., Mackay, E.B., Spears, B.M., 
# Campbell, R., Carvalho, L., Dudley, B., Gunn, I.D.M., Johns, D.G., Maberly, 
# S.C., May, L., Newell, M.A., Wanless, S., Winfield, I.J., Thackeray, S.J. 
# and Daunt, F. (2016), Do early warning indicators consistently predict 
# nonlinear change in long-term ecological data?. J Appl Ecol, 53: 666-676. 
# https://doi.org/10.1111/1365-2664.12519

curve.direction <- function(x=0, L.CI, U.CI){
  pos <- ifelse(x<U.CI, 1, -1) 
  negs<-ifelse(x>L.CI, -1, 1)
  return(pos+negs)}