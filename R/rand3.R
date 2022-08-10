rand3 <-  function(x, na.rm = FALSE){x <- ifelse(x<0.001,0.001,x)
round(runif(n=1, x*0.97,x*1.03),4)}