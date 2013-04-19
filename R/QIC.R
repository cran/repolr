QIC <-
function(x){
 qlike <- x$gee$y*x$gee$linear.predictor+log(1-x$gee$fitted.values)
 QIC <- -2*sum(qlike)+2*length(x$gee$coeff)
 QIC
}
