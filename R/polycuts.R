polycuts <-
function(x){
 mat1 <- x$corr$x[1:(x$corr$categories-1),1:(x$corr$poly+1)]
 mat2 <- x$gee$coeff[1:(x$corr$poly+1)]
 cuts <- mat1%*%mat2
 cuts
}
