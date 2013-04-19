summary.gee <-
function (object, working.correlation = FALSE, ...) 
{
    coef <- object$coefficients
    resid <- object$residuals
    n <- length(resid)
    p <- object$rank
    if (is.null(p)) 
        p <- sum(!is.na(coef))
    if (!p) {
        warning("This model has zero rank --- no summary is provided")
        return(object)
    }
    nas <- is.na(coef)
    cnames <- names(coef[!nas])
    coef <- matrix(rep(coef[!nas], 5), ncol = 5)
    dimnames(coef) <- list(cnames, c("Estimate", "Naive S.E.", 
        "Naive z", "Robust S.E.", "Robust z"))
    rse <- sqrt(diag(object$robust.variance))
    nse <- sqrt(diag(object$naive.variance))
    coef[, 2] <- nse
    coef[, 3] <- coef[, 1]/coef[, 2]
    coef[, 4] <- rse
    coef[, 5] <- coef[, 1]/coef[, 4]
    summary <- list()
    summary$call <- object$call
    summary$version <- object$version
    summary$nobs <- object$nobs
    summary$residual.summary <- quantile(as.vector(object$residuals))
    names(summary$residual.summary) <- c("Min", "1Q", "Median", 
        "3Q", "Max")
    summary$model <- object$model
    summary$title <- object$title
    summary$coefficients <- coef
    summary$scale <- object$scale
    summary$error <- paste("Error code was", object$error)
    summary$iterations <- object$iterations
    if(working.correlation==TRUE){summary$working.correlation <- object$working.correlation}
    attr(summary, "class") <- "summary.gee"
    summary
}
