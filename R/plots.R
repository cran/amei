`plot.epiman` <-
function (x, type = c("epi", "costs", "params", "fracs", "stops"), 
    showd = FALSE, main = NULL, true = NULL, ...) 
{
    type <- match.arg(type)
    if (type == "epi") {
        if (is.null(main)) 
            main <- "Evolution of Epidemic"
        PlotEpi(x$soln, showd = showd, main = main)
    }
    else if (type == "costs") {
        if (is.null(main)) 
            main <- "Evolution of Costs"
        PlotCosts(x$soln, main = main)
    }
    else if (type == "params") {
        if (is.null(main)) 
            main <- "MCMC Inference"
        PlotParams(x$samp, NULL, true)
    }
    else if (type == "fracs") {
        ylab <- "vaccination fraction"
        TimeSeriesOfDensities(x$vachist$fracs, x$vactimes, 
            c(-0.1, 1.1), ylab)
        if (!is.null(main)) 
            title(main)
    }
    else {
        ylim <- "stop number"
        TimeSeriesOfDensities(x$vachist$stops, x$vactimes, 
            c(0, x$soln$S[1]), ylim)
        if (!is.null(main)) 
            title(main)
    }
}



`plot.optvac` <-
function (x, main = NULL, ...) 
{
    if (is.null(main)) 
        main <- "Optimal vaccination policy surface"
    image(x$vacgrid$fracs, x$vacgrid$stops, x$C, main = main, 
        xlab = "fraction", ylab = "stop number", ...)
    grid(length(x$vacgrid$fracs), length(x$vacgrid$stops), lty = 1, 
        col = "black")
    best <- getpolicy(x)
    worst <- getpolicy(x, "worst")
    text(best$frac, best$stop, best$cost)
    text(worst$frac, worst$stop, worst$cost)
}


`plot.MCepi` <-
function (x, type = c("epi", "costs", "fracs", "stops"), showd = FALSE, 
    showv = FALSE, main = NULL, ...) 
{
  type <- match.arg(type)
  if (type == "epi") {
    if (is.null(main)) main <- "Monte Carlo Epidemics"
    PlotEpi(x$Median, showd = showd, showv = showv, main = main, 
            ...)
    PlotEpi(x$Q1, add = TRUE, showd = showd, showv = showv)
    PlotEpi(x$Q3, add = TRUE, showd = showd, showv = showv)
  }
  else if (type == "costs") {
    if (is.null(main)) main <- "Monte Carlo Costs"
    PlotCosts(x$Median, ylim = c(min(x$Q1$C), max(x$Q3$C)), 
              main = main, ...)
    PlotCosts(x$Q1, add = TRUE)
    PlotCosts(x$Q3, add = TRUE)
  }
  else if (type == "fracs") {
    if (is.null(main)) main <- "Monte Carlo Fraction Vaccinated"
    plot(x$Median$frac, type = "l", lty = 1, lwd = 2, 
         xlab = "time", ylab = "fraction",
         ylim = c(min(x$Q1$frac), max(x$Q3$frac)), main = main, ...)
    lines(x$Q1$frac, lty = 2, lwd = 2)
    lines(x$Q3$frac, lty = 2, lwd = 2)
  }
  else {
    if (is.null(main)) main <- "Monte Carlo Stopping Threshold"
    plot(x$Median$stop, type = "l", lty = 1, lwd = 2, 
         xlab = "time", ylab = "stop time",
         ylim = c(min(x$Q1$stop), max(x$Q3$stop)), main=main, ...)
    lines(x$Q1$stop, lty = 2, lwd = 2)
    lines(x$Q3$stop, lty = 2, lwd = 2)
  }
}
