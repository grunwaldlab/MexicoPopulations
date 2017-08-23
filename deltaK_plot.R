deltaK_plot <- function (sr, plot = TRUE) {
  
  k.tbl <- table(sapply(sr, function(x) x$summary["k"]))
  
  sr.smry <- t(sapply(sr, function(x) x$summary))
  ln.k <- tapply(sr.smry[, "est.ln.prob"], sr.smry[, "k"], 
                 mean)
  sd.ln.k <- tapply(sr.smry[, "est.ln.prob"], sr.smry[, "k"], 
                    sd)
  ln.pk <- diff(ln.k)
  ln.ppk <- abs(diff(ln.pk))
  delta.k <- sapply(2:(length(ln.k) - 1), function(i) {
    abs(ln.k[i + 1] - (2 * ln.k[i]) + ln.k[i - 1])/sd.ln.k[i]
  })
  result <- data.frame(k = as.numeric(names(ln.k)), reps = as.numeric(table(sr.smry[, 
                                                                                    "k"])), mean.ln.k = as.numeric(ln.k), sd.ln.k = as.numeric(sd.ln.k), 
                       ln.pk = c(NA, ln.pk), ln.ppk = c(NA, ln.ppk, NA), delta.k = c(NA, 
                                                                                     delta.k, NA))

  rownames(result) <- NULL
  if (plot) {
    xlim <- range(c(0, result$k))
    tick.spacing <- 1
    xticks <- seq(0, max(result$k), tick.spacing)
   
    #layout(matrix(1, nrow = 1, byrow = TRUE))
    

    plot.func <- function(y, sd = NULL) {
      ylim <- if (is.null(sd)) {
        range(y, na.rm = TRUE)
      }
      else {
        sd[is.na(sd)] <- 0
        range(c(y, y + sd, y - sd), na.rm = TRUE)
      }
      plot(result$k, y, xlim = xlim, ylim = ylim, type = "b", 
           xlab = "K", ylab = "Delta(K)", axes = T, bty = "o", 
           pch = 19, bg = "black", cex.lab = 1.5, cex.axis = 1.5, font.lab = 2, 
           font = 2, family = "Microsoft Sans Serif", ps =50)
      
    }
    
    plot.func(result$delta.k)
    #tiff(filename = "evanno.tiff", width = 3.2, height = 2.5, units = "in", res = 600) 

  }
  result
  

  
}


