plotroc <- function(A, est1, est2, others = NULL, xmax = NA, mytitle = "", ylimit = c(0.5,1)) {
  glasso = ROC(A, est1)
  bar = ROC(A, est2)
  glasso$mycolors = 1
  bar$mycolors = 2
  combined = rbind(glasso, bar)
  
  if (is.na(xmax)) {
    xrange = c(0, max(combined$falsePos, na.rm = T))
  } else {
    xrange = c(0, xmax)
    bar = bar[bar$falsePos <= xmax+.02, ]
    combined = combined[combined$falsePos <= xmax+.02, ]
  }
  
  combined = combined[combined$truePos > ylimit[1],]
  bar = bar[bar$truePos > ylimit[1],]
  #png("./pngs/roc.png", width = 720, height = 720, units = "px", pointsize = 18)
  #plot(combined$falsePos[combined$mycolors == 1], combined$truePos[combined$mycolors == 1], 
  plot(glasso$falsePos, glasso$truePos, 
       xlim = xrange, ylim = ylimit, type = "l", col = "black", lwd = 3, 
       xlab="False Positive Rate", ylab="True Positive Rate", main = mytitle)
  lines(bar$falsePos,bar$truePos,
        col = "red", lwd = 3, ylim = ylimit, lty = 1, xlim = xrange)
  if(!is.null(others)) {
    for (i in 1:length(others)) {
      temp = others[[i]]
      temp = ROC(A, temp$network)
      lines(temp$falsePos,temp$truePos,
            col = "black", lwd = 3, lty = i+1, ylim = ylimit, xlim = xrange)
    }
  }
  # legend((xrange[2] - xrange[1])*.49, 0.71, 
  #        legend=c("BAR", "LPGM th=0", "LPGM th=0.005", "LPGM th=0.001", "LPGM th=0.0001"),
  #         col=c("red", rep("black", 4)), lty=c(1, 1:4), cex=0.8)
  # legend(c(0.35, 0.6), c(0.4, 0.6), legend=c("BAR", "LPGM th=0", "LPGM th=0.005", "LPGM th=0.001", "LPGM th=0.0001"),
  #        col=c("red", rep("black", 4)), lwd = 2, lty = 1:4, cex=0.9, inset = 0.03, adj = c(0), bty = "n")
  # legend("topright", legend=c("BAR", "LPGM th=0", "LPGM th=0.005", "LPGM th=0.001", "LPGM th=0.0001"),
  #        col=c("red", rep("black", 4)), lwd = 2, lty = 1:4, cex=0.9, inset = 0, bty = "n")
  legend("bottomright", legend=c("BAR", "LPGM th=0", "LPGM th=0.005", "LPGM th=0.001", "LPGM th=0.0001"),
         col=c("red", rep("black", 4)), lty=c(1, 1:4), cex=0.8, inset = 0.03)
}
