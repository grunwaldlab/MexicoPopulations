scatterdapc <- function (x, xax = 1, yax = 2, grp = x$grp, col = seasun(length(levels(grp))), 
                               pch = 20, bg = "white", solid = 0.7, scree.da = TRUE, scree.pca = FALSE, 
                               posi.da = "bottomright", posi.pca = "bottomleft", bg.inset = "white", 
                               ratio.da = 0.25, ratio.pca = 0.3, inset.da = 0.02, inset.pca = 0.01, 
                               inset.solid = 0.5, onedim.filled = TRUE, mstree = FALSE, 
                               lwd = 1, lty = 1, segcol = "black", legend = FALSE, posi.leg = "topright", 
                               cleg = 1, txt.leg = levels(grp), cstar = 1, cellipse = 1.5, 
                               axesell = FALSE, label = levels(grp), clabel = 1, xlim = NULL, 
                               ylim = NULL, grid = F, addaxes = T, origin = c(0, 
                                                                              0), include.origin = T, sub = "", csub = 1, possub = "bottomleft", 
                               cgrid = 1, pixmap = NULL, contour = NULL, area = NULL, label.inds = NULL, 
                               ...) 
{
  ONEDIM <- xax == yax | ncol(x$ind.coord) == 1
  col <- rep(col, length(levels(grp)))
  pch <- rep(pch, length(levels(grp)))
  col <- transp(col, solid)
  bg.inset <- transp(bg.inset, inset.solid)
  if (is.null(grp)) {
    grp <- x$grp
  }
  if (is.null(xax) || is.null(yax)) {
    xax <- 1L
    yax <- ifelse(ncol(x$ind.coord) == 1L, 1L, 2L)
    ONEDIM <- TRUE
  }
  if (!ONEDIM) {
    opar <- par(mar = par("mar"))
    par(mar = c(0.1, 0.1, 0.1, 0.1), bg = bg)
    on.exit(par(opar))
    axes <- c(xax, yax)
    s.class(x$ind.coord[, axes], fac = grp, col = col, cpoint = 0, 
            cstar = cstar, cellipse = cellipse, axesell = axesell, 
            label = label, clabel = clabel, xlim = xlim, ylim = ylim, 
            grid = grid, addaxes = addaxes, origin = origin, 
            include.origin = include.origin, sub = sub, csub = csub, 
            possub = possub, cgrid = cgrid, pixmap = pixmap, 
            contour = contour, area = area)
    colfac <- pchfac <- grp
    levels(colfac) <- col
    levels(pchfac) <- pch
    colfac <- as.character(colfac)
    pchfac <- as.character(pchfac)
    if (is.numeric(col)) 
      colfac <- as.numeric(colfac)
    if (is.numeric(pch)) 
      pchfac <- as.numeric(pchfac)
    points(x$ind.coord[, xax], x$ind.coord[, yax], col = colfac, 
           pch = pchfac, ...)
    s.class(x$ind.coord[, axes], fac = grp, col = col, cpoint = 0, 
            add.plot = TRUE, cstar = cstar, cellipse = cellipse, 
            axesell = axesell, label = label, clabel = clabel, 
            xlim = xlim, ylim = ylim, grid = grid, addaxes = addaxes, 
            origin = origin, include.origin = include.origin, 
            sub = sub, csub = csub, possub = possub, cgrid = cgrid, 
            pixmap = pixmap, contour = contour, area = area)
    if (!is.null(label.inds) & is.list(label.inds)) {
      appendList <- function(x, val) {
        stopifnot(is.list(x), is.list(val))
        xnames <- names(x)
        for (v in names(val)) {
          x[[v]] <- if (v %in% xnames && is.list(x[[v]]) && 
                        is.list(val[[v]])) 
            appendList(x[[v]], val[[v]])
          else c(x[[v]], val[[v]])
        }
        x
      }
      do.call("orditorp", c(appendList(list(x = x$ind.coord[, 
                                                            c(xax, yax)], display = "species"), label.inds)))
    }
    if (mstree) {
      meanposi <- apply(x$tab, 2, tapply, grp, mean)
      D <- dist(meanposi)^2
      tre <- ade4::mstree(D)
      x0 <- x$grp.coord[tre[, 1], axes[1]]
      y0 <- x$grp.coord[tre[, 1], axes[2]]
      x1 <- x$grp.coord[tre[, 2], axes[1]]
      y1 <- x$grp.coord[tre[, 2], axes[2]]
      segments(x0, y0, x1, y1, lwd = lwd, lty = lty, col = segcol)
    }
  }
  else {
    scree.da <- FALSE
    if (ncol(x$ind.coord) == 1) {
      pcLab <- 1
    }
    else {
      pcLab <- xax
    }
    ldens <- tapply(x$ind.coord[, pcLab], grp, density)
    allx <- unlist(lapply(ldens, function(e) e$x))
    ally <- unlist(lapply(ldens, function(e) e$y))
    par(bg = bg)
    plot(allx, ally, type = "n", xlab = paste("Discriminant function", 
                                              pcLab), ylab = "Density")
    for (i in 1:length(ldens)) {
      if (!onedim.filled) {
        lines(ldens[[i]]$x, ldens[[i]]$y, col = col[i], 
              lwd = 2)
      }
      else {
        polygon(c(ldens[[i]]$x, rev(ldens[[i]]$x)), c(ldens[[i]]$y, 
                                                      rep(0, length(ldens[[i]]$x))), col = col[i], 
                lwd = 2, border = col[i])
      }
      points(x = x$ind.coord[grp == levels(grp)[i], pcLab], 
             y = rep(0, sum(grp == levels(grp)[i])), pch = "|", 
             col = col[i])
    }
  }
  if (legend) {
    temp <- list(...)$cex
    if (is.null(temp)) 
      temp <- 1
    if (ONEDIM | temp < 0.5 | all(pch == "")) {
      legend(posi.leg, fill = col, legend = txt.leg, cex = cleg, 
             bg = bg.inset)
    }
    else {
      legend(posi.leg, col = col, legend = txt.leg, cex = cleg, 
             bg = bg.inset, pch = pch, pt.cex = temp)
    }
  }
  if (scree.da && ratio.da > 0.01) {
    inset <- function() {
      myCol <- rep("white", length(x$eig))
      myCol[1:x$n.da] <- "grey"
      myCol[c(xax, yax)] <- "black"
      myCol <- transp(myCol, inset.solid)
      barplot(x$eig, col = myCol, xaxt = "n", yaxt = "n", 
              ylim = c(0, x$eig[1] * 1.1))
      mtext(side = 3, "DA eigenvalues", line = -1.2, adj = 0.8)
      box()
    }
    add.scatter(inset(), posi = posi.da, ratio = ratio.da, 
                bg.col = bg.inset, inset = inset.da)
  }
  if (scree.pca && !is.null(x$pca.eig) && ratio.pca > 0.01) {
    inset <- function() {
      temp <- 100 * cumsum(x$pca.eig)/sum(x$pca.eig)
      myCol <- rep(c("black", "grey"), c(x$n.pca, length(x$pca.eig)))
      myCol <- transp(myCol, inset.solid)
      plot(temp, col = myCol, ylim = c(0, 115), type = "h", 
           xaxt = "n", yaxt = "n", xlab = "", ylab = "", 
           lwd = 2)
      mtext(side = 3, "Eigenvalues", line = -1.2, adj = 0.1, cex = 0.7)
    }
    add.scatter(inset(), posi = posi.pca, ratio = ratio.pca, 
                bg.col = bg.inset, inset = inset.pca)
  }
  return(invisible(match.call()))
}

# scatter.dapc.test(P.infx$DAPC, col=myCol, clabel = 0.75, pch=15:19, scree.pca = TRUE, scree.da = FALSE, 
#                   posi.pca = "bottomright", posi.leg = "topright", legend = F, 
#                   cleg = 0.9, inset.solid = 1, xax = 1, yax = 2, solid = 1, cstar = 0)
