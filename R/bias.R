#
# bias.R
# created 2007-12-05
# Aron Eklund
#
#
# 


getBiasMetrics <- function(x.batch, x.rma = rma(x.batch)) {
  x.pm <- log2(pm(x.batch))
  message('Calculating degradation scores...')
  degradation <- affy::AffyRNAdeg(x.batch)$slope
  data.frame(
          pm.median   = apply(x.pm, 2, median),
          pm.IQR      = apply(x.pm, 2, IQR),
          rma.IQR     = apply(exprs(x.rma), 2, IQR),
          degradation = degradation
  )
}


getBiasMetrics2 <- function(x.batch, x.rma = rma(x.batch)) {
  ## first do the four original bias metrics
  x.pm <- log2(pm(x.batch))
  pm.median   <- apply(x.pm, 2, median)
  pm.IQR      <- apply(x.pm, 2, IQR)
  rma.IQR     <- apply(exprs(x.rma), 2, IQR)
  message('Calculating degradation scores...')
  degradation <- affy::AffyRNAdeg(x.batch)$slope
  
  ## now the "extended" bias metrics
  message('Calculating MAS5 calls...')
  present.calls <- try(apply(exprs(mas5calls(x.batch)) == 'P', 2, mean), silent = TRUE)
  if(class(present.calls) == "try-error") {
    present.calls <- rep(as.numeric(NA), length(pm.median))
  }  
  
  ## probe-specific bias metrics 
  ## (all probe sets may not be present on all arrays -- so we use safety functions)  
  availableProbes <- featureNames(x.rma)
  getMeans <- function(probes) {
    ok <- intersect(probes, availableProbes)
    colMeans(exprs(x.rma)[ok, , drop=FALSE])
  }

  pr.spikes.mRNA <- c("AFFX-LysX-3_at", "AFFX-LysX-5_at", "AFFX-LysX-M_at",             
                   "AFFX-PheX-3_at", "AFFX-PheX-5_at", "AFFX-PheX-M_at", 
                   "AFFX-ThrX-3_at", "AFFX-ThrX-5_at", "AFFX-ThrX-M_at",
                   "AFFX-DapX-3_at", "AFFX-DapX-5_at", "AFFX-DapX-M_at"  )
  pr.spikes.cRNA <- c("AFFX-BioB-3_at", "AFFX-BioB-5_at", "AFFX-BioB-M_at",
                   "AFFX-BioC-3_at", "AFFX-BioC-5_at", 
                   "AFFX-BioDn-3_at", "AFFX-BioDn-5_at",
                   "AFFX-CreX-3_at", "AFFX-CreX-5_at"   )
  pr.rRNA <- c("AFFX-HUMRGE/M10098_3_at", "AFFX-HUMRGE/M10098_5_at", 
            "AFFX-HUMRGE/M10098_M_at", "AFFX-M27830_5_at", "AFFX-M27830_M_at")
  pr.alu <- "AFFX-hum_alu_at"
  
  rma.spikes.mRNA <- getMeans(pr.spikes.mRNA)
  rma.spikes.cRNA <- getMeans(pr.spikes.cRNA)
  rma.rRNA <- getMeans(pr.rRNA)
  rma.alu <- getMeans(pr.alu)
  
  border.plus  <- sapply(borders(x.batch), function(x) median(unlist(x$plus)))
  border.minus <- sapply(borders(x.batch), function(x) median(unlist(x$minus)))

  data.frame( 
    pm.median, pm.IQR, rma.IQR, degradation,
    present.calls, 
    rma.spikes.mRNA, rma.spikes.cRNA, 
    rma.rRNA, rma.alu,
    border.plus, border.minus
  )
}


borders <- function(x.batch) {
  nr <- nrow(x.batch)
  nc <- ncol(x.batch)
  xy2i <- function(x, y) x + (nr * (y - 1))
  sb <- (nr %% 2)  # sign on bottom
  sr <- (nc %% 2)  # sign on right
  index.plus <- list( bottom = xy2i(seq(1 + sb, nc, by=2), nr),
                        left = xy2i(1, seq(1, nr, by=2)),
                         top = xy2i(seq(1, nc, by=2), 1), 
                       right = xy2i(nc, seq(1 + sr, nr, by=2))  )
  index.minus <- list(bottom = xy2i(seq(2 - sb, nc, by=2), nr),
                        left = xy2i(1, seq(2, nr, by=2)),
                         top = xy2i(seq(2, nc, by=2), 1), 
                       right = xy2i(nc, seq(2 - sr, nr, by=2))  )
  out <- lapply(1:ncol(exprs(x.batch)), function(j) {
    plus  <- lapply(index.plus,  function(i) exprs(x.batch)[i, j])    
    minus <- lapply(index.minus, function(i) exprs(x.batch)[i, j])    
    list(plus = plus, minus = minus)
  })
  names(out) <- sampleNames(x.batch)
  out
}


biasCorrection <- function(x, metrics) {
  if( class(x) %in% c('exprSet', 'ExpressionSet')) {
    in.mat <- Biobase::exprs(x)
  } else {
    in.mat <- x
  }
  fit <- lm(t(in.mat) ~ ., data = metrics)
  res <- t(residuals(fit))
  fixed <- sweep(res, 1, rowMeans(in.mat), FUN = '+')
  if( class(x) %in% c('exprSet', 'ExpressionSet')) {
    out <- x
    Biobase::exprs(out) <- fixed
  } else {
    out <- fixed
  }
  return(out)
}



getScanDate <- function(filenames = list.celfiles()) {
  DatHeaders <- rep('', length(filenames))
  for (i in seq(along = filenames)) {
    con <- gzfile(filenames[i])
    DatHeaders[i] <- grep('DatHeader', readLines(con, n = 20), value = TRUE)
    close(con)
  }
  dateInUSAformat <- sub("^.*(../../..).*$", "\\1", DatHeaders, perl = TRUE)
  dateInDateFormat <- as.Date(dateInUSAformat, "%m/%d/%y")
  names(dateInDateFormat) <- filenames
  dateInDateFormat
}


calcCM <- function(x, n = 50, 
       use = 'all', method = 'pearson', 
       FUN = median, ...) {
  if( class(x) %in% c('exprSet', 'ExpressionSet')) 
    x <- Biobase::exprs(x)
  x <- t(as.matrix(x))
  b <- round(seq(0, ncol(x), length = n + 1))
  z <- matrix(NA, nrow = n, ncol = n)
  count <- matrix(NA, nrow = n, ncol = n)
  for (i in 1:n) {
    which.i <- (b[i] + 1):b[i + 1]
    for (j in i:n) {
      which.j <- (b[j] + 1):b[j + 1]
      myCor <- cor( x[,which.i], x[, which.j], 
                    use = use, method = method )
      if( i == j ) { 
        myCor <- myCor[lower.tri(myCor)]
        z[i,j] <- FUN(myCor, ...) 
        count[i,j] <- length(myCor) 
      } else {
        z[i,j] <- z[j,i] <- FUN(myCor, ...)
        count[i,j] <- count[j,i] <- length(myCor)
      }
    }
  }
  return(list(x = b, y = b, z = z, count = count, 
    call = paste( as.character(match.call()), collapse=" " ) ))
}



quantileNormalize <- function(x, margin = 2) {
  v <- rowMeans(apply(x, margin, sort))
  x.qn <- apply(x, margin, function(y) v[rank(y, ties.method = 'random')] )
  rownames(x.qn) <- rownames(x)
  return(x.qn)
}


batchCorrection <- function(x, b) {
    stopifnot( all( ! is.na(b) ) )
    stopifnot(length(b) == ncol(x))
    if (class(x) %in% c("exprSet", "ExpressionSet")) {
        in.mat <- Biobase::exprs(x)
    }
    else {
        in.mat <- x
    }
    b <- factor(b)
    batchMeans <- sapply(levels(b), function(a) 
                     rowMeans(in.mat[,b == a, drop = FALSE]) )
    adj <- batchMeans - rowMeans(in.mat)
    fixed <- in.mat - adj[,b]
    if (class(x) %in% c("exprSet", "ExpressionSet")) {
        out <- x
        Biobase::exprs(out) <- fixed
    }
    else {
        out <- fixed
    }
    out
}
