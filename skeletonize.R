thin_zhang_suen <- function(bin) {
  bin <- bin > 0
  pad0 <- function(m) rbind(FALSE, cbind(FALSE, m, FALSE), FALSE)
  sh <- function(m, dy, dx) {
    M <- pad0(m); nr <- nrow(M); nc <- ncol(M)
    yy <- (2 + dy):(nr - 1 + dy); xx <- (2 + dx):(nc - 1 + dx)
    M[yy, xx, drop = FALSE]
  }
  changed <- TRUE; m <- bin
  while (changed) {
    changed <- FALSE
    p2 <- sh(m,-1,0); p3 <- sh(m,-1,1); p4 <- sh(m,0,1); p5 <- sh(m,1,1)
    p6 <- sh(m,1,0);  p7 <- sh(m,1,-1); p8 <- sh(m,0,-1); p9 <- sh(m,-1,-1)
    N  <- p2+p3+p4+p5+p6+p7+p8+p9
    trans <- (!p2&p3)+(!p3&p4)+(!p4&p5)+(!p5&p6)+(!p6&p7)+(!p7&p8)+(!p8&p9)+(!p9&p2)
    c1 <- m & (N>=2)&(N<=6)&(trans==1)&!(p2&p4&p6)&!(p4&p6&p8)
    if (any(c1)) { m[c1] <- FALSE; changed <- TRUE }
    p2 <- sh(m,-1,0); p3 <- sh(m,-1,1); p4 <- sh(m,0,1); p5 <- sh(m,1,1)
    p6 <- sh(m,1,0);  p7 <- sh(m,1,-1); p8 <- sh(m,0,-1); p9 <- sh(m,-1,-1)
    N  <- p2+p3+p4+p5+p6+p7+p8+p9
    trans <- (!p2&p3)+(!p3&p4)+(!p4&p5)+(!p5&p6)+(!p6&p7)+(!p7&p8)+(!p8&p9)+(!p9&p2)
    c2 <- m & (N>=2)&(N<=6)&(trans==1)&!(p2&p4&p8)&!(p2&p6&p8)
    if (any(c2)) { m[c2] <- FALSE; changed <- changed|any(c2) }
  }
  m
}
