#------------------------------------------------------
# Momentary time sampling 
#------------------------------------------------------

MTS_info <- function(phi, zeta, c_M) {
  rho <- zeta / (phi * (1 - phi))
  erc <- exp(-rho * c_M)
  p0 <- p_0(c_M, phi, zeta)
  p1 <- p_1(c_M, phi, zeta)
  eq_pi <- rep(c(1 - phi, phi), each = 2) / c(1 - p0, p0, 1 - p1, p1)
  dpi_dphi <- c(erc - 1, 1 - erc, erc - 1, 1 - erc)
  dpi_drho <- c(-phi * c_M * erc, phi * c_M * erc, (1 - phi) * c_M * erc, -(1 - phi) * c_M * erc)
  info <- rbind(dpi_dphi, dpi_drho) %*% (eq_pi * cbind(dpi_dphi, dpi_drho))
  rownames(info) <- colnames(info) <- c("phi","rho")
  info
}



AIR_info <- function(phi, zeta, c_A, d_A) {
  rho <- zeta / (phi * (1 - phi))
  p0d <- p_0(c_A, phi, zeta)
  p1d <- p_0(c_A, phi, zeta)
  p0cd <- p_0(c_A + d_A, phi, zeta)
  p1cd <- p_1(c_A + d_A, phi, zeta)
  erpc <- exp(-rho * phi * c_A)
  erqc <- exp(-rho * (1 - phi) * c_A)
  pi_trans <- c((1 - p0d) * erpc, p0d * erpc, 1 - p0cd - (1 - p0d) * erpc, p0cd - p0d * erpc,
                1 - p1cd - (1 - p1d) * erqc, p1cd - p1d * erqc, (1 - p1d) * erqc, p1d * erqc)
  eq_pi <- rep(c(1 - phi, phi), each = 4) / pi_trans
  
  dpi_dphi <- c()
  dpi_drho <- c()
  info <- rbind(dpi_dphi, dpi_drho) %*% (eq_pi * cbind(dpi_dphi, dpi_drho))
  rownames(info) <- colnames(info) <- c("phi","rho")
  info
}


