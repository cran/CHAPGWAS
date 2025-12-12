
#' RANDOM: Mixed model and haplotype random effect test
#'
#' Internal core function used in SEL.HAP() to fit mixed models and
#' test the random effect of a haplotype design matrix \code{z}.
#'
#' @param z Haplotype design matrix (n x r). When \code{PAR} is NULL,
#'   \code{z} can be NULL and the function only estimates variance
#'   components in the mixed model.
#' @param YFIX A data.frame or matrix with phenotype in the first column
#'   and fixed-effect covariates in the remaining columns.
#' @param KIN A list of kinship matrices.
#' @param PAR Optional vector of (log) variance components. If NULL, they
#'   are estimated by the internal mixed model.
#'
#' @return If \code{PAR} is NULL, returns the estimated parameters (vector).
#'   Otherwise returns a list with likelihood ratio test statistics and
#'   parameter estimates.
#' @importFrom MASS ginv
#' @importFrom stats optim pchisq
#' @export
RANDOM <- function(z, YFIX, KIN, PAR) {

  MIXED <- function(YFIX, KIN) {

    loglike <- function(theta) {
      if (n.var > 1) {
        lambda <- as.list(exp(theta))
        Kin <- Map("*", KIN, lambda)
        KK <- Reduce("+", Kin)
        qq <- eigen(KK)
        yu <- t(qq[[2]]) %*% y
        xu <- t(qq[[2]]) %*% x
      } else {
        qq <- QQ
        qq[[1]] <- exp(theta) * qq[[1]]
      }
      delta <- qq[[1]]
      logdt <- sum(log(delta + 1))
      h <- 1 / (delta + 1)
      yy <- sum(yu * h * yu)
      lx <- t(xu) %*% diag(h)
      yx <- lx %*% yu
      xx <- lx %*% xu

      CLxx <- tryCatch(chol(xx), error = function(e) NULL)
      if (is.null(CLxx)) {
        ev <- eigen(xx)[[1]]
        INVxx <- MASS::ginv(xx)
      } else {
        ev <- diag(CLxx)^2
        INVxx <- chol2inv(CLxx)
      }
      DETxx <- prod(ev[ev > 0])

      loglike <- -0.5 * logdt -
        0.5 * (n - q) * log(yy - t(yx) %*% INVxx %*% yx) -
        0.5 * log(DETxx)
      if (!is.finite(loglike)) loglike <- -1e+10
      return(-loglike)
    }

    x <- as.matrix(YFIX[, -1])
    y <- as.matrix(YFIX[, 1])
    q <- ncol(x)
    n <- nrow(x)
    n.var <- length(KIN)
    if (n.var == 1) {
      QQ <- eigen(KIN[[1]])
      uu <- QQ[[2]]
      yu <- t(uu) %*% y
      xu <- t(uu) %*% x
    }
    theta <- rep(0, n.var)
    Parm <- optim(
      par = theta, fn = loglike, hessian = TRUE,
      method = "L-BFGS-B",
      lower = rep(-50, length(KIN)), upper = rep(10, length(KIN))
    )
    return(Parm)
  }

  Loglike <- function(theta) {
    xi <- exp(theta)
    tt <- zz * xi + diag(r)

    CLtt <- tryCatch(chol(tt), error = function(e) NULL)
    if (is.null(CLtt)) {
      ev <- eigen(tt)[[1]]
      INVtt <- MASS::ginv(tt)
    } else {
      ev <- diag(CLtt)^2
      INVtt <- chol2inv(CLtt)
    }
    DETtt <- prod(ev[ev > 0])

    logdt2 <- log(DETtt)
    zm <- xi * zx %*% INVtt
    my <- xi * INVtt %*% zy
    yHy <- yy - t(zy) %*% my
    yHx <- yx - zm %*% zy
    xHx <- xx - zm %*% t(zx)
    CLxHx <- tryCatch(chol(xHx), error = function(e) NULL)
    if (is.null(CLxHx)) {
      ev <- eigen(xHx)[[1]]
      INVxHx <- MASS::ginv(xHx)
    } else {
      ev <- diag(CLxHx)^2
      INVxHx <- chol2inv(CLxHx)
    }
    DETxHx <- prod(ev[ev > 0])
    loglike <- -0.5 * logdt2 -
      0.5 * (n - s) * log(yHy - t(yHx) %*% INVxHx %*% yHx) -
      0.5 * log(DETxHx)
    if (!is.finite(loglike)) loglike <- -1e+10
    return(-loglike)
  }

  fixed <- function(xi) {
    tmp0 <- zz * xi + diag(r)
    CLtmp0 <- tryCatch(chol(tmp0), error = function(e) NULL)
    if (is.null(CLtmp0)) {
      INVtmp0 <- MASS::ginv(tmp0)
    } else {
      INVtmp0 <- chol2inv(CLtmp0)
    }
    tmp <- xi * INVtmp0
    zxm <- zx %*% tmp
    zzm <- zz %*% tmp
    mzy <- tmp %*% zy

    xyzmzyx <- Map("%*%", list(zx, tmp, zz), list(tmp, zy, tmp))
    BAL <- list(
      xyzmzyx[[1]], zx, t(zy),
      xyzmzyx[[1]], xyzmzyx[[3]], xyzmzyx[[3]]
    )
    BAR <- list(
      t(zx), xyzmzyx[[2]], xyzmzyx[[2]],
      zz, zy, zz
    )
    XYZHZYX <- Map("-", list(xx, yx, yy, zx, zy, zz),
                   Map("%*%", BAL, BAR))

    yHy <- yy - t(zy) %*% mzy
    yHx <- yx - zx %*% mzy
    xHx <- xx - zxm %*% t(zx)
    zHx <- zx - zxm %*% zz
    zHy <- zy - zzm %*% zy
    zHz <- zz - zzm %*% zz
    CLxHx <- tryCatch(chol(xHx), error = function(e) NULL)
    if (is.null(CLxHx)) {
      INVxHx <- MASS::ginv(xHx)
    } else {
      INVxHx <- chol2inv(CLxHx)
    }
    beta <- INVxHx %*% yHx
    sigma2 <- (yHy - t(yHx) %*% INVxHx %*% yHx) / (n - s)
    gamma <- xi * zHy - xi * t(zHx) %*% INVxHx %*% yHx
    var <- abs((xi * diag(r) - xi * zHz * xi) * as.numeric(sigma2))
    stderr <- var
    result <- list(gamma, stderr, beta, sigma2)
    return(result)
  }

  if (is.null(PAR)) {
    PAR <- MIXED(YFIX, KIN)$par
    return(PAR)
  } else {
    x <- as.matrix(YFIX[, -1])
    y <- as.matrix(YFIX[, 1])
    n <- nrow(y)
    s <- ncol(x)
    n.var <- length(KIN)

    lambda <- exp(PAR)
    if (n.var > 1) {
      lambda <- as.list(lambda)
      Kin <- Map("*", KIN, lambda)
      KK <- Reduce("+", Kin)
      qq <- eigen(KK)
    } else {
      qq <- eigen(KIN[[1]])
      qq[[1]] <- lambda * qq[[1]]
    }

    yu <- t(qq[[2]]) %*% y
    xu <- t(qq[[2]]) %*% x
    h <- 1 / (qq[[1]] + 1)
    yy <- sum(yu * h * yu)
    lx <- t(xu) %*% diag(h)
    yx <- lx %*% yu
    xx <- lx %*% xu

    r <- ncol(z)
    zu <- t(qq[[2]]) %*% z
    lz <- t(zu) %*% diag(h)
    zy <- lz %*% yu
    zz <- lz %*% zu
    zx <- t(lz %*% xu)

    Low <- rep(-10, 1)
    Up <- -Low

    theta <- rep(0, 1)
    par <- optim(
      par = theta, fn = Loglike, hessian = TRUE,
      method = "L-BFGS-B", lower = Low, upper = Up
    )
    xi <- exp(par$par)
    parmfix <- fixed(xi)

    gamma <- parmfix[[1]]
    stderr <- parmfix[[2]]
    beta <- parmfix[[3]]
    sigma2 <- parmfix[[4]]
    lambda <- xi
    sigma2g <- lambda * sigma2
    fn0 <- Loglike(-Inf)
    fn1 <- par$value
    lrt <- 2 * (fn0 - fn1)
    if (lrt <= 1e-03) {
      p_lrt <- 0.5 + pchisq(lrt, df = 1, lower.tail = FALSE) * 0.5
    } else {
      p_lrt <- 1 - (0.5 + pchisq(lrt, df = 1) * 0.5)
    }
    p.l <- p_lrt

    parm <- c(beta, sigma2, sigma2g, gamma, stderr, lrt, p_lrt)
    RESULT <- list(c(lrt, p.l, p_lrt), parm)
    return(RESULT)
  }
}

#' test.HAP: Test haplotype effects for a local region
#'
#' @param HAP.X Haplotype genotype matrix for a candidate genomic region.
#'   Each column corresponds to a SNP within the region. Every two
#'   consecutive rows correspond to the two homologous haplotypes of
#'   the same individual (rows 1–2 = individual 1, rows 3–4 = individual 2, etc.).
#'
#' @param YFIX Phenotype and fixed-effect covariates.
#' @param KIN List of kinship matrices.
#' @param PAR Optional vector of variance components; if NULL, estimated internally.
#'
#' @return A list containing:
#'   * my.scan – scan statistics
#'   * z – haplotype design matrix (n x number of haplotypes)
#'   * ELEMENT – labels of distinct haplotype categories
#'
#' @importFrom plyr mapvalues
#' @export
test.HAP <- function(HAP.X, YFIX, KIN, PAR) {
  n <- nrow(YFIX)

  NAME <- apply(HAP.X, 1, paste, collapse = "")
  HAP.CAT <- unique(NAME)
  NAME.hap <- plyr::mapvalues(NAME, from = HAP.CAT, to = seq_along(HAP.CAT))

  code <- 1:n
  NAME.HAP <- t(matrix(NAME.hap[c(2 * code - 1, 2 * code)], n, 2))
  NAME.HAPLOTYPE <- paste(apply(NAME.HAP, 2, min),
                          apply(NAME.HAP, 2, max),
                          sep = "")
  ELEMENT <- unique(NAME.HAPLOTYPE)
  NUM.ELEMENT <- length(ELEMENT)
  my.scan <- c(0, 1, 1, NUM.ELEMENT)
  z <- NULL

  if (NUM.ELEMENT > 1) {
    Element <- diag(NUM.ELEMENT)
    z <- matrix(-1, n, NUM.ELEMENT)
    for (ii in 1:NUM.ELEMENT) {
      pos.ii <- which(NAME.HAPLOTYPE == ELEMENT[ii])
      z[pos.ii, ] <- rep(Element[ii, 1:NUM.ELEMENT],
                         each = length(pos.ii))
    }
    my.scan[1:3] <- RANDOM(z, YFIX, KIN, PAR)[[1]]
  }

  return(list(my.scan, z, ELEMENT))
}

#' SEL.HAP: Haplotype selection and extension along the genome
#'
#' Perform genome-wide haplotype selection and extension using the
#' CHAP-GWAS framework. The function scans along each chromosome,
#' builds local haplotype segments, and adaptively extends them based
#' on association evidence with the phenotype.
#'
#' @param GEN Genotype matrix with rows corresponding to markers and
#'   columns corresponding to individuals. The first two columns give
#'   chromosome (\code{chr}) and physical position (\code{pos}); the
#'   remaining columns contain alleles for each individual (e.g.
#'   \code{"A"}, \code{"C"}, \code{"G"}, \code{"T"}), one allele per
#'   haplotype copy.
#' @param YFIX A matrix or data.frame with phenotype in the first column
#'   and fixed-effect covariates (e.g. intercept, PCs) in the remaining
#'   columns, one row per individual.
#' @param KIN A list of kinship matrices, each of dimension
#'   \eqn{n \times n}, where \eqn{n} is the number of individuals.
#' @param nHap Initial haplotype window size (number of consecutive
#'   markers).
#' @param p.threshold P-value threshold for haplotype extension.
#' @param PAR Optional variance component parameters passed to
#'   \code{RANDOM()}. If \code{NULL}, they are estimated internally.
#'
#' @return A list of three matrices summarizing:
#'   \itemize{
#'     \item \strong{FINAL[[1]]}: initial haplotype segments
#'     \item \strong{FINAL[[2]]}: extended haplotype segments
#'     \item \strong{FINAL[[3]]}: final selected segments after extension
#'   }
#'
#' @examples
#' ## Minimal example with small simulated data (alleles encoded as A/C/G/T)
#' set.seed(1)
#'
#' ## Number of individuals and markers
#' n_ind  <- 200
#' n_mark <- 50
#'
#' ## Construct a simple GEN matrix:
#' ## first two columns: chromosome and position
#' ## remaining columns: alleles "A","C","G","T"
#' chr <- rep(1, n_mark)
#' pos <- seq_len(n_mark) * 100
#' alleles <- c("A", "C", "G", "T")
#' geno <- matrix(sample(alleles, n_mark * n_ind, replace = TRUE),
#'                nrow = n_mark, ncol = n_ind)
#' GEN <- cbind(chr, pos, geno)
#'
#' ## Phenotype + intercept as fixed effect
#' y <- rnorm(n_ind)
#' X <- cbind(1, rnorm(n_ind))  # intercept + one covariate
#' YFIX <- cbind(y, X)
#'
#' ## Simple kinship: identity matrix
#' KIN <- list(diag(n_ind))
#'
#' ## Run SEL.HAP with a small initial window and mild threshold
#' res <- SEL.HAP(GEN, YFIX, KIN,
#'                nHap = 2,
#'                p.threshold = 0.05,
#'                PAR = NULL)
#'
#' ## Inspect the structure of the result (three matrices)
#' str(res)
#'
#' @importFrom utils tail
#' @export
SEL.HAP <- function(GEN, YFIX, KIN, nHap, p.threshold, PAR) {
  MAP <- matrix(as.numeric(GEN[, 1:2]), ncol = 2)
  chr <- unique(MAP[, 1])

  PoS <- list()
  Pos.Chr <- list()
  for (i in chr) {
    Pos.Chr[[i]] <- which(MAP[, 1] == i)
    PoS[[i]] <- Pos.Chr[[i]][-length(Pos.Chr[[i]])]
  }

  Extension <- function(POS, POS.CHR, GEN, YFIX, KIN, p.threshold, PAR) {
    Former <- c(
      MAP[POS[1], 1], range(POS), MAP[range(POS), 2],
      test.HAP(
        matrix(t(GEN[POS, -(1:2)]), ncol = length(POS)),
        YFIX, KIN, PAR
      )[[1]]
    )

    WL <- NULL
    FORMER <- Former
    LATER <- FORMER
    if (FORMER[8] < p.threshold) {
      iter <- 0
      RE.WL <- list()
      while (((LATER[6] > FORMER[6] &&
               LATER[7] <= FORMER[7] &&
               LATER[9] > FORMER[9]) || iter == 0) &&
             length(POS) > 1) {
        if (iter > 0) {
          FORMER <- LATER
          POS <- Pos
        }
        POS.ALL <- c(
          max(POS[1] - 1, POS.CHR[1]),
          min(tail(POS, 1) + 1, tail(POS.CHR, 1))
        )
        POS.T <- list(
          unique(c(POS.ALL[1], POS)),
          unique(c(POS, POS.ALL[2]))
        )
        POS.T <- POS.T[order(lengths(POS.T))]

        Later <- NULL
        for (i in 1:2) {
          T.LATER <- c(
            MAP[POS.T[[i]][1], 1],
            range(POS.T[[i]]),
            MAP[range(POS.T[[i]]), 2],
            test.HAP(
              matrix(t(GEN[POS.T[[i]], -(1:2)]),
                     ncol = length(POS.T[[i]])),
              YFIX, KIN, PAR
            )[[1]]
          )
          RE.WL[[2 * iter + i]] <- T.LATER
          Later <- rbind(Later, T.LATER)
        }

        POS.PT <- c(which.min(Later[, 7]), which.max(Later[, 6]))
        Pos.LR <- max(
          POS.PT[1],
          POS.PT[2 - abs(diff(lengths(POS.T)))]
        )
        Pos <- POS.T[[Pos.LR]]
        LATER <- Later[Pos.LR, ]

        iter <- iter + 1
      }
      WL <- do.call(rbind, RE.WL)
    }
    return(list(Former, FORMER, WL))
  }

  if (is.null(PAR)) {
    PAR <- RANDOM(NULL, YFIX, KIN, NULL)
  }

  Final <- list(NULL, NULL, NULL)
  for (i in chr) {
    for (j in PoS[[i]]) {
      Final[[1]][[j]] <- Extension(
        j, Pos.Chr[[i]], GEN, YFIX, KIN,
        p.threshold, PAR
      )[[1]]
      Re <- Extension(
        j:(j + nHap - 1), Pos.Chr[[i]], GEN, YFIX, KIN,
        p.threshold, PAR
      )
      Final[[3]][[j]] <- Re[[2]]
      Final[[2]][[j]] <- Re[[1]]
    }
    Final[[1]][[max(PoS[[i]]) + 1]] <-
      Extension(
        max(PoS[[i]]) + 1, Pos.Chr[[i]], GEN, YFIX, KIN,
        p.threshold, PAR
      )[[1]]
  }

  FINAL <- lapply(Final, do.call, what = rbind)
  return(FINAL)
}

