#
#     Two-fold recalibration of prediction models under selective predictor measurement
#     Delord M. et al. -- King's College London -- marc.delord@kcl.ac.uk
#
#     Generic R code for the test-based recalibration procedure, illustrated on a
#     simulated dataset. Applies to any published model of the form
#         risk = 1 - S0(t)^exp(LP)
#     in which one predictor Z is measured in some patients only (S = 1) and imputed
#     from the predictors available in all patients otherwise (S = 0).
#
#     Input data (one row per patient):
#       time    follow-up time
#       status  0 = censored, 1 = event of interest, 2 = competing event
#       S       1 if Z is measured, 0 otherwise
#       lp      linear predictor of the published model, computed with the observed
#               Z (S = 1) or with the imputed Z (S = 0)
#
#     Steps
#       1. pseudo-observations of the cumulative incidence at horizon H
#       2. calibration model (van Geloven et al. 2022): cloglog(F1) = a + b cloglog(p)
#       3. Test 1: do the two patterns need the same correction? (H0: a_S = b_S = 0)
#       4. Test 2: closed testing (Vergouwe et al. 2017) within the patients concerned
#       5. recalibrated equations: S0_new and coefficients multiplied by b
#
#     Requires: survival, prodlim, geepack
#

library(survival)
library(prodlim)
library(geepack)

clamp   <- function(p, eps = 1e-6) pmin(pmax(p, eps), 1 - eps)
cloglog <- function(p) log(-log(1 - p))
icll    <- function(eta) 1 - exp(-exp(eta))

# ==============================================================================
# 1. Imputation of the selectively measured predictor (fitted where S = 1)
# ==============================================================================
# formula: e.g. z ~ x1 + x2 ; returns the fitted model; predict() gives the imputed z
impute_model <- function(formula, d) lm(formula, data = d, subset = S == 1)

# ==============================================================================
# 2. Observed outcome: pseudo-observations of the cumulative incidence at H
# ==============================================================================
pseudo_obs <- function(d, H)
  as.vector(jackknife(prodlim(Hist(time, status) ~ 1, data = d), times = H, cause = 1))

# Aalen-Johansen cumulative incidence at H (observed risk)
aj <- function(d, H) {
  st  <- factor(d$status, levels = 0:2, labels = c("censor", "event", "competing"))
  fit <- survfit(Surv(d$time, st) ~ 1)
  s   <- summary(fit, times = H, extend = TRUE)
  unname(s$pstate[1, which(fit$states == "event")])
}

# ==============================================================================
# 3. Calibration model: GEE on pseudo-observations, cloglog link, binomial variance
# ==============================================================================
gee_fit <- function(formula, d)
  geese(formula, id = id, data = d, family = gaussian, mean.link = "cloglog",
        variance = "binomial", scale.fix = TRUE, corstr = "independence", jack = TRUE)

# calibration summary: O/E, intercept (slope fixed to 1), slope
calib <- function(d, p, H) {
  d$x <- cloglog(clamp(p))
  m1  <- gee_fit(pv ~ 1 + offset(x), d)
  m2  <- gee_fit(pv ~ x, d)
  a   <- unname(m1$beta[1]); se_a <- sqrt(m1$vbeta.ajs[1, 1])
  b   <- unname(m2$beta[2]); se_b <- sqrt(m2$vbeta.ajs[2, 2])
  o   <- aj(d, H)
  c(n = nrow(d), observed = o, predicted = mean(p), OE = o / mean(p),
    intercept = a, int_lo = a - 1.96 * se_a, int_hi = a + 1.96 * se_a,
    slope = b, slope_lo = b - 1.96 * se_b, slope_hi = b + 1.96 * se_b)
}

# ==============================================================================
# 4. Test 2: closed testing of the required correction
#    level 0 = none, 1 = intercept, 2 = intercept + slope
# ==============================================================================
closed_test <- function(d, alpha = 0.05) {
  m1 <- gee_fit(pv ~ 1 + offset(cll), d)
  m2 <- gee_fit(pv ~ cll, d)
  ab <- unname(m2$beta); V <- m2$vbeta.ajs
  dl <- c(ab[1], ab[2] - 1)
  p_any   <- pchisq(drop(t(dl) %*% solve(V) %*% dl), 2, lower.tail = FALSE)
  p_slope <- pchisq((ab[2] - 1)^2 / V[2, 2], 1, lower.tail = FALSE)
  level   <- if (p_any >= alpha) 0L else if (p_slope >= alpha) 1L else 2L
  list(par   = rbind(`0` = c(a = 0, b = 1),
                     `1` = c(a = unname(m1$beta[1]), b = 1),
                     `2` = c(a = ab[1], b = ab[2])),
       p_any = p_any, p_slope = p_slope, level = level, theta = ab, V = V)
}

# ==============================================================================
# 5. Full procedure
# ==============================================================================
twofold_recalibration <- function(d, S0, H, alpha = 0.05) {
  d$p   <- clamp(1 - S0^exp(d$lp))
  d$cll <- cloglog(d$p)
  d$id  <- seq_len(nrow(d))
  d$pv  <- pseudo_obs(d, H)

  before <- rbind(`S = 0` = calib(d[d$S == 0, ], d$p[d$S == 0], H),
                  `S = 1` = calib(d[d$S == 1, ], d$p[d$S == 1], H))

  ct <- list(S0 = closed_test(d[d$S == 0, ], alpha), S1 = closed_test(d[d$S == 1, ], alpha),
             pooled = closed_test(d, alpha))

  # Test 1: H0 a_S = b_S = 0 (separate fits = fully interacted model under independence GEE)
  dth   <- ct$S1$theta - ct$S0$theta
  W     <- drop(t(dth) %*% solve(ct$S0$V + ct$S1$V) %*% dth)
  p_het <- pchisq(W, 2, lower.tail = FALSE)
  twofold <- p_het < alpha

  # corrections selected by Test 2 within the patients concerned
  sel <- function(x) x$par[as.character(x$level), ]
  par <- if (twofold) list(`0` = sel(ct$S0), `1` = sel(ct$S1))
         else         list(`0` = sel(ct$pooled), `1` = sel(ct$pooled))
  d$p_new <- ifelse(d$S == 0, icll(par$`0`[["a"]] + par$`0`[["b"]] * d$cll),
                              icll(par$`1`[["a"]] + par$`1`[["b"]] * d$cll))

  after <- rbind(`S = 0` = calib(d[d$S == 0, ], d$p_new[d$S == 0], H),
                 `S = 1` = calib(d[d$S == 1, ], d$p_new[d$S == 1], H))

  # recalibrated equations: risk = 1 - S0_new^exp(b * LP), log(-log S0_new) = a + b log(-log S0)
  eq <- t(sapply(par, function(x) c(a = x[["a"]], b = x[["b"]],
                                     S0_new = exp(-exp(x[["a"]] + x[["b"]] * log(-log(S0)))))))
  rownames(eq) <- c("S = 0", "S = 1")

  list(test1 = data.frame(a_S = dth[1], b_S = dth[2], W = W, df = 2, p = p_het,
                          decision = if (twofold) "two-fold" else "pooled"),
       test2 = data.frame(group = c("S = 0", "S = 1", "pooled"),
                          p_any_correction = sapply(ct, `[[`, "p_any"),
                          p_slope = sapply(ct, `[[`, "p_slope"),
                          level = c("none", "intercept", "intercept + slope")[sapply(ct, `[[`, "level") + 1],
                          row.names = NULL),
       before = before, after = after, equations = eq, data = d)
}

# ==============================================================================
# Example on simulated data (no real patient data)
# ==============================================================================
set.seed(2026)
n  <- 10000
x1 <- rnorm(n)                                   # predictors available in all patients
x2 <- rbinom(n, 1, 0.5)
w  <- rnorm(n)                                   # information driving measurement, not in the model
z  <- 0.5 * x1 + 0.5 * w + rnorm(n, sd = 0.8)    # selectively measured predictor
S  <- rbinom(n, 1, plogis(-1 + 0.5 * x1 + 1.2 * w))

# event of interest and competing event (exponential), uniform censoring
eta    <- -4 + 0.8 * x1 + 0.3 * x2 + 0.6 * z
t_ev   <- rexp(n, exp(eta))
t_comp <- rexp(n, 0.04)
t_cens <- runif(n, 0, 10)
time   <- pmin(t_ev, t_comp, t_cens)
status <- ifelse(time == t_ev, 1L, ifelse(time == t_comp, 2L, 0L))

d <- data.frame(time, status, S, x1, x2, z = ifelse(S == 1, z, NA))

# published model (assumed): LP = 0.8 x1 + 0.3 x2 + 0.6 z, S0(5) = 0.97
imp   <- impute_model(z ~ x1 + x2, d)
z_use <- ifelse(d$S == 1, d$z, predict(imp, newdata = d))
d$lp  <- 0.8 * d$x1 + 0.3 * d$x2 + 0.6 * z_use

res <- twofold_recalibration(d, S0 = 0.97, H = 5)
res$test1
res$test2
round(res$before, 3)
round(res$after, 3)
round(res$equations, 4)
