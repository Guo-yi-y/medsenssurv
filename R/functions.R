#' generate conditional probability of U

#' @param data   data frame including time, outcome, m, t, and covariates
#' @param inter  	a logical value indicating the existence of exposure-mediator interaction
#' @param p.U    marginal probability of unmeasured confounder
#' @param outcome survival outcome
#' @param time    survival time
#' @param m       mediator
#' @param t       exposure
#' @param covariates the vector of covariates
#' @param b.y,b.m,b.t the coefficients between U and exposure, mediator, and survival outcome
#' @param scale.m,scale.t,scale.u scale of mediator, exposure and U, binary or continuous
#' @param sigma.u The standard deviation of the prior distribution of U if it is continuous. The default is 1.
#' @param iter_num number of iteration
#' @param lower_u,upper_u  lower and upper bound of the integration interval for U if it is continuous.
#' @return U and its conditional probability (if U is binary)
#' @importFrom stats rbinom dnorm gaussian binomial glm as.formula model.matrix
#' @importFrom stats pnorm sd var sigma setNames
#' @importFrom survival coxph basehaz
#' @importFrom dplyr bind_cols bind_rows mutate
#' @importFrom magrittr %>%
#' @importFrom glue glue
#' @importFrom stats coef integrate rnorm
#' @importFrom distr AbscontDistribution r
#' @importFrom CMAverse cmest
#' @importFrom parallel makeCluster detectCores clusterExport clusterEvalQ parLapply
#' @importFrom plotly plot_ly add_trace layout
#' @export

gen_p <- function(data, inter = FALSE, p.U = 0.5, outcome, time, m, t, covariates, b.y, b.m, b.t, scale.m, scale.t, scale.u, sigma.u = 1, iter_num = 20, lower_u = -10, upper_u = 10) {

  covariates_fm <- if (length(covariates) == 0) "" else paste(covariates, collapse = " + ")

  add_rhs <- function(...) {
    terms <- c(...)
    terms <- terms[terms != ""]
    paste(terms, collapse = " + ")
  }

  make_init_formulas <- function(inter, outcome, time, m, t, covariates_fm, u_name) {
    out_rhs <- if (!inter) {
      add_rhs(t, m, covariates_fm, u_name)
    } else {
      add_rhs(t, m, paste0(t, ":", m), covariates_fm, u_name)
    }

    med_rhs <- add_rhs(t, covariates_fm, u_name)
    exp_rhs <- add_rhs(covariates_fm, u_name)

    list(
      out = as.formula(glue("Surv({time}, {outcome}) ~ {out_rhs}")),
      med = as.formula(glue("{m} ~ {med_rhs}")),
      exp = as.formula(glue("{t} ~ {exp_rhs}"))
    )
  }

  make_refit_formulas <- function(inter, outcome, time, m, t, covariates_fm, scale.u) {
    if (scale.u == "binary") {
      out_rhs <- if (!inter) {
        add_rhs(t, m, covariates_fm, "offset(log(exp(b.y) * p + (1 - p)))")
      } else {
        add_rhs(t, m, paste0(t, ":", m), covariates_fm, "offset(log(exp(b.y) * p + (1 - p)))")
      }
    } else {
      out_rhs <- if (!inter) {
        add_rhs(t, m, covariates_fm, "offset(b.y * U_simu)")
      } else {
        add_rhs(t, m, paste0(t, ":", m), covariates_fm, "offset(b.y * U_simu)")
      }
    }

    med_rhs <- add_rhs(t, covariates_fm, "offset(b.m * U_simu)")
    exp_rhs <- add_rhs(covariates_fm, "offset(b.t * U_simu)")

    list(
      out = as.formula(glue("Surv({time}, {outcome}) ~ {out_rhs}")),
      med = as.formula(glue("{m} ~ {med_rhs}")),
      exp = as.formula(glue("{t} ~ {exp_rhs}"))
    )
  }

  fit_glm_scale <- function(formula, data, scale_type) {
    if (scale_type == "continuous") {
      glm(formula, data = data, family = gaussian())
    } else {
      glm(formula, data = data, family = binomial(link = "logit"))
    }
  }

  add_u_coef <- function(fit, b_u, u_name) {
    cc <- coef(fit)
    cc[is.na(cc)] <- 0
    cc[u_name] <- b_u
    cc
  }

  fit_init_models <- function(data, inter, outcome, time, m, t, covariates_fm,
                              scale.m, scale.t, u_name, b.y, b.m, b.t) {
    fms <- make_init_formulas(inter, outcome, time, m, t, covariates_fm, u_name)

    o_fit <- coxph(fms$out, data = data)
    m_fit <- fit_glm_scale(fms$med, data, scale.m)
    x_fit <- fit_glm_scale(fms$exp, data, scale.t)

    list(
      o_fit = o_fit,
      m_fit = m_fit,
      x_fit = x_fit,
      o_coef = add_u_coef(o_fit, b.y, u_name),
      m_coef = add_u_coef(m_fit, b.m, u_name),
      x_coef = add_u_coef(x_fit, b.t, u_name)
    )
  }

  fit_refit_models <- function(data, inter, outcome, time, m, t, covariates_fm,
                               scale.m, scale.t, scale.u, b.y, b.m, b.t) {
    fms <- make_refit_formulas(inter, outcome, time, m, t, covariates_fm, scale.u)

    o_fit <- coxph(fms$out, data = data)
    m_fit <- fit_glm_scale(fms$med, data, scale.m)
    x_fit <- fit_glm_scale(fms$exp, data, scale.t)

    list(
      o_fit = o_fit,
      m_fit = m_fit,
      x_fit = x_fit,
      o_coef = add_u_coef(o_fit, b.y, "U_simu"),
      m_coef = add_u_coef(m_fit, b.m, "U_simu"),
      x_coef = add_u_coef(x_fit, b.t, "U_simu")
    )
  }

  get_basehaz_info <- function(o_fit, data, time) {
    bh <- basehaz(o_fit, centered = FALSE)
    bh$time <- round(bh$time, 3)
    time_obs <- round(data[, time], 3)
    idx <- match(time_obs, bh$time)
    H0 <- bh$hazard[idx]
    list(bh = bh, H0 = H0, idx = idx)
  }

  ## 模型矩阵去掉 u_name 列，并按名字对齐系数
  get_mm_coef_without_u <- function(fit_obj, coef_vec, u_name) {
    mm <- model.matrix(fit_obj)

    if (u_name %in% colnames(mm)) {
      mm0 <- mm[, colnames(mm) != u_name, drop = FALSE]
    } else {
      mm0 <- mm
    }

    coef0 <- coef_vec[names(coef_vec) != u_name]
    coef0 <- coef0[colnames(mm0)]

    list(mm0 = mm0, coef0 = coef0)
  }

  ## 二元 U：所有样本在固定 u = 0/1 下的 eta
  eta_binary <- function(fit_obj, coef_vec, u_val, u_name) {
    obj <- get_mm_coef_without_u(fit_obj, coef_vec, u_name)
    mm0 <- obj$mm0
    coef0 <- obj$coef0
    eta0 <- as.numeric(mm0 %*% coef0)
    bu <- unname(coef_vec[u_name])
    eta0 + bu * u_val
  }

  ## 连续 U：单个样本拆成 eta0 + b*u
  split_eta_single <- function(fit_obj, coef_vec, i, u_name) {
    obj <- get_mm_coef_without_u(fit_obj, coef_vec, u_name)
    Xi <- obj$mm0[i, , drop = FALSE]
    coef0 <- obj$coef0
    eta0 <- as.numeric(Xi %*% coef0)
    bu <- unname(coef_vec[u_name])
    list(eta0 = eta0, bu = bu)
  }

  ## m/t 条件核
  obs_kernel <- function(y_i, eta, fit_obj, scale_type) {
    if (scale_type == "binary") {
      p1 <- 1 / (1 + exp(-eta))
      p1^y_i * (1 - p1)^(1 - y_i)
    } else {
      dnorm(y_i, mean = eta, sd = sigma(fit_obj))
    }
  }

  ## binary U 的 Cox 核：保留原修正项
  surv_kernel_binary <- function(delta, eta, H0, p_vec, b.y) {
    corr <- exp(mean(log(exp(b.y) * p_vec + (1 - p_vec))))
    exp(eta)^delta * exp(-(H0 / corr) * exp(eta))
  }

  ## continuous U 的 Cox 核
  surv_kernel_cont <- function(delta, eta, H0) {
    exp(delta * eta - H0 * exp(eta))
  }

  ## -------------------------
  ## binary U update
  ## -------------------------
  update_binary_U <- function(data, fits, theta, p_current, outcome, time, m, t,
                              scale.m, scale.t, b.y, u_name) {
    o_fit  <- fits$o_fit
    m_fit  <- fits$m_fit
    x_fit  <- fits$x_fit
    o_coef <- fits$o_coef
    m_coef <- fits$m_coef
    x_coef <- fits$x_coef

    bh_info <- get_basehaz_info(o_fit, data, time)
    H0 <- bh_info$H0

    eta_o1 <- eta_binary(o_fit, o_coef, 1, u_name)
    eta_o0 <- eta_binary(o_fit, o_coef, 0, u_name)

    eta_m1 <- eta_binary(m_fit, m_coef, 1, u_name)
    eta_m0 <- eta_binary(m_fit, m_coef, 0, u_name)

    eta_x1 <- eta_binary(x_fit, x_coef, 1, u_name)
    eta_x0 <- eta_binary(x_fit, x_coef, 0, u_name)

    pyu1 <- surv_kernel_binary(data[, outcome], eta_o1, H0, p_current, b.y)
    pyu0 <- surv_kernel_binary(data[, outcome], eta_o0, H0, p_current, b.y)

    if (scale.m == "binary") {
      pm1u1 <- 1 / (1 + exp(-eta_m1))
      pm1u0 <- 1 / (1 + exp(-eta_m0))
      pmu1 <- pm1u1^data[, m] * (1 - pm1u1)^(1 - data[, m])
      pmu0 <- pm1u0^data[, m] * (1 - pm1u0)^(1 - data[, m])
    } else {
      pmu1 <- dnorm(data[, m], mean = eta_m1, sd = sigma(m_fit))
      pmu0 <- dnorm(data[, m], mean = eta_m0, sd = sigma(m_fit))
    }

    if (scale.t == "binary") {
      pt1u1 <- 1 / (1 + exp(-eta_x1))
      pt1u0 <- 1 / (1 + exp(-eta_x0))
      pxu1 <- pt1u1^data[, t] * (1 - pt1u1)^(1 - data[, t])
      pxu0 <- pt1u0^data[, t] * (1 - pt1u0)^(1 - data[, t])
    } else {
      pxu1 <- dnorm(data[, t], mean = eta_x1, sd = sigma(x_fit))
      pxu0 <- dnorm(data[, t], mean = eta_x0, sd = sigma(x_fit))
    }

    pu1 <- pmu1 * theta * pyu1 * pxu1
    pu0 <- pmu0 * (1 - theta) * pyu0 * pxu0

    p <- as.vector(pu1 / (pu1 + pu0))
    p[pu1 == 0 & pu0 == 0] <- 0

    U_simu <- rbinom(nrow(data), 1, p)

    list(p = p, U_simu = U_simu)
  }

  ## -------------------------
  ## continuous U update
  ## -------------------------
  update_continuous_U <- function(data, fits, outcome, time, m, t,
                                  scale.m, scale.t, sigma.u,
                                  lower_u, upper_u, u_name) {
    o_fit  <- fits$o_fit
    m_fit  <- fits$m_fit
    x_fit  <- fits$x_fit
    o_coef <- fits$o_coef
    m_coef <- fits$m_coef
    x_coef <- fits$x_coef

    bh_info <- get_basehaz_info(o_fit, data, time)
    H0_all <- bh_info$H0

    n <- nrow(data)

    res_inte <- sapply(seq_len(n), function(i) {
      eo <- split_eta_single(o_fit, o_coef, i, u_name)
      em <- split_eta_single(m_fit, m_coef, i, u_name)
      ex <- split_eta_single(x_fit, x_coef, i, u_name)

      H0_i <- H0_all[i]

      integrand <- function(u) {
        eta_o <- eo$eta0 + eo$bu * u
        eta_m <- em$eta0 + em$bu * u
        eta_x <- ex$eta0 + ex$bu * u

        pyu <- surv_kernel_cont(data[i, outcome], eta_o, H0_i)
        pmu <- obs_kernel(data[i, m], eta_m, m_fit, scale.m)
        pxu <- obs_kernel(data[i, t], eta_x, x_fit, scale.t)

        pyu * pmu * pxu * dnorm(u, mean = 0, sd = sigma.u)
      }

      integrate(integrand, lower = lower_u, upper = upper_u)$value
    })

    if (any(!is.finite(res_inte)) || any(res_inte <= 0)) {
      stop("Some normalization integrals are non-finite or non-positive in continuous U update.")
    }

    U_simu <- sapply(seq_len(n), function(i) {
      eo <- split_eta_single(o_fit, o_coef, i, u_name)
      em <- split_eta_single(m_fit, m_coef, i, u_name)
      ex <- split_eta_single(x_fit, x_coef, i, u_name)

      H0_i <- H0_all[i]

      conditional.u <- function(u) {
        eta_o <- eo$eta0 + eo$bu * u
        eta_m <- em$eta0 + em$bu * u
        eta_x <- ex$eta0 + ex$bu * u

        pyu <- surv_kernel_cont(data[i, outcome], eta_o, H0_i)
        pmu <- obs_kernel(data[i, m], eta_m, m_fit, scale.m)
        pxu <- obs_kernel(data[i, t], eta_x, x_fit, scale.t)

        pyu * pmu * pxu * dnorm(u, mean = 0, sd = sigma.u) / res_inte[i]
      }

      r(AbscontDistribution(d = conditional.u))(1)
    })

    list(U_simu = U_simu, res_inte = res_inte)
  }

  ## -------------------------
  ## initialization
  ## -------------------------
  n <- nrow(data)

  if (scale.u == "binary") {
    data$U_ini <- rbinom(n, 1, p.U)
    theta <- p.U
    data$p <- rep(p.U, n)
  } else {
    data$U_ini <- rnorm(n, mean = 0, sd = sigma.u)
    theta <- NA_real_
    data$p <- NA_real_
  }

  ## 初始拟合：把 U_ini 当普通协变量
  fits <- fit_init_models(
    data = data, inter = inter, outcome = outcome, time = time,
    m = m, t = t, covariates_fm = covariates_fm,
    scale.m = scale.m, scale.t = scale.t,
    u_name = "U_ini", b.y = b.y, b.m = b.m, b.t = b.t
  )

  ## 初始更新 U
  if (scale.u == "binary") {
    upd <- update_binary_U(
      data = data, fits = fits, theta = theta, p_current = data$p,
      outcome = outcome, time = time, m = m, t = t,
      scale.m = scale.m, scale.t = scale.t, b.y = b.y,
      u_name = "U_ini"
    )
    data$p <- upd$p
    data$U_simu <- upd$U_simu
  } else {
    upd <- update_continuous_U(
      data = data, fits = fits,
      outcome = outcome, time = time, m = m, t = t,
      scale.m = scale.m, scale.t = scale.t,
      sigma.u = sigma.u, lower_u = lower_u, upper_u = upper_u,
      u_name = "U_ini"
    )
    data$U_simu <- upd$U_simu
  }

  ## 用 offset 重新拟合
  fits <- fit_refit_models(
    data = data, inter = inter, outcome = outcome, time = time,
    m = m, t = t, covariates_fm = covariates_fm,
    scale.m = scale.m, scale.t = scale.t,
    scale.u = scale.u, b.y = b.y, b.m = b.m, b.t = b.t
  )

  ## -------------------------
  ## iterations
  ## -------------------------
  for (iter in seq_len(iter_num)) {
    if (scale.u == "binary") {
      upd <- update_binary_U(
        data = data, fits = fits, theta = theta, p_current = data$p,
        outcome = outcome, time = time, m = m, t = t,
        scale.m = scale.m, scale.t = scale.t, b.y = b.y,
        u_name = "U_simu"
      )
      data$p <- upd$p
      data$U_simu <- upd$U_simu
    } else {
      upd <- update_continuous_U(
        data = data, fits = fits,
        outcome = outcome, time = time, m = m, t = t,
        scale.m = scale.m, scale.t = scale.t,
        sigma.u = sigma.u, lower_u = lower_u, upper_u = upper_u,
        u_name = "U_simu"
      )
      data$U_simu <- upd$U_simu
    }

    fits <- fit_refit_models(
      data = data, inter = inter, outcome = outcome, time = time,
      m = m, t = t, covariates_fm = covariates_fm,
      scale.m = scale.m, scale.t = scale.t,
      scale.u = scale.u, b.y = b.y, b.m = b.m, b.t = b.t
    )
  }

  ## -------------------------
  ## return
  ## -------------------------
  if (scale.u == "binary") {
    return(list(
      p = data$p,
      U_simu = data$U_simu
    ))
  } else {
    return(list(
      p = NULL,
      U_simu = data$U_simu
    ))
  }
}

#' estimate the NIE and NDE under each sensitivity parameter
#' @param data   data frame including time, outcome, m, t, and covariates
#' @param inter  	a logical value indicating the existence of exposure-mediator interaction
#' @param p.U    marginal probability of unmeasured confounder
#' @param outcome survival outcome
#' @param time    survival time
#' @param m       mediator
#' @param t       exposure
#' @param covariates the vector of covariates
#' @param scale.m,scale.t,scale.u scale of mediator, exposure and U, binary or continuous
#' @param sigma.u The standard deviation of the prior distribution of U if it is continuous. The default is 1.
#' @param lower_u,upper_u  lower and upper bound of the integration interval for U if it is continuous.
#' @param range.b.y,range.b.m,range.b.t the range of coefficients between U and exposure, mediator, and survival outcome
#' @param ngrid defines how many intervals are used to divide range.b.y, range.b.m, and range.b.t into grid points.
#' @param iter_num number of iteration
#' @param multiple_draw whether to draw multiple samples of U
#' @param R the number of draws for multiple samples of U
#' @param B the number of bootstrap
#' @param parallel whether to use parallel computation
#' @return data frame containing sensitivity parameters, the corresponding mediation effects, and p-values.
#' @importFrom stats rbinom dnorm gaussian binomial glm as.formula model.matrix
#' @importFrom stats pnorm sd var sigma setNames
#' @importFrom survival coxph basehaz
#' @importFrom dplyr bind_cols bind_rows mutate
#' @importFrom stats coef integrate rnorm
#' @importFrom distr AbscontDistribution r
#' @importFrom magrittr %>%
#' @importFrom glue glue
#' @importFrom CMAverse cmest
#' @importFrom parallel makeCluster detectCores clusterExport clusterEvalQ parLapply
#' @importFrom plotly plot_ly add_trace layout
#' @export

medsurvsens = function(data, inter = FALSE, p.U = 0.5, outcome, time, m, t, covariates, scale.m, scale.t, scale.u, sigma.u = 1, lower_u = -10, upper_u = 10, range.b.y, range.b.m, range.b.t, ngrid = 10, iter_num = 20, multiple_draw = F, R = 1, B = 100, parallel = F){


  par_dat = expand.grid(b.y = seq(range.b.y[1], range.b.y[2], length.out = ngrid),
                        b.m = seq(range.b.m[1], range.b.m[2], length.out = ngrid),
                        b.t = seq(range.b.t[1], range.b.t[2], length.out = ngrid)


                        )



  if(parallel == F){
    all_res = lapply(1:nrow(par_dat), function(i){

      res = gen_p(data, inter, p.U, outcome, time, m, t, covariates, b.y = par_dat$b.y[i], b.m = par_dat$b.m[i], b.t = par_dat$b.t[i], scale.m, scale.t, scale.u, sigma.u,iter_num, lower_u,upper_u)

      if(scale.u == "binary" & multiple_draw == T){

        data <- bind_cols(data, as.data.frame(
          setNames(
            replicate(R, rbinom(nrow(data), 1, res$p), simplify = FALSE),
            paste0("U_simu", 1:R)
          )
        ))


        me_giu = summary(cmest(data = data, model = "gformula", outcome = time,
                               exposure = t, mediator = m,
                               basec = c(covariates),
                               EMint = inter, event = outcome,
                               mreg = list("logistic"),
                               yreg = "coxph",
                               astar = 0, a = 1, mval = list( 1), multimp = F, yval = 0,
                               estimation = "imputation", inference = "bootstrap", nboot = B,
                               boot.ci.type = "per"))$summarydf


        me_gsu_df = lapply(1:R, function(j){



          temp_df = lapply(1:B, function(B){

            B_res = summary(cmest(data = data, model = "gformula", outcome = time,
                                  exposure = t, mediator = m,
                                  basec = c(covariates, glue("U_simu{j}") ),
                                  EMint = T, event = outcome,
                                  mreg = list("logistic"),
                                  yreg = "coxph",
                                  astar = 0, a = 1, mval = list( 1), multimp = F, yval = 0,
                                  estimation = "imputation", inference = "bootstrap", nboot = B,
                                  boot.ci.type = "per"))$summarydf
            B_df = data.frame(log_nde = log(B_res[2,1]), log_nie = log(B_res[5,1]))
            return(B_df)


          }) %>% bind_rows()

          temp_res = data.frame(log_nie = mean(temp_df$log_nie), sig_log_nie = sd(temp_df$log_nie), log_nde = mean(temp_df$log_nde), sig_log_nde = sd(temp_df$log_nde))


          return(temp_res)

        }) %>% bind_rows()

        nie_sens = data.frame(medef = mean(me_gsu_df$log_nie), medef_std = (mean(me_gsu_df$sig_log_nie^2)+(1+1/R)*var(me_gsu_df$log_nie)   )^0.5 ) %>%
          mutate(p_ef = 2 * (1 - pnorm(abs(medef/medef_std))),
                 b_o = par_dat$b.y[i], b_m = par_dat$b.m[i], b_x = par_dat$b.t[i])


        nde_sens = data.frame(medef = mean(me_gsu_df$log_nde), medef_std = (mean(me_gsu_df$sig_log_nde^2)+(1+1/R)*var(me_gsu_df$log_nde)   )^0.5 ) %>%
          mutate(p_ef = 2 * (1 - pnorm(abs(medef/medef_std))),
                 b_o = par_dat$b.y[i], b_m = par_dat$b.m[i], b_x = par_dat$b.t[i])




      }else{


        data$U_simu = res$U_simu

        temp_df = summary(cmest(data = data, model = "gformula", outcome = time,
                                exposure = t, mediator = m,
                                basec = c(covariates, glue("U_simu") ),
                                EMint = T, event = outcome,
                                mreg = list("logistic"),
                                yreg = "coxph",
                                astar = 0, a = 1, mval = list( 1), multimp = F, yval = 0,
                                estimation = "imputation", inference = "bootstrap", nboot = B,
                                boot.ci.type = "per"))$summarydf

        me_gsu_df = data.frame(log_nie = log(temp_df$Estimate[5]), p_nie = temp_df$P.val[5], log_nde = log(temp_df$Estimate[2]), p_nde = temp_df$P.val[2] )

        nie_sens = data.frame(medef = me_gsu_df$log_nie, p_ef = me_gsu_df$p_nie,
                              b_o = par_dat$b.y[i], b_m = par_dat$b.m[i], b_x = par_dat$b.t[i] )


        nde_sens = data.frame(medef = me_gsu_df$log_nde, p_ef = me_gsu_df$p_nde,
                              b_o = par_dat$b.y[i], b_m = par_dat$b.m[i], b_x = par_dat$b.t[i] )


      }









      return(list(nie_sens = nie_sens, nde_sens = nde_sens))
    })
  } else{

    cl = makeCluster(detectCores())

    clusterExport(cl, c("gen_p", "data", "par_dat", "p.U", "outcome", "m", "t", "time", "covariates", "scale.m", "scale.t", "scale.u", "sigma.u", "iter_num", "multiple_draw", "inter", "B", "R", "lower_u", "upper_u"),
                  envir=environment())
    clusterEvalQ(cl, c(library(dplyr), library(CMAverse), library(glue), library(survival), library(distr)))

    all_res = parLapply(cl, 1:nrow(par_dat), function(i){

      res = gen_p(data, inter, p.U, outcome, time, m, t, covariates, b.y = par_dat$b.y[i], b.m = par_dat$b.m[i], b.t = par_dat$b.t[i], scale.m, scale.t, scale.u, sigma.u,iter_num, lower_u,upper_u)


      if(scale.u == "binary" & multiple_draw == T){

        data <- bind_cols(data, as.data.frame(
          setNames(
            replicate(R, rbinom(nrow(data), 1, res$p), simplify = FALSE),
            paste0("U_simu", 1:R)
          )
        ))


        me_giu = summary(cmest(data = data, model = "gformula", outcome = time,
                               exposure = t, mediator = m,
                               basec = c(covariates),
                               EMint = inter, event = outcome,
                               mreg = list("logistic"),
                               yreg = "coxph",
                               astar = 0, a = 1, mval = list( 1), multimp = F, yval = 0,
                               estimation = "imputation", inference = "bootstrap", nboot = B,
                               boot.ci.type = "per"))$summarydf


        me_gsu_df = lapply(1:R, function(j){



          temp_df = lapply(1:B, function(B){

            B_res = summary(cmest(data = data, model = "gformula", outcome = time,
                                  exposure = t, mediator = m,
                                  basec = c(covariates, glue("U_simu{j}") ),
                                  EMint = T, event = outcome,
                                  mreg = list("logistic"),
                                  yreg = "coxph",
                                  astar = 0, a = 1, mval = list( 1), multimp = F, yval = 0,
                                  estimation = "imputation", inference = "bootstrap", nboot = B,
                                  boot.ci.type = "per"))$summarydf
            B_df = data.frame(log_nde = log(B_res[2,1]), log_nie = log(B_res[5,1]))
            return(B_df)


          }) %>% bind_rows()

          temp_res = data.frame(log_nie = mean(temp_df$log_nie), sig_log_nie = sd(temp_df$log_nie), log_nde = mean(temp_df$log_nde), sig_log_nde = sd(temp_df$log_nde))


          return(temp_res)

        }) %>% bind_rows()

        nie_sens = data.frame(medef = mean(me_gsu_df$log_nie), medef_std = (mean(me_gsu_df$sig_log_nie^2)+(1+1/R)*var(me_gsu_df$log_nie)   )^0.5 ) %>%
          mutate(p_ef = 2 * (1 - pnorm(abs(medef/medef_std))),
                 b_o = par_dat$b.y[i], b_m = par_dat$b.m[i], b_x = par_dat$b.t[i])


        nde_sens = data.frame(medef = mean(me_gsu_df$log_nde), medef_std = (mean(me_gsu_df$sig_log_nde^2)+(1+1/R)*var(me_gsu_df$log_nde)   )^0.5 ) %>%
          mutate(p_ef = 2 * (1 - pnorm(abs(medef/medef_std))),
                 b_o = par_dat$b.y[i], b_m = par_dat$b.m[i], b_x = par_dat$b.t[i])




      }else{


        data$U_simu = res$U_simu

        temp_df = summary(cmest(data = data, model = "gformula", outcome = time,
                                exposure = t, mediator = m,
                                basec = c(covariates, glue("U_simu") ),
                                EMint = T, event = outcome,
                                mreg = list("logistic"),
                                yreg = "coxph",
                                astar = 0, a = 1, mval = list( 1), multimp = F, yval = 0,
                                estimation = "imputation", inference = "bootstrap", nboot = B,
                                boot.ci.type = "per"))$summarydf

        me_gsu_df = data.frame(log_nie = log(temp_df$Estimate[5]), p_nie = temp_df$P.val[5], log_nde = log(temp_df$Estimate[2]), p_nde = temp_df$P.val[2] )

        nie_sens = data.frame(medef = me_gsu_df$log_nie, p_ef = me_gsu_df$p_nie,
                              b_o = par_dat$b.y[i], b_m = par_dat$b.m[i], b_x = par_dat$b.t[i] )


        nde_sens = data.frame(medef = me_gsu_df$log_nde, p_ef = me_gsu_df$p_nde,
                              b_o = par_dat$b.y[i], b_m = par_dat$b.m[i], b_x = par_dat$b.t[i] )


      }









      return(list(nie_sens = nie_sens, nde_sens = nde_sens))
    })
  }

  nie_sens_df = lapply(1:length(all_res), function(k){
    all_res[[k]]$nie_sens

  }) %>% bind_rows()

  nde_sens_df = lapply(1:length(all_res), function(k){
    all_res[[k]]$nde_sens

  }) %>% bind_rows()


  nie_rev = nie_sens_df[nie_sens_df$p_ef>0.05,]

  nde_rev = nde_sens_df[nde_sens_df$p_ef>0.05,]

  edges <- list(
    # —— bottom）——
    list(x = range.b.t,        y = rep(range.b.y[1], 2), z = rep(range.b.m[1], 2)),
    list(x = range.b.t,        y = rep(range.b.y[2], 2), z = rep(range.b.m[1], 2)),
    list(x = rep(range.b.t[1], 2), y = range.b.y,        z = rep(range.b.m[1], 2)),
    list(x = rep(range.b.t[2], 2), y = range.b.y,        z = rep(range.b.m[1], 2)),

    # —— top——
    list(x = range.b.t,        y = rep(range.b.y[1], 2), z = rep(range.b.m[2], 2)),
    list(x = range.b.t,        y = rep(range.b.y[2], 2), z = rep(range.b.m[2], 2)),
    list(x = rep(range.b.t[1], 2), y = range.b.y,        z = rep(range.b.m[2], 2)),
    list(x = rep(range.b.t[2], 2), y = range.b.y,        z = rep(range.b.m[2], 2)),

    # —— side）——
    list(x = rep(range.b.t[1], 2), y = rep(range.b.y[1], 2), z = range.b.m),
    list(x = rep(range.b.t[2], 2), y = rep(range.b.y[1], 2), z = range.b.m),
    list(x = rep(range.b.t[1], 2), y = rep(range.b.y[2], 2), z = range.b.m),
    list(x = rep(range.b.t[2], 2), y = rep(range.b.y[2], 2), z = range.b.m)
  )

  # plot 3D scatter
  p_nie <- plot_ly() %>%
    add_trace(
      data = nie_rev,
      x = ~b_m, y = ~b_o, z = ~b_x,
      type = "scatter3d",
      mode = "markers",
      marker = list(size = 2, color = 'red'),
      name = "points"
    )

  p_nde <- plot_ly() %>%
    add_trace(
      data = nde_rev,
      x = ~b_m, y = ~b_o, z = ~b_x,
      type = "scatter3d",
      mode = "markers",
      marker = list(size = 2, color = 'red'),
      name = "points"
    )

  # edges, lines
  for (e in edges) {
    p_nie <- p_nie %>%
      add_trace(
        x = e$x, y = e$y, z = e$z,
        type = "scatter3d",
        mode = "lines",
        line = list(color = 'black', width = 2),
        showlegend = FALSE
      )

    p_nde <- p_nde %>%
      add_trace(
        x = e$x, y = e$y, z = e$z,
        type = "scatter3d",
        mode = "lines",
        line = list(color = 'black', width = 2),
        showlegend = FALSE
      )

  }

  # 5. 等比例立方体、范围一致、刻度一致，且去掉 0 处的零线
  p_nie <- p_nie %>%
    layout(
      scene = list(
        aspectmode = 'cube',
        xaxis = list(
          title   = "Beta of U-mediator",
          range   = range.b.m,
          tick0   = range.b.m[1],
          dtick   = (range.b.m[2]-range.b.m[1])/5,
          zeroline = FALSE     # remove zero line
        ),
        yaxis = list(
          title   = "Beta of U-outcome",
          range   = range.b.y,
          tick0   = range.b.y[1],
          dtick   = (range.b.y[2]-range.b.y[1])/5,
          zeroline = FALSE
        ),
        zaxis = list(
          title   = "Beta of U-exposure",
          range   = range.b.t,
          tick0   = range.b.t[1],
          dtick   = (range.b.t[2]-range.b.t[1])/5,
          zeroline = FALSE
        )
      )
    )

  p_nde <- p_nde %>%
    layout(
      scene = list(
        aspectmode = 'cube',
        xaxis = list(
          title   = "Beta of U-mediator",
          range   = range.b.m,
          tick0   = range.b.m[1],
          dtick   = (range.b.m[2]-range.b.m[1])/5,
          zeroline = FALSE     # remove zero line
        ),
        yaxis = list(
          title   = "Beta of U-outcome",
          range   = range.b.y,
          tick0   = range.b.y[1],
          dtick   = (range.b.y[2]-range.b.y[1])/5,
          zeroline = FALSE
        ),
        zaxis = list(
          title   = "Beta of U-exposure",
          range   = range.b.t,
          tick0   = range.b.t[1],
          dtick   = (range.b.t[2]-range.b.t[1])/5,
          zeroline = FALSE
        )
      )
    )


  return(list(nie_sens_df = nie_sens_df, nde_sens_df = nde_sens_df, p_nie = p_nie, p_nde = p_nde))

}
