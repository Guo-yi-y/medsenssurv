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
#' @param scale.m,scale.t scale of mediator and exposure, binary or continuous
#' @param iter_num number of iteration
#' @return conditional probability of U
#' @importFrom stats rbinom dnorm gaussian binomial glm as.formula model.matrix
#' @importFrom stats pnorm sd var sigma setNames
#' @importFrom survival coxph basehaz
#' @importFrom dplyr bind_cols bind_rows mutate
#' @importFrom magrittr %>%
#' @importFrom glue glue
#' @importFrom CMAverse cmest
#' @importFrom parallel makeCluster detectCores clusterExport clusterEvalQ parLapply
#' @importFrom plotly plot_ly add_trace layout
#' @export

gen_p = function(data, inter = FALSE, p.U, outcome, time, m, t, covariates, b.y, b.m, b.t, scale.m, scale.t, iter_num = 20){
  data$U_ini = rbinom(nrow(data), 1, p.U)
  theta = p.U
  data$p = p.U
  covariates_fm = paste({covariates}, collapse = '+')
  if(inter == FALSE){

    out_formula = glue("Surv({time}, {outcome})~{t}+{m}+{covariates_fm}+U_ini")
    med_formula = glue("{m}~{t}+{covariates_fm}+U_ini+{covariates_fm}")
    exp_formula = glue("{t}~{covariates_fm}+U_ini")

    #fit outcome
    o_fit <- coxph(as.formula(out_formula) , data = data)
    #o_coef <- c(o_fit$coef, b.y)
    o_coef <- c(o_fit$coef)
    o_coef['U_ini'] = b.y


    #fit m
    if (scale.m == "continuous"){
      m_fit <- glm(as.formula(med_formula), data = data, family = gaussian)
    } else{
      m_fit <- glm(as.formula(med_formula), data = data, family = binomial(link = "logit"))
    }

    #m_coef <- c(m_fit$coef, b.m)
    m_coef <- c(m_fit$coef)
    m_coef['U_ini'] = b.m


    #fit t
    if (scale.t == "continuous"){
      x_fit <- glm(as.formula(exp_formula) , data = data, family = gaussian)
    } else{
      x_fit <- glm(as.formula(exp_formula) , data = data, family = binomial(link = "logit"))
    }

    #x_coef <- c(x_fit$coef, b.t)
    x_coef <- c(x_fit$coef)
    x_coef['U_ini'] = b.t


    o_coef[is.na(o_coef)] <- 0
    m_coef[is.na(m_coef)] <- 0
    x_coef[is.na(x_coef)] <- 0

    p <- rep(0,nrow(data))

    bh1 = basehaz(o_fit, centered=F) #cumulative baseline hazard for event with mean(offset)
    bh1$time = round(bh1$time, 3)
    index1 = match(data[, time], bh1$time)

    if (scale.t == "continuous"){
      mean.xu1 = (cbind(model.matrix(x_fit)%>% .[, -ncol(.)], U_simu = 1) %*% matrix(x_coef, ncol=1) )
      sd.xu1 = sigma(x_fit)
      mean.xu0 = (cbind(model.matrix(x_fit)%>% .[, -ncol(.)], U_simu = 0) %*% matrix(x_coef, ncol=1)  )
      sd.xu0 = sigma(x_fit)
      pxu1 = dnorm(data[, t], mean = mean.xu1, sd = sd.xu1)
      pxu0 = dnorm(data[, t], mean = mean.xu0, sd = sd.xu0)
    } else{
      px1u1 = 1/(1 + exp(-(cbind(model.matrix(x_fit)%>% .[, -ncol(.)], U_simu = 1) %*% matrix(x_coef, ncol=1) )))
      px1u0 = 1/(1 + exp(-(cbind(model.matrix(x_fit)%>% .[, -ncol(.)], U_simu = 0) %*% matrix(x_coef, ncol=1)  )))
      pxu1 = px1u1^data[, t] * (1 - px1u1)^(1 - data[, t])
      pxu0 = px1u0^data[, t] * (1 - px1u0)^(1 - data[, t])
    }


    if (scale.m == "continuous"){
      mean.mu1 = (cbind(model.matrix(m_fit)%>% .[, -ncol(.)], U_simu = 1) %*% matrix(m_coef, ncol=1) )
      sd.mu1 = sigma(m_fit)
      mean.mu0 = (cbind(model.matrix(m_fit)%>% .[, -ncol(.)], U_simu = 0) %*% matrix(m_coef, ncol=1)  )
      sd.mu0 = sigma(m_fit)
      pmu1 = dnorm(data[, m], mean = mean.mu1, sd = sd.mu1)
      pmu0 = dnorm(data[, m], mean = mean.mu0, sd = sd.mu0)
    } else {
      pm1u1 = 1/(1 + exp(-(cbind(model.matrix(m_fit)%>% .[, -ncol(.)], U_simu = 1) %*% matrix(m_coef, ncol=1) )))
      pm1u0 = 1/(1 + exp(-(cbind(model.matrix(m_fit)%>% .[, -ncol(.)], U_simu = 0) %*% matrix(m_coef, ncol=1)  )))
      pmu1 = pm1u1^data[, m] * (1 - pm1u1)^(1 - data[, m])
      pmu0 = pm1u0^data[, m] * (1 - pm1u0)^(1 - data[, m])
    }



    pyu1 = exp( cbind(model.matrix(o_fit)%>% .[, -ncol(.)], U_simu = 1) %*% o_coef )^data[, outcome] *
      exp(-bh1[index1,1]/exp(mean(log(exp(b.y)*p+(1-p)))) * exp( cbind(model.matrix(o_fit)%>% .[, -ncol(.)], U_simu = 1) %*% o_coef ) )


    pyu0 = exp( cbind(model.matrix(o_fit)%>% .[, -ncol(.)], U_simu = 0) %*% o_coef )^data[, outcome] *
      exp(-bh1[index1,1]/exp(mean(log(exp(b.y)*p+(1-p)))) * exp( cbind(model.matrix(o_fit)%>% .[, -ncol(.)], U_simu = 0) %*% o_coef ) )


    pu1 = pmu1*theta*pyu1*pxu1
    pu0 = pmu0*(1-theta)*pyu0*pxu0



    p = as.vector(pu1/(pu1 + pu0)) #posterior dist of U

    p[pu1==0 & pu0==0] = 0


    U_simu = rbinom(nrow(data), 1, p)

    data$p = p

    #fit time to event

    o_fit = coxph(as.formula(glue("Surv({time}, {outcome})~{t}+{m}+{covariates_fm}+ offset(log(exp(b.y)*p+(1-p)))")), data = data)

    o_coef = c(o_fit$coef)

    o_coef["U_simu"] = b.y

    #fit t
    if (scale.t == "continuous"){
      x_fit = glm(as.formula(glue("{t}~{covariates_fm}+ offset(b.t*U_simu)")), data = data, family = gaussian)
    } else{
      x_fit = glm(as.formula(glue("{t}~{covariates_fm}+ offset(b.t*U_simu)")), data = data, family = binomial(link = "logit"))
    }

    x_coef = c(x_fit$coefficients)

    x_coef["U_simu"] = b.t

    if (scale.m == "continuous"){
      m_fit = glm(as.formula(glue("{m}~{t}+{covariates_fm}+ offset(b.m*U_simu)")), data = data, family = gaussian)
    } else{
      m_fit = glm(as.formula(glue("{m}~{t}+{covariates_fm}+ offset(b.m*U_simu)")), data = data, family = binomial(link = "logit"))
    }

    m_coef = c(m_fit$coefficients)

    m_coef["U_simu"] = b.m

    o_coef[is.na(o_coef)] = 0
    x_coef[is.na(x_coef)] = 0
    m_coef[is.na(m_coef)] = 0


    data$U_simu <- rbinom(nrow(data), 1, p)



    n_iter = 0



    while( n_iter<=iter_num) {
      bh1 = basehaz(o_fit, centered=F) #cumulative baseline hazard for event with mean(offset)
      bh1$time = round(bh1$time, 3)
      index1 = match(data[, time], bh1$time)

      if (scale.t == "continuous"){
        mean.xu1 = (cbind(model.matrix(x_fit), U_simu = 1) %*% matrix(x_coef, ncol=1) )
        sd.xu1 = sigma(x_fit)
        mean.xu0 = (cbind(model.matrix(x_fit), U_simu = 0) %*% matrix(x_coef, ncol=1)  )
        sd.xu0 = sigma(x_fit)
        pxu1 = dnorm(data[, t], mean = mean.xu1, sd = sd.xu1)
        pxu0 = dnorm(data[, t], mean = mean.xu0, sd = sd.xu0)
      } else{
        px1u1 = 1/(1 + exp(-(cbind(model.matrix(x_fit), U_simu = 1) %*% matrix(x_coef, ncol=1) )))
        px1u0 = 1/(1 + exp(-(cbind(model.matrix(x_fit), U_simu = 0) %*% matrix(x_coef, ncol=1)  )))
        pxu1 = px1u1^data[, t] * (1 - px1u1)^(1 - data[, t])
        pxu0 = px1u0^data[, t] * (1 - px1u0)^(1 - data[, t])
      }


      if (scale.m == "continuous"){
        mean.mu1 = (cbind(model.matrix(m_fit), U_simu = 1) %*% matrix(m_coef, ncol=1) )
        sd.mu1 = sigma(m_fit)
        mean.mu0 = (cbind(model.matrix(m_fit), U_simu = 0) %*% matrix(m_coef, ncol=1)  )
        sd.mu0 = sigma(m_fit)
        pmu1 = dnorm(data[, m], mean = mean.mu1, sd = sd.mu1)
        pmu0 = dnorm(data[, m], mean = mean.mu0, sd = sd.mu0)
      } else {
        pm1u1 = 1/(1 + exp(-(cbind(model.matrix(m_fit), U_simu = 1) %*% matrix(m_coef, ncol=1) )))
        pm1u0 = 1/(1 + exp(-(cbind(model.matrix(m_fit), U_simu = 0) %*% matrix(m_coef, ncol=1)  )))
        pmu1 = pm1u1^data[, m] * (1 - pm1u1)^(1 - data[, m])
        pmu0 = pm1u0^data[, m] * (1 - pm1u0)^(1 - data[, m])
      }


      pyu1 = exp( cbind(model.matrix(o_fit), U_simu = 1) %*% o_coef )^data[, outcome] *
        exp(-bh1[index1,1]/exp(mean(log(exp(b.y)*p+(1-p)))) * exp( cbind(model.matrix(o_fit), U_simu = 1) %*% o_coef ) )


      pyu0 = exp( cbind(model.matrix(o_fit), U_simu = 0) %*% o_coef )^data[, outcome] *
        exp(-bh1[index1,1]/exp(mean(log(exp(b.y)*p+(1-p)))) * exp( cbind(model.matrix(o_fit), U_simu = 0) %*% o_coef ) )


      pu1 = pmu1*theta*pyu1*pxu1
      pu0 = pmu0*(1-theta)*pyu0*pxu0



      p = as.vector(pu1/(pu1 + pu0)) #posterior dist of U
      p[pu1==0 & pu0==0] = 0


      data$U_simu = rbinom(nrow(data), 1, p)

      data$p = p


      o_fit = coxph(as.formula(glue("Surv({time}, {outcome})~{t}+{m}+{covariates_fm}+ offset(log(exp(b.y)*p+(1-p)))")), data = data)

      o_coef = c(o_fit$coef)

      o_coef["U_simu"] = b.y

      #fit t
      if (scale.t == "continuous"){
        x_fit = glm(as.formula(glue("{t}~{covariates_fm}+ offset(b.t*U_simu)")), data = data, family = gaussian)
      } else {
        x_fit = glm(as.formula(glue("{t}~{covariates_fm}+ offset(b.t*U_simu)")), data = data, family = binomial(link = "logit"))
      }

      x_coef = c(x_fit$coefficients)

      x_coef["U_simu"] = b.t

      if (scale.m == "continuous"){
        m_fit = glm(as.formula(glue("{m}~{t}+{covariates_fm}+ offset(b.m*U_simu)")), data = data, family = gaussian)
      } else{
        m_fit = glm(as.formula(glue("{m}~{t}+{covariates_fm}+ offset(b.m*U_simu)")), data = data, family = binomial(link = "logit"))
      }

      m_coef = c(m_fit$coefficients)

      m_coef["U_simu"] = b.m

      o_coef[is.na(o_coef)] = 0
      x_coef[is.na(x_coef)] = 0
      m_coef[is.na(m_coef)] = 0


      data$U_simu <- rbinom(nrow(data), 1, p)



      n_iter= n_iter+1
    }



  } else {  #fit outcome

    out_formula = glue("Surv({time}, {outcome})~{t}+{m}+{t}:{m}+{covariates_fm}+U_ini")
    med_formula = glue("{m}~{t}+{covariates_fm}+U_ini")
    exp_formula = glue("{t}~{covariates_fm}+U_ini")

    #fit outcome
    o_fit <- coxph(as.formula(out_formula) , data = data)
    #o_coef <- c(o_fit$coef, b.y)
    o_coef <- c(o_fit$coef)
    o_coef['U_ini'] = b.y


    #fit m
    if (scale.m == "continuous"){
      m_fit <- glm(as.formula(med_formula), data = data, family = gaussian)
    } else{
      m_fit <- glm(as.formula(med_formula), data = data, family = binomial(link = "logit"))
    }

    #m_coef <- c(m_fit$coef, b.m)
    m_coef <- c(m_fit$coef)
    m_coef['U_ini'] = b.m


    #fit t
    if (scale.t == "continuous"){
      x_fit <- glm(as.formula(exp_formula) , data = data, family = gaussian)
    } else{
      x_fit <- glm(as.formula(exp_formula) , data = data, family = binomial(link = "logit"))
    }
    #x_coef <- c(x_fit$coef, b.t)
    x_coef <- c(x_fit$coef)
    x_coef['U_ini'] = b.t


    o_coef[is.na(o_coef)] <- 0
    m_coef[is.na(m_coef)] <- 0
    x_coef[is.na(x_coef)] <- 0

    p <- rep(0,nrow(data))



    bh1 = basehaz(o_fit, centered=F) #cumulative baseline hazard for event with mean(offset)
    bh1$time = round(bh1$time, 3)
    index1 = match(data[, time], bh1$time)

    if (scale.t == "continuous"){
      mean.xu1 = (cbind(model.matrix(x_fit)%>% .[, -ncol(.)], U_simu = 1) %*% matrix(x_coef, ncol=1) )
      sd.xu1 = sigma(x_fit)
      mean.xu0 = (cbind(model.matrix(x_fit)%>% .[, -ncol(.)], U_simu = 0) %*% matrix(x_coef, ncol=1)  )
      sd.xu0 = sigma(x_fit)
      pxu1 = dnorm(data[, t], mean = mean.xu1, sd = sd.xu1)
      pxu0 = dnorm(data[, t], mean = mean.xu0, sd = sd.xu0)
    } else{
      px1u1 = 1/(1 + exp(-(cbind(model.matrix(x_fit)%>% .[, -ncol(.)], U_simu = 1) %*% matrix(x_coef, ncol=1) )))
      px1u0 = 1/(1 + exp(-(cbind(model.matrix(x_fit)%>% .[, -ncol(.)], U_simu = 0) %*% matrix(x_coef, ncol=1)  )))
      pxu1 = px1u1^data[, t] * (1 - px1u1)^(1 - data[, t])
      pxu0 = px1u0^data[, t] * (1 - px1u0)^(1 - data[, t])
    }


    if (scale.m == "continuous"){
      mean.mu1 = (cbind(model.matrix(m_fit)%>% .[, -ncol(.)], U_simu = 1) %*% matrix(m_coef, ncol=1) )
      sd.mu1 = sigma(m_fit)
      mean.mu0 = (cbind(model.matrix(m_fit)%>% .[, -ncol(.)], U_simu = 0) %*% matrix(m_coef, ncol=1)  )
      sd.mu0 = sigma(m_fit)
      pmu1 = dnorm(data[, m], mean = mean.mu1, sd = sd.mu1)
      pmu0 = dnorm(data[, m], mean = mean.mu0, sd = sd.mu0)
    } else {
      pm1u1 = 1/(1 + exp(-(cbind(model.matrix(m_fit)%>% .[, -ncol(.)], U_simu = 1) %*% matrix(m_coef, ncol=1) )))
      pm1u0 = 1/(1 + exp(-(cbind(model.matrix(m_fit)%>% .[, -ncol(.)], U_simu = 0) %*% matrix(m_coef, ncol=1)  )))
      pmu1 = pm1u1^data[, m] * (1 - pm1u1)^(1 - data[, m])
      pmu0 = pm1u0^data[, m] * (1 - pm1u0)^(1 - data[, m])
    }


    pyu1 = exp( cbind(model.matrix(o_fit)%>% .[, -ncol(.)], U_simu = 1) %*% o_coef )^data[, outcome] *
      exp(-bh1[index1,1]/exp(mean(log(exp(b.y)*p+(1-p)))) * exp( cbind(model.matrix(o_fit)%>% .[, -ncol(.)], U_simu = 1) %*% o_coef ) )


    pyu0 = exp( cbind(model.matrix(o_fit)%>% .[, -ncol(.)], U_simu = 0) %*% o_coef )^data[, outcome] *
      exp(-bh1[index1,1]/exp(mean(log(exp(b.y)*p+(1-p)))) * exp( cbind(model.matrix(o_fit)%>% .[, -ncol(.)], U_simu = 0) %*% o_coef ) )


    pu1 = pmu1*theta*pyu1*pxu1
    pu0 = pmu0*(1-theta)*pyu0*pxu0



    p = as.vector(pu1/(pu1 + pu0)) #posterior dist of U
    p[pu1==0 & pu0==0] = 0


    data$U_simu = rbinom(nrow(data), 1, p)

    data$p = p

    #fit time to event

    o_fit = coxph(as.formula(glue("Surv({time}, {outcome})~{t}+{m}+{t}:{m}+{covariates_fm}++ offset(log(exp(b.y)*p+(1-p)))")), data = data)

    o_coef = c(o_fit$coef)

    o_coef["U_simu"] = b.y

    #fit t
    if (scale.t == "continuous"){
      x_fit = glm(as.formula(glue("{t}~{covariates_fm}+ offset(b.t*U_simu)")), data = data, family = gaussian)
    } else{
      x_fit = glm(as.formula(glue("{t}~{covariates_fm}+ offset(b.t*U_simu)")), data = data, family = binomial(link = "logit"))
    }

    x_coef = c(x_fit$coefficients)

    x_coef["U_simu"] = b.t

    if (scale.m == "continuous"){
      m_fit = glm(as.formula(glue("{m}~{t}+{covariates_fm}+ offset(b.m*U_simu)")), data = data, family = gaussian)
    } else{
      m_fit = glm(as.formula(glue("{m}~{t}+{covariates_fm}+ offset(b.m*U_simu)")), data = data, family = binomial(link = "logit"))
    }

    m_coef = c(m_fit$coefficients)

    m_coef["U_simu"] = b.m

    o_coef[is.na(o_coef)] = 0
    x_coef[is.na(x_coef)] = 0
    m_coef[is.na(m_coef)] = 0




    data$U_simu <- rbinom(nrow(data), 1, p)

    n_iter = 0

    while(  n_iter<=iter_num) {
      bh1 = basehaz(o_fit, centered=F) #cumulative baseline hazard for event with mean(offset)
      bh1$time = round(bh1$time, 3)
      index1 = match(data[, time], bh1$time)

      if (scale.t == "continuous"){
        mean.xu1 = (cbind(model.matrix(x_fit), U_simu = 1) %*% matrix(x_coef, ncol=1) )
        sd.xu1 = sigma(x_fit)
        mean.xu0 = (cbind(model.matrix(x_fit), U_simu = 0) %*% matrix(x_coef, ncol=1)  )
        sd.xu0 = sigma(x_fit)
        pxu1 = dnorm(data[, t], mean = mean.xu1, sd = sd.xu1)
        pxu0 = dnorm(data[, t], mean = mean.xu0, sd = sd.xu0)
      } else{
        px1u1 = 1/(1 + exp(-(cbind(model.matrix(x_fit), U_simu = 1) %*% matrix(x_coef, ncol=1) )))
        px1u0 = 1/(1 + exp(-(cbind(model.matrix(x_fit), U_simu = 0) %*% matrix(x_coef, ncol=1)  )))
        pxu1 = px1u1^data[, t] * (1 - px1u1)^(1 - data[, t])
        pxu0 = px1u0^data[, t] * (1 - px1u0)^(1 - data[, t])
      }


      if (scale.m == "continuous"){
        mean.mu1 = (cbind(model.matrix(m_fit), U_simu = 1) %*% matrix(m_coef, ncol=1) )
        sd.mu1 = sigma(m_fit)
        mean.mu0 = (cbind(model.matrix(m_fit), U_simu = 0) %*% matrix(m_coef, ncol=1)  )
        sd.mu0 = sigma(m_fit)
        pmu1 = dnorm(data[, m], mean = mean.mu1, sd = sd.mu1)
        pmu0 = dnorm(data[, m], mean = mean.mu0, sd = sd.mu0)
      } else {
        pm1u1 = 1/(1 + exp(-(cbind(model.matrix(m_fit), U_simu = 1) %*% matrix(m_coef, ncol=1) )))
        pm1u0 = 1/(1 + exp(-(cbind(model.matrix(m_fit), U_simu = 0) %*% matrix(m_coef, ncol=1)  )))
        pmu1 = pm1u1^data[, m] * (1 - pm1u1)^(1 - data[, m])
        pmu0 = pm1u0^data[, m] * (1 - pm1u0)^(1 - data[, m])
      }


      pyu1 = exp( cbind(model.matrix(o_fit), U_simu = 1) %*% o_coef )^data[, outcome] *
        exp(-bh1[index1,1]/exp(mean(log(exp(b.y)*p+(1-p)))) * exp( cbind(model.matrix(o_fit), U_simu = 1) %*% o_coef ) )


      pyu0 = exp( cbind(model.matrix(o_fit), U_simu = 0) %*% o_coef )^data[, outcome] *
        exp(-bh1[index1,1]/exp(mean(log(exp(b.y)*p+(1-p)))) * exp( cbind(model.matrix(o_fit), U_simu = 0) %*% o_coef ) )


      pu1 = pmu1*theta*pyu1*pxu1
      pu0 = pmu0*(1-theta)*pyu0*pxu0



      p = as.vector(pu1/(pu1 + pu0)) #posterior dist of U
      p[pu1==0 & pu0==0] = 0


      data$U_simu = rbinom(nrow(data), 1, p)

      data$p = p

      #fit time to event

      o_fit = coxph(as.formula(glue("Surv({time}, {outcome})~{t}+{m}+{t}:{m}+{covariates_fm}++ offset(log(exp(b.y)*p+(1-p)))")), data = data)

      o_coef = c(o_fit$coef)

      o_coef["U_simu"] = b.y

      #fit t
      if (scale.t == "continuous"){
        x_fit = glm(as.formula(glue("{t}~{covariates_fm}+ offset(b.t*U_simu)")), data = data, family = gaussian)
      } else {
        x_fit = glm(as.formula(glue("{t}~{covariates_fm}+ offset(b.t*U_simu)")), data = data, family = binomial(link = "logit"))
      }

      x_coef = c(x_fit$coefficients)

      x_coef["U_simu"] = b.t

      if (scale.m == "continuous"){
        m_fit = glm(as.formula(glue("{m}~{t}+{covariates_fm}+ offset(b.m*U_simu)")), data = data, family = gaussian)
      } else {
        m_fit = glm(as.formula(glue("{m}~{t}+{covariates_fm}+ offset(b.m*U_simu)")), data = data, family = binomial(link = "logit"))
      }

      m_coef = c(m_fit$coefficients)

      m_coef["U_simu"] = b.m

      o_coef[is.na(o_coef)] = 0
      x_coef[is.na(x_coef)] = 0
      m_coef[is.na(m_coef)] = 0




      data$U_simu <- rbinom(nrow(data), 1, p)



      n_iter = n_iter+1

    }







  }



  return(p)



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
#' @param scale.m,scale.t scale of mediator and exposure, binary or continuous
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
#' @importFrom magrittr %>%
#' @importFrom glue glue
#' @importFrom CMAverse cmest
#' @importFrom parallel makeCluster detectCores clusterExport clusterEvalQ parLapply
#' @importFrom plotly plot_ly add_trace layout
#' @export

medsurvsens = function(data, inter = FALSE, p.U = 0.5, outcome, time, m, t, covariates, scale.m, scale.t, range.b.y, range.b.m, range.b.t, ngrid = 10, iter_num = 20, multiple_draw = F, R = 5, B = 100, parallel = F){


  par_dat = expand.grid(b.y = seq(range.b.y[1], range.b.y[2], length.out = ngrid),
                        b.m = seq(range.b.m[1], range.b.m[2], length.out = ngrid),
                        b.t = seq(range.b.t[1], range.b.t[2], length.out = ngrid)


                        )



  if(parallel == F){
    all_res = lapply(1:nrow(par_dat), function(i){

      p = gen_p(data, inter, p.U, outcome, time, m, t, covariates, b.y = par_dat$b.y[i], b.m = par_dat$b.m[i], b.t = par_dat$b.t[i], scale.m, scale.t, iter_num)



      data <- bind_cols(data, as.data.frame(
        setNames(
          replicate(5, rbinom(nrow(data), 1, p), simplify = FALSE),
          paste0("U_simu", 1:5)
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

      if(multiple_draw == T){
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

      } else{
        temp_df = summary(cmest(data = data, model = "gformula", outcome = time,
                                exposure = t, mediator = m,
                                basec = c(covariates, glue("U_simu{j}") ),
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

    clusterExport(cl, c("gen_p", "data", "par_dat", "p.U", "outcome", "m", "t", "time", "covariates", "scale.m", "scale.t", "iter_num", "multiple_draw", "inter", "B", "R"),
                  envir=environment())
    clusterEvalQ(cl, c(library(dplyr), library(CMAverse), library(glue), library(survival)))

    all_res = parLapply(cl, 1:nrow(par_dat), function(i){

      p = gen_p(data, inter, p.U, outcome, time, m, t, covariates, b.y = par_dat$b.y[i], b.m = par_dat$b.m[i], b.t = par_dat$b.t[i], scale.m, scale.t, iter_num)



      data <- bind_cols(data, as.data.frame(
        setNames(
          replicate(5, rbinom(nrow(data), 1, p), simplify = FALSE),
          paste0("U_simu", 1:5)
        )
      ))



      if(multiple_draw == T){
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

      } else{
        data$U_simu = rbinom(nrow(data), 1, p)
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
