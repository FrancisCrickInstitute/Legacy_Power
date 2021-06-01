library(simr)
library(lme4)
library(tidyverse)
library(rms)
library(emmeans)

## Sample sizes

# unspecified Ns in the case cohort
Rheum1 <- 100 
Rheum2 <- 100
Lymph1 <- 100

# Ns in the treatment disease cohort

n_per_drug_per_cohort <- list(
  case=list(
    Cancer=c(700),
    MS=c(125,153,175),
    Rheumatology=c(180,Rheum1, Rheum2),
    IBD=c(500,197,19,163),
    Lymphoma=c(Lymph1)
  )
)

# Ns in the control cohort
cc_prop  = 0.2; cc_const = 20 # multiply by prop, add const, to get nunber of controls
# 0.2 and 0 above give us a total of 500 controls - as would 0 and 100
#cc_prop = 0; cc_const = 100
n_per_drug_per_cohort$control <- lapply(
  n_per_drug_per_cohort$case, function(x) sum(x)*cc_prop + cc_const
)

# scale back to 2000 in total
n_per_drug_per_cohort <- lapply(n_per_drug_per_cohort, function(x) {lapply(x, function(y) y * 2000/3114)})
## Effect sizes


# for controls with all doses, response in the vaccinated-for variant
# odds scale. So e.g. 19 = 95% chance success, 5% failure
# final_vaccine_odds <- percent / (100-percent)
prob2odd <- function(p) p/(1-p)
odd2prob <- function(o) o/(1+o)

outcomes <- expand.grid(strain=c("vaccine", "voc"),
                       cohort=c("control","case"))
baseline=0.05
voc <- 0.25
case <- 0.25
outcomes$absolute <- c(baseline, case, voc, case+voc-baseline)
outcomes$relative <- c(baseline, case, voc, case * voc / baseline)
outcomes$odds <- c(baseline, case, voc, odd2prob(prob2odd(voc) * prob2odd(case)/prob2odd(baseline) ))
  

final_vaccine_odds <- prob2odd(0.90) # 

# Relative for the other variants
variant_multiplier <- list(
  vaccine_target = 1,
  voc = prob2odd(.7)/final_vaccine_odds
)

# Overall case effect
case_effect <- list(
  control = 1,
  case    = prob2odd(.7)/final_vaccine_odds
)

# effect of interest
eoi <- list(mult = 0.5, # 0.5 is about half the log-odds effect size of 0.7
           target = quote(variant=="voc" & cohort=="case"))

# sd of effect across diseases
# either side of .70. 95% confit is approx double stderr, so +/- 10 percentage points
# half the width, and then get average of extent
# (prob2odd(.75)/final_vaccine_odds - prob2odd(.65)/final_vaccine_odds )/2
# is 0.06.  

disease_sd <- 0.06

# sd of treatment effect within disease

drug_sd <- 0.02



## Build data's covariates

case_control <- imap(n_per_drug_per_cohort,
                    function(cohort, i) {
                      drug_names <- 
                        paste(i,
                              rep(names(cohort), times=sapply(cohort, length)),
                              do.call(c, lapply(cohort, function(x) seq(along=x))),
                              sep="_")
                      data.frame(
                          disease = rep(names(cohort),
                                        times=sapply(cohort, sum)),
                          drug = rep(drug_names,
                                     times=unlist(cohort)),
                          cohort=i
                        )
                      }
                    )
case_control$control$drug <- "none"
case_control_frame <- do.call(rbind, case_control)


  
measure_frame <- expand.grid(
  variant = names(variant_multiplier)
)

ind <- expand.grid(
  i_cc = 1:nrow(case_control_frame),
  i_measure = 1:nrow(measure_frame)
)

case_control_frame <- cbind(
  case_control_frame[ind$i_cc,],
  measure_frame[ind$i_measure,,drop=FALSE]
  )


## Fill in fixed effects multipliers
case_control_frame <- case_control_frame %>%
  dplyr::mutate(
    odds = final_vaccine_odds *
      unlist(variant_multiplier[as.character(variant)])
  )


bootRes <- data.frame(i=1:100,
                     p=0.0,
                     effect=0.0)
diseaseRes <- lapply(n_per_drug_per_cohort$case,
                    function(x) bootRes)


for (i in 1:nrow(bootRes)) { 
  cat(i, " ")
  ## Get instance of random effects:
  diseases <- names(n_per_drug_per_cohort$case)
  offset <- rnorm(length(diseases), 0, disease_sd)
  disease_effect <- lapply(
    case_effect,
    function(mult) {setNames(offset+mult, diseases)}
  )
    
  drug_effect <- imap(disease_effect, function(diseases, cc) {
    imap(diseases, function(mult, disease) {
      if (cc=="case") {
        x <- rnorm(n=length(n_per_drug_per_cohort[[cc]][[disease]]),
                  mult,
                  drug_sd)
        setNames(x, paste(cc, disease, seq(along=n_per_drug_per_cohort[[cc]][[disease]]), sep="_"))
      } else {
        x <- c("none"=mult)
      }
    }
    )
  }
  )

  drug_effect_frame <-
    map_dfr(drug_effect,
            function(cohort) {
              map_dfr(cohort,
                      function(disease){
                        data.frame(mult=disease, drug=names(disease))
                      },
                      .id="disease")
            },
            .id="cohort")
  dat <- dplyr::inner_join(case_control_frame, drug_effect_frame)
  dat <- mutate(
    dat,
    mult = mult * ifelse(eval(eoi$target), eoi$mult, 1)
  )
  dat <- mutate(dat,
               response = rlogis(n()) < log(odds*mult))
  if (alternative_fits <- FALSE) {
    dat_agg <- dat %>%
      dplyr::select(-odds, -mult) %>%
      group_by(cohort, disease, drug, variant) %>%
      summarise(respond=sum(response),
                n=length(response),
                .groups="drop"
                ) %>%
      mutate(prob_response=respond/n)
    # ggplot(dat_agg, aes(x=variant, y=prob_response, colour=cohort)) + geom_point() + facet_wrap(~drug)
    fit <- lrm(response ~ cohort* variant, data=dat, x=TRUE, y=TRUE)
    g <- robcov(fit)
    h <- bootcov(fit, cluster=dat$drug)
    fit <- ols(prob_response~cohort * variant, data=dat_agg, x=TRUE, y=TRUE)
    h <- bootcov(fit, cluster=dat_agg$drug)
  }
  fit <- glmer(response ~ cohort * variant +   (1|disease), data=dat, family=binomial)
  bootRes[i,-1] <- summary(fit)$coef["cohortcontrol:variantvoc", c(4,1)]
}



additive_prob <- function(start, end1, end2, link="logit") {
  fns <- make.link(link)
  odd_effect_1 <- fns$linkfun(end1)-fns$linkfun(start)
  odd_effect_2 <- fns$linkfun(end2)-fns$linkfun(start)
  overall_odd_effect <- odd_effect_1 + odd_effect_2
  final_odd <- fns$linkfun(start) + overall_odd_effect
  fns$linkinv(final_odd)
}



### Individual responsiveness
p_wuhan <- list(control = 0.8,
               case     = 0.3)
p_voc  <- list("wuhan+" = 0.4,
              "wuhan-"  = 0.05)

patient_frame <- do.call(rbind, case_control)
patient_frame <- mutate(
  patient_frame,
  drug=ifelse(cohort=="control","none",drug),
  cohort=relevel(factor(cohort), "control"))


boot_individ <- data.frame(i=1:1000, coef=0.0)
diseases <- unique(patient_frame$disease)
drugs <- unique(patient_frame$drug[patient_frame$cohort=="case"])
drug_offset <- sapply( unique(patient_frame$drug), function(x) 0.0)
#boot_mat <- matrix(0,4,nrow(boot_individ))

for (i in 1:nrow(boot_individ)) {
  disease_offset <- setNames(runif(length(diseases), -0.05, .05), diseases)
  drug_offset[drugs] <- runif(length(drugs), -0.01, .01)
  patient_frame <- mutate(
    patient_frame,
    wuhan=ifelse(runif(n())<unlist(p_wuhan[as.character(cohort)]) + disease_offset[disease] + drug_offset[drug],"wuhan+","wuhan-"),
    voc = ifelse(runif(n()) < unlist(p_voc[as.character(wuhan)]), "voc+", "voc-"),
    wuhan=relevel(factor(wuhan), "wuhan+"),
    voc=relevel(factor(voc), "voc-")
  )
  dat_agg <- group_by(patient_frame , disease, drug, cohort, wuhan, voc) %>%
    summarise(n=n(), .groups="drop")
#  fit <- glm(n~cohort*wuhan*voc + disease, data=dat_agg, family=poisson)
#  boot_individ$coef[i] <- coef(fit)["cohortcase:wuhanwuhan-:vocvoc+"]
#  next
  fit <- rms::lrm(voc ~ wuhan * cohort + disease, data=patient_frame,x=TRUE, y=TRUE)
#  fit <- robcov(fit, cluster=patient_frame$drug)
#  em <- emmeans(fit, spec=~cohort+wuhan,  tran="logit")
#  boot_mat[,i] <- summary(em, type="response")$response
  boot_individ$coef[i] <- coef(fit)["wuhan=wuhan- * cohort=case"]
  #boot_individ$coef[i] <- as.data.frame(em$contrasts)$estimate
}
quantile(exp(na.omit(boot_individ$coef)),c(0.025,0.975))
quantile(na.omit(boot_individ$coef),c(0.025,0.975))

