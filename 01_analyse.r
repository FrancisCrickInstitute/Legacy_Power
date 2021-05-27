library(simr)
library(lme4)
library(tidyverse)

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
cc_prop  = 1 # multiply corresponding case numbers by this for controls, then 
cc_const = 0 # add this amount
n_per_drug_per_cohort$control <- lapply(
  n_per_drug_per_cohort$case, function(x) sum(x)*cc_prop + cc_const
)

## Effect sizes

# for controls with all doses, response in the vaccinated-for variant
# odds scale. So e.g. 19 = 95% chance success, 5% failure
# final_vaccine_odds <- percent / (100-percent)
final_vaccine_odds <- 19

# Relative for the other variants
variant_multiplier <- list(
  vaccine_target = 1,
  voc = 0.9)

# Relative for other doses 
dose_multiplier <- list(
  first = 0.95,
  second = 1
)

# Overall case effect
case_effect <- list(
  control = 1,
  case    = 0.8
)

# effect of interest
eoi <- list(mult = 0.25,
           target = quote(variant=="voc" & cohort=="case"))

# sd of effect across diseases
disease_sd <- list(
  control = .1,
  case = 0.1
)

# sd of treatment effect within disease
drug_sd <- list(
  control= 0.1,
  case = 0.1
)



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
case_control_frame <- do.call(rbind, case_control)


  
measure_frame <- expand.grid(
  variant = names(variant_multiplier),
  dose = names(dose_multiplier))

ind <- expand.grid(
  i_cc = 1:nrow(case_control_frame),
  i_measure = 1:nrow(measure_frame)
)

case_control_frame <- cbind(
  case_control_frame[ind$i_cc,],
  measure_frame[ind$i_measure,]
  )


## Fill in fixed effects multipliers
case_control_frame <- case_control_frame %>%
  dplyr::mutate(
    odds = final_vaccine_odds *
      unlist(variant_multiplier[variant]) *
      unlist(dose_multiplier[dose])
  )


bootRes <- data.frame(i=1:100,
                     p=0.0,
                     effect=0.0)

for (i in 1:nrow(bootRes)) { 
  cat(i, " ")
  
## Get instance of random effects:
disease_effect <- imap(case_effect, function(mult, cc) {
  is_cc <- case_control_frame$cohort==cc
  diseases <- unique(case_control_frame$disease[is_cc])
  setNames(
    rnorm(length(diseases),
          case_effect[[cc]],disease_sd[[cc]]
          ),
    diseases
  )
}
)

drug_effect <- imap(case_effect, function(mult_cc, cc) {
  diseases <- disease_effect[[cc]]
  imap(diseases, function(mult, disease) {
    x <- rnorm(n=length(n_per_drug_per_cohort[[cc]][[disease]]),
          mult,
          drug_sd[[cc]])
    setNames(x, paste(cc, disease, seq(along=n_per_drug_per_cohort[[cc]][[disease]]), sep="_"))
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
                            


fit <- glmer(response ~ cohort * variant + (1|disease/drug), data=dat, family=binomial)
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
  
                  
