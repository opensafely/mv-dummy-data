library('tidyverse')

censor <- function(date, censor_date, na.censor=TRUE){
  if (na.censor)
    if_else(date>censor_date, as.Date(NA_character_, origin = "1970-01-01"), as.Date(date, origin = "1970-01-01"))
  else
    if_else(date>censor_date, as.Date(censor_date, origin = "1970-01-01"), as.Date(date, origin = "1970-01-01"))
}

missing = function(x, p){
  n <- length(x)
  stopifnot("p must be length 1 or length of x"=length(p) %in% c(1, n))
  NA_type_ <- NA
  mode(NA_type_) <- typeof(x)
  if_else(runif(n)>p, x, NA_type_)
}


pop_size <- 1000
index_date <- as.Date("2020-03-01")
end_date <- as.Date("2023-01-01")

simulated_data <- 
  tibble(
  patient_id=seq_len(pop_size)
  ) %>%
  # simulate raw variables
  mutate(
    age = rnorm(n=n(), mean=60, sd=15),
    sex = rcat(n=n(), levels = c("F", "M"), p = c(0.51, 0.49)),
    diabetes = purrr::rbernoulli(n=n(), p = plogis(-1 + age*0.002 + I(sex=='F')*-0.2)),
    hosp_admission_count = rpois(n=n(), lambda = exp(-2.5 + age*0.03 + I(sex=='F')*-0.2 +diabetes*1)),
    death_date = {
        rexp(
          n=n(), 
          rate = exp(-5 + age*0.01 + I(age^2)*0.0001 + diabetes*1.5 + hosp_admission_count*1)/365
        ) + index_date
      }
  ) %>%
  # add missing values
  mutate(
    age = missing(age, 0.05),
    sex = missing(sex, 0.01),
    death_date = censor(death_date, end_date, na.censor=TRUE),
  ) %>%
  # opensafelify
  mutate(
    age = floor(age),
  )

write_csv(simulated_data, here::here("simulated_data.csv"))
