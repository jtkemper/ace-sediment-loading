# //////////////////////////////////////////////////////////////////////////////
# //////////////////////////////////////////////////////////////////////////////

# This script does WHAAAAAAAAAAAAAT

## ** Inputs **

## ** Outputs 

# //////////////////////////////////////////////////////////////////////////////
# //////////////////////////////////////////////////////////////////////////////

# Housekeeping ----

## Load packages ----

### Data processing 

require(tidyverse)
require(tsibble)

### File & environment mgmt

require(usethis)
require(here)
require(fs)

### Modeling
require(ranger)

## Set locations ----

### Where are we?

here::i_am("01_data_processing.R")

## Set some global vars ----

input_data_dir <- here("input", "cleaned-data")

output_data_dir <- here("output", "processed-data")

output_fig_dir <- here("output", "plots")

# //////////////////////////////////////////////////////////////////////////////
# //////////////////////////////////////////////////////////////////////////////

# Data download ----

input_sens_files <- dir_ls(input_data_dir)

# ******************************

## Load in sensor data from file ----

#i <- 1

sens <- list()

for(i in 1:length(input_sens_files)) {
  
  # Track progress
  cat(crayon::green("Reading ", input_sens_files[i]))
  
  # Read file
  # Pivot longer and wider to get var names independent of site
  # Create a new column with site 
  # And then add new columns that say which reservoir and site type
  # (upstream or downstream)

  sens[[i]] <- read_csv(input_sens_files[i]) %>%
    rename_all(tolower) %>%
    rename(obs = 1) %>%
    pivot_longer(
      cols = -c(obs, datetime_utc),
      names_to = c("parameter", "site"),
      names_sep = "_",
      values_to = "value"
    ) %>%
    pivot_wider(values_from = "value",
                names_from = "parameter") %>%
    mutate(
      site = str_replace(site, "(up|down)$", "_\\1")) %>%
    separate(
      site,
      into = c("site_reservoir", "site_type"),
      sep = "_",
      remove = FALSE
    )
  
  
}

zz <- read_csv(input_sens_files[1]) %>%
  # mutate(DateTime_UTC = as.POSIXct(DateTime_UTC,
  #                                  tz = "UTC")) %>%
  mutate(day = as_date(DateTime_UTC)) %>%
  filter(day == "2025-06-05")

xx <- sens[[1]] %>%
  # mutate(DateTime_UTC = as.POSIXct(DateTime_UTC,
  #                                  tz = "UTC")) %>%
  mutate(day = as_date(datetime_utc)) %>%
  filter(day == "2025-06-05")


### Bind all the data together in one dataframe 
sens_data_all <- bind_rows(sens) %>%
  dplyr::select(!obs) 

# ******************************

## Load in TSS-Turb coefficient data from file ----

coeff <- read_csv(dir_ls(here("input", 
                              "grab-samples"), 
                         glob = "*.csv"))

coeff <- coeff %>%
  rename_all(tolower) %>%
  mutate(site = tolower(site)) %>%
  mutate(site = str_replace(site, 
                            "(up|down)$", 
                            "_\\1"))

# ******************************

# Data processing & calculation

## Transform turb to TSS at each site ----

transform_turb <- function(sensor_df, coeff_df, smear = 1.1) {
  
  sensor_df <- sens[[2]]
  
  tss_coeff <- coeff %>%
    filter(site == unique(sensor_df$site)) %>%
    .$coefficient
  
  sensor_df <- sensor_df %>%
    mutate(tss = turb*tss_coeff*smear)
  
  return(sensor_df)
  
}

# ******************************

## Convert to hourly ----

#### The transformation to hourly data for turbidity is more complicated
#### Because we do not want to simply take a arithmetic mean of turbidity
#### But instead take a "flow-weighted" average
#### Meaning we have to convert turbidity measures to "flux" using discharge,
#### integrate under the flux curve over each hour,
#### and then divide by the water "load" for that hour 
#### (which we get by integrating under the discharge curve).
#### This gives us a mean turbidity that is more reflective of the central 
#### tendency of turbidity in the real world than the arithmetic mean,
#### which might be overly skewed by high or low values

#### To do this, we must do the following:
#### 1) Multiply discharge and turbidity to get a flux
#### 2) Multiply the flux by the unit of time reflected in each measurement 
#### (importantly, this is one hour if there is one measure per hour,
#### 30 minutes if there are two per hour
#### 20 minutes if there are three per hour
#### and 15 minutes if there are four per hour)
#### 3) Do the same for discharge to get water load
#### 4) Divide turbidity "load" by water "load"
#### 5) And finally fill in gaps that are 24 hours are less
#### With linear interpolation

# wq_and_q_df = sens[[1]]
# constit = "turb"

make_hourly <- function(wq_and_q_df = NULL,
                        constit = "turb") {
  
  
  ### Convert to hourly data
  
  hourly <- wq_and_q_df %>%
    mutate(date = as_date(datetime_utc),
           hour = hour(datetime_utc)) %>%
    dplyr::group_by(date, hour) %>%
    mutate(wq_flux = discharge*get(constit)) %>%
    mutate(obs_per_hour = as.numeric(n())) %>%
    dplyr::ungroup() %>%
    mutate(wq_load = case_when(obs_per_hour == 4 ~ 
                                 wq_flux*60*15,
                               obs_per_hour == 2  ~ 
                                 wq_flux*60*30,
                               obs_per_hour == 1 ~ 
                                 wq_flux*60*60,
                               obs_per_hour == 3 ~ 
                                 wq_flux*60*20)) %>%
    mutate(q_load = case_when(obs_per_hour == 4 ~ 
                                discharge*60*15,
                              obs_per_hour == 2 ~ 
                                discharge*60*30,
                              obs_per_hour == 1 ~ 
                                discharge*60*60,
                              obs_per_hour == 3 ~ 
                                discharge*60*20)) %>%
    dplyr::group_by(date, hour,
                    site, site_reservoir, site_type) %>%
    summarise(datetime_utc = floor_date(datetime_utc[1], 
                                        unit = "hour"),
              !!constit := sum(wq_load)/
                sum(q_load),
              discharge = mean(discharge)
              ) %>%
    dplyr::ungroup()

  
  return(hourly)
  
}

hourly_abiq_down <- make_hourly(sens[[1]],
                                "turb")

  
  
  


# ******************************

# ******************************

# //////////////////////////////////////////////////////////////////////////////

## Fill gaps in WQ record ----

sensor_df = hourly_abiq_down
constit = "turb"
gap_length = 12
gap_info = TRUE

replace_gaps <- function(sensor_df = NULL,
                         constit = "turb",
                         gap_length = 12,
                         gap_info = TRUE){
  
  ### Store site
  
  sens_site <- sensor_df$site[1]
  
  ### Make empty list for returning things
  
  ### Check for duplicate observations at identical timestamps
  ### And extract only one observation per timestamp
  ### if those duplicates exist
  
  dups <- is_duplicated(sensor_df, index = datetime_utc)

  
  if(dups == TRUE){
    
    sensor_df <- sensor_df %>%
      dplyr::group_by(datetime_utc) %>%
      dplyr::slice(1) %>%
      dplyr::ungroup()
    
  }
  
  ### Make into tsibble
  
  sensor_df_ts <- sensor_df %>%
    as_tsibble(index = datetime_utc)
  
  ### Check for literal gaps (as in missing entries, not just NAs)
  ### And then fill those gaps with additional NAs
  ### This will give us a row for every fifteen minute timestamp 
  ### in the period of record
  
  
  if(has_gaps(sensor_df_ts) == TRUE) {
    
    ### Fill gaps
    sensor_df_ts <- sensor_df_ts %>%
      fill_gaps()
    
  }
  
  ### Now, determine the number of NAs
  
  #### Variable naming conventions
  constit_gap <- paste0("days_w_", constit, "_gaps")
  
  #### For turbidity 
  
  #### Find the gaps
  #### And calcuate the number of days with at least twelve hour gaps
  
  wq_nas_per_day <- sensor_df_ts %>%
    as_tibble() %>%
    dplyr::group_by(date) %>%
    summarise(na_per_day = sum(is.na(get(constit)))) 
  
  wq_days_w_gaps <- wq_nas_per_day %>%
    mutate(year = year(date)) %>%
    filter(na_per_day >= gap_length) %>%
    dplyr::group_by(year) %>%
    summarise(!!constit_gap := n()) %>%
    mutate(site = sens_site)
  
  
  #### And for discharge
  q_nas_per_day <- sensor_df_ts %>%
    as_tibble() %>%
    mutate(day = as_date(datetime_utc)) %>%
    dplyr::group_by(day) %>%
    summarise(na_per_day = sum(is.na(discharge))) 
  
  q_days_w_gaps <- q_nas_per_day %>%
    mutate(year = year(day)) %>%
    filter(na_per_day >= gap_length) %>%
    dplyr::group_by(year) %>%
    summarise(days_with_q_gaps = n())
  
  ### Return information about the amount of gaps
  ### (if we so wish) 
  
  if(gap_info == TRUE) {
    
    
    filled_sensor_df[[2]] <- wq_days_w_gaps
    
    filled_sensor_df[[3]] <- q_days_w_gaps
    
    
    
  }
  
  ### And finally fill the gaps
  
  #### Using linear interpolation (gaps < 12 hrs)
  
  ### Now interpolate gaps that are less than 24 hours
  
  sensor_df <- sensor_df %>%
    as_tsibble(index = "datetime_utc") %>%
    fill_gaps(.full = FALSE) %>%
    mutate(!!constit := 
             zoo::na.approx(get(constit), 
                            maxgap = gap_length)) %>%
    as_tibble() %>%
    dplyr::ungroup()

  return(days_with_twelve_turb)

  
}

turb_gaps <- map_dfr(sens, replace_gaps)

# //////////////////////////////////////////////////////////////////////////////

# //////////////////////////////////////////////////////////////////////////////

# Calculate loads ----

sensor_df = hourly_abiq_down
interval = "hourly"
type = "annual"


## Annual ----


calculate_loads <- function(sensor_df = NULL, 
                            interval = "hourly",
                            type = "annual") {
  
  
  ### Determine flux
  flux_df <- sensor_df %>%
    mutate(flux = turb*discharge) %>%
    mutate(hourly_load = flux*60*60)
  
  ### Calculate loads
  if(type == "daily") {
    
    load_df <- flux_df %>%
      dplyr::group_by(date, 
                      site,
                      site_reservoir,
                      site_type) %>%
      dplyr::summarise(load = sum(hourly_load,
                                  na.rm = TRUE)) %>%
      mutate(type = type)
    
  } else if(type == "monthly") {
    
    load_df <- flux_df %>%
      mutate(month = month(date,
                           label = TRUE,
                           abbr = FALSE),
             year = year(date)) %>%
      dplyr::group_by(month,
                      year,
                      site,
                      site_reservoir,
                      site_type) %>%
      dplyr::summarise(load = sum(hourly_load,
                                  na.rm = TRUE)) %>%
      mutate(type = type)
    
  } else if(type == "annual") {
    
    load_df <- flux_df %>%
      mutate(year = year(date)) %>%
      dplyr::group_by(
                      year,
                      site,
                      site_reservoir,
                      site_type) %>%
      dplyr::summarise(load = sum(hourly_load,
                                  na.rm = TRUE)) %>%
      mutate(type = type)
    
    
  } else if(type == "event"){
    
    ### ID events
    
    ### Calc. loads
    
    load_df <- flux_df %>%
      dplyr::group_by(
        event_id,
        site,
        site_reservoir,
        site_type) %>%
      dplyr::summarise(load = sum(hourly_load,
                                  na.rm = TRUE)) %>%
      mutate(type = type)
    
    
    
    
  }
  
  
  else {
    
    stop("Error - type must be one of annual, monthly, daily, or event ")
  }
  
  
  
  
  
  
  
  
}




  
sensor_df_clean <- sensor_df_ts %>%
  as_tibble() %>%
  filter(!is.na(turb)) %>%
  filter(!is.na(discharge)) %>%
  filter(turb > 0) %>%
  mutate(log_turb = log10(turb),
         log_q = log10(discharge))





z <- sensor_df_clean %>%
  mutate(date = as_date(datetime_utc),
         hour = hour(datetime_utc)) %>%
  dplyr::group_by(date, hour) %>%
  summarise(mean_log_q = mean(log_q),
            mean_log_turb = mean(log_turb)) #%>%
  ggplot() +
    geom_point(aes(x = log_q, y = log_turb), color = "tan4") +
    theme_bw()
  
mod <- lm(log_turb ~ log_q, data = sensor_df_clean)
  
summary(mod)

mod2 <- lm(mean_log_turb ~ mean_log_q, data = z)

summary(mod2)
