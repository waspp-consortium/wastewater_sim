# Reshape linelist --> number shedding and wastewater detections per day

## This function takes this line list + a shedding profile -->
## number shedding per day and effective number shedding per day

#' Compute total and effective shedding per day (with optional isolation filtering)
#'
#' @param df A linelist dataframe containing: onset_date, asymptomatic, lab_confirm_date, isolated_from_date
#' @param shedding_dist Normalised discrete shedding distribution with p = 1 on peak day
#' @param maxX Maximum number of shedding days (default 16)
#' @param only_after_isolation Logical; if TRUE, shedding before isolation is removed (default FALSE)
#'
#' @return A tibble with one row per day containing:
#'   - n_shedding (number of people shedding that day)
#'   - effective_shedding (sum of weighted probabilities)
#'
#' @export
compute_daily_shedding <- function(df,
                                   shedding_dist,
                                   only_after_isolation = FALSE) {

  # prepare input data
  df <- df %>%
    mutate(across(c(onset_date, lab_confirm_date, isolated_from_date), as.Date)) %>%
    mutate(day0 = if_else(asymptomatic == "No", onset_date, lab_confirm_date))

  # expand individuals across their shedding periods
  shedding_df <- df %>%
    crossing(shedding_dist) %>%
    mutate(shedding_date = day0 + days(rel_day))

  # keep only shedding on/after isolation if specified
  if (only_after_isolation) {
    shedding_df <- shedding_df %>%
      filter(shedding_date >= isolated_from_date)
  }

  # summarise by calendar day
  daily_shedding <- shedding_df %>%
    group_by(shedding_date) %>%
    summarise(
      n_shedding = n(),
      effective_shedding = sum(prob, na.rm = TRUE),
      .groups = "drop"
    )
}




