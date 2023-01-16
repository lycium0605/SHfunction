#' Create a data frame with year-long intervals prior to specific target dates
#'
#' @param target_df A data frame that includes the columns sname, sex, grp, and date
#' @param babase A DBI connection to the babase database
#' @param members_l A subset of members table produced by the function 'subset_members'
#' @param window_length Length in years of the time window for the social index
#' @param window_length_day Length in days in addition to years of the time window
#' @param .by_grp Logical indicating whether to separate by group. Default is TRUE
#' @param .adults_only Logical indicating whether to include adults only. Default is TRUE
#' @param series logical indicating whether to create a series instead of a single value
#' @param direction "backward" or "forward", direction to create window, default to backward from obs_date
#' @param overlap in term of days, default to zero
#' @param min_pres_day minimum number of days in a group, default to 0
#'
#' @return A tibble with one row per animal (and optionally, per group) and target date, with contextual data
#' @export
#'
#' @examples
#'

# target_df<-ge_meta
# .adults_only<-T

make_target_date_df_sh <- function(target_df, babase, members_l,
                                window_length = 1, window_length_day = 0,
                                min_pres_day = 0,
                                .by_grp = TRUE,
                                .adults_only = TRUE,
                                series=F,
                                overlap=0,
                                direction="backward") {

  if (class(babase) != "PostgreSQLConnection") {
    stop("Invalid connection to babase.")
  }

  # Return an empty tibble if the subset is empty
  if (is.null(target_df) |
      !all(c("sname", "sex", "date") %in% names(target_df))) {
    stop("Problem with input data. Target data frame must include rows 'sname', 'sex', and 'date'.")
  }

  # babase-tables -----------------------------------------------------------

  message("Creating connections to babase tables...")

  # Database connections
  biograph <- dplyr::tbl(babase, "biograph")
  maturedates <- dplyr::tbl(babase, "maturedates")
  rankdates <- dplyr::tbl(babase, "rankdates")

  # Local
  biograph_l <- dplyr::collect(biograph)

  md_females <- maturedates %>%
    dplyr::semi_join(dplyr::filter(biograph, sex == "F"), by = "sname") %>%
    collect()

  rd_males <- rankdates %>%
    dplyr::semi_join(dplyr::filter(biograph, sex == "M"), by = "sname") %>%
    collect()

  # Find last date
  last_date <- max(members_l$date)

  message("Creating target-date data set...")

  # target_df_0<-ge_meta
  target_df <- target_df %>% #target_df_0 %>%
    dplyr::left_join(biograph_l, by = c("sname", "sex")) %>%
    dplyr::left_join(dplyr::select(md_females, sname, matured), by = "sname") %>%
    dplyr::left_join(dplyr::select(rd_males, sname, ranked), by = "sname") %>%
    dplyr::select(sname, obs_date = date, sex, birth, statdate, matured, ranked)

  if (.adults_only) {
    target_df <- target_df %>%
      dplyr::mutate(first_start_date = dplyr::case_when(
        sex == "F" ~ matured,
        sex == "M" ~ ranked
      )) %>%
      tidyr::drop_na(first_start_date) %>%
      dplyr::select(sname, obs_date, sex, birth, first_start_date, statdate, -ranked, -matured)
  } else {
    target_df <- target_df %>%
      dplyr::mutate(first_start_date = dplyr::case_when(
        sex == "F" ~ birth,
        sex == "M" ~ birth
      )) %>%
      tidyr::drop_na(first_start_date) %>%
      dplyr::select(sname, obs_date, sex, birth, first_start_date, statdate, -ranked, -matured)
  }

  # target_df <- target_df %>%
  #   dplyr::mutate(first_start_date = dplyr::case_when(
  #     sex == "F" ~ matured,
  #     sex == "M" ~ ranked
  #   )) %>%
  #   dplyr::select(sname, obs_date, sex, birth, first_start_date, statdate, -ranked, -matured)
  # window_length_day<-6
  # window_length<-1

# Backward calculation ----------------------------------------------------

  if(direction=="backward"){
    target_df_0<-target_df
## Create the first time window --------------------------------------------
    target_df <- target_df_0 %>%
      dplyr::mutate(start = dplyr::case_when(
        first_start_date >= obs_date - lubridate::years(window_length) + days(1) ~ first_start_date,
        TRUE ~ obs_date - lubridate::years(window_length) - days(window_length_day) + days(1)))

    target_df <- target_df %>%
      dplyr::mutate(end = obs_date) %>%
      select(sname, sex, birth, obs_date, first_start_date, statdate, start, end)
    target_df <- target_df %>%
      dplyr::filter(start < end) %>%
      arrange(sname, obs_date)

    if(series){
      target_df<-distinct(target_df)
      tmp<-target_df%>%
        dplyr::filter(start != first_start_date)%>%
        dplyr::mutate(end=end - lubridate::years(window_length) - lubridate::days(window_length_day) + lubridate::days(overlap))%>%
        dplyr::mutate(start = dplyr::case_when(
          first_start_date >= end - lubridate::years(window_length) - lubridate::days(window_length_day) + lubridate::days(1) ~ first_start_date,
          TRUE ~ end - lubridate::years(window_length) - lubridate::days(window_length_day) + lubridate::days(1)))%>%
        dplyr::filter(start < end) %>%
        arrange(sname, obs_date)
      tmp2<-tmp
      while(nrow(tmp)!=0){
        tmp<-tmp%>%
          dplyr::filter(start != first_start_date)%>%
          dplyr::mutate(end=end - lubridate::years(window_length) - lubridate::days(window_length_day) + lubridate::days(overlap))%>%
          dplyr::mutate(start = dplyr::case_when(
            first_start_date >= end - lubridate::years(window_length) + lubridate::days(1) ~ first_start_date,
            TRUE ~ end - lubridate::years(window_length) - lubridate::days(window_length_day) + days(1)))%>%
          dplyr::filter(start < end) %>%
          arrange(sname, obs_date)
        tmp2<-rbind(tmp2,tmp)
        #message(nrow(tmp))
      }
      target_df<-rbind(target_df,tmp2)%>%
        dplyr::filter(start < end) %>%
        arrange(sname, obs_date)
    }
# Forward calculation ----------------------------------------------------
  }else if(direction=="forward"){
    target_df_0<-target_df
    ## Create the first time window --------------------------------------------
    target_df <- target_df_0 %>%
      dplyr::mutate(end = dplyr::case_when(
        obs_date <=  birth + lubridate::years(window_length) + lubridate::days(window_length_day) - days(1) ~ obs_date,
        TRUE ~ birth + lubridate::years(window_length) + days(window_length_day) - days(1)))

    target_df <- target_df %>%
      dplyr::mutate(start = birth) %>%
      select(sname, sex, birth, obs_date, first_start_date, statdate, start, end)
    target_df <- target_df %>%
      dplyr::filter(start < end) %>%
      arrange(sname, obs_date)

      target_df<-distinct(target_df)
      tmp<-target_df%>%
        dplyr::filter(end!=obs_date)%>%
        dplyr::mutate(start=start + lubridate::years(window_length) + days(window_length_day)-days(overlap))%>%
        dplyr::mutate(end = dplyr::case_when(
          obs_date <=  start + lubridate::years(window_length) + days(window_length_day) - days(1) ~ obs_date,
          TRUE ~ start + lubridate::years(window_length) + days(window_length_day) - days(1)))%>%
        dplyr::filter(start < end) %>%
        arrange(sname, obs_date)
      tmp2<-tmp
      while(nrow(tmp)!=0){
        tmp<-tmp%>%
          dplyr::filter(end!=obs_date)%>%
          dplyr::mutate(start=start + lubridate::years(window_length) + days(window_length_day)-days(overlap))%>%
          dplyr::mutate(end = dplyr::case_when(
            obs_date <=  start + lubridate::years(window_length) + days(window_length_day) - days(1) ~ obs_date,
            TRUE ~ start + lubridate::years(window_length) + days(window_length_day) - days(1)))%>%
          dplyr::filter(start < end) %>%
          arrange(sname, obs_date)
        tmp2<-rbind(tmp2,tmp)
      }
      target_df<-rbind(target_df,tmp2)%>%
        dplyr::filter(start < end) %>%
        dplyr::filter(end > first_start_date) %>%
        dplyr::mutate(start = dplyr::case_when(
          first_start_date >= start ~ first_start_date,
          TRUE ~ start))%>%
        dplyr::filter(start < end) %>%
        arrange(sname, obs_date)

    }else(stop("direction has to be specified as forward or backward"))


  ## Check in which groups the individual was present in the focal year
  ## and create one row per focal year per group
  temp <- target_df %>%
    dplyr::left_join(dplyr::select(members_l, sname, date, grp), by = c("sname")) %>%
    dplyr::filter(date >= start & date <= end) %>%
    dplyr::distinct(sname, start, end, grp)

  zdata <- target_df %>%
    dplyr::inner_join(temp, by = c("sname", "start", "end")) %>%
    tibble::rownames_to_column()

  ## And check how many days the focal was present in the group in a focal year
  zdata <- zdata %>%
    dplyr::inner_join(dplyr::select(members_l, sname, grp, date), by = c("sname", "grp")) %>%
    dplyr::filter(date >= start & date <= end) %>%
    dplyr::group_by(sname, grp, start, end, rowname) %>%
    dplyr::summarise(days_present = n()) %>%
    dplyr::arrange(sname, grp, start, end)

  target_df <- zdata %>%
    dplyr::inner_join(target_df, by = c("sname", "start", "end")) %>%
    dplyr::arrange(sname, grp, start, end) %>%
    dplyr::select(-rowname)


  # Calculate date variables
  target_df <- target_df %>%
    dplyr::mutate(midpoint = start + floor((end - start) / 2),
                  age_start_yrs = as.numeric(start - birth) / 365.25,
                  age_class = floor(plyr::round_any(age_start_yrs, 0.005)) + 1)

  target_df <- dplyr::ungroup(target_df) %>%
    distinct()%>%
    dplyr::filter(days_present>=min_pres_day)

  return(target_df)
}
