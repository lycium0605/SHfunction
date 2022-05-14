# Social index -------------------------------------------------
# from https://gist.github.com/camposfa/50e138a48c89dcd31f5a41d5177d66d4

get_centrality <- function(my_focal, my_grp, my_subset) {

  my_network <- my_subset %>%
    dplyr::ungroup() %>%
    dplyr::rename(from = sname, to = partner) %>%
    dplyr::mutate(grp_fct = factor(grp),
           from = paste(from, grp, sep = "_"),
           to = paste(to, grp, sep = "_")) %>%
    dplyr::filter(res_i_adj > -10)

  grps <- dplyr::bind_rows(
    dplyr::select(my_network, name = from, grp),
    dplyr::select(my_network, name = to, grp)) %>%
    dplyr::distinct(name, grp)

  sxs <- dplyr::bind_rows(dplyr::select(my_network, name = from, sex = sname_sex),
                          dplyr::select(my_network, name = to, sex = partner_sex)) %>%
    dplyr::distinct(name, sex)

  f_graph <- tidygraph::as_tbl_graph(my_network, directed = FALSE, my_network) %>%
    tidygraph::activate(nodes) %>%
    dplyr::left_join(grps, by = "name") %>%
    dplyr::left_join(sxs, by = "name") %>%
    dplyr::mutate(sname = str_sub(name, 1, 3))

  # Calculate centrality metrics
  res <- f_graph %>%
    tidygraph::morph(to_split, grp) %>%
    dplyr::mutate(eigen_wt = centrality_eigen(weights = res_i_adj)
           # eigen = centrality_eigen(),
           # betweenness = centrality_betweenness(normalized = TRUE),
           # degree = centrality_degree(normalized = TRUE)
    ) %>%
    tidygraph::unmorph()

  return(res)
}

#' merge_sr
#'
#' @description Merging sci, agi, dsi_pop and dsi_pop_summary, calculate eigen_wt and sumbond
#' @param sci SCI data frame
#' @param agi agi data frame
#' @param dsi_pop dsi_pop data frame
#' @param dsi_pop_summary dsi_pop_summary data frame
#'
#' @return a merged social info data frame
#' @export
#'
merge_sr<-function(sci,agi,dsi_pop,dsi_pop_summary){
  # calculate-eigen_wt
  dsi<-dsi_pop %>%
    dplyr::mutate(centrality = suppressMessages(purrr::pmap(list(sname, grp, subset), get_centrality)))%>%
    dplyr::select(-di, -subset) %>%
    dplyr::mutate(network = purrr::map(centrality, as_tibble)) %>%
    dplyr::select(-centrality) %>%
    tidyr::unnest() %>%
    dplyr::mutate(name = str_sub(name, 1, 3)) %>%
    dplyr::filter(sname == name & grp == grp1) %>%
    dplyr::select(sname, grp, start, end, sex, age_class, eigen_wt)

  # Merge data
    social_info<-sci%>%
      dplyr::select(-subset)%>%
      dplyr::left_join(agi)%>%
      dplyr::left_join(dsi_pop_summary)%>%
      dplyr::left_join(dsi)%>%
      dplyr::select(-subset)%>%
      dplyr::mutate(SumBond_F=StronglyBonded_F+VeryStronglyBonded_F+WeaklyBonded_F)%>%
      dplyr::mutate(SumBond_M=StronglyBonded_M+VeryStronglyBonded_M+WeaklyBonded_M)

    return(social_info)
}

# Make target date df custom ----------------------------------------------
#' Create a data frame with year-long intervals given start and end date, and backwards.
#'
#' @param target_df A data frame that includes the columns sname, sex, grp, and date
#' @param babase A DBI connection to the babase database
#' @param members_l A subset of members table produced by the function 'subset_members'
#' @param window_length Length in years of the time window for the social index
#' @param .by_grp Logical indicating whether to separate by group. Default is TRUE
#' @param .adults_only Logical indicating whether to include adults only. Default is TRUE
#' @param window_shift Length in years of the time shift for the social index. i.e. 1: 2-1 years prior
#' @param .early_life Logical indicating whether to calculate early life: birth - rnk/mature date, default is FALSE
#' @param .cumulative_adult_life Logical indicating whether to calculate adult life: rnk/mature date - darting date, default is FALSE
#'
#' @return A tibble with one row per animal (and optionally, per group) and target date, with contextual data

make_target_date_custom_df <- function(target_df, babase, members_l, window_shift=0,window_length = 1, .by_grp = TRUE,
                                       .adults_only = TRUE, .early_life=FALSE,.cumulative_adult_life=F) {

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

  target_df <- target_df %>%
    dplyr::left_join(biograph_l, by = c("sname", "sex")) %>%
    dplyr::left_join(dplyr::select(md_females, sname, matured), by = "sname") %>%
    dplyr::left_join(dplyr::select(rd_males, sname, ranked), by = "sname") %>%
    dplyr::select(sname, obs_date = date, sex, birth, statdate, matured, ranked)

  if (.early_life){
    target_df <- target_df %>%
      dplyr::mutate(first_start_date = dplyr::case_when(
        sex == "F" ~ birth,
        sex == "M" ~ birth
      )) %>%
      tidyr::drop_na(first_start_date) %>%
      dplyr::mutate(start = first_start_date)%>%
      dplyr::mutate(end = dplyr::case_when(
        sex == "F" ~ matured,
        sex == "M" ~ ranked
      )) %>%
      tidyr::drop_na(end)%>%
      select(sname, sex, birth, obs_date, first_start_date, statdate, start, end)

  } else if(.cumulative_adult_life){
    target_df <- target_df %>%
      dplyr::mutate(first_start_date = dplyr::case_when(
        sex == "F" ~ matured,
        sex == "M" ~ ranked
      )) %>%
      tidyr::drop_na(first_start_date) %>%
      dplyr::mutate(start = first_start_date)%>%
      dplyr::mutate(end = obs_date) %>%
      tidyr::drop_na(end)%>%
      select(sname, sex, birth, obs_date, first_start_date, statdate, start, end)

  }else {
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

    target_df <- target_df %>%
      dplyr::mutate(start = dplyr::case_when(
        first_start_date >= obs_date - lubridate::years(window_length+window_shift) + lubridate::days(1) ~ first_start_date,
        TRUE ~ obs_date - lubridate::years(window_length+window_shift) + lubridate::days(1)))

    target_df <- target_df %>%
      dplyr::mutate(end = obs_date - lubridate::years(window_shift)) %>%
      select(sname, sex, birth, obs_date, first_start_date, statdate, start, end)
  }


  target_df <- target_df %>%
    dplyr::filter(start <= end) %>%
    arrange(sname, obs_date)

  # .by_grp <- TRUE

  if (.by_grp) {
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
  } else {
    ## Check how many days the focal was present in ANY group in a focal year
    # temp <- target_df %>%
    #   dplyr::inner_join(dplyr::select(members_l, sname, date), by = c("sname")) %>%
    #   dplyr::filter(date >= start & date <= end) %>%
    #   dplyr::group_by(sname, start, end) %>%
    #   dplyr::summarise(days_present = n()) %>%
    #   dplyr::arrange(sname, start, end)
    #
    # target_df <- temp %>%
    #   dplyr::inner_join(target_df, by = c("sname", "start", "end")) %>%
    #   dplyr::arrange(sname, start, end)

    stop("Not Yet Completed.")
  }

  # Calculate date variables
  target_df <- target_df %>%
    dplyr::mutate(midpoint = start + floor((end - start) / 2),
                  age_start_yrs = as.numeric(start - birth) / 365.25,
                  age_class = floor(plyr::round_any(age_start_yrs, 0.005)) + 1)

  target_df <- dplyr::ungroup(target_df) %>%
    distinct()

  return(target_df)
}
