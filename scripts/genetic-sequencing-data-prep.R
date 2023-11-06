library(tidyverse)

sewershed_metadata <- snakemake@input[["sewershed"]]
voc_metadata <- snakemake@input[["variants_of_concern"]]
lineages_metadata <- snakemake@input[["lineage_map"]]
freyja_output <- snakemake@input[["freyja"]]
concentration_data <- snakemake@input[["concentration"]]
output_file <- snakemake@output[["rds_data"]]

import.data = function(file.name) {
    readr::read_delim(file.name, delim = ',', show_col_types = FALSE)
}

# sewer shed metadata
meta = import.data(sewershed_metadata)

# variants of concern
voc = import.data(voc_metadata) |> 
  janitor::clean_names() |> 
  mutate(label = na_if(label, ''))|>
  filter(!is.na(label))|>
  mutate(variant = gsub('(.*)\\(.*\\)','\\1', variant))|>
  filter(!duplicated(variant))

# lineage map
lineage.map = import.data(lineages_metadata) |> 
  rename(hex = hex_code,
         variant = lineage,
         lineage = callout_group)

# genetic sequencing variant data
var.data = import.data(freyja_output) |> 
  separate(sample_id, into = c('date', 'cdc_id'), sep = 8) |> 
  mutate(date = ymd(date)) |>
  left_join(meta, by = 'cdc_id') |> 
  filter(!is.na(variant)) |>
  filter(variant != '') |> 
  left_join(lineage.map, by = 'variant') |> 
  filter(!is.na(lineage))

# add monitored or not
var.data$monitored = ifelse(var.data$lineage %in% voc$variant,
                            'Monitored',
                            'Not monitored')

# VOC'S DETECTED

# aggregate by week for every site, every week

# object with each week for the loop
weeks <- unique(floor_date(var.data$date, unit = 'week'))

# empty list to store df
voc_data_list <- list()

for(i in unique(weeks)){
  # create date window for summarizing data
  date_window <- floor_date(as.Date(i, origin = "1970-01-01"), unit='week')-40
  
  # filter for the current date range
  gs_data_recent <- var.data |>
    # remove any variants that do not have a max prev of at least 5%
    filter(variant_pct >= 0.05) %>%
    mutate(year_week = floor_date(date, 'weeks')) %>%
    filter(year_week  >= date_window) |>
    group_by(cdc_id) |>
    filter(year_week == max(year_week, na.rm = TRUE)) %>%
    ungroup() 
  
  # LIST THE VOC'S IN THE SEWERSHED
  var_sewer.voc <- gs_data_recent %>%
    group_by(sw_id) %>%
    filter(!is.na(variant)) %>%
    filter(monitored == "Monitored") %>%
    filter(!duplicated(variant)) %>%
    ungroup()
  
  # paste together the variants if they have multiple
  var_sewer.voc <- aggregate(variant ~ sw_id, unique(var_sewer.voc), paste, collapse = ", ")
  colnames(var_sewer.voc) <- c("sw_id", "vocs_detected")
  var_sewer.voc$year_week <- as.Date(i, origin = "1970-01-01")
  
  # store in list
  voc_data_list[[i]] <- var_sewer.voc
  
}

# transform into df
var_sewer.voc <- do.call(rbind, voc_data_list)

# ----------------------------------------------------------------------------------------------------------------------------------------------

# SUMMARY FILE CREATION # 

# collapse across variants into parent lineages
var.data = var.data |>
  # calculate total prevalence of each lineage (lineage) per day
  group_by(date, cdc_id, lineage) |>
  mutate(variant_pct_sewershed = sum(variant_pct, na.rm = TRUE)) |>
  ungroup() |>
  distinct(date,
           cdc_id,
           sw_id,
           county,
           region,
           epaid,
           wwtp_name,
           population_served,
           lineage,
           monitored,
           variant_pct_sewershed)

# # if variant pct < 5%, set to other, then resummarize again
var.data$lineage <- ifelse(var.data$variant_pct_sewershed < 0.05, "Other", var.data$lineage)
var.data$monitored <- ifelse(var.data$lineage == "Other", "Not monitored", var.data$monitored)
var.data = var.data |>
  # calculate total prevalence of each lineage (lineage) per day
  group_by(date, cdc_id, lineage) |>
  mutate(variant_pct_sewershed = sum(variant_pct_sewershed, na.rm = TRUE)) |>
  ungroup()

# add the week to the data
var.data = var.data |>
  mutate(year_week = floor_date(date, 'weeks'))

# for each week, we need sewer, county, region prevalence of each variant/lineage detected
var.data = var.data |>
  # average the prevalence across lineages for each week
  group_by(year_week, cdc_id, sw_id, county, lineage, region, epaid,wwtp_name, population_served, monitored) |>
  summarize(mean_variant_pct_sewershed = mean(variant_pct_sewershed,
                                              na.rm = TRUE)) |>
  ungroup()

# sewer shed percentages
var.data_sewershed = var.data |>
  # add the total prevalence at the sampling site
  group_by(year_week, cdc_id) |>
  mutate(total_variant_pct_sewershed = sum(mean_variant_pct_sewershed,
                                           na.rm = TRUE))|>
  mutate(percent_sewershed = mean_variant_pct_sewershed / total_variant_pct_sewershed) |>
  ungroup() |>
  distinct(year_week,
           cdc_id,
           sw_id,
           county,
           region,
           epaid,
           population_served,
           lineage,
           monitored,
           percent_sewershed)

# wwtp population-weighted mean of prevalence (for the upstream sampling wwtps, e.g., Binghamton)
var.data_wwtp = var.data |>
  group_by(year_week, epaid, lineage) |>
  mutate(mean_variant_pct_wwtp = weighted.mean(mean_variant_pct_sewershed,
                                               population_served,
                                               na.rm = TRUE))|>
  ungroup() |>
  distinct(year_week,
           epaid,
           county,
           region,
           wwtp_name,
           # population_served,
           lineage,
           mean_variant_pct_wwtp)

# sum the total prev by wwtp
# and divide the mean_variant_pct_wwtp by the total to get wwtp variant_pct
var.data_wwtp = var.data_wwtp |>
  group_by(year_week, epaid) |>
  mutate(total_variant_pct_wwtp = sum(mean_variant_pct_wwtp,
                                      na.rm = TRUE))|>
  mutate(percent_wwtp = mean_variant_pct_wwtp / total_variant_pct_wwtp) |>
  ungroup() |>
  distinct(year_week,
           epaid,
           county,
           region,
           wwtp_name,
           lineage,
           percent_wwtp)

# pop weighted mean for the county
var.data_county = var.data |>
  # county pop weighted average for each lineage
  group_by(year_week, county, lineage) |>
  mutate(mean_variant_pct_county = weighted.mean(mean_variant_pct_sewershed,
                                                 population_served,
                                                 na.rm = TRUE)) |>
  ungroup() |>
  distinct(year_week,
           county,
           region,
           lineage,
           mean_variant_pct_county)

# sum the total prevalence across the county
# then divide the mean_variant_pct_county by the total county prevalence
# to get the percent prevalence for each variant
var.data_county = var.data_county |>
  # county total prevalence (denominator for county percent)
  group_by(year_week, county) |>
  mutate(total_variant_pct_county = sum(mean_variant_pct_county,
                                        na.rm = TRUE)) |>
  mutate(percent_county = mean_variant_pct_county / total_variant_pct_county) |>
  ungroup() |>
  distinct(year_week,
           county,
           region,
           lineage, 
           percent_county)

# pop weighted mean for the region
var.data_region = var.data |>
  # region pop weighted average
  group_by(year_week, region, lineage) |>
  mutate(mean_variant_pct_region = weighted.mean(mean_variant_pct_sewershed,
                                                 population_served, na.rm = TRUE)) |>
  ungroup() |>
  distinct(year_week,region, lineage, mean_variant_pct_region)

# sum the total prevalence across the region
# then divide the mean_variant_pct_region by the total region prevalence
# to get the percent prevalence for each variant
var.data_region = var.data_region |>
  # region total prevalence (denominator for regional percent)
  group_by(year_week, region) |>
  mutate(total_variant_pct_region = sum(mean_variant_pct_region, na.rm = TRUE))|>
  mutate(percent_region = mean_variant_pct_region / total_variant_pct_region) |>
  ungroup() |>
  distinct(year_week,region, lineage, percent_region)

# combine data sets
var.data_summary = var.data_sewershed |> 
  left_join(var.data_county, by = c('year_week', 'county', 'region', 'lineage')) |> 
  left_join(var.data_region, by = c('year_week', 'region', 'lineage')) |> 
  left_join(var.data_wwtp, by = c('year_week', 'county', 'region', 'epaid', 'lineage'))

# add hex codes

# remove duplicate lineages
lineage.map_brief <- lineage.map %>% select(lineage, hex) %>% distinct()

# add hex codes for the three new categories
na <- as.data.frame(cbind(c("lineage","No sequencing data available"), c("hex", "#3B3B3B")))
colnames(na) <- c("lineage", "hex")
nocovid <- as.data.frame(cbind(c("lineage","No SARS-CoV-2 detected in sample"), c("hex", "white")))
colnames(nocovid) <- c("lineage", "hex")
na <- tail(na,1)
nocovid <- tail(nocovid,1)
lineage.map_brief <- bind_rows(lineage.map_brief, nocovid, na)

# -------------------------------------------------------------------------------------------------

# LINK WITH QUANTIFICATION DATA

# load quantification data to identify the categories of : no sample, no sars-cov-2 detected, and no seq data available
quant.data = import.data(concentration_data) |> 
  separate(sample_id, into = c('date', 'cdc_id'), sep = 8) |> 
  mutate(date = ymd(date))

# Calculating copies of SARS2 RNA #
quant.data$copies <- 3.5
quant.data$copies <- ifelse(!is.na(quant.data$sars2_copies_ml), as.numeric(as.character(quant.data$sars2_copies_ml)), quant.data$copies)
quant.data$copies <- ifelse(quant.data$sars_pos==0, 1, quant.data$copies)
quant.data$copies <- ifelse(is.na(quant.data$copies), as.numeric(as.character(quant.data$sars2_copies_ml)), quant.data$copies)
quant.data$copies <- ifelse(is.na(quant.data$copies), 1, quant.data$copies)
quant.data$copies[quant.data$sars2_copies_ml == 0] <- 1
quant.data$copies <- ifelse(quant.data$copies == 0, 1, quant.data$copies)
quant.data$copies[quant.data$copies < 1] <- 1

# summarize sars data to the weekly level picking highest copies value
quant.data_weekly <- quant.data %>%
  group_by(year_week = floor_date(date, 'weeks'), cdc_id)%>%
  slice(which.max(copies)
        )%>%
  ungroup()%>%
  select(year_week, cdc_id, sw_id, copies)

quant.data_weekly <- left_join(quant.data_weekly, meta, by = c("sw_id", "cdc_id")) %>%
  select(year_week, cdc_id, sw_id, county, region, epaid, population_served, copies)

# merge quant data and var.data
# add max date
max_week <- max(var.data_summary$year_week, na.rm = TRUE)
var.data_summary$max_week <- max(var.data_summary$year_week, na.rm = TRUE)
var.data_summary <- full_join(var.data_summary, quant.data_weekly, by = c("year_week", "cdc_id", "sw_id", "epaid", "county", "population_served", "region"))
var.data_summary$max_week <- ifelse(is.na(var.data_summary$max_week), max_week, var.data_summary$max_week)
var.data_summary$max_week <- as.Date(var.data_summary$max_week, origin = "1970-01-01")

# remove data before seq data started and after the max gen seq date
var.data_summary <- var.data_summary %>%
  filter(year_week >= "2022-12-28") %>%
  filter(year_week <= max_week)

# add no sample, no sars-cov-2 detected, or no seq data available for those sites
var.data_summary$lineage <- ifelse(var.data_summary$copies == 1 & is.na(var.data_summary$lineage), "No SARS-CoV-2 detected in sample", var.data_summary$lineage)
var.data_summary$lineage <- ifelse(is.na(var.data_summary$lineage) & !is.na(var.data_summary$copies) , "No sequencing data available", var.data_summary$lineage)
var.data_summary$lineage <- ifelse(is.na(var.data_summary$copies) & is.na(var.data_summary$lineage), "No sequencing data available", var.data_summary$lineage)
# 
# # set the value of the percentages to 1 for all these categories
# missing_cat <- c("No SARS-CoV-2 detected in sample", "No sequencing data available", "No sample")
# var.data_summary$variant_pct_sewershed <- ifelse(var.data_summary$lineage %in% missing_cat, 1, var.data_summary$variant_pct_sewershed)
# var.data_summary$monitored <- ifelse(var.data_summary$lineage %in% missing_cat, "Not monitored", var.data_summary$monitored)

# merge lineage map to get hex codes
var.data_summary <- left_join(var.data_summary, lineage.map_brief, by = c("lineage"))

# -------------------------------------------------------------------------------------------------

var.data_summary <- left_join(var.data_summary, var_sewer.voc, by = c("sw_id", "year_week")) 

# save to RDS format
saveRDS(var.data_summary, file = output_file)
