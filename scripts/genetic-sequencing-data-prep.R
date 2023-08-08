library(tidyverse)

sewershed_metadata <- snakemake@input[["sewershed"]]
voc_metadata <- snakemake@input[["variants_of_concern"]]
lineages_metadata <- snakemake@input[["lineage_map"]]
freyja_output <- snakemake@input[["freyja"]]
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

# merge lineage map to get hex codes
var.data_summary <- left_join(var.data_summary, lineage.map_brief, by = c("lineage"))

# add character string of detected voc's
# isolate the recent data week range
max_week <- max(var.data_summary$year_week)

# 40 day window for most recent data, otherwise grey for no data
# what about historical spatial data? slider will need data at regular intervals, will need to maybe pick weekly or monthly or two week rolling windows

date_window <- floor_date(as.Date(max_week), unit='week')-40

# filter for the current date range
gs_data_recent <- var.data_summary |>
  # remove any variants that do not have a max prev of at least 5%
  filter(percent_sewershed >= 0.05) %>%
  filter(year_week >= date_window) |>
  group_by(cdc_id) |>
  filter(year_week == max(year_week, na.rm = TRUE)) %>%
  ungroup()

# LIST THE VOC'S IN THE SEWERSHED
var_sewer.voc <- gs_data_recent %>%
  filter(!is.na(lineage)) %>%
  group_by(sw_id) %>%
  filter(monitored != "Not monitored") %>%
  filter(!duplicated(lineage)) %>%
  ungroup()

# paste together the variants if they have multiple
var_sewer.voc <- aggregate(lineage ~ sw_id, unique(var_sewer.voc), paste, collapse = ", ")
colnames(var_sewer.voc) <- c("sw_id", "vocs_detected")

# merge to sewer list
var.data_summary <- left_join(var.data_summary, var_sewer.voc, by = c("sw_id"))

# save to RDS format
saveRDS(var.data_summary, file = output_file)
