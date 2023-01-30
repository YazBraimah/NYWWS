#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

freyja_parse <- args[1]
coverageReport <- args[2]
sampling_locations <- args[3]
coveragePointPlotDate <- args[4]
coveragePointCovClass <- args[5]
varFreqBarPlotsAll <- args[6]
varFreqBarPlotsMin20X <- args[7]


### Load packages:
library("tidyverse")
# library(gtsummary)
library("gt")

### create a custom colour pallette:
my_20_cols <- c("#ffb7af",
                "#01c53a",
                "#e02fe6",
                "#91c800",
                "#e972ff",
                "#f2e300",
                "#4b0078",
                "#72ff9a",
                "#f6009d",
                "#01d186",
                "#f9003a",
                "#01f4f0",
                "#e34300",
                "#787dff",
                "#cfff8c",
                "#04004a",
                "#ffee82",
                "#004ca7",
                "#b39100",
                "#0177ba",
                "#c97200",
                "#00a4db",
                "#ca0029",
                "#02c4b6",
                "#e90072",
                "#00773c",
                "#ff6abe",
                "#7c8700",
                "#b9b7ff",
                "#a94700",
                "#7bbcff",
                "#ffad54",
                "#001c24",
                "#d9ffab",
                "#8a004f",
                "#e4ffe1",
                "#1a0011",
                "#ffdda3",
                "#560036",
                "#8be5ff",
                "#630600",
                "#d6d4ff",
                "#00320f",
                "#ff7296",
                "#008975",
                "#ff94ac",
                "#465300",
                "#005669",
                "#553600")


### Load data files:
fre <- read.csv(freyja_parse, header = F, sep = ",")
cov <- read.csv(coverageReport, header = F, sep = "\t") %>%
  rename(sample_id = "V1", site = "V2", coverage = "V3")
man <- read.csv(sampling_locations, header = T, sep = ",")



### Generate coverage table
cov %>%
  separate(sample_id, into = c("sampling_date", "cdc_id"), remove = F, sep = "NY") %>%
  mutate(cdc_id = gsub("^", "NY", cdc_id),
         sampling_date = gsub("^(.{4})(.*)$", "\\1-\\2", sampling_date),
         sampling_date = gsub("^(.{7})(.*)$", "\\1-\\2", sampling_date),
         sampling_date = as.Date(sampling_date)) %>%
  merge(y = man, by.x = "cdc_id", by.y = "cdc_id", all.x = T) %>%
  group_by(cdc_id, sample_id, sampling_date, county, sampling_location, Seq_site) %>%
  summarize(mean_coverage = mean(coverage), n = n(), sd = sd(coverage), se = sd(coverage)/sqrt(n)) -> coverage_table


### produced an abridged version with sample ID and coverage class only.
coverage_table %>%
  filter(!is.na(Seq_site)) %>%
  mutate(status = ifelse(mean_coverage < 50 & mean_coverage > 20, "medium coverage", ifelse(mean_coverage < 20, "low coverage", "high coverage"))) %>%
  ungroup() %>%
  select(sample_id, mean_coverage, sd, se, status) -> coverage_status

### Output sample coverage point plots ordered by date
coverage_table %>%
  filter(!is.na(Seq_site)) %>%
  ggplot(aes(reorder(sample_id, sampling_date), log2(mean_coverage+1), colour = sampling_date)) +
  geom_point() +
  geom_hline(yintercept = log2(20+1), linetype = "dashed", colour = "#918c00") +
  geom_errorbar(aes(ymin = log2((mean_coverage + 1) - se), ymax = log2((mean_coverage +1) + se))) +
  # scale_colour_manual(values = c("#ab40c6",
  #                                "#bc002f",
  #                                "#84abff")) +
  facet_grid(.~Seq_site, scales = "free_x", space = "free_x") +
  theme_minimal(base_size = 16) +
  theme(axis.text.x = element_blank()) +
  labs(y = "average coverage\n(log scale)", x = "Sample", title = "Sample coverage", subtitle = "Sample coverage should be bimodal, with most samples falling in the high coverage range\nand a few samples in the low coverage range", caption = "Dashed line represents 20X coverage\nError bars represent standard error of coverage mean")
ggsave(coveragePointPlotDate, height = 5, width = 10)


### Output sample coverage point plots ordered by coverage
coverage_table %>%
  filter(!is.na(Seq_site)) %>%
  mutate(status = ifelse(mean_coverage < 50 & mean_coverage > 20, "medium coverage", ifelse(mean_coverage < 20, "low coverage", "high coverage"))) %>%
  ggplot(aes(reorder(sample_id, mean_coverage), log2(mean_coverage+1), colour = status)) +
  geom_point() +
  geom_hline(yintercept = log2(20+1), linetype = "dashed", colour = "#918c00") +
  geom_errorbar(aes(ymin = log2((mean_coverage + 1) - se), ymax = log2((mean_coverage +1) + se))) +
  scale_colour_manual(values = c("#ab40c6",
                                 "#bc002f",
                                 "#84abff")) +
  facet_grid(.~Seq_site, scales = "free_x", space = "free_x") +
  theme_bw(base_size = 16) +
  theme(axis.text.x = element_blank()) +
  labs(y = "average coverage\n(log scale)", x = "Sample", title = "Sample coverage", caption = "Dashed line represents 20X coverage\nError bars represent standard error of coverage mean")
ggsave(coveragePointCovClass, height = 5, width = 10)

### flag samples with coverage less than 20X
coverage_table %>%
  filter(mean_coverage < 20) %>%
  pull(sample_id) -> low_coverage_samples

### Generate variant frequency table:
fre %>%
  separate(V1, into = c("sampling_date", "cdc_id"), remove = F, sep = "NY") %>%
  mutate(cdc_id = gsub("^", "NY", cdc_id),
         sampling_date = gsub("^(.{4})(.*)$", "\\1-\\2", sampling_date),
         sampling_date = gsub("^(.{7})(.*)$", "\\1-\\2", sampling_date),
         sampling_date = as.Date(sampling_date)) %>%
  merge(y = man, by.x = "cdc_id", by.y = "cdc_id", all.x = T) %>%
  merge(y = coverage_status, by.x = "V1", by.y = "sample_id", all.x = T) %>%
  rename(sample_id = "V1", variant = "V2", frequency = "V3") -> freq_table

### Make a list of counties:
counties <- freq_table %>% filter(!is.na(county)) %>% pull(county) %>% unique

### Make variant frequency plotting function
freq_plot <- function(table, cty){
  table %>%
    filter(county == cty) %>%
    mutate(variant = ifelse(frequency < 0.1, "other", variant)) %>%
    ggplot(aes(sampling_date, frequency, fill = variant)) +
    geom_col(position = "fill") +
    facet_wrap(~sampling_location) +
    scale_fill_manual(values = my_20_cols) +
    # geom_text(aes(label = variant)) +
    labs(title = paste(cty, "County", sep = " "),
         subtitle = paste("Sequencing site: ", table$Seq_site, sep = ""),
         y = "Frequency",
         x = "Sample date",
         caption = "variants under 5% are designated \"other\"") +
    theme_bw(base_size = 16) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  # ggsave(filename = paste("figures/", cty, "_County_variant_freq_barplots.pdf", sep = ""), width = 15, height = 5)
}

### Make variant frequency table function
freq_table_output <- function(cty){
  freq_table %>%
    filter(frequency >= 0.01 & county == cty) %>%
    select(sampling_location, Seq_site, sampling_date, variant, frequency) %>%
    rename("Sampling location" = "sampling_location", "Sequencing site" = "Seq_site", "Sample date" = "sampling_date", Variant = "variant", Frequency = "frequency") %>%
    gt() %>%
    tab_header(
      title = paste(cty, "County", sep = " "),
      subtitle = "Wastewater SARS-CoV2 variants (only variants with ≥1% frequency are shown)") %>%
    gtsave(paste(cty, "_County_variant_table.pdf", sep = ""))
}


### Output PDFs of variant frequncies by county
lapply(counties, freq_table_output)


### Output a single PDF of all variant barplots by county and sampling location
pdf(varFreqBarPlotsAll, width = 15, height = 5)
lapply(counties, freq_plot, table = freq_table)
dev.off()


### Output a single PDF of variant barplots by county and sampling location, but only those above 20X coverage
pdf(varFreqBarPlotsMin20X, width = 15, height = 5)
lapply(counties, freq_plot, table = filter(freq_table, status != "low_coverage"))
dev.off()
