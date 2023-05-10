
library("tidyverse")

freyja_parse <- snakemake@input[["parse"]]
coverageReport <- snakemake@input[["cov"]]
sampling_locations <- snakemake@input[["samLoc"]]
sample_Info <- snakemake@input[["samInfo"]]
Lineages_Info <- snakemake@input[["linInfo"]]
coveragePointPlotDate <- snakemake@output[["coveragePointPlot"]]
varFreqBarPlotsMin20X <- snakemake@output[["varFreqBarPlotsMin20X"]]
comprehensive_results_table <- snakemake@output[["compResultsTable"]]


### Load data files:
fre <- read.csv(freyja_parse, header = T, sep = ",")
cov <- read.csv(coverageReport, header = T, sep = "\t") %>%
  rename(site = "X.Position..bp.", coverage = "Coverage")
man <- read.csv(sampling_locations, header = T, sep = ",")
sampInf <- read.csv(sample_Info, header = T, sep = "\t")
LinInf <- read.csv(Lineages_Info, header = T, sep = ",")



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

### flag samples with coverage less than 20X at 50% of genome
sampInf %>%
  filter(enough_coverage == "no") %>%
  pull(sample_id) -> low_coverage_samples

### Generate variant frequency table:
fre %>%
  separate(sample_id, into = c("sampling_date", "cdc_id"), remove = F, sep = "NY") %>%
  mutate(cdc_id = gsub("^", "NY", cdc_id),
         sampling_date = gsub("^(.{4})(.*)$", "\\1-\\2", sampling_date),
         sampling_date = gsub("^(.{7})(.*)$", "\\1-\\2", sampling_date),
         sampling_date = as.Date(sampling_date)) %>%
  merge(y = man, by.x = "cdc_id", by.y = "cdc_id", all.x = T) %>%
  merge(y = coverage_table, by.x = "sample_id", by.y = "sample_id", all.x = T) %>%
  merge(y = sampInf, by.x = "sample_id", by.y = "sample_id", all.x = T) %>%
  merge(y = LinInf, by.x = "variant", by.y = "lineage", all.x = T) -> freq_table

### Make a list of counties:
counties <- freq_table %>% filter(!is.na(`county.x`)) %>% pull(`county.x`) %>% unique

### Set the CDC HEX color code for lineage callouts
col <- as.character(freq_table$hex_code)
names(col) <- as.character(freq_table$callout_group)

### Make variant frequency plotting function
freq_plot <- function(table, cty){
  seq_site <- unique(filter(table, `county.x` == cty)$`Seq_site.x`)
  table %>%
    filter(`county.x` == cty) %>%
    select(`sampling_date.x`, variant_pct, `county.x`, `sampling_location.x`, `Seq_site.x`, callout_group, hex_code) %>%
    group_by(`sampling_date.x`, `county.x`, `sampling_location.x`, `Seq_site.x`, callout_group, hex_code) %>%
    summarize(callout_freq = sum(variant_pct)) %>%
    ggplot(aes(`sampling_date.x`, callout_freq, fill = callout_group)) +
    geom_col(position = "fill")  +
    facet_wrap(~`sampling_location.x`) +
    scale_fill_manual(values = col) +
    labs(title = paste(cty, "County", sep = " "),
       subtitle = paste("Sequencing site: ", seq_site, sep = ""),
       y = "Frequency",
       x = "Sample date",
       # caption = "variants under 5% are designated \"other\"",
       ) +
    theme_bw(base_size = 16) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

### Output a single PDF of variant barplots by county and sampling location, but only those above 20X coverage
pdf(varFreqBarPlotsMin20X, width = 15, height = 5)
lapply(counties, freq_plot, table = filter(freq_table, enough_coverage == "yes"))
dev.off()


## output the freq_table results to table:
write.table(freq_table, comprehensive_results_table, quote = F, sep = "\t", row.names = F)
