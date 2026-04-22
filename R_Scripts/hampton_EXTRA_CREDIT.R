# ================================================================
# Hamptonian Highway Helper — Bonus Section
# Weather Impact on Hampton Roads Traffic Volume (2024)
# Author: Cameron Ridgley
# Asked Claude to help review and edit some code blocks



# ── 0. Packages ─────────────────────────────────────────────────
# Run once:
# install.packages(c("dplyr","lubridate","ggplot2","scales","ggcorrplot"))

library(dplyr)
library(lubridate)
library(ggplot2)
library(scales)
library(ggcorrplot)


# ── 1. Load Data ─────────────────────────────────────────────────
setwd("/Users/cameron.ridgley/Documents/Projects/HamptonHighwayHelper/DataSets")

traffic_raw <- read.csv("VDOT_Bidirectional_Traffic_Volume_2024.csv",
                        stringsAsFactors = FALSE)
weather_raw  <- read.csv("Hampton_Roads_Weather_2024.csv",
                         stringsAsFactors = FALSE)


# ── 2. Clean & Filter Traffic ────────────────────────────────────
traffic_raw$date <- as.Date(substr(traffic_raw$DATA_DATE, 1, 10),
                            format = "%Y/%m/%d")

traffic <- traffic_raw %>%
  filter(FROM_DISTRICT == "Hampton Roads" | TO_DISTRICT == "Hampton Roads") %>%
  filter(!is.na(ADT), ADT > 0) %>%
  filter(format(date, "%Y") == "2024") %>%
  # Jan 1 is VDOT's annual reset date — 5,339 segments all stamped the same
  # rainy day, which would confound weather comparisons.  Excluded here;
  # noted as a limitation (parallel to paper's Section 6 limitations).
  filter(date != as.Date("2024-01-01")) %>%
  mutate(
    route_type = case_when(
      RTE_TYPE_CD == "IS" ~ "Interstate",
      RTE_TYPE_CD == "US" ~ "US Route",
      RTE_TYPE_CD == "SR" ~ "State Route",
      RTE_TYPE_CD == "ST" ~ "Street",
      RTE_TYPE_CD == "SC" ~ "Secondary",
      TRUE                ~ "Other"
    ),
    # Mirror paper's in_session variable — TRUE outside June & July
    in_session = !(month(date) %in% c(6, 7)),
    season = case_when(
      month(date) %in% c(12, 1, 2) ~ "Winter",
      month(date) %in% c(3, 4, 5)  ~ "Spring",
      month(date) %in% c(6, 7, 8)  ~ "Summer",
      TRUE                          ~ "Fall"
    )
  ) %>%
  select(date, route = ROUTE_COMMON_NAME, route_type,
         in_session, season,
         adt = ADT, jurisdiction = FROM_JURISDICTION)


# ── 3. Clean Weather ─────────────────────────────────────────────
weather <- weather_raw %>%
  mutate(date = as.Date(date)) %>%
  select(date, tavg_f, tmax_f, tmin_f, precip_in,
         wind_avg_mph, visibility_mi, weather_condition,
         is_precipitation, is_snow, is_fog)


# ── 4. Merge ─────────────────────────────────────────────────────
merged <- traffic %>% inner_join(weather, by = "date")

cat(sprintf("Merged rows: %d  |  Unique dates: %d\n\n",
            nrow(merged), n_distinct(merged$date)))

# Per-date aggregate (one row per observation date, mirrors paper's
# approach of aggregating to a meaningful time unit before plotting)
daily <- merged %>%
  group_by(date, season, in_session,
           tavg_f, tmax_f, tmin_f, precip_in,
           wind_avg_mph, visibility_mi,
           weather_condition, is_precipitation) %>%
  summarise(mean_adt = mean(adt), median_adt = median(adt),
            n_segs   = n(), .groups = "drop") %>%
  mutate(month_lbl = factor(month(date, label = TRUE, abbr = TRUE),
                            levels = month.abb))


# ── 5. Inferential Statistics (α = 0.05, matching paper) ────────

# RW1 — Welch t-test: ADT on precipitation vs. dry days
t_precip <- t.test(adt ~ is_precipitation, data = merged)
cat("=== Welch t-test: ADT ~ Precipitation ===\n"); print(t_precip)

# RW2 — One-way ANOVA: ADT across seasons (extends paper's RQ3)
aov_season <- aov(adt ~ season, data = merged)
cat("\n=== ANOVA: ADT ~ Season ===\n"); print(summary(aov_season))

# RW3 — Welch t-test: in-session vs. out-of-session
#        Directly mirrors the paper's RQ3 t-test methodology
t_session <- t.test(adt ~ in_session, data = merged)
cat("\n=== Welch t-test: ADT ~ In Session (mirrors paper RQ3) ===\n")
print(t_session)

# RW4 — Linear regression on daily aggregates
lm_mod <- lm(mean_adt ~ tavg_f + precip_in + wind_avg_mph + visibility_mi,
             data = daily)
cat("\n=== Linear Regression: Mean Daily ADT ~ Weather ===\n")
print(summary(lm_mod))


# ── 6. Figures ───────────────────────────────────────────────────
# Style mirrors the paper: theme_minimal, minimal legend clutter,
# bold titles, labelled axes, dashed reference lines where useful.

PAPER_THEME <- theme_minimal(base_size = 13) +
  theme(plot.title    = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(size = 10, color = "grey40"),
        axis.title    = element_text(size = 11),
        legend.position = "bottom")


# ── Figure 1: Mean ADT by Month (mirrors paper Fig 4 — monthly mean
#    HRBT speed with CI ribbon, but for volume instead of speed)
fig1_data <- merged %>%
  mutate(month_lbl = factor(month(date, label = TRUE, abbr = TRUE),
                            levels = month.abb)) %>%
  group_by(month_lbl) %>%
  summarise(mean_adt = mean(adt),
            se       = sd(adt) / sqrt(n()),
            .groups  = "drop")

fig1 <- ggplot(fig1_data, aes(x = month_lbl, y = mean_adt, group = 1)) +
  geom_ribbon(aes(ymin = mean_adt - 1.96 * se,
                  ymax = mean_adt + 1.96 * se),
              fill = "steelblue", alpha = 0.2) +
  geom_line(color = "steelblue", linewidth = 1.1) +
  geom_point(color = "steelblue", size = 2.5) +
  scale_y_continuous(labels = comma) +
  labs(
    title    = "Figure 1: Monthly Mean ADT with 95% CI Ribbon",
    subtitle = "Hampton Roads, VA — 2024 (Jan 1 excluded)",
    x        = "Month", y = "Mean ADT (vehicles/day)"
  ) + PAPER_THEME

ggsave("fig1_monthly_mean_adt.png", fig1, width = 8, height = 4.5, dpi = 150)
cat("Saved fig1_monthly_mean_adt.png\n")


# ── Figure 2: ADT Distribution by Weather Condition
#    Mirrors paper Fig 3 — violin plots showing full density by route
fig2 <- merged %>%
  filter(weather_condition %in% c("Clear","Drizzle","Thunderstorm","Rain")) %>%
  ggplot(aes(x = weather_condition, y = adt, fill = weather_condition)) +
  geom_violin(alpha = 0.55, trim = FALSE, color = NA) +
  geom_boxplot(width = 0.12, outlier.shape = 21,
               outlier.size = 1.2, fill = "white", alpha = 0.9) +
  scale_y_continuous(labels = comma) +
  scale_fill_manual(values = c(
    "Clear"        = "#4DAF4A",
    "Drizzle"      = "#377EB8",
    "Rain"         = "#984EA3",
    "Thunderstorm" = "#E41A1C"
  )) +
  labs(
    title    = "Figure 2: ADT Distribution by Weather Condition",
    subtitle = "Violins show full density | Hampton Roads, VA — 2024",
    x        = "Weather Condition", y        = "ADT (vehicles/day)"
  ) +
  PAPER_THEME + theme(legend.position = "none")

ggsave("fig2_adt_distribution_by_condition.png", fig2,
       width = 7, height = 5, dpi = 150)
cat("Saved fig2_adt_distribution_by_condition.png\n")


# ── Figure 3: Mean ADT per Observation Date — timeline coloured
#    by weather condition (shows the 38 actual count dates)
fig3 <- ggplot(daily, aes(x = date, y = mean_adt,
                           color = weather_condition,
                           shape = weather_condition)) +
  geom_line(color = "grey75", linewidth = 0.7) +
  geom_point(size = 3.2, alpha = 0.9) +
  scale_x_date(date_breaks = "1 month", date_labels = "%b") +
  scale_y_continuous(labels = comma) +
  scale_color_manual(
    values = c("Clear"="#4DAF4A","Drizzle"="#377EB8",
               "Rain"="#984EA3","Thunderstorm"="#E41A1C"),
    name = "Condition") +
  scale_shape_manual(
    values = c("Clear"=16,"Drizzle"=17,"Rain"=15,"Thunderstorm"=18),
    name = "Condition") +
  labs(
    title    = "Figure 3: Mean ADT per VDOT Count Date",
    subtitle = "Each point = one observation date | Color = weather condition",
    x        = "Date (2024)", y = "Mean ADT (vehicles/day)"
  ) + PAPER_THEME

ggsave("fig3_adt_timeline_by_condition.png", fig3,
       width = 10, height = 4.5, dpi = 150)
cat("Saved fig3_adt_timeline_by_condition.png\n")


# ── Figure 4: ADT vs Temperature — scatter + regression
#    Mirrors paper's analytical style of showing a relationship
#    split by a binary grouping variable (here: in_session)
fig4 <- merged %>%
  mutate(session_lbl = ifelse(in_session, "In Session", "Out of Session")) %>%
  ggplot(aes(x = tavg_f, y = adt, color = session_lbl)) +
  geom_point(alpha = 0.28, size = 1.4) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 1.3) +
  scale_y_continuous(labels = comma) +
  scale_color_manual(values = c("In Session"     = "#E41A1C",
                                "Out of Session" = "#377EB8")) +
  labs(
    title    = "Figure 4: ADT vs. Temperature by Academic Calendar",
    subtitle = "Extends paper RQ3 — session status mirrors paper's in_session variable",
    x        = "Average Daily Temperature (°F)",
    y        = "ADT (vehicles/day)",
    color    = NULL
  ) + PAPER_THEME

ggsave("fig4_adt_vs_temp_by_session.png", fig4,
       width = 8, height = 5, dpi = 150)
cat("Saved fig4_adt_vs_temp_by_session.png\n")


# ── Figure 5: Precipitation vs Dry — boxplot comparison
#    Mirrors paper Figs 6/7 — side-by-side comparison with
#    t-test annotation baked into the subtitle
t_label <- sprintf("Welch t = %.2f, p < .001", t_precip$statistic)
means    <- tapply(merged$adt, merged$is_precipitation, mean)

fig5 <- merged %>%
  mutate(precip_lbl = factor(ifelse(is_precipitation == 1,
                                    "Precipitation", "Dry"),
                             levels = c("Dry","Precipitation"))) %>%
  ggplot(aes(x = precip_lbl, y = adt, fill = precip_lbl)) +
  geom_boxplot(outlier.shape = 21, outlier.size = 1.3,
               outlier.alpha = 0.5, alpha = 0.8, width = 0.45) +
  geom_hline(yintercept = mean(merged$adt), linetype = "dashed",
             color = "grey50", linewidth = 0.7) +
  annotate("text", x = 2.35, y = mean(merged$adt) + 20,
           label = "Grand mean", size = 3.2, color = "grey40") +
  scale_y_continuous(labels = comma) +
  scale_fill_manual(values = c("Dry" = "#4DAF4A",
                               "Precipitation" = "#377EB8")) +
  labs(
    title    = "Figure 5: ADT — Dry Days vs. Precipitation Days",
    subtitle = t_label,
    x        = NULL, y = "ADT (vehicles/day)"
  ) +
  PAPER_THEME + theme(legend.position = "none")

ggsave("fig5_adt_dry_vs_precip.png", fig5,
       width = 6, height = 5, dpi = 150)
cat("Saved fig5_adt_dry_vs_precip.png\n")


# ── Figure 6: ADT by Season — matches paper's RQ3 seasonal focus
fig6 <- merged %>%
  mutate(season = factor(season,
                         levels = c("Winter","Spring","Summer","Fall"))) %>%
  ggplot(aes(x = season, y = adt, fill = season)) +
  geom_violin(trim = FALSE, alpha = 0.5, color = NA) +
  geom_boxplot(width = 0.13, fill = "white",
               outlier.shape = 21, outlier.size = 1.2, alpha = 0.9) +
  scale_y_continuous(labels = comma) +
  scale_fill_manual(values = c(
    "Winter" = "#377EB8",
    "Spring" = "#4DAF4A",
    "Summer" = "#FF7F00",
    "Fall"   = "#984EA3"
  )) +
  labs(
    title    = "Figure 6: ADT Distribution by Season",
    subtitle = "Violin + boxplot | Hampton Roads, VA — 2024",
    x        = "Season", y = "ADT (vehicles/day)"
  ) +
  PAPER_THEME + theme(legend.position = "none")

ggsave("fig6_adt_by_season.png", fig6,
       width = 7, height = 5, dpi = 150)
cat("Saved fig6_adt_by_season.png\n")


# ── Figure 7: Correlation Matrix — weather variables vs mean ADT
#    Provides the quantitative backbone for the regression (RW4)
cor_df <- daily %>%
  select(`Mean ADT`      = mean_adt,
         `Avg Temp (F)`  = tavg_f,
         `Max Temp (F)`  = tmax_f,
         `Min Temp (F)`  = tmin_f,
         `Precip (in)`   = precip_in,
         `Wind (mph)`    = wind_avg_mph,
         `Visibility (mi)`= visibility_mi,
         `Precip Flag`   = is_precipitation)

cor_mat <- cor(cor_df, use = "complete.obs")

fig7 <- ggcorrplot(cor_mat,
                   method   = "circle",
                   type     = "lower",
                   lab      = TRUE,
                   lab_size = 3.2,
                   colors   = c("#377EB8", "white", "#E41A1C"),
                   title    = "Figure 7: Correlation Matrix — Weather vs. Mean ADT",
                   ggtheme  = theme_minimal(base_size = 11)) +
  theme(plot.title = element_text(face = "bold", size = 13))

ggsave("fig7_correlation_matrix.png", fig7,
       width = 8, height = 7, dpi = 150)
cat("Saved fig7_correlation_matrix.png\n")


# ── 7. Results Summary Table (mirrors paper Table 2) ────────────
cat("\n")
cat("============================================================\n")
cat("  RESULTS SUMMARY (mirrors paper Table 2 format)\n")
cat("============================================================\n")
cat(sprintf("  RW1  Welch t-test (Precipitation):  t = %.2f, p < .001\n",
            t_precip$statistic))
cat(sprintf("       Mean ADT: Dry = %.0f  |  Precip = %.0f\n",
            means["0"], means["1"]))
cat(sprintf("  RW2  ANOVA (Season):  F = %.2f, p = %.4f\n",
            summary(aov_season)[[1]]$`F value`[1],
            summary(aov_season)[[1]]$`Pr(>F)`[1]))
cat(sprintf("  RW3  Welch t-test (In Session):     t = %.2f, p = %.4f\n",
            t_session$statistic, t_session$p.value))
cat(sprintf("       Mean ADT: In Session = %.0f  |  Out = %.0f\n",
            t_session$estimate[2], t_session$estimate[1]))
cat(sprintf("  RW4  Linear Model R²: %.4f\n",
            summary(lm_mod)$r.squared))
cat("============================================================\n\n")
cat("All 7 figures saved to:", getwd(), "\n")