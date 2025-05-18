library(tidyverse)
library(scales)
library(stats)
library(corrplot)
library(car)
library(cluster)
library(plotly)

# Load TB dataset
tb_data <- read.csv("TB_Burden_Country.csv", stringsAsFactors = FALSE)

# Filter for latest available year in dataset
latest_year <- max(tb_data$Year, na.rm = TRUE)
tb_latest <- tb_data %>% filter(Year == latest_year)

# Select and clean relevant columns
tb_clean <- tb_latest %>%
  select(
    country = `Country.or.territory.name`,
    population = `Estimated.total.population.number`,
    prevalence_per_100k = `Estimated.prevalence.of.TB..all.forms..per.100.000.population`,
    mortality_per_100k = `Estimated.mortality.of.TB.cases..all.forms..excluding.HIV..per.100.000.population`,
    incidence_per_100k = `Estimated.incidence..all.forms..per.100.000.population`,
    hiv_percent = `Estimated.HIV.in.incident.TB..percent.`,
    case_detection_rate = `Case.detection.rate..all.forms...percent`
  ) %>%
  drop_na()

tb_clean[38, "country"] <- 'Korea'

# Correlation Analysis
numeric_vars <- tb_clean %>% 
  select(prevalence_per_100k, mortality_per_100k, incidence_per_100k, hiv_percent, case_detection_rate)

correlation_matrix <- cor(numeric_vars, use = "complete.obs")
correlation_plot <- corrplot(correlation_matrix, method = "color", 
                             type = "upper", order = "hclust",
                             addCoef.col = "black", number.cex = 0.7)

# Linear Modelling
lm_model1 <- lm(mortality_per_100k ~ incidence_per_100k + hiv_percent + prevalence_per_100k, data = tb_clean)
cooks_d <- cooks.distance(lm_model1)
outlier <- which(cooks_d > 1)
tb_clean_no_outlier <- tb_clean[-outlier, ]

lm_model2 <- lm(mortality_per_100k ~ incidence_per_100k + hiv_percent + prevalence_per_100k, data = tb_clean_no_outlier)
model_summary <- summary(lm_model2)
conf_intervals <- confint(lm_model2)

#Regression plot
regression_plot <- ggplot(tb_clean_no_outlier, aes(x = incidence_per_100k, y = mortality_per_100k)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, color = "darkred") +
  theme_minimal() +
  labs(title = "TB Mortality vs. Incidence with Confidence Interval",
       subtitle = paste("R² =", round(model_summary$r.squared, 3)),
       x = "Incidence per 100k", y = "Mortality per 100k")

#HIV prevalence in TB cases and Mortality plot 
hiv_mortality_plot <- ggplot(tb_clean, aes(x = hiv_percent, y = mortality_per_100k)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, color = "purple") +
  theme_minimal() +
  labs(
    title = "Relationship Between HIV Prevalence in TB cases and TB Mortality",
    x = "HIV % in Incident TB",
    y = "TB Mortality per 100k"
  )

#Bar Plot: Top 10 TB Incidence by Country
tb_top_10 <- tb_clean %>% arrange(desc(incidence_per_100k)) %>% slice(1:10)  
incidence_plot <- ggplot(tb_top_10, aes(x = reorder(country, -incidence_per_100k), y = incidence_per_100k)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  labs(title = "Top 10 TB Incidence Rate by Country", x = "Country", y = "Incidence Rate per 100,000")

#K-Means Clustering
scaled_data <- scale(numeric_vars)
kmeans_result <- kmeans(scaled_data, centers = 3, nstart = 25)

tb_clustered <- tb_clean %>%
  mutate(cluster = as.factor(kmeans_result$cluster))

#Cluster Summary
cluster_stats <- tb_clustered %>%
  group_by(cluster) %>%
  summarize(
    count = n(),
    mean_incidence = mean(incidence_per_100k),
    mean_mortality = mean(mortality_per_100k),
    mean_hiv_percent = mean(hiv_percent),
    mean_detection = mean(case_detection_rate)
  )

#Cluster Plot
tb_clustered <- tb_clustered %>%
  mutate(cluster_label = case_when(
    cluster == 3 ~ "High Burden",
    cluster == 1 ~ "Moderate Burden",
    cluster == 2 ~ "Low Burden"
  ))

plotly_cluster <- ggplot(tb_clustered, aes(x = incidence_per_100k, y = mortality_per_100k, color = cluster_label, text = country)) +
  geom_point(size = 3) +
  labs(title = 'Interactive Country Cluster', x = "Incidence Rate per 100k", y = 'Mortality Rate per 100k')
  theme_minimal() 

plotly_cluster <- ggplotly(plotly_cluster, tooltip = "text")

