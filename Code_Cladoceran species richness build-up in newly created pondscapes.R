#-------------------------------------------------------------------------------

# Read packages

#-------------------------------------------------------------------------------

library(lme4)
library(effects)
library(car)
library(ggplot2)
library(tidyverse)
library(RColorBrewer)
library(vegan)
library(export)
library(tidyverse)
library(brms)
library(patchwork)

#-------------------------------------------------------------------------------

# Load data

#-------------------------------------------------------------------------------

fysicochemie <- read.csv("FysicochemTNTP_oktober_2024.csv")
cumsumspec <- read.csv("Cumsumspec.csv")
inSR <- read.csv("insSR.csv")
datahydromacro <- read.csv("Hydromacrophyte.csv") %>%
  rename(Hydroperiod = percentage)

#-------------------------------------------------------------------------------

# Data wrangling fysicochemie

#-------------------------------------------------------------------------------

fysicochemie <- fysicochemie %>%
  mutate(
    Date = as.Date(Date, format = '%d/%m/%Y'),
    day = yday(Date),
    month = format(Date, '%m'),
    year = format(Date, '%Y'),
    meanCHLA = rowMeans(select(., 14:16), na.rm = TRUE),
    meanPC = rowMeans(select(., 17:19), na.rm = TRUE)
  ) %>%
  #calibration curves
  mutate(meanCHLA = case_when(
    Session %in% 1:2 ~ meanCHLA / 19.39,
    Session %in% 3:4 ~ meanCHLA / 20.25,
    Session %in% 5:11 ~ meanCHLA / 18.77
  )) %>%
  mutate(
    DSC = as.numeric(case_when(
      PondCode %in% c("RUL_1", "RUL_2") ~ difftime(Date, as.Date("01/06/2023", format = "%d/%m/%Y"), units = "days"),
      PondCode == "RUL_16" ~ difftime(Date, as.Date("01/05/2022", format = "%d/%m/%Y"), units = "days"),
      TRUE ~ difftime(Date, as.Date("01/11/2021", format = "%d/%m/%Y"), units = "days")
    ))
  ) %>%
  filter(DSC >= 0, Max.depth != 0) %>%
  select(-c(5, 12:19, 24)) %>%
  rename(Visit = Session, POND = PondCode) %>%
  mutate(
    across(c(1:3, 18:19), as.factor),
    across(c(5:11, 16), as.numeric),
    across(12:14, ~ . / 100)
  )


#-------------------------------------------------------------------------------

# PCA physicochemistry

#-------------------------------------------------------------------------------

# Create Macrophyte variable #
fysicochemie <- fysicochemie %>%
  mutate(Macrophyte = if_else(Emerg != 0 | Subm != 0 | Float != 0, "Present", "Absent"))

# PCA Preparation #
fysico_pca_data <- fysicochemie %>%
  select(POND, Pondscape, Visit, year, DSC, month, Macrophyte, Max.depth, Conductivity, Oxygen, Oxygen.., pH, Water.temperature, meanCHLA, meanPC, TN, TP) %>%
  na.omit()

# PCA Analysis (Global) #
pca_result <- prcomp(
  fysico_pca_data %>% select(-POND, -Pondscape, -Visit, -year, -month, -Macrophyte, -DSC, -Water.temperature, -Max.depth, -Oxygen..),
  center = TRUE, scale. = TRUE
)

pca_ind <- as.data.frame(pca_result$x) %>%
  mutate(
    Pondscape = fysico_pca_data$Pondscape,
    Pond = fysico_pca_data$POND,
    Visit = fysico_pca_data$Visit,
    DSC = fysico_pca_data$DSC,
    year = fysico_pca_data$year,
    month = fysico_pca_data$month,
    Hydroperiod = if_else(Pond %in% c("RUL1", "RUL2", "RUL16", "BOF3", "BOF4", "BOF5", "BOF6"), "Permanent", "Temporary"),
    Macrophyte = fysico_pca_data$Macrophyte
  )

pca_var <- as.data.frame(pca_result$rotation) %>%
  mutate(Variable = rownames(.))

# PCA Plot Global #
ggplot() +
  geom_point(data = pca_ind, aes(x = PC1, y = PC2, color = as.factor(year), shape = Pondscape), size = 3, alpha = 0.8) +
  stat_ellipse(data = pca_ind, aes(x = PC1, y = PC2, color = as.factor(year)), level = 0.95) +
  geom_segment(data = pca_var, aes(x = 0, y = 0, xend = PC1 * 5, yend = PC2 * 5),
               arrow = arrow(length = unit(0.2, "cm")), color = "black") +
  geom_text(data = pca_var, aes(x = PC1 * 5, y = PC2 * 5, label = Variable),
            vjust = 1.5, hjust = 1.2, size = 4) +
  scale_color_manual(values = c("#A7C7E7", "#74C69D", "#40916C")) +
  theme_minimal() +
  labs(
    x = paste0("PC1 (", round(summary(pca_result)$importance[2, 1] * 100, 1), "%)"),
    y = paste0("PC2 (", round(summary(pca_result)$importance[2, 2] * 100, 1), "%)"),
    color = "Year",
    shape = "Pondscape"
  ) +
  theme(legend.position = "right")
ggsave("PCARullenBoffereth.png", width = 14, height = 14, units = "cm", dpi = 600)


#-------------------------------------------------------------------------------

# PCA physicochemistry seperatly (Boffereth and Rullen)

#-------------------------------------------------------------------------------

# Boffereth #
fysico_pca_data_BOF <- fysico_pca_data %>% filter(Pondscape == "Boffereth")
pca_resultBOF <- prcomp(
  fysico_pca_data_BOF %>% select(-POND, -Pondscape, -Visit, -year, -month, -Macrophyte, -DSC, -Water.temperature, -Max.depth, -Oxygen..),
  center = TRUE, scale. = TRUE
)
pca_indBOF <- as.data.frame(pca_resultBOF$x) %>%
  mutate(year = fysico_pca_data_BOF$year)
pca_varBOF <- as.data.frame(pca_resultBOF$rotation) %>%
  mutate(Variable = rownames(.))

ggplot() +
  geom_point(data = pca_indBOF, aes(x = PC1, y = PC2, color = as.factor(year)), size = 3, alpha = 0.8) +
  stat_ellipse(data = pca_indBOF, aes(x = PC1, y = PC2, color = as.factor(year)), level = 0.95) +
  geom_segment(data = pca_varBOF, aes(x = 0, y = 0, xend = PC1 * 5, yend = PC2 * 5),
               arrow = arrow(length = unit(0.2, "cm")), color = "black") +
  geom_text(data = pca_varBOF, aes(x = PC1 * 5, y = PC2 * 5, label = Variable),
            vjust = 1.5, hjust = 1.2, size = 4) +
  scale_color_manual(values = c("#A7C7E7", "#74C69D", "#40916C")) +
  theme_minimal() +
  labs(
    x = paste0("PC1 (", round(summary(pca_resultBOF)$importance[2, 1] * 100, 1), "%)"),
    y = paste0("PC2 (", round(summary(pca_resultBOF)$importance[2, 2] * 100, 1), "%)"),
    color = "Year"
  ) +
  theme(legend.position = "right")
ggsave("PCABoffereth.png", width = 14, height = 14, units = "cm", dpi = 600)


# Rullen #
fysico_pca_data_RUL <- fysico_pca_data %>% filter(Pondscape == "Rullen")
pca_resultRUL <- prcomp(
  fysico_pca_data_RUL %>% select(-POND, -Pondscape, -Visit, -year, -month, -Macrophyte, -DSC, -Water.temperature, -Max.depth, -Oxygen..),
  center = TRUE, scale. = TRUE
)
pca_indRUL <- as.data.frame(pca_resultRUL$x) %>%
  mutate(year = fysico_pca_data_RUL$year)
pca_varRUL <- as.data.frame(pca_resultRUL$rotation) %>%
  mutate(Variable = rownames(.))


ggplot() +
  geom_point(data = pca_indRUL, aes(x = PC1, y = PC2, color = as.factor(year)), size = 3, alpha = 0.8) +
  stat_ellipse(data = pca_indRUL, aes(x = PC1, y = PC2, color = as.factor(year)), level = 0.95) +
  geom_segment(data = pca_varRUL, aes(x = 0, y = 0, xend = PC1 * 5, yend = PC2 * 5),
               arrow = arrow(length = unit(0.2, "cm")), color = "black") +
  geom_text(data = pca_varRUL, aes(x = PC1 * 5, y = PC2 * 5, label = Variable),
            vjust = 1.5, hjust = 1.2, size = 4) +
  scale_color_manual(values = c("#A7C7E7", "#74C69D", "#40916C")) +
  theme_minimal() +
  labs(
    x = paste0("PC1 (", round(summary(pca_resultRUL)$importance[2, 1] * 100, 1), "%)"),
    y = paste0("PC2 (", round(summary(pca_resultRUL)$importance[2, 2] * 100, 1), "%)"),
    color = "Year"
  ) +
  theme(legend.position = "right")
ggsave("PCARullen.png", width = 14, height = 14, units = "cm", dpi = 600)


#-------------------------------------------------------------------------------

# Mean cumulative species richness plot

#-------------------------------------------------------------------------------


cumsumspec$age <- as.factor(cumsumspec$age)
cumsumspec_long <- cumsumspec %>%
  pivot_longer(cols = -age, names_to = "pond", values_to = "species_richness") %>%
  group_by(age) %>%
  summarise(
    mean_richness = mean(species_richness, na.rm = TRUE),
    se_richness = sd(species_richness, na.rm = TRUE) / sqrt(sum(!is.na(species_richness)))
  )
cumsumspec_rullen <- cumsumspec %>%
  select(age, starts_with("RUL")) %>%
  pivot_longer(cols = -age, names_to = "pond", values_to = "species_richness") %>%
  group_by(age) %>%
  summarise(
    mean_richness = mean(species_richness, na.rm = TRUE),
    se_richness = sd(species_richness, na.rm = TRUE) / sqrt(sum(!is.na(species_richness)))
  )

cumsumspec_boffereth <- cumsumspec %>%
  select(age, starts_with("BOF")) %>%
  pivot_longer(cols = -age, names_to = "pond", values_to = "species_richness") %>%
  group_by(age) %>%
  summarise(
    mean_richness = mean(species_richness, na.rm = TRUE),
    se_richness = sd(species_richness, na.rm = TRUE) / sqrt(sum(!is.na(species_richness)))
  )

ggplot() +
  geom_point(data = cumsumspec_rullen, aes(x = as.numeric(as.character(age)), y = mean_richness, color = "Rullen"), size = 3) +
  geom_errorbar(data = cumsumspec_rullen, aes(x = as.numeric(as.character(age)), ymin = mean_richness - se_richness, ymax = mean_richness + se_richness, color = "Rullen"), width = 0.2) +
  geom_point(data = cumsumspec_boffereth, aes(x = as.numeric(as.character(age)), y = mean_richness, color = "Boffereth"), size = 3, shape = 17) +  # shape = 17 for triangles
  geom_errorbar(data = cumsumspec_boffereth, aes(x = as.numeric(as.character(age)), ymin = mean_richness - se_richness, ymax = mean_richness + se_richness, color = "Boffereth"), width = 0.2) +
  labs(x = "Age (month)", y = "Mean cumulative species richness", color = "Pondscape") +
  scale_color_manual(values = c("Rullen" = "#A7C7E7", 
                                "Boffereth" = "#74C69D")) +  
  theme_minimal()
ggsave("Mean cumulative species richness.png", width = 16, height = 10, units = "cm", dpi = 600)


#-------------------------------------------------------------------------------

# Put together dataset to use in models

#-------------------------------------------------------------------------------

# Define ponds with macrophytes #
poelen_met_macrofyten <- c("BOF1", "BOF2", "BOF3", "BOF4", "BOF5", "BOF6",
                           "RUL16", "RUL4", "RUL5", "RUL6")

# Make long format of cumulative richness data #
long_data <- cumsumspec %>%
  gather(key = "pond", value = "cumulative_individuals", matches("BOF|RUL")) %>%
  filter(!is.na(cumulative_individuals)) %>%
  mutate(
    pondscape = ifelse(grepl("^RUL", pond), "Rullen", "Boffereth"),
    Time = as.numeric(as.character(age)) / 12
  ) %>%
  mutate(
    Visit = case_when(
      pond %in% c("RUL1", "RUL2") ~ match(Time * 12, c(4, 5, 7, 11, 17)) + 6,
      pond == "RUL16" ~ match(Time * 12, c(4, 5, 7, 11, 17, 19, 22, 24, 28)) + 2,
      TRUE ~ match(Time * 12, c(4, 5, 7, 11, 17, 19, 22, 24, 28, 32, 34))
    )
  ) %>%
  select(-age)

# Get ponds with emergent vegetation cover #
ponds_met_emerg_df <- fysicochemie %>%
  filter(Emerg >= 0.25) %>%
  distinct(POND) %>%
  arrange(POND)

# Add macrophytes and hydroperiod information #
long_data <- long_data %>%
  mutate(macrophytes = if_else(pond %in% poelen_met_macrofyten, 1, 0)) %>%
  left_join(data %>% select(POND, Hydroperiod), by = c("pond" = "POND"))

# Reshape inSR (Instantaneous species richness) to long format #
inSR_long <- inSR %>%
  mutate(across(-Visit, as.character)) %>%
  pivot_longer(
    cols = -Visit,
    names_to = "pond",
    values_to = "value"
  ) %>%
  mutate(
    Visit = as.factor(Visit),
    inSR_value = as.numeric(value)
  ) %>%
  select(-value)

# Join pca scores and inSR values #
long_data <- long_data %>%
  mutate(Visit = as.factor(Visit)) %>%
  left_join(
    pca_ind %>%
      mutate(Visit = as.factor(Visit)) %>%
      select(Visit, Pond, PC1, PC2),
    by = c("Visit", "pond" = "Pond")
  ) %>%
  left_join(inSR_long, by = c("Visit", "pond"))

# Add emergent vegetation #
long_data <- long_data %>%
  left_join(fysicochemie %>% select(Visit, POND, Emerg),
            by = c("Visit", "pond" = "POND"))

# Rename variables #
long_data <- long_data %>%
  rename(
    PondID = pond,
    Pondscape = pondscape,
    Age = Time,
    Hydroperiod = Hydroperiod,
    Macrophytes = macrophytes,
    cumSR = cumulative_individuals,
    insSR = inSR_value
  ) %>%
  select(PondID, Pondscape, Age, Hydroperiod, Macrophytes, PC1, PC2, cumSR, insSR, everything())

# Add Year variable #
long_data <- long_data %>%
  mutate(Year = ceiling(Age)) %>%
  relocate(Year, .after = Age)

# Calculate overall PC1 and PC2 means per pond #
long_data <- long_data %>%
  group_by(PondID) %>%
  mutate(
    PCmeanoverall1 = mean(PC1, na.rm = TRUE),
    PCmeanoverall2 = mean(PC2, na.rm = TRUE)
  ) %>%
  ungroup()

# Calculate PC means per pond per year (years 1-3) #
pc_means_years <- long_data %>%
  filter(Year %in% 1:3) %>%
  group_by(PondID, Year) %>%
  summarise(
    PC1 = mean(PC1, na.rm = TRUE),
    PC2 = mean(PC2, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = Year,
    values_from = c(PC1, PC2),
    names_glue = "PCmeanyear{Year}.{.value}"
  )

# Add year-specific PC means to long_data #
long_data <- long_data %>%
  left_join(pc_means_years, by = "PondID")

# Save to file #
write.csv(long_data, "data_modelvoeren.csv", row.names = FALSE)


#-------------------------------------------------------------------------------

# Model: Read and preprocess data

#-------------------------------------------------------------------------------

data <- read.csv("data_modelvoeren.csv") %>%
  select(PondID, Pondscape, Age, Hydroperiod, Macrophytes, PCmeanoverall1, PCmeanoverall2, cumSR) %>%
  mutate(Hydroperiod = Hydroperiod/100,
         PCmeanoverall1 = c(scale(PCmeanoverall1)),
         PCmeanoverall2 = c(scale(PCmeanoverall2)))

# Augment the data with the known species richness (0) at age 0

data <- rbind(data %>%
                group_by(PondID) %>%
                filter(row_number() == 1) %>%
                mutate(Age = 0,
                       cumSR = 0),
              data)

#-------------------------------------------------------------------------------

# Fit the Bayesian LMM

#-------------------------------------------------------------------------------

fit <- brm(cumSR ~ 0 + Age + Age:Pondscape + Age:Hydroperiod + Age:Macrophytes + Age:PCmeanoverall1 + Age:PCmeanoverall2 + (0 + Age | PondID),
           data = data,
           prior = c(prior(normal(0, 1), class = "b"),
                     prior(normal(0, 1), class = "sd"),
                     prior(normal(0, 1), class = "sigma")))

#-------------------------------------------------------------------------------

# First exploration of model outcomes, incl. convergence and goodness-of-fit checks

#-------------------------------------------------------------------------------

summary(fit)
plot(fit)
pp_check(fit)
conditional_effects(fit)

#-------------------------------------------------------------------------------

# Visualization of accumulation across conditions

#-------------------------------------------------------------------------------

conditional_effects(fit, effects = "Age", conditions = data.frame(Pondscape = c("Boffereth","Rullen"), Hydroperiod = mean(distinct(data, PondID, Hydroperiod)$Hydroperiod), Macrophytes = mean(distinct(data, PondID, Macrophytes)$Macrophytes)))[[1]] %>%
  group_by(Age) %>%
  summarize(estimate__ = mean(estimate__),
            lower__ = mean(lower__),
            upper__ = mean(upper__)) %>%
  ggplot(aes(x = Age)) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__), alpha = 0.3, color = NA) +
  geom_point(shape = 16, size = 0.25, data = data, aes(y = cumSR), position = position_jitter(width = 0.05, height = 0.1, seed = 0), alpha = 0.3) +
  geom_path(linewidth = 0.2, data = data, aes(y = cumSR, group = PondID), position = position_jitter(width = 0.05, height = 0.1, seed = 0), alpha = 0.3) +
  geom_line(aes(y = estimate__)) +
  scale_x_continuous("Age (years)", expand = c(0,0), limits = c(0,3)) +
  scale_y_continuous("Cumulative species richness", expand = c(0,0), limits = c(0,12)) +
  coord_cartesian(clip = "off") +
  ggtitle("(a) Overall") +
  theme(plot.title = element_text(face = "bold", size = 8),
        panel.background = element_blank(),
        panel.grid.major = element_line(color = "grey93", linewidth = 0.3),
        panel.grid.minor = element_line(color = "grey93", linewidth = 0.15),
        axis.line.x = element_line(color = "black", linewidth = 0.3),
        axis.title = element_text(face = "bold", size = 8),
        axis.text = element_text(size = 7),
        axis.ticks = element_line(linewidth = 0.3),
        legend.title = element_blank(),
        legend.text = element_text(size = 7),
        legend.position.inside = c(0.37, 0.87),
        legend.direction = "horizontal",
        legend.key = element_blank(),
        legend.key.size = unit(0.4, "cm")) +
  
  conditional_effects(fit, effects = "Age:Pondscape", conditions = data.frame(Hydroperiod = mean(distinct(data, PondID, Hydroperiod)$Hydroperiod), Macrophytes = mean(distinct(data, PondID, Macrophytes)$Macrophytes)))[[1]] %>%
  ggplot(aes(x = Age, color = Pondscape, fill = Pondscape, group = Pondscape)) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__), alpha = 0.3, color = NA) +
  geom_point(shape = 16, size = 0.25, data = data, aes(y = cumSR), position = position_jitter(width = 0.05, height = 0.1, seed = 0), alpha = 0.3) +
  geom_path(linewidth = 0.2, data = data, aes(y = cumSR, group = PondID), position = position_jitter(width = 0.05, height = 0.1, seed = 0), alpha = 0.3) +
  geom_line(aes(y = estimate__)) +
  scale_x_continuous("Age (years)", expand = c(0,0), limits = c(0,3)) +
  scale_y_continuous("Cumulative species richness", expand = c(0,0), limits = c(0,12)) +
  scale_color_manual(values = c("#3288bd","#d53e4f"), guide = guide_legend(position = "inside")) +
  scale_fill_manual(values = c("#3288bd","#d53e4f"), guide = guide_legend(position = "inside")) +
  coord_cartesian(clip = "off") +
  ggtitle("(b) Pondscape") +
  theme(plot.title = element_text(face = "bold", size = 8),
        panel.background = element_blank(),
        panel.grid.major = element_line(color = "grey93", linewidth = 0.3),
        panel.grid.minor = element_line(color = "grey93", linewidth = 0.15),
        axis.line.x = element_line(color = "black", linewidth = 0.3),
        axis.title = element_text(face = "bold", size = 8),
        axis.text = element_text(size = 7),
        axis.ticks = element_line(linewidth = 0.3),
        legend.title = element_blank(),
        legend.text = element_text(size = 7),
        legend.position.inside = c(0.37, 0.87),
        legend.direction = "horizontal",
        legend.key = element_blank(),
        legend.key.size = unit(0.3, "cm")) +
  
  conditional_effects(fit, effects = "Age:Hydroperiod", conditions = data.frame(Pondscape = c("Boffereth","Rullen"), Macrophytes = mean(distinct(data, PondID, Macrophytes)$Macrophytes)), int_conditions = list(Hydroperiod = c(min(data$Hydroperiod), max(data$Hydroperiod))))[[1]] %>%
  group_by(Age, Hydroperiod) %>%
  summarize(estimate__ = mean(estimate__),
            lower__ = mean(lower__),
            upper__ = mean(upper__)) %>%
  ggplot(aes(x = Age, color = Hydroperiod*100, fill = Hydroperiod*100, group = Hydroperiod*100)) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__), alpha = 0.3, color = NA) +
  geom_point(shape = 16, size = 0.25, data = data, aes(y = cumSR), position = position_jitter(width = 0.05, height = 0.1, seed = 0), alpha = 0.3) +
  geom_path(linewidth = 0.2, data = data, aes(y = cumSR, group = PondID), position = position_jitter(width = 0.05, height = 0.1, seed = 0), alpha = 0.3) +
  geom_line(aes(y = estimate__)) +
  scale_x_continuous("Age (years)", expand = c(0,0), limits = c(0,3)) +
  scale_y_continuous("Cumulative species richness", expand = c(0,0), limits = c(0,12)) +
  scale_color_distiller(palette = "Spectral", guide = guide_colorbar(position = "inside", barheight = unit(0.1, "cm")), breaks = c(0.5,0.75,1)*100, labels = c("50%", "75%", "100%")) +
  scale_fill_distiller(palette = "Spectral", guide = guide_colorbar(position = "inside", barheight = unit(0.1, "cm")), breaks = c(0.5,0.75,1)*100, labels = c("50%", "75%", "100%")) +
  coord_cartesian(clip = "off") +
  ggtitle("(c) Hydroperiod") +
  theme(
    plot.title = element_text(face = "bold", size = 8),
    panel.background = element_blank(),
    panel.grid.major = element_line(color = "grey93", linewidth = 0.3),
    panel.grid.minor = element_line(color = "grey93", linewidth = 0.15),
    axis.line.x = element_line(color = "black", linewidth = 0.3),
    axis.title = element_text(face = "bold", size = 8),
    axis.text = element_text(size = 7),
    axis.ticks = element_line(linewidth = 0.3),
    legend.title = element_blank(),
    legend.text = element_text(size = 7),
    legend.position.inside = c(0.37, 0.87),
    legend.direction = "horizontal",
    legend.key = element_blank(),
    legend.key.size = unit(0.4, "cm")) +

  conditional_effects(fit, effects = "Age:Macrophytes", conditions = data.frame(Pondscape = c("Boffereth","Rullen"), Hydroperiod = mean(distinct(data, PondID, Hydroperiod)$Hydroperiod)), int_conditions = list(Macrophytes = c(0,1)))[[1]] %>%
  group_by(Age, Macrophytes) %>%
  summarize(estimate__ = mean(estimate__), lower__ = mean(lower__), upper__ = mean(upper__)) %>%
  ggplot(aes(x = Age, color = factor(Macrophytes), fill = factor(Macrophytes), group = factor(Macrophytes))) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__, color = factor(Macrophytes), fill = factor(Macrophytes)), alpha = 0.3, color = NA) +
  geom_point(shape = 16, size = 0.25, data = data, aes(y = cumSR, color = factor(Macrophytes), fill = factor(Macrophytes)), position = position_jitter(width = 0.05, height = 0.1, seed = 0), alpha = 0.3) +
  geom_path(linewidth = 0.2, data = data, aes(y = cumSR, group = PondID, color = factor(Macrophytes), fill = factor(Macrophytes)), position = position_jitter(width = 0.05, height = 0.1, seed = 0), alpha = 0.3) +
  geom_line(aes(y = estimate__, color = factor(Macrophytes), fill = factor(Macrophytes))) +
  scale_x_continuous("Age (years)", expand = c(0,0), limits = c(0,3)) +
  scale_y_continuous("Cumulative species richness", expand = c(0,0), limits = c(0,12)) +
  scale_color_manual(values = c("#3288bd","#d53e4f"), guide = guide_legend(position = "inside"), limits = c(0,1), breaks = c(0,1), labels = c("<25%", ">25%")) +
  scale_fill_manual(values = c("#3288bd","#d53e4f"), guide = guide_legend(position = "inside"), limits = c(0,1), breaks = c(0,1), labels = c("<25%", ">25%")) +
  coord_cartesian(clip = "off") +
  ggtitle("(d) Macrophyte coverage") +
  theme(
    plot.title = element_text(face = "bold", size = 8),
    panel.background = element_blank(),
    panel.grid.major = element_line(color = "grey93", linewidth = 0.3),
    panel.grid.minor = element_line(color = "grey93", linewidth = 0.15),
    axis.line.x = element_line(color = "black", linewidth = 0.3),
    axis.title = element_text(face = "bold", size = 8),
    axis.text = element_text(size = 7),
    axis.ticks = element_line(linewidth = 0.3),
    legend.title = element_blank(),
    legend.text = element_text(size = 7),
    legend.position.inside = c(0.37, 0.87),
    legend.direction = "horizontal",
    legend.key = element_blank(),
    legend.key.size = unit(0.3, "cm")) + 
  
  conditional_effects(fit, effects = "Age:PCmeanoverall1", conditions = data.frame(Pondscape = c("Boffereth","Rullen"), Hydroperiod = mean(distinct(data, PondID, Hydroperiod)$Hydroperiod), Macrophytes = mean(distinct(data, PondID, Macrophytes)$Macrophytes)), int_conditions = list(PCmeanoverall1 = c(min(data$PCmeanoverall1), max(data$PCmeanoverall1))))[[1]] %>%
  group_by(Age, PCmeanoverall1) %>%
  summarize(estimate__ = mean(estimate__),
            lower__ = mean(lower__),
            upper__ = mean(upper__)) %>%
  ggplot(aes(x = Age, color = PCmeanoverall1, fill = PCmeanoverall1, group = PCmeanoverall1)) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__), alpha = 0.3, color = NA) +
  geom_point(shape = 16, size = 0.25, data = data, aes(y = cumSR), position = position_jitter(width = 0.05, height = 0.1, seed = 0), alpha = 0.3) +
  geom_path(linewidth = 0.2, data = data, aes(y = cumSR, group = PondID), position = position_jitter(width = 0.05, height = 0.1, seed = 0), alpha = 0.3) +
  geom_line(aes(y = estimate__)) +
  scale_x_continuous("Age (years)", expand = c(0,0), limits = c(0,3)) +
  scale_y_continuous("Cumulative species richness", expand = c(0,0), limits = c(0,12)) +
  scale_color_distiller(palette = "Spectral", limits = c(-max(abs(data$PCmeanoverall1)), max(abs(data$PCmeanoverall1))), guide = guide_colorbar(position = "inside", barheight = unit(0.1, "cm"))) +
  scale_fill_distiller(palette = "Spectral", limits = c(-max(abs(data$PCmeanoverall1)), max(abs(data$PCmeanoverall1))), guide = guide_colorbar(position = "inside", barheight = unit(0.1, "cm"))) +
  coord_cartesian(clip = "off") +
  ggtitle("(e) Environmental PC1") +
  theme(plot.title = element_text(face = "bold", size = 8),
        panel.background = element_blank(),
        panel.grid.major = element_line(color = "grey93", linewidth = 0.3),
        panel.grid.minor = element_line(color = "grey93", linewidth = 0.15),
        axis.line.x = element_line(color = "black", linewidth = 0.3),
        axis.title = element_text(face = "bold", size = 8),
        axis.text = element_text(size = 7),
        axis.ticks = element_line(linewidth = 0.3),
        legend.title = element_blank(),
        legend.text = element_text(size = 7),
        legend.position.inside = c(0.37, 0.87),
        legend.direction = "horizontal",
        legend.key = element_blank(),
        legend.key.size = unit(0.4, "cm")) +
  
  conditional_effects(fit, effects = "Age:PCmeanoverall2", conditions = data.frame(Pondscape = c("Boffereth","Rullen"), Hydroperiod = mean(distinct(data, PondID, Hydroperiod)$Hydroperiod), Macrophytes = mean(distinct(data, PondID, Macrophytes)$Macrophytes)), int_conditions = list(PCmeanoverall2 = c(min(data$PCmeanoverall2), max(data$PCmeanoverall2))))[[1]] %>%
  group_by(Age, PCmeanoverall2) %>%
  summarize(estimate__ = mean(estimate__),
            lower__ = mean(lower__),
            upper__ = mean(upper__)) %>%
  ggplot(aes(x = Age, color = PCmeanoverall2, fill = PCmeanoverall2, group = PCmeanoverall2)) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__), alpha = 0.3, color = NA) +
  geom_point(shape = 16, size = 0.25, data = data, aes(y = cumSR), position = position_jitter(width = 0.05, height = 0.1, seed = 0), alpha = 0.3) +
  geom_path(linewidth = 0.2, data = data, aes(y = cumSR, group = PondID), position = position_jitter(width = 0.05, height = 0.1, seed = 0), alpha = 0.3) +
  geom_line(aes(y = estimate__)) +
  scale_x_continuous("Age (years)", expand = c(0,0), limits = c(0,3)) +
  scale_y_continuous("Cumulative species richness", expand = c(0,0), limits = c(0,12)) +
  scale_color_distiller(palette = "Spectral", limits = c(-max(abs(data$PCmeanoverall2)), max(abs(data$PCmeanoverall2))), guide = guide_colorbar(position = "inside", barheight = unit(0.1, "cm"))) +
  scale_fill_distiller(palette = "Spectral", limits = c(-max(abs(data$PCmeanoverall2)), max(abs(data$PCmeanoverall2))), guide = guide_colorbar(position = "inside", barheight = unit(0.1, "cm"))) +
  coord_cartesian(clip = "off") +
  ggtitle("(f) Environmental PC2") +
  theme(plot.title = element_text(face = "bold", size = 8),
        panel.background = element_blank(),
        panel.grid.major = element_line(color = "grey93", linewidth = 0.3),
        panel.grid.minor = element_line(color = "grey93", linewidth = 0.15),
        axis.line.x = element_line(color = "black", linewidth = 0.3),
        axis.title = element_text(face = "bold", size = 8),
        axis.text = element_text(size = 7),
        axis.ticks = element_line(linewidth = 0.3),
        legend.title = element_blank(),
        legend.text = element_text(size = 7),
        legend.position.inside = c(0.37, 0.87),
        legend.direction = "horizontal",
        legend.key = element_blank(),
        legend.key.size = unit(0.4, "cm")) +
  
  plot_layout(ncol = 3, axis_titles = "collect", axes = "collect")

ggsave("Slope_plots.png", width = 16, height = 10, units = "cm", dpi = 600)


#-------------------------------------------------------------------------------

# Quantitative model output

#-------------------------------------------------------------------------------


# Extract all regression coefficients

fixedeffects <- fixef(fit, summary = F)

# Derive meaningful slope estimates for each setting (ceteris paribus)

slope_mean <- fixedeffects[,1] + mean(distinct(data, PondID, Hydroperiod)$Hydroperiod)*fixedeffects[,4] + mean(distinct(data, PondID, Macrophytes)$Macrophytes)*fixedeffects[,5]
slope_boffereth <- fixedeffects[,1] + fixedeffects[,2] + mean(distinct(data, PondID, Hydroperiod)$Hydroperiod)*fixedeffects[,4] + mean(distinct(data, PondID, Macrophytes)$Macrophytes)*fixedeffects[,5]
slope_rullen <- fixedeffects[,1] + fixedeffects[,3] + mean(distinct(data, PondID, Hydroperiod)$Hydroperiod)*fixedeffects[,4] + mean(distinct(data, PondID, Macrophytes)$Macrophytes)*fixedeffects[,5]
slope_permanent <- fixedeffects[,1] + fixedeffects[,4] + mean(distinct(data, PondID, Macrophytes)$Macrophytes)*fixedeffects[,5]
slope_temporary <- fixedeffects[,1] + 0.5*fixedeffects[,4] + mean(distinct(data, PondID, Macrophytes)$Macrophytes)*fixedeffects[,5]
slope_macrophytes <- fixedeffects[,1] + mean(distinct(data, PondID, Hydroperiod)$Hydroperiod)*fixedeffects[,4] + fixedeffects[,5]
slope_nomacrophytes <- fixedeffects[,1] + mean(distinct(data, PondID, Hydroperiod)$Hydroperiod)*fixedeffects[,4]

# Average slope for the average condition

mean(slope_mean);quantile(slope_mean, 0.025);quantile(slope_mean, 0.975)

# Average slope for a Boffereth pond wtih an average condition

mean(slope_boffereth);quantile(slope_boffereth, 0.025);quantile(slope_boffereth, 0.975)

# Average slope for a Rullen pond with an average condition

mean(slope_rullen);quantile(slope_rullen, 0.025);quantile(slope_rullen, 0.975)

# Average slope for a permanent pond with an average condition

mean(slope_permanent);quantile(slope_permanent, 0.025);quantile(slope_permanent, 0.975)

# Average slope for a temporary pond with an average condition

mean(slope_temporary);quantile(slope_temporary, 0.025);quantile(slope_temporary, 0.975)

# Average slope for a macrophyte pond with an average condition

mean(slope_macrophytes);quantile(slope_macrophytes, 0.025);quantile(slope_macrophytes, 0.975)

# Average slope for a no-macrophyte pond with an average condition

mean(slope_nomacrophytes);quantile(slope_nomacrophytes, 0.025);quantile(slope_nomacrophytes, 0.975)



# Effect of Boffereth vs Rullen

mean(slope_boffereth - slope_rullen);quantile(slope_boffereth - slope_rullen, 0.025);quantile(slope_boffereth - slope_rullen, 0.975); mean(slope_boffereth - slope_rullen > 0)

# Effect of permanent vs temporary

mean(slope_permanent - slope_temporary);quantile(slope_permanent - slope_temporary, 0.025);quantile(slope_permanent - slope_temporary, 0.975); mean(slope_permanent - slope_temporary > 0)

# Effect of macrophytes vs no macrophytes

mean(slope_macrophytes - slope_nomacrophytes);quantile(slope_macrophytes - slope_nomacrophytes, 0.025);quantile(slope_macrophytes - slope_nomacrophytes, 0.975); mean(slope_macrophytes - slope_nomacrophytes > 0)

# Effect of a one standard deviation increase in PC1

mean(fixedeffects[,6]);quantile(fixedeffects[,6], 0.025);quantile(fixedeffects[,6], 0.975);mean(fixedeffects[,6] > 0)

# Effect of a one standard deviation increase in PC2

mean(fixedeffects[,7]);quantile(fixedeffects[,7], 0.025);quantile(fixedeffects[,7], 0.975);mean(fixedeffects[,7] > 0)



# Explore correlated joint posterior effects of hydroperiod and macrophytes

plot(fixedeffects[,4], fixedeffects[,5])
cor.test(fixedeffects[,4], fixedeffects[,5])

mean((slope_permanent - slope_temporary) > 0 & (slope_macrophytes - slope_nomacrophytes) > 0)
mean((slope_permanent - slope_temporary) < 0 & (slope_macrophytes - slope_nomacrophytes) > 0)
mean((slope_permanent - slope_temporary) > 0 & (slope_macrophytes - slope_nomacrophytes) < 0)
mean((slope_permanent - slope_temporary) < 0 & (slope_macrophytes - slope_nomacrophytes) < 0)

data.frame(Hydroperiod = slope_permanent - slope_temporary,
           Macrophytes = slope_macrophytes - slope_nomacrophytes) %>%
  ggplot() +
  geom_point(aes(x = Hydroperiod, y = Macrophytes), shape = 16, size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  geom_text(data = data.frame(Hydroperiod = c(2,-2,2,-2),
                              Macrophytes = c(2,2,-2,-2),
                              label = c(paste0(sprintf("%0.1f", 100*mean((slope_permanent - slope_temporary) > 0 & (slope_macrophytes - slope_nomacrophytes) > 0)), "% posterior\nprobability\n of a positive effect\nof hydroperiod\nand a positive effect\nof macrophyte status"),
                                        paste0(sprintf("%0.1f", 100*mean((slope_permanent - slope_temporary) < 0 & (slope_macrophytes - slope_nomacrophytes) > 0)), "% posterior\nprobability\n of a negative effect\nof hydroperiod\nand a positive effect\nof macrophyte status"),
                                        paste0(sprintf("%0.1f", 100*mean((slope_permanent - slope_temporary) > 0 & (slope_macrophytes - slope_nomacrophytes) < 0)), "% posterior\nprobability\n of a positive effect\nof hydroperiod\nand a negative effect\nof macrophyte status"),
                                        paste0(sprintf("%0.1f", 100*mean((slope_permanent - slope_temporary) < 0 & (slope_macrophytes - slope_nomacrophytes) < 0)), "% posterior\nprobability\n of a negative effect\nof hydroperiod\nand a negative effect\nof macrophyte status"))),
            aes(x = Hydroperiod, y = Macrophytes, label = label), size = 3) +
  scale_x_continuous("Posterior effect of hydroperiod") +
  scale_y_continuous("Posterior effect of macrophyte status") +
  coord_equal(xlim = c(-2.5,2.5), ylim = c(-2.5, 2.5)) +
  theme(panel.background = element_blank(),
        panel.grid = element_line(color = "grey93"),
        axis.title = element_text(face = "bold", size = 8),
        axis.text = element_text(size = 7))

ggsave("JointPD_plot.png", width = 14, height = 14, units = "cm", dpi = 600)

#-------------------------------------------------------------------------------

# plot hydroperiod and macrophytes

#-------------------------------------------------------------------------------

# Calculate mean percentage per pond group
data_mod <- datahydromacro %>%
  mutate(POND = ifelse(row_number() <= 11, "Group 1", POND)) %>%  # First 11 rows assigned to Group 1
  group_by(POND) %>%
  summarise(Hydroperiod = mean(Hydroperiod, na.rm = TRUE)) %>%  # Calculate mean percentage
  ungroup()

# Ensure factor levels are correctly ordered
data_mod$POND <- factor(data_mod$POND, levels = c("Group 1", datahydromacro$POND[12:nrow(data)]))

# Compute error metrics for first 11 ponds
data_modmac <- datahydromacro %>%
  mutate(POND = ifelse(row_number() <= 11, "Group 1", POND)) %>%
  group_by(POND)

data_modmac$POND <- factor(data_modmac$POND, levels = c("Group 1", datahydromacro$POND[12:nrow(data)]))

data_modmac <- data_modmac %>%
  mutate(macrophyte.percentage = aantal.maanden.met.macrophyten / pond.age * 100)  # Convert to percentage

# Prepare data for visualization
data_long <- data_mod %>%
  select(POND, Hydroperiod) %>%
  rename(Value = Hydroperiod) %>%
  mutate(Variable = "Hydroperiode") %>%
  bind_rows(
    data_modmac %>%
      select(POND, macrophyte.percentage) %>%
      rename(Value = macrophyte.percentage) %>%
      mutate(Variable = "Macrophyten")
  )

# Plot the data
ggplot(data_long) +
  geom_col(aes(x = POND, y = Value, group = Variable, fill = Variable), 
           position = position_dodge(width = 0.7), width = 0.6, color = "black") +  
  scale_fill_manual(values = c("Hydroperiode" = "#AEC6CF", "Macrophyten" = "#B5EAD7")) +
  labs(x = "Pond", y = "Percentage of time holding water", fill = "Variable", color = "Trend") +
  scale_y_continuous(
    name = "Percentage",
    sec.axis = sec_axis(~ (.), 
      name = "Months containing macrophytes",
      breaks = seq(0, 100, by = 25), # same spacing as left axis
      labels = function(x) round(x * 36 / 100))) +
  theme_classic(base_size = 14) +  
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = "#5A5A5A"),
    axis.text.y.left = element_text(color = "#AEC6CF"),
    axis.title.y.left = element_text(color = "#AEC6CF", face = "bold"),
    axis.text.y.right = element_text(color = "#B5EAD7"),
    axis.title.y.right = element_text(color = "#B5EAD7", face = "bold"),
    axis.line = element_line(color = "black", linewidth = 0.8),
    legend.position = "none")

# Export plot to PowerPoint
ggsave("Hydroperiod_macrophyteplot.png", width = 14, height = 14, units = "cm", dpi = 600)

