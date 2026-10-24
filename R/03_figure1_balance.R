#!/usr/bin/env Rscript
# =============================================================================
# 03_figure1_balance.R -- Figure 1 and Figures S2-S4
#
# (a,b) case/control locations before and after matching  -> Figure 1.pdf
# (c,d) propensity-score density before and after matching -> Figure 1.pdf
#       jitter plot and Love plot of covariate balance     -> Figures S2-S3
#       correlation matrix of the 12 greenness metrics      -> Figure S4
#
# Requires : R/00_setup.R, R/02_psm_table1.R (match.it, matchdatan),
#            china.geojson in the working directory (see data/README.md)
# Outputs  : Figure 1.pdf
# =============================================================================

if (!exists("matchdatan")) source(file.path("R", "02_psm_table1.R"))

#######Figure 1 S2-4######
#install.packages(c("sf","ggspatial"))
library(sf)
library(ggplot2)
library(ggspatial)
library(patchwork) 

#China GS（2024）0650
#install.packages(c("rnaturalearth", "rnaturalearthdata"))
library(rnaturalearth)
library(rnaturalearthdata)
china_prov <- st_read("china.geojson")
Location <- ggplot() +
  geom_sf(data = china_prov, fill = "lightblue", color = "darkgrey", linewidth = 0.3)+
  theme_void()

points_before <- st_as_sf(data, coords = c("lon_wgs", "lat_wgs"), crs = 4326)
points_after <- st_as_sf(matchdatan, coords = c("lon_wgs", "lat_wgs"), crs = 4326)
points_before$group <- factor(points_before$group_01, levels = c(0, 1), labels = c("Control", "Case"))
points_after$group <- factor(points_after$group_01, levels = c(0, 1), labels = c("Control", "Case"))

control_before <- subset(points_before, group == "Control")
case_before     <- subset(points_before, group == "Case")
control_after  <- subset(points_after,  group == "Control")
case_after     <- subset(points_after,  group == "Case")

p_before <- ggplot() +
  geom_sf(data = china_prov, fill = "#F5F5F5", color = "#2C3E50", linewidth = 0.4) +
  geom_sf(data = control_before, aes(fill = group), color = "black", shape = 21, size = 2.5, stroke = 0.7, alpha = 1) +
  geom_sf(data = case_before,     aes(fill = group), color = "black", shape = 21, size = 2.5, stroke = 0.7, alpha = 0.75) +
  scale_fill_manual(values = c("Control" = "#63baab", "Case" = "#D55E00"), name = "Group") +
  annotation_scale(location = "bl", width_hint = 0.2) +
  coord_sf() +
  theme_void() +
  labs(title = "Before Matching") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "plain"),
    legend.position = c(0.05, 0.05),
    legend.justification = c(0, 0),
    legend.background = element_blank(),
    legend.key = element_blank()
  )

p_after <- ggplot() +
  geom_sf(data = china_prov, fill = "#F5F5F5", color = "#2C3E50", linewidth = 0.4) +
  geom_sf(data = control_after, aes(fill = group), color = "black", shape = 21, size = 2.5, stroke = 0.7, alpha = 1) +
  geom_sf(data = case_after,     aes(fill = group), color = "black", shape = 21, size = 2.5, stroke = 0.7, alpha = 0.75) +
  scale_fill_manual(values = c("Control" = "#63baab", "Case" = "#D55E00"), name = "Group") +
  annotation_scale(location = "bl", width_hint = 0.2) +  
  coord_sf() +
  theme_void() +
  labs(title = "After Matching") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "plain"),
    legend.position = c(0.05, 0.05),
    legend.justification = c(0, 0),
    legend.background = element_blank(),
    legend.key = element_blank()
  )

##PD-PS Figure 1c,d
distance_scores <- match.it$distance
pre_data <- data
pre_data$distance <- distance_scores
dens0_pre <- density(pre_data$distance[pre_data$group_01 == 0])
dens1_pre <- density(pre_data$distance[pre_data$group_01 == 1])
back_pre <- rbind(
  data.frame(x = dens0_pre$x, y = dens0_pre$y, group = "0", stage = "Before Matching"),
  data.frame(x = dens1_pre$x, y = -dens1_pre$y, group = "1", stage = "Before Matching")
)
dens0 <- density(matchdatan$distance[matchdatan$group_01 == 0])
dens1 <- density(matchdatan$distance[matchdatan$group_01 == 1])
back_post <- rbind(
  data.frame(x = dens0$x, y = dens0$y, group = "0", stage = "After Matching"),
  data.frame(x = dens1$x, y = -dens1$y, group = "1", stage = "After Matching")
)
back_all <- rbind(back_pre, back_post)
back_all$group <- factor(back_all$group, levels = c("0", "1"), labels = c("Control", "Case"))
back_before <- subset(back_all, stage == "Before Matching")
back_after  <- subset(back_all, stage == "After Matching")

p_C <- ggplot(back_before, aes(x = x, y = y, fill = group)) +
  geom_area(alpha = 0.6, position = "identity") +
  geom_hline(yintercept = 0) +
  labs(x = "Propensity Scores", y = "Propensity Density", title = "Before Matching") +
  scale_y_continuous(labels = abs) +
  scale_fill_manual(values = c("Control" = "#63baab", "Case" = "#D55E00"), name = "Group") +
  coord_cartesian(xlim = c(0, 0.6)) +
  theme_minimal() +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 16, face = "plain"),
    legend.position = "none" 
  )

p_D <- ggplot(back_after, aes(x = x, y = y, fill = group)) +
  geom_area(alpha = 0.6, position = "identity") +
  geom_hline(yintercept = 0) +
  labs(x = "Propensity Scores", y = "Propensity Density", title = "After Matching") +
  scale_y_continuous(labels = abs) +
  scale_fill_manual(values = c("Control" = "#63baab", "Case" = "#D55E00"), name = "Group") +
  coord_cartesian(xlim = c(0, 0.6)) +
  theme_minimal() +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 16, face = "plain"),
    legend.position = "none"
  )

p_legend <- ggplot(back_before, aes(x = x, y = y, fill = group)) +
  geom_area() +
  scale_fill_manual(values = c("Control" = "#63baab", "Case" = "#D55E00"), name = "Group") +
  theme_void() +
  theme(legend.position = "bottom", legend.justification = "center")
legend_grob <- cowplot::get_legend(p_legend)

layout <- (p_before | p_after) / (p_C | p_D) / wrap_elements(legend_grob) +
  plot_layout(heights = c(2, 1, 0.05)) + 
  plot_annotation(tag_levels = list(c('a', 'b', 'c', 'd', ''))) & 
  theme(
    plot.tag = element_text(face = 'bold', size = 18) 
  )

pdf("Figure 1.pdf", width = 12, height = 12) 
print(layout)
dev.off()

###Figure S2-4
par(mfrow=c(1,1))
plot(match.it, type = "jitter", interactive = FALSE)#S3 10*8
cobalt::love.plot(match.it, binary = "std", abs = TRUE, thresholds = c(m = .1),#S2 10*8
                  var.order = "unadjusted", colors = c("#D55E00", "#63baab"),, title = "Covariate Balance (PSM)")

#install.packages("corrplot")
library(corrplot)
library(RColorBrewer)
selected_data <- data[, c("NDVI250m_mean","NDVI250m_max","NDVI500m_mean","NDVI500m_max","NDVI1000m_mean","NDVI1000m_max",
                          "EVI250m_mean","EVI250m_max","EVI500m_mean","EVI500m_max","EVI1000m_mean","EVI1000m_max")]
colnames(selected_data) <- c("NDVI mean 250m","NDVI max 250m","NDVI mean 500m","NDVI max 500m","NDVI mean 1000m","NDVI max 1000m",
                             "EVI mean 250m","EVI max 250m","EVI mean 500m","EVI max 500m","EVI mean 1000m","EVI max 1000m")
cor_matrix <- cor(selected_data, method = "pearson", use = "complete.obs")
print(cor_matrix)

cor_P <- cor.mtest(selected_data, conf.level = 0.95)
print(cor_P)

custom_col <- colorRampPalette(c( "#B2182B","#D6604D","#FDDBC7", "#4393C3","#2166AC"))(200)

corrplot.mixed(cor_matrix, upper = "circle", lower = "number", upper.col = custom_col, lower.col = custom_col,
              tl.pos = "lt",tl.col = "black",tl.cex = 0.8,cl.cex = 0.8,number.cex = 0.8,
              p.mat = cor_P$p,  sig.level = 0.05, insig = "blank")
