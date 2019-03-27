#Midterm

# ---- load dependencies ----
library(tidyverse)
library(sf)
library(car)
library(plotROC)
library(caret)
library(stargazer)
library(gridExtra)

# ---- set global params ----

setwd("/Users/johnmichaellasalle/Dropbox/Classes/CPLN\ 675\ Land\ Use\ and\ Environmental\ Monitoring/midterm/data")

options(scipen = 3)
epsg <- 26911
cell_size <- 100

map_theme <- function(base_size = 12) {
  theme(
    text = element_text( color = "black", family = "Helvetica"),
    plot.title = element_text(size = 14,colour = "black"),
    plot.subtitle=element_text(face="italic"),
    plot.caption=element_text(hjust=0),
    axis.ticks = element_blank(),
    panel.background = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()
  )
}

# ---- get initial data ----

#calgary boundary
clg_boundary <- st_read("https://data.calgary.ca/api/geospatial/7t9h-2z9s?method=export&format=GeoJSON") %>% 
  st_transform(., crs=epsg)
#st_write(clg_boundary, "clg_boundary.shp", delete_layer = TRUE)

#Edmonton Boundary
edm_boundary <- st_read("https://data.edmonton.ca/api/geospatial/m45c-6may?method=export&format=GeoJSON") %>% 
  st_transform(., crs=epsg)

#calgary fishnet
clg_fnet <- st_make_grid(clg_boundary, cellsize=cell_size) %>% 
  st_sf() %>%
  mutate(cell_id = row_number())
#st_write(clg_fnet, "clg_fnet.shp", delete_layer = TRUE)

clg_ctrnds <- clg_fnet %>% 
  st_centroid() %>%
  mutate(x = st_coordinates(.)[,1],
         y = st_coordinates(.)[,2]) %>%
  as.data.frame() %>%
  select(cell_id, x, y)

#Edmonton Fishnet
edm_fnet <- st_make_grid(edm_boundary, cellsize=cell_size) %>% 
  st_sf() %>%
  mutate(cell_id = row_number())
#st_write(edm_fnet, "edm_fnet.shp", delete_layer = TRUE)

edm_ctrnds <- edm_fnet %>% 
  st_centroid() %>%
  mutate(x = st_coordinates(.)[,1],
         y = st_coordinates(.)[,2]) %>%
  as.data.frame() %>%
  select(cell_id, x, y)

# ---- download and classify land cover ----
#calgary land cover
clg_lc <- st_read("https://data.calgary.ca/api/geospatial/as2i-6z3n?method=export&format=GeoJSON") %>% 
  st_transform(., crs=epsg) %>%
  mutate(permeable = ifelse(lc_cat %in% c("Building/Paved", 
                                          "Construction", 
                                          "PavedTrails", 
                                          "Roads/Raillines", 
                                          "Stream/Rivers", 
                                          "Reservoir"), 
                            0,1))

#Edmonton Land cover
edm_lc <- st_read("https://data.edmonton.ca/api/geospatial/5x9p-z4dg?method=export&format=GeoJSON") %>% 
  st_transform(., crs=epsg)

edm_lc <- edm_lc %>% 
  mutate(permeable = ifelse(stype1 %in% c("Established Residental Community", 
                                         "Aggregates and/or fill site", 
                                         "Anthropogenic Water Body", 
                                         "Commercial/Industrial Development", 
                                         "Natural Water Body", 
                                         "Residential Development",
                                         "Transportation Surface",
                                         "Building/Parking Complex",
                                         "Established Commercial / Industrial"), 
                            0, 1))

# ---- pull out categories of land cover ----

#calgary hydrology
clg_water <- clg_lc %>%
  filter(lc_cat %in% c("Stream/Rivers", "Reservoir"))
clg_wetlands <- clg_lc %>%
  filter(lc_cat %in% c("Natural Wetland", "Storm ponds/Modified Wetlands"))

#edmonton hydrology
edm_water <- edm_lc %>%
  filter(stype1 %in% c("Natural Water Body", "Anthropogenic Water Body"))
edm_wetlands <- edm_lc %>%
  filter(wetland1 %in% levels(wetland1) | wetland2 %in% levels(wetland2) | wetland3 %in% levels(wetland3))

#calgary permeable
clg_permeable <- clg_lc %>% 
  filter(permeable == 1)

#edmonton permeable
edm_permeable <- edm_lc %>% 
  filter(permeable == 1)

# ---- summarize land cover by fishnet cell ----
#calgary permeable
clg_perm_sum <- st_intersection(clg_permeable, clg_fnet) %>% 
  mutate(area = st_area(.)) %>% 
  group_by(cell_id) %>% 
  summarise(area_permeable = sum(area)) %>% 
  as.data.frame(.) %>%
  select(-geometry)

#calgary water
clg_water_sum <- st_intersection(clg_water, clg_fnet) %>% 
  mutate(area = st_area(.)) %>% 
  group_by(cell_id) %>% 
  summarise(area_water = sum(area)) %>% 
  as.data.frame(.) %>%
  select(-geometry)

#calgary wetlands
clg_wet_sum <- st_intersection(clg_wetlands, clg_fnet) %>% 
  mutate(area = st_area(.)) %>% 
  group_by(cell_id) %>% 
  summarise(area_wetlands = sum(area)) %>% 
  as.data.frame(.) %>%
  select(-geometry)

#join all calgary
clg_lc_sum <- full_join(clg_perm_sum, clg_water_sum) %>% 
  left_join(clg_wet_sum)


#edmonton permeable
edm_perm_sum <- st_intersection(edm_permeable, edm_fnet) %>% 
  mutate(area = st_area(.)) %>% 
  group_by(cell_id) %>% 
  summarise(area_permeable = sum(area)) %>% 
  as.data.frame(.) %>%
  select(-geometry)

#edmonton water
edm_water_sum <- st_intersection(edm_water, edm_fnet) %>% 
  mutate(area = st_area(.)) %>% 
  group_by(cell_id) %>% 
  summarise(area_water = sum(area)) %>% 
  as.data.frame(.) %>%
  select(-geometry)

#join all calgary
edm_lc_sum <- full_join(edm_perm_sum, edm_water_sum) %>% 
  left_join(edm_wet_sum)

#edmonton wetlands
edm_wet_sum <- st_intersection(edm_wetlands, edm_fnet) %>% 
  mutate(area = st_area(.)) %>% 
  group_by(cell_id) %>% 
  summarise(area_wetlands = sum(area)) %>% 
  as.data.frame(.) %>%
  select(-geometry)

#join all edmonton
edm_lc_sum <- full_join(edm_perm_sum, edm_water_sum) %>% 
  left_join(edm_wet_sum)

# ---- download and clip roads ----
al_rd <- st_read("alberta_roads.shp") %>% 
  st_transform(., crs=epsg)
clg_rd <- al_rd[clg_boundary,]
edm_rd <- al_rd[edm_boundary,]
rm(al_rd)

# ---- summarize roads by fishnet cell ----
clg_road_sum <- st_intersection(clg_rd, clg_fnet) %>% 
  mutate(length = st_length(.)) %>% 
  group_by(cell_id) %>% 
  summarise(road_length = sum(length)) %>% 
  as.data.frame(.) %>%
  select(-geometry)

edm_road_sum <- st_intersection(edm_rd, edm_fnet) %>% 
  mutate(length = st_length(.)) %>% 
  group_by(cell_id) %>% 
  summarise(road_length = sum(length)) %>% 
  as.data.frame(.) %>%
  select(-geometry)

# ---- read arc features ----
clg_innundation <- read.csv("clg_max_inundation.txt") %>%
  rename(cell_id = CELL_ID,
         flood = MAX) %>%
  select(cell_id, flood) %>%
  mutate(flood = factor(flood))

clg_str_dist <- read.csv("clg_max_strdist.txt") %>%
  rename(cell_id = CELL_ID,
         strdist = MAX) %>%
  select(cell_id, strdist)

edm_str_dist <- read.csv("edm_minmaxmean_strdist.txt") %>%
  rename(cell_id = CELL_ID,
         strdist = MAX) %>%
  select(cell_id, strdist)

clg_elev <- read.csv("clg_minmax_elevation.txt") %>%
  rename(cell_id = CELL_ID,
         elev = MIN) %>%
  select(cell_id, elev)

edm_elev <- read.csv("edm_minmaxmean_elevation.txt") %>%
  rename(cell_id = CELL_ID,
         elev = MIN) %>%
  select(cell_id, elev)

clg_slope <- read.csv("clg_minmaxmean_slope.txt") %>%
  rename(cell_id = CELL_ID,
         slope = MEAN) %>%
  select(cell_id, slope)

edm_slope <- read.csv("edm_minmaxmean_slope.txt") %>%
  rename(cell_id = CELL_ID,
         slope = MEAN) %>%
  select(cell_id, slope)

clg_str_order <- read.csv("clg_max_strorder.txt") %>%
  rename(cell_id = CELL_ID,
         strorder = MAX) %>%
  select(cell_id, strorder) %>%
  mutate(strorder = ifelse(is.na(strorder), 0, strorder))

edm_str_order <- read.csv("edm_max_strorder.txt") %>%
  rename(cell_id = CELL_ID,
         strorder = MAX) %>%
  select(cell_id, strorder) %>%
  mutate(strorder = ifelse(is.na(strorder), 0, strorder))

clg_flow_dist <- read.csv("clg_max_flowlen.txt") %>%
  rename(cell_id = CELL_ID,
         flow_dist = MAX) %>%
  select(cell_id, flow_dist)

edm_flow_dist <- read.csv("edm_max_flowlen.txt") %>%
  rename(cell_id = CELL_ID,
         flow_dist = MAX) %>%
  select(cell_id, flow_dist)

clg_flow_acc <- read.csv("clg_max_flowaccu.txt") %>%
  rename(cell_id = CELL_ID,
         flow_acc = MAX) %>%
  select(cell_id, flow_acc)

edm_flow_acc <- read.csv("edm_max_flowaccu.txt") %>%
  rename(cell_id = CELL_ID,
         flow_acc = MAX) %>%
  select(cell_id, flow_acc)

clg_arc_features <- left_join(clg_innundation, clg_elev) %>%
  left_join(clg_str_dist) %>%
  left_join(clg_slope) %>%
  left_join(clg_str_order) %>%
  left_join(clg_flow_dist) %>%
  left_join(clg_flow_acc) %>%
  mutate(strorder = ifelse(is.na(strorder), 0, strorder))

edm_arc_features <- left_join(edm_elev, edm_str_dist) %>%
  left_join(edm_slope) %>%
  left_join(edm_str_order) %>%
  left_join(edm_flow_dist) %>%
  left_join(edm_flow_acc) %>%
  mutate(strorder = ifelse(is.na(strorder), 0, strorder))

# ---- join features to fishnet ----

clg_fnet <- left_join(clg_fnet, clg_lc_sum) %>% 
  left_join(clg_road_sum) %>%
  inner_join(clg_arc_features) %>%
  left_join(clg_ctrnds) %>%
  mutate(area_permeable = ifelse(is.na(area_permeable), 0, area_permeable),
         area_water = ifelse(is.na(area_water), 0, area_water),
         area_wetlands = ifelse(is.na(area_wetlands), 0, area_wetlands),
         road_length = ifelse(is.na(road_length), 0, road_length),
         strorder = ifelse(is.na(strorder), 0, strorder))
st_write(clg_fnet, "clg_fnet_backup.shp", delete_layer = TRUE)

edm_fnet <- left_join(edm_fnet, edm_lc_sum) %>% 
  left_join(edm_road_sum) %>%
  inner_join(edm_arc_features) %>%
  left_join(edm_ctrnds) %>%
  mutate(area_permeable = ifelse(is.na(area_permeable), 0, area_permeable),
         area_water = ifelse(is.na(area_water), 0, area_water),
         area_wetlands = ifelse(is.na(area_wetlands), 0, area_wetlands),
         road_length = ifelse(is.na(road_length), 0, road_length),
         strorder = ifelse(is.na(strorder), 0, strorder)
         )
st_write(edm_fnet, "edm_fnet_backup.shp", delete_layer = TRUE)

# flow_acc_poly <- clg_fnet %>% select(flow_acc)
# st_write(clg_fnet, "flow_acc.shp", delete_layer = TRUE)

clg_fnet <- clg_fnet[clg_boundary,]
edm_fnet <- edm_fnet[edm_boundary,]
edm_fnet <- edm_fnet[edm_lc,]

# ---- re-engineering features ----

clg_fnet <- clg_fnet %>%
  mutate(permeable = as.factor(ifelse(area_permeable > 5000, 1, 0)),
         water = as.factor(ifelse(area_water > 0, 1, 0)),
         wetland = as.factor(ifelse(area_wetlands > 0, 1, 0)),
         elev_norm = elev-min(elev))

edm_fnet <- edm_fnet %>%
  mutate(permeable = as.factor(ifelse(area_permeable > 5000, 1, 0)),
         water = as.factor(ifelse(area_water > 0, 1, 0)),
         wetland = as.factor(ifelse(area_wetlands > 0, 1, 0)),
         elev_norm = elev-min(elev))


# ---- exploratory plots ----

ggplot() + 
  geom_point(data = clg_fnet, mapping = aes(x,y,col=log(flow_acc))) 

# ---- regression prep ----
#create train and test sets
set.seed(336)
trainIndex <- createDataPartition(clg_fnet$flood, p = .75,
                                  list = FALSE,
                                  times = 1)
clg_train <- clg_fnet[ trainIndex,]
clg_test  <- clg_fnet[-trainIndex,]

#create data frame for results
results_all <- clg_fnet %>%
  select(cell_id, flood, x, y)
results_test <- clg_test %>%
  select(cell_id, flood, x, y)

# ---- run regressions ----
formula <- flood ~ area_permeable + area_water + area_wetlands + elev_norm + road_length + strdist + slope + strorder + flow_dist + flow_acc
reg1 <- train(formula, data=clg_train,
              method="glm", 
              family="binomial",
              trControl = trainControl(method = "repeatedcv", 
                                       number = 100, 
                                       savePredictions = TRUE))

# ---- make predictions -----

#test set predictions
results_test$prob <- plogis(predict(reg1$finalModel, clg_test))
results_test <- results_test %>%
  mutate(pred = factor(ifelse(prob>=0.25, 1,0)),
         accuracy = factor(ifelse(pred == flood, 1, 0)),
         confusion = factor(ifelse(pred==1 & flood==1, "True Positive",
                                   ifelse(pred==0 & flood == 0, "True Negative",
                                          ifelse(pred==1 & flood == 0, "False Positive", "False Negative")
                                          )
                                   )
                            )
         )
#full dataset predictions
results_all$prob <- plogis(predict(reg1$finalModel, clg_fnet))
results_all <- results_all %>%
  mutate(pred = factor(ifelse(prob>=0.25, 1,0)),
         accuracy = factor(ifelse(pred == flood, 1, 0)),
         confusion = factor(ifelse(pred==1 & flood==1, "True Positive",
                                   ifelse(pred==0 & flood == 0, "True Negative",
                                          ifelse(pred==1 & flood == 0, "False Positive", "False Negative")
                                   )
         )
         )
  )
st_write(results_all, "clg_results.shp", delete_layer = TRUE)

results_edm <- edm_fnet %>%
  select(cell_id, x,y)

results_edm$prob <- plogis(predict(reg1$finalModel, edm_fnet))
results_edm <- results_edm %>%
  mutate(pred = factor(ifelse(prob>=0.25, 1,0)))
st_write(results_edm, "edm_results.shp", delete_layer = TRUE)

# ---- model analysis ----

confusionMatrix(results_test$pred, results_test$flood)

#cross validation plot
reg1$resample %>%
  select(-Resample, 
         -Kappa) %>%
  gather()%>%
  ggplot(., aes(value)) + 
  geom_histogram(bins=100, fill="#4d0054") +
  scale_x_continuous(limits = c(0, 1)) +
  labs(x="Accuracy Rate",
       y="Count") +
  theme_minimal()


#ROC curve
ggplot(results_test, 
         aes(d = as.numeric(flood), m = prob)) + 
  geom_roc(n.cuts = 50, labels = FALSE) + 
  style_roc(theme = theme_minimal)+
  geom_abline(slope = 1, intercept = 0, size = 1.5, color = 'grey') +
  labs(title = "ROC Curve for Test Set")



ggplot(results_test, mapping = aes(x=reg4, col=flood, fill=flood)) + geom_density(alpha=.25) + theme_minimal()

#confusion matrix plot
ggplot(results_test, mapping = aes(x=confusion, fill=confusion)) + geom_bar() + scale_fill_viridis_d() + theme_minimal()

#test set confusion map
ggplot() +
  geom_point(data=results_test, mapping = aes(x,y,col=confusion)) + scale_color_viridis_d() + theme_minimal()

#innundation prediction for calgary
ggplot(data=results_all, mapping = aes(x=x,y=y,col=pred)) + geom_point()+ scale_color_viridis_d(labels=c("Not Flooded", "Flooded")) + coord_equal() + theme_minimal()

#innundation prediction for edmonton
ggplot(data=results_edm, mapping = aes(x=x,y=y,col=pred)) + geom_point() + scale_color_viridis_d(labels=c("Not Flooded", "Flooded")) + coord_equal() + theme_minimal()


# ---- deliverables ----

# cross validation
reg1$resample %>%
  select(-Resample, 
         -Kappa) %>%
  gather()%>%
  ggplot(., aes(value)) + 
  geom_histogram(bins=100, fill="#4d0054") +
  scale_x_continuous(limits = c(0, 1)) +
  labs(title = "100 Subset Variation in Accuracy",
       x="Accuracy Rate",
       y="Count") +
  theme_minimal()
ggsave("cv.pdf",
       path = "/Users/johnmichaellasalle/Desktop", 
       units = "in", width = 3.5, height = 2.5,
       useDingbats = FALSE)

#ROC
ggplot(results_test, 
       aes(d = as.numeric(flood), m = prob)) + 
  geom_roc(n.cuts = 50, labels = FALSE, col="#4d0054") + 
  style_roc(theme = theme_minimal)+
  geom_abline(slope = 1, intercept = 0, size = 1.5, color = 'grey') +
  labs(title = "ROC Curve for Test Set")
ggsave("roc.pdf",
       path = "/Users/johnmichaellasalle/Desktop", 
       units = "in", width = 3.5, height = 2.5,
       useDingbats = FALSE)

#confusion matrix
ggplot(results_test, mapping = aes(x=confusion, fill=confusion)) + 
  geom_bar() + 
  scale_fill_viridis_d() + 
  labs(title = "Prediction Results",
       fill = "",
       y="Count",
       x="Result") +
  theme_minimal() +
  theme(legend.justification = c(0,1), 
        legend.position = c(0.1, 1))
ggsave("confusion.pdf",
       path = "/Users/johnmichaellasalle/Desktop", 
       units = "in", width = 4, height = 3.8,
       useDingbats = FALSE)
