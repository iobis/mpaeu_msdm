## code to prepare `sdmdata` dataset goes here
set.seed(2023)

sp_data <- robis::occurrence("Solea solea")
sp_data <- sp_data[sample(1:nrow(sp_data), 1000),]

env_data <- sdmpredictors::load_layers(c("BO22_tempmean_ss", "BO22_salinitymean_ss"))
env_data <- terra::rast(env_data)
# env_data <- terra::as.data.frame(env_data, xy = T, na.rm = F)

usethis::use_data(sp_data, env_data, overwrite = TRUE)


# my_data <- mp_prepare_data(sp_data, species_id = "teste", env_layers = env_data)
# 
# my_data <- mp_prepare_blocks(my_data)
# library(arrow);library(obissdm);library(terra);library(blockCV)
# 
# sp_data <- read_parquet("~/Research/mpa_europe/mpaeu_sdm/data/virtual_species/key=101/date=20230828/ftype=vsp/occurrences_high_rep1.parquet")
# 
# library(sdmpredictors)
# layers_1 <- c("BO22_tempmean_ss",
#               "BO22_salinitymin_ss",
#               "BO22_phosphatemean_ss",
#               "BO22_parmean")
# 
# env <- load_layers(layers_1,
#                    datadir = "data/raw/env_layers")
# 
# names(env) <- c("sst", "sal", "pho", "par")
# 
# env <- rast(env)
# 
# env <- mask(env, terra::app(env, prod))
# 
# sdm_data <- mp_prepare_data(sp_data, species_id = "vsp1", env_layers = env, quad_number = 30000)
# 
# sdm_data
# 
# sdm_data <- mp_prepare_blocks(sdm_data, block_types = "spatial_grid")
# 
# sdm_data
# 
# original_suit <- rast("~/Research/mpa_europe/mpaeu_sdm/data/virtual_species/key=101/date=20230828/ftype=vsp/suitability_key101.tif")
# 
# plot(original_suit)
# plot(res_pred)
