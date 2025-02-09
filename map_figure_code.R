#Author: Carter Boyd
#Date: 04/04/2024
#Purpose: Making maps for manuscript

library(tidyverse)
theme_set(theme_classic())
library(ggspatial)
library(sf)
library(raster)
library(patchwork)
library(rnaturalearth)
library(rnaturalearthdata)
library(tmap)
library(terra)
library(cowplot)
library(viridis)
library(scales)
library(paletteer)
library(ggthemes)
library(RColorBrewer)
library(ggimage)
# library(extrafont)
# font_import(path = "C:/Users/carterboyd/OneDrive - Virginia Tech/Desktop/Projects/2024_01_22_Code_for_Publication/Figure Making/Draft 3 Figures/arial.ttf", prompt = F)


# #Fig 1: Water Mask Figure
# tmap_mode("plot")
# tmap_options(check.and.fix = TRUE)
# 
# #Plot Panel A, S2 water mask over S2 true color image
# s2_tc <- rast("C:/Users/Carter/OneDrive - Virginia Tech/Desktop/Projects/2024_01_22_Code_for_Publication/Figure Making/True Color Images/s2_truecolor_new.tif")
# crs(s2_tc) <- '+init=EPSG:26914'
# s2_tc[s2_tc == 0] <- NA
# s2_tc_s <- stretch(s2_tc, minv = 0, maxv = 188, minq = 0.02, maxq = 0.98)
# s2_extent <- ext(s2_tc)
# s2_wm <- rast("C:/Users/Carter/OneDrive - Virginia Tech/Desktop/Projects/2023_05_04_Platte_Take_5_Labeling/RivWidth/NK14/product/NK14_chan_class.img")
# crs(s2_wm) <- '+init=EPSG:26914'
# s2_extent <- align(s2_extent, s2_wm, snap = "in")
# s2_wm <- crop(s2_wm, s2_extent)
# s2_wm[s2_wm == 0] <- NA
# # s2_wm[is.na(s2_tc_s)] <- NA
#   
# tm_shape(s2_tc_s, raster.downsample = F)+
#   tm_rgb(r=1, g=2, b=3, max.value = 188)+
#   tm_shape(s2_wm, raster.downsample = F)+
#   tm_raster(style = "cat", labels = "Water Mask", palette = "#9a031e", alpha = 0.7, legend.show = F)+
#   tm_scale_bar(breaks = c(0,0.5),
#                bg.alpha = 0,
#                color.dark = "white",
#                text.color = "white",
#                bg.color = "white", position = c(0.65, 0.04), text.size = 1)+
#   # tm_compass(text.color = "white",
#   #            color.dark = "white",
#   #            color.light = "black",
#   #            lwd=2,
#   #   #bg.color = "white", 
#   #            position = c(0.9, 0.05))+
#   tm_layout(frame = F, inner.margins = 0)
#   
# 
# #Plot Panel B, NAIP water mask over NAIP true color image
# naip_tc <- rast("C:/Users/Carter/OneDrive - Virginia Tech/Desktop/Projects/2024_01_22_Code_for_Publication/Figure Making/True Color Images/naip_truecolor_smaller.tif")
# crs(naip_tc) <- '+init=EPSG:26913'
# naip_tc[naip_tc == 0] <- NA
# naip_extent <- ext(naip_tc)
# naip_wm <- rast("C:/Users/Carter/OneDrive - Virginia Tech/Desktop/Projects/2024_01_22_Code_for_Publication/Figure Making/Water Masks/SVC_wm.tif")
# crs(naip_wm) <- '+init=EPSG:26913'
# naip_wm <- crop(naip_wm, naip_extent)
# naip_wm[naip_wm == 0] <- NA
# 
# tm_shape(naip_tc, raster.downsample = F)+
#   tm_rgb(r=1, g=2, b=3)+
#   tm_shape(naip_wm, raster.downsample = F)+
#   tm_raster(style = "cat", labels = "Water Mask", palette = "#e36414", alpha = 0.7, legend.show = F)+
#   tm_scale_bar(breaks = c(0,0.25),
#                bg.alpha = 0,
#                color.dark = "white",
#                text.color = "white",
#                bg.color = "white", position = c(0.6, 0.04), text.size = 1)+
#   # tm_compass(text.color = "white",
#   #            color.dark = "white",
#   #            color.light = "black",
#   #            lwd=2,
#   #            #bg.color = "white", 
#   #            position = c(0.9, 0.05))+
#   tm_layout(frame = F, inner.margins = 0)


#Fig 1: Width map
# coastline <- ne_coastline(scale = "medium", returnclass = "sf")
countries <- ne_countries(country = c("united states of america", "canada", "mexico"),scale = "medium", returnclass = "sf")

mississippi_geom <- st_read("C:/Users/Carter/OneDrive - Virginia Tech/Desktop/Projects/2024_01_22_Code_for_Publication/Figure Making/Basin Shapefiles/mississippi_river_basin.shp")
platte_geom <- st_read("C:/Users/Carter/OneDrive - Virginia Tech/Desktop/Projects/2024_01_22_Code_for_Publication/Figure Making/Basin Shapefiles/platte_geom.shp")
svc_geom <- st_read("C:/Users/Carter/OneDrive - Virginia Tech/Desktop/Projects/2024_01_22_Code_for_Publication/Figure Making/Basin Shapefiles/stvrain_geom.shp")
#SVC geom with SSVC basin cut out, to use for SVC and SSVC scale mapping
svc_geom_union <- st_read("C:/Users/Carter/OneDrive - Virginia Tech/Desktop/Projects/2024_01_22_Code_for_Publication/Figure Making/Basin Shapefiles/stvrain_geom_union.shp")
ssvc_geom <- st_read("C:/Users/Carter/OneDrive - Virginia Tech/Desktop/Projects/2024_01_22_Code_for_Publication/Figure Making/Basin Shapefiles/ssvc_geom.shp")

mississippi_lakes <- st_read("C:/Users/Carter/OneDrive - Virginia Tech/Desktop/Projects/2024_01_22_Code_for_Publication/Figure Making/Lakes Shapefile/mississippi_lakes_gt5km2.shp")
platte_lakes <- st_read("C:/Users/Carter/OneDrive - Virginia Tech/Desktop/Projects/2024_01_22_Code_for_Publication/Figure Making/Lakes Shapefile/platte_lakes_gt2km2.shp")
svc_lakes <- st_read("C:/Users/Carter/OneDrive - Virginia Tech/Desktop/Projects/2024_01_22_Code_for_Publication/Figure Making/Lakes Shapefile/svc_lakes.shp")
ssvc_lakes <- st_read("C:/Users/Carter/OneDrive - Virginia Tech/Desktop/Projects/2024_01_22_Code_for_Publication/Figure Making/Lakes Shapefile/ssvc_lakes.shp")

mississippi_lakes_wm <- st_read("C:/Users/Carter/OneDrive - Virginia Tech/Desktop/Projects/2024_01_22_Code_for_Publication/Data/lakes from wm/mississippi_lakes_wm.shp")
platte_lakes_wm <- st_read("C:/Users/Carter/OneDrive - Virginia Tech/Desktop/Projects/2024_01_22_Code_for_Publication/Data/lakes from wm/platte_lakes_wm.shp")
svc_lakes_wm <- st_read("C:/Users/Carter/OneDrive - Virginia Tech/Desktop/Projects/2024_01_22_Code_for_Publication/Data/lakes from wm/svc_lakes_wm_new.shp")
ssvc_lakes_wm <- st_read("C:/Users/Carter/OneDrive - Virginia Tech/Desktop/Projects/2024_01_22_Code_for_Publication/Data/lakes from wm/ssvc_lakes_wm.shp")

all_centerlines <- st_read("C:/Users/Carter/OneDrive - Virginia Tech/Desktop/Projects/2024_01_22_Code_for_Publication/Figure Making/Centerlines/all_centerlines.shp")
centerlines1 <- dplyr::filter(all_centerlines, SO == 1)
centerlines2 <- dplyr::filter(all_centerlines, SO == 2)
centerlines3 <- dplyr::filter(all_centerlines, SO == 3)
centerlines4 <- dplyr::filter(all_centerlines, SO == 4)
centerlines5 <- dplyr::filter(all_centerlines, SO == 5)
centerlines6 <- dplyr::filter(all_centerlines, SO == 6)
centerlines7 <- dplyr::filter(all_centerlines, SO == 7)
centerlines8 <- dplyr::filter(all_centerlines, SO == 8)
centerlines9 <- dplyr::filter(all_centerlines, SO == 9)
centerlines10 <- dplyr::filter(all_centerlines, SO == 10)
centerlines11 <- dplyr::filter(all_centerlines, SO == 11)
centerlines12 <- dplyr::filter(all_centerlines, SO == 12)
centerlines13 <- dplyr::filter(all_centerlines, SO == 13)

all_centerlines_platte <- st_intersection(all_centerlines, platte_geom)
centerlines1_platte <- dplyr::filter(all_centerlines_platte, SO == 1)
centerlines2_platte <- dplyr::filter(all_centerlines_platte, SO == 2)
centerlines3_platte <- dplyr::filter(all_centerlines_platte, SO == 3)
centerlines4_platte <- dplyr::filter(all_centerlines_platte, SO == 4)
centerlines5_platte <- dplyr::filter(all_centerlines_platte, SO == 5)
centerlines6_platte <- dplyr::filter(all_centerlines_platte, SO == 6)
centerlines7_platte <- dplyr::filter(all_centerlines_platte, SO == 7)
centerlines8_platte <- dplyr::filter(all_centerlines_platte, SO == 8)
centerlines9_platte <- dplyr::filter(all_centerlines_platte, SO == 9)
centerlines10_platte <- dplyr::filter(all_centerlines_platte, SO == 10)
centerlines11_platte <- dplyr::filter(all_centerlines_platte, SO == 11)
centerlines12_platte <- dplyr::filter(all_centerlines_platte, SO == 12)
centerlines13_platte <- dplyr::filter(all_centerlines_platte, SO == 13)

all_centerlines_svc <- st_intersection(all_centerlines_platte, svc_geom)
centerlines1_svc <- dplyr::filter(all_centerlines_svc, SO == 1)
centerlines2_svc <- dplyr::filter(all_centerlines_svc, SO == 2)
centerlines3_svc <- dplyr::filter(all_centerlines_svc, SO == 3)
centerlines4_svc <- dplyr::filter(all_centerlines_svc, SO == 4)
centerlines5_svc <- dplyr::filter(all_centerlines_svc, SO == 5)
centerlines6_svc <- dplyr::filter(all_centerlines_svc, SO == 6)
centerlines7_svc <- dplyr::filter(all_centerlines_svc, SO == 7)
centerlines8_svc <- dplyr::filter(all_centerlines_svc, SO == 8)
centerlines9_svc <- dplyr::filter(all_centerlines_svc, SO == 9)
centerlines10_svc <- dplyr::filter(all_centerlines_svc, SO == 10)
centerlines11_svc <- dplyr::filter(all_centerlines_svc, SO == 11)
centerlines12_svc <- dplyr::filter(all_centerlines_svc, SO == 12)
centerlines13_svc <- dplyr::filter(all_centerlines_svc, SO == 13)

all_centerlines_ssvc <- st_intersection(all_centerlines_svc, ssvc_geom)
centerlines1_ssvc <- dplyr::filter(all_centerlines_ssvc, SO == 1)
centerlines2_ssvc <- dplyr::filter(all_centerlines_ssvc, SO == 2)
centerlines3_ssvc <- dplyr::filter(all_centerlines_ssvc, SO == 3)
centerlines4_ssvc <- dplyr::filter(all_centerlines_ssvc, SO == 4)
centerlines5_ssvc <- dplyr::filter(all_centerlines_ssvc, SO == 5)
centerlines6_ssvc <- dplyr::filter(all_centerlines_ssvc, SO == 6)
centerlines7_ssvc <- dplyr::filter(all_centerlines_ssvc, SO == 7)
centerlines8_ssvc <- dplyr::filter(all_centerlines_ssvc, SO == 8)
centerlines9_ssvc <- dplyr::filter(all_centerlines_ssvc, SO == 9)
centerlines10_ssvc <- dplyr::filter(all_centerlines_ssvc, SO == 10)
centerlines11_ssvc <- dplyr::filter(all_centerlines_ssvc, SO == 11)
centerlines12_ssvc <- dplyr::filter(all_centerlines_ssvc, SO == 12)
centerlines13_ssvc <- dplyr::filter(all_centerlines_ssvc, SO == 13)



my_breaks <- c(0,10,200,4000)
bottom_color <- "grey65"
bg_color <- "grey30"
fg_color <- "grey5"
nested_color <- "grey45"
linewidth_miss <- 0.35
linewidth <- 0.2
centerline_palette <- "Spectral" #"YlGnBu"

#####Version 1
#mississippi basin (done)
mississippi_width_map <- ggplot()+
  # geom_sf(data = coastline, color = bottom_color)+
  geom_sf(data = countries, color = bottom_color, fill = bottom_color)+
  geom_sf(data = mississippi_geom, color = fg_color, fill = fg_color)+
  geom_sf(data = platte_geom, color = nested_color, fill = nested_color)+
  # geom_sf(data = svc_geom, color = "white", fill = fg_color)+
  # geom_sf(data = ssvc_geom, color = "white", fill = fg_color)+
  # scale_color_viridis_c(direction = -1, breaks = my_breaks, labels = my_breaks,
  #                       trans = scales::pseudo_log_trans())+
  scale_color_distiller(name = "Width (m)\n",
                        palette = centerline_palette,
                        direction = -1,
                        breaks = my_breaks,
                        labels = my_breaks,
                        trans = scales::pseudo_log_trans(),
                        limits = c(0,4003))+
  # scale_color_gradientn(name = "Width (m)\n",
  #   colors = c("00ffc8", "00f0d0", "00e2d8", "00d3e0", "00c5e7", "00b6ef", "00a8f7", "0099ff"),
  #   breaks = my_breaks,
  #   guide = "colourbar",
  #   limits = c(0,4003),
  #   labels = my_breaks,
  #   transform = "pseudo_log"
  # )+
  geom_sf(data = centerlines1, aes(color = width_m), linewidth = linewidth)+
  geom_sf(data = centerlines2, aes(color = width_m), linewidth = linewidth)+
  geom_sf(data = centerlines3, aes(color = width_m), linewidth = linewidth)+
  geom_sf(data = centerlines4, aes(color = width_m), linewidth = linewidth)+
  geom_sf(data = centerlines5, aes(color = width_m), linewidth = linewidth)+
  geom_sf(data = centerlines6, aes(color = width_m), linewidth = linewidth)+
  geom_sf(data = centerlines7, aes(color = width_m), linewidth = linewidth)+
  geom_sf(data = centerlines8, aes(color = width_m), linewidth = linewidth)+
  geom_sf(data = centerlines9, aes(color = width_m), linewidth = linewidth)+
  geom_sf(data = centerlines10, aes(color = width_m), linewidth = linewidth)+
  geom_sf(data = centerlines11, aes(color = width_m), linewidth = linewidth)+
  geom_sf(data = centerlines12, aes(color = width_m), linewidth = linewidth)+
  geom_sf(data = centerlines13, aes(color = width_m), linewidth = linewidth)+
  scale_fill_manual(values = c("Lakes and River Masks"="#99ccff"), name = element_blank())+
  # scale_color_manual(values = c("Lakes"="#7BD9F6"))+
  geom_sf(data = mississippi_lakes_wm, aes(fill = "Lakes and River Masks"), color = "transparent", show.legend = T)+
  
  labs(color = "Width (m)\n")+
  scale_x_continuous(limits = c(-119,-73.77))+
  scale_y_continuous(limits = c(28,51 #51
  ))+
  geom_segment(aes(x=-109.5,y=37.51,xend=-119,yend=28),
               color = "white", linewidth = 0.4)+
  annotation_scale(location = "bl", width_hint = 0.3, style = "bar",
                   pad_x = unit(0.06, "in"),
                   pad_y = unit(0.06, "in"),
                   text_cex = (5/6),
                   height = unit(0.1, "in"),
                   bar_cols = "white",
                   line_col = "white",
                   text_col = "white")+
  geom_rect(aes(xmin = -109.5, xmax = -95.5, ymin = 37.51, ymax = 44.49), color = "white", fill = NA, linewidth = 2*linewidth)+
  coord_sf(expand = F)+
  theme_void()+
  theme(text = element_text(family = "sans", size = 8),
        # axis.text = element_text(color = "white"),
        # axis.line = element_line(color="white"),
        # axis.ticks = element_line(color = "white"),
        legend.position = c(0.267,0.357 #0.85,0.80
        ),
        legend.background = element_blank(), #element_rect(color = "black", fill = "white"),
        legend.text = element_text(color="black", size = 8),
        legend.title = element_text(color = "black", size = 8,
                                    margin = margin(b=0,unit='in'),
                                    vjust = -1),
        # legend.key = element_rect(margin(t=0,unit='in')),
        legend.key.size = unit(0.09, 'in'),
        legend.margin = margin(t=-0.05,b=0.1,l=0.05,r=0.05,unit = 'in'),
        # legend.spacing.x = unit(0, "in"),
        legend.spacing.y = unit(0.02, "in"),
        panel.background = element_rect(fill=bg_color, color = NA_character_),
        plot.background = element_rect(fill='transparent', color=NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # legend.box.background = element_rect(fill='transparent')
  # )+
  # geom_image(
  #   x = -75.1, y = 46.1, size = 0.35, aes(image = "C:/Users/Carter/OneDrive - Virginia Tech/Desktop/Projects/2024_01_22_Code_for_Publication/Figure Making/Draft 7 Figures/Fig 1 Maps and Hists/mississippi_hist.png")
  # )+
  # geom_image(
  #   x = -73.82, y = 32.49, size = 0.35, aes(image = "C:/Users/Carter/OneDrive - Virginia Tech/Desktop/Projects/2024_01_22_Code_for_Publication/Figure Making/Draft 7 Figures/Fig 1 Maps and Hists/mississippi_watermask_map.png")
  )+
  annotate("text", x = -81, y = 49.5, label = "N = 1,384,979", size = 3)
# guides(color = guide_legend(override.aes = list(size = 0.5)))
# miss_hist <- system.file("extdata","C:/Users/Carter/OneDrive - Virginia Tech/Desktop/Projects/2024_01_22_Code_for_Publication/Figure Making/Draft 7 Figures/Fig 1 Maps and Hists/mississippi_hist.png", package = "cowplot")
# 
# ggdraw(mississippi_width_map)+
#   draw_image(miss_hist, x = 0.9, y = 0.9, scale = 1)
ggsave(filename = "C:/Users/Carter/OneDrive - Virginia Tech/Desktop/Projects/2024_01_22_Code_for_Publication/Figure Making/Draft 7 Figures/Fig 1 Maps and Hists/.pdf/mississippi_width_map.pdf",
       plot = mississippi_width_map,
       width = 2.65,
       height = 1.75,
       units = "in",
       dpi = 1200,
       bg = "white")



#platte basin (done)
platte_width_map <- ggplot()+
  # geom_sf(data = coastline, color = bottom_color)+
  geom_sf(data = countries, color = bg_color, fill = bg_color)+
  geom_sf(data = mississippi_geom, color = bottom_color, fill = bottom_color)+
  geom_sf(data = platte_geom, color = fg_color, fill = fg_color)+
  geom_sf(data = svc_geom, color = nested_color, fill = nested_color)+
  # geom_sf(data = ssvc_geom, color = "white", fill = bg_color)+
  # scale_color_viridis_c(direction = -1, breaks = my_breaks, labels = my_breaks,
  #                       trans = scales::pseudo_log_trans())+
  scale_color_distiller(name = "Width (m)\n",
                        palette = centerline_palette,
                        direction = -1,
                        breaks = my_breaks,
                        labels = my_breaks,
                        trans = scales::pseudo_log_trans(),
                        limits = c(0,4003))+
  geom_sf(data = centerlines1_platte, aes(color = width_m), linewidth = linewidth)+
  geom_sf(data = centerlines2_platte, aes(color = width_m), linewidth = linewidth)+
  geom_sf(data = centerlines3_platte, aes(color = width_m), linewidth = linewidth)+
  geom_sf(data = centerlines4_platte, aes(color = width_m), linewidth = linewidth)+
  geom_sf(data = centerlines5_platte, aes(color = width_m), linewidth = linewidth)+
  geom_sf(data = centerlines6_platte, aes(color = width_m), linewidth = linewidth)+
  geom_sf(data = centerlines7_platte, aes(color = width_m), linewidth = linewidth)+
  geom_sf(data = centerlines8_platte, aes(color = width_m), linewidth = linewidth)+
  geom_sf(data = centerlines9_platte, aes(color = width_m), linewidth = linewidth)+
  geom_sf(data = centerlines10_platte, aes(color = width_m), linewidth = linewidth)+
  geom_sf(data = centerlines11_platte, aes(color = width_m), linewidth = linewidth)+
  geom_sf(data = centerlines12_platte, aes(color = width_m), linewidth = linewidth)+
  geom_sf(data = centerlines13_platte, aes(color = width_m), linewidth = linewidth)+
  geom_sf(data = platte_lakes_wm, color = "transparent", fill = "#99ccff")+
  labs(color = "Width (m)\n")+
  scale_x_continuous(limits = c(-109.5,-95.5))+
  scale_y_continuous(limits = c(37.51,44.49 #37.51,44.49
  ))+
  geom_segment(aes(x=-105.844,y=39.8,xend=-109.5,yend=37.51),
               color = "white", linewidth = 0.4)+
  geom_rect(aes(xmin = -105.844, xmax = -104.756, ymin = 39.8, ymax = 40.35), color = "white", fill = NA, linewidth = 2*linewidth)+
  coord_sf(expand = F)+
  annotation_scale(location = "br", width_hint = 0.3, style = "bar",
                   pad_x = unit(0.05, "in"),
                   pad_y = unit(0.05, "in"),
                   text_cex = (5/6),
                   height = unit(0.1, "in"),
                   bar_cols = "black",
                   line_col = "black",
                   text_col = "black")+
  theme_void()+
  theme(text = element_text(family = "sans", size = 12),
        legend.position = "none",
        legend.background = element_rect(color = NA, fill = NA),
        legend.text = element_text(color="white"),
        legend.title = element_text(color = "white"),
        panel.background = element_rect(fill='transparent', color = NA_character_),
        plot.background = element_rect(fill='transparent', color=NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.box.background = element_rect(fill='transparent'))+
  annotate("text", x = -97.5, y = 44, label = "N = 378,179", size = 3)
# geom_image(
#   x = -97, y = 39, size = 0.3, aes(image = "C:/Users/Carter/OneDrive - Virginia Tech/Desktop/Projects/2024_01_22_Code_for_Publication/Figure Making/Draft 7 Figures/Fig 1 Maps and Hists/platte_hist.png")
# )

ggsave(filename = "C:/Users/Carter/OneDrive - Virginia Tech/Desktop/Projects/2024_01_22_Code_for_Publication/Figure Making/Draft 7 Figures/Fig 1 Maps and Hists/.pdf/platte_width_map.pdf",
       plot = platte_width_map,
       width = 2.65,
       height = 1.75,
       units = "in",
       dpi = 1200,
       bg = "white")



#st vrain basin (done)
svc_width_map <- ggplot()+
  # geom_sf(data = coastline, color = bottom_color)+
  geom_sf(data = countries, color = bg_color, fill = bg_color)+
  # geom_sf(data = mississippi_geom, color = "white", fill = bg_color)+
  geom_sf(data = platte_geom, color = bottom_color, fill = bottom_color)+
  geom_sf(data = svc_geom_union, color = fg_color, fill = fg_color)+
  geom_sf(data = ssvc_geom, color = nested_color, fill = nested_color)+
  # scale_color_viridis_c(direction = -1, breaks = my_breaks, labels = my_breaks,
  #                       trans = scales::pseudo_log_trans())+
  scale_color_distiller(name = "Width (m)\n",
                        palette = centerline_palette,
                        direction = -1,
                        breaks = my_breaks,
                        labels = my_breaks,
                        trans = scales::pseudo_log_trans(),
                        limits = c(0,4003))+
  geom_sf(data = centerlines1_svc, aes(color = width_m), linewidth = linewidth)+
  geom_sf(data = centerlines2_svc, aes(color = width_m), linewidth = linewidth)+
  geom_sf(data = centerlines3_svc, aes(color = width_m), linewidth = linewidth)+
  geom_sf(data = centerlines4_svc, aes(color = width_m), linewidth = linewidth)+
  geom_sf(data = centerlines5_svc, aes(color = width_m), linewidth = linewidth)+
  geom_sf(data = centerlines6_svc, aes(color = width_m), linewidth = linewidth)+
  geom_sf(data = centerlines7_svc, aes(color = width_m), linewidth = linewidth)+
  geom_sf(data = centerlines8_svc, aes(color = width_m), linewidth = linewidth)+
  geom_sf(data = centerlines9_svc, aes(color = width_m), linewidth = linewidth)+
  geom_sf(data = centerlines10_svc, aes(color = width_m), linewidth = linewidth)+
  geom_sf(data = centerlines11_svc, aes(color = width_m), linewidth = linewidth)+
  geom_sf(data = centerlines12_svc, aes(color = width_m), linewidth = linewidth)+
  geom_sf(data = centerlines13_svc, aes(color = width_m), linewidth = linewidth)+
  geom_sf(data = svc_lakes_wm, color = "transparent", fill = "#7BD9F6")+
  labs(color = "Width (m)\n")+
  scale_x_continuous(limits = c(-105.844,-104.756), breaks = c(-105.7,-105.5,-105.3,-105.1,-104.9))+
  scale_y_continuous(limits = c(39.8,40.35 #51
  ))+
  geom_segment(aes(x=-105.66,y=40.0448,xend=-105.844,yend=39.8),
               color = "white", linewidth = 0.4)+
  geom_rect(aes(xmin = -105.66, xmax = -105.57, ymin = 40.0448, ymax = 40.0903), color = "white", fill = NA, linewidth = 2*linewidth)+
  coord_sf(expand = F)+
  annotation_scale(location = "br", width_hint = 0.25, style = "bar",
                   pad_x = unit(0.05, "in"),
                   pad_y = unit(0.05, "in"),
                   text_cex = (5/6),
                   height = unit(0.1, "in"),
                   bar_cols = "black",
                   line_col = "black",
                   text_col = "black")+
  theme_void()+
  theme(text = element_text(family = "sans", size = 12),
        legend.position = "none",
        legend.background = element_rect(color = NA, fill = NA),
        legend.text = element_text(color="white"),
        legend.title = element_text(color = "white"),
        panel.background = element_rect(fill='transparent', color = NA_character_),
        plot.background = element_rect(fill='transparent', color=NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.box.background = element_rect(fill='transparent'))+
  annotate("text", x = -104.91, y = 40.315, label = "N = 436,208", size = 3)

ggsave(filename = "C:/Users/Carter/OneDrive - Virginia Tech/Desktop/Projects/2024_01_22_Code_for_Publication/Figure Making/Draft 7 Figures/Fig 1 Maps and Hists/.pdf/svc_width_map.pdf",
       plot = svc_width_map,
       width = 2.65,
       height = 1.75,
       units = "in",
       dpi = 1200,
       bg = "white")



#south st vrain basin (done)
ssvc_width_map <- ggplot()+
  # geom_sf(data = coastline, color = bottom_color)+
  geom_sf(data = countries, color = bg_color, fill = bg_color)+
  # geom_sf(data = mississippi_geom, color = "white", fill = bg_color)+
  # geom_sf(data = platte_geom, color = "white", fill = bg_color)+
  geom_sf(data = svc_geom_union, color = bottom_color, fill = bottom_color)+
  geom_sf(data = ssvc_geom, color = fg_color, fill = fg_color)+
  # scale_color_viridis_c(direction = -1, breaks = my_breaks, labels = my_breaks,
  #                       trans = scales::pseudo_log_trans())+
  scale_color_distiller(name = "Width (m)\n",
                        palette = centerline_palette,
                        direction = -1,
                        breaks = my_breaks,
                        labels = my_breaks,
                        trans = scales::pseudo_log_trans(),
                        limits = c(0,4003))+
  geom_sf(data = centerlines1_ssvc, aes(color = width_m), linewidth = linewidth)+
  geom_sf(data = centerlines2_ssvc, aes(color = width_m), linewidth = linewidth)+
  geom_sf(data = centerlines3_ssvc, aes(color = width_m), linewidth = linewidth)+
  geom_sf(data = centerlines4_ssvc, aes(color = width_m), linewidth = linewidth)+
  geom_sf(data = centerlines5_ssvc, aes(color = width_m), linewidth = linewidth)+
  geom_sf(data = centerlines6_ssvc, aes(color = width_m), linewidth = linewidth)+
  geom_sf(data = centerlines7_ssvc, aes(color = width_m), linewidth = linewidth)+
  geom_sf(data = centerlines8_ssvc, aes(color = width_m), linewidth = linewidth)+
  geom_sf(data = centerlines9_ssvc, aes(color = width_m), linewidth = linewidth)+
  geom_sf(data = centerlines10_ssvc, aes(color = width_m), linewidth = linewidth)+
  geom_sf(data = centerlines11_ssvc, aes(color = width_m), linewidth = linewidth)+
  geom_sf(data = centerlines12_ssvc, aes(color = width_m), linewidth = linewidth)+
  geom_sf(data = centerlines13_ssvc, aes(color = width_m), linewidth = linewidth)+
  geom_sf(data = ssvc_lakes_wm, color = "transparent", fill = "#7BD9F6")+
  labs(color = "Width (m)\n")+
  scale_x_continuous(limits = c(-105.66,-105.57)#, breaks = c(-105.7,-105.5,-105.3,-105.1,-104.9)
  )+
  scale_y_continuous(limits = c(40.0448,40.0903 #51
  ))+
  # geom_rect(aes(xmin = -105.66, xmax = -105.57, ymin = 40.045, ymax = 40.09), color = "white", fill = NA, linewidth = 1)+
  coord_sf(expand = F)+
  annotation_scale(location = "br", width_hint = 0.3, style = "bar",
                   pad_x = unit(0.06, "in"),
                   pad_y = unit(0.06, "in"),
                   text_cex = (5/6),
                   height = unit(0.1, "in"),
                   bar_cols = "black",
                   line_col = "black",
                   text_col = "black")+
  theme_void()+
  theme(text = element_text(family = "sans", size = 8),
        legend.position = "none",
        legend.background = element_rect(color = NA, fill = NA),
        legend.text = element_text(color="white"),
        legend.title = element_text(color = "white"),
        panel.background = element_rect(fill='transparent', color = NA_character_),
        plot.background = element_rect(fill='transparent', color=NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.box.background = element_rect(fill='transparent'))+
  annotate("text", x = -105.581, y = 40.0872, label = "N = 3,443", size = 3)


ggsave(filename = "C:/Users/Carter/OneDrive - Virginia Tech/Desktop/Projects/2024_01_22_Code_for_Publication/Figure Making/Draft 7 Figures/Fig 1 Maps and Hists/.pdf/ssvc_width_map.pdf",
       plot = ssvc_width_map,
       width = 2.65,
       height = 1.75,
       units = "in",
       dpi = 1200,
       bg = "white")

#####Version 2
#mississippi basin (done)
mississippi_width_map_v2 <- ggplot()+
  # geom_sf(data = coastline, color = bottom_color)+
  geom_sf(data = countries, color = bottom_color, fill = bottom_color)+
  geom_sf(data = mississippi_geom, color = fg_color, fill = fg_color)+
  geom_sf(data = platte_geom, color = nested_color, fill = nested_color)+
  # geom_sf(data = svc_geom, color = "white", fill = fg_color)+
  # geom_sf(data = ssvc_geom, color = "white", fill = fg_color)+
  # scale_color_viridis_c(direction = -1, breaks = my_breaks, labels = my_breaks,
  #                       trans = scales::pseudo_log_trans())+
  scale_color_distiller(name = "Width (m)\n",
                        palette = centerline_palette,
                        direction = -1,
                        breaks = my_breaks,
                        labels = my_breaks,
                        trans = scales::pseudo_log_trans(),
                        limits = c(0,4003))+
  # scale_color_gradientn(name = "Width (m)\n",
  #   colors = c("00ffc8", "00f0d0", "00e2d8", "00d3e0", "00c5e7", "00b6ef", "00a8f7", "0099ff"),
  #   breaks = my_breaks,
  #   guide = "colourbar",
  #   limits = c(0,4003),
  #   labels = my_breaks,
  #   transform = "pseudo_log"
  # )+
  geom_sf(data = centerlines1, aes(color = width_m), linewidth = linewidth_miss)+
  geom_sf(data = centerlines2, aes(color = width_m), linewidth = linewidth_miss)+
  geom_sf(data = centerlines3, aes(color = width_m), linewidth = linewidth_miss)+
  geom_sf(data = centerlines4, aes(color = width_m), linewidth = linewidth_miss)+
  geom_sf(data = centerlines5, aes(color = width_m), linewidth = linewidth_miss)+
  geom_sf(data = centerlines6, aes(color = width_m), linewidth = linewidth_miss)+
  geom_sf(data = centerlines7, aes(color = width_m), linewidth = linewidth_miss)+
  geom_sf(data = centerlines8, aes(color = width_m), linewidth = linewidth_miss)+
  geom_sf(data = centerlines9, aes(color = width_m), linewidth = linewidth_miss)+
  geom_sf(data = centerlines10, aes(color = width_m), linewidth = linewidth_miss)+
  geom_sf(data = centerlines11, aes(color = width_m), linewidth = linewidth_miss)+
  geom_sf(data = centerlines12, aes(color = width_m), linewidth = linewidth_miss)+
  geom_sf(data = centerlines13, aes(color = width_m), linewidth = linewidth_miss)+
  scale_fill_manual(values = c("Lakes and River Masks"="#99ccff"), name = element_blank())+
  # scale_color_manual(values = c("Lakes"="#7BD9F6"))+
  geom_sf(data = mississippi_lakes_wm, aes(fill = "Lakes and River Masks"), color = "transparent", show.legend = T)+

  labs(color = "Width (m)\n")+
  scale_x_continuous(limits = c(-114.25,-69))+
  scale_y_continuous(limits = c(28.76,50 #51
                                ))+
  # geom_segment(aes(x=-109.5,y=37.51,xend=-118,yend=28),
  #              color = "white", linewidth = 0.4)+
  annotation_scale(location = "bl", width_hint = 0.3, style = "bar",
                   pad_x = unit(0.06, "in"),
                   pad_y = unit(0.06, "in"),
                   text_cex = (6/6),
                   height = unit(0.1, "in"),
                   bar_cols = "white",
                   line_col = "white",
                   text_col = "white")+
  # geom_rect(aes(xmin = -109.5, xmax = -95.5, ymin = 37.51, ymax = 44.49), color = "white", fill = NA, linewidth = 2*linewidth)+
  coord_sf(expand = F)+
  theme_void()+
  theme(text = element_text(family = "sans", size = 12),
        # axis.text = element_text(color = "white"),
        # axis.line = element_line(color="white"),
        # axis.ticks = element_line(color = "white"),
        legend.position = c(0.165,0.3 #0.85,0.80
                            ),
        legend.background = element_blank(), #element_rect(color = "black", fill = "white"),
        legend.text = element_text(color="black", size = 10),
        legend.title = element_text(color = "black", size = 10,
                                    margin = margin(b=0,unit='in'),
                                    vjust = -1),
        # legend.key = element_rect(margin(t=0,unit='in')),
        legend.key.size = unit(0.2, 'in'),
        legend.margin = margin(t=-0.05,b=0.1,l=0.05,r=0.05,unit = 'in'),
        # legend.spacing.x = unit(0, "in"),
        legend.spacing.y = unit(0.05, "in"),
        panel.background = element_rect(fill=bg_color, color = NA_character_),
        plot.background = element_rect(fill='transparent', color=NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # legend.box.background = element_rect(fill='transparent')
        )+
  geom_image(
    x = -75.1, y = 46.1, size = 0.35, aes(image = "C:/Users/Carter/OneDrive - Virginia Tech/Desktop/Projects/2024_01_22_Code_for_Publication/Figure Making/Draft 7 Figures/Fig 1 Maps and Hists/mississippi_hist.png")
  )+
  geom_image(
    x = -73.82, y = 32.49, size = 0.35, aes(image = "C:/Users/Carter/OneDrive - Virginia Tech/Desktop/Projects/2024_01_22_Code_for_Publication/Figure Making/Draft 7 Figures/Fig 1 Maps and Hists/mississippi_watermask_map.png")
  )
  # guides(color = guide_legend(override.aes = list(size = 0.5)))
# miss_hist <- system.file("extdata","C:/Users/Carter/OneDrive - Virginia Tech/Desktop/Projects/2024_01_22_Code_for_Publication/Figure Making/Draft 7 Figures/Fig 1 Maps and Hists/mississippi_hist.png", package = "cowplot")
# 
# ggdraw(mississippi_width_map)+
#   draw_image(miss_hist, x = 0.9, y = 0.9, scale = 1)
ggsave(filename = "C:/Users/Carter/OneDrive - Virginia Tech/Desktop/Projects/2024_01_22_Code_for_Publication/Figure Making/Draft 7 Figures/Fig 1 Maps and Hists_v2/mississippi_width_map_v2.png",
       plot = mississippi_width_map_v2,
       width = 5.6,
       height = 3.4,
       units = "in",
       dpi = 1200,
       bg = "white")



#platte basin (done)
platte_width_map_v2 <- ggplot()+
  # geom_sf(data = coastline, color = bottom_color)+
  # geom_sf(data = countries, color = bottom_color, fill = bottom_color)+
  # geom_sf(data = mississippi_geom, color = bg_color, fill = bg_color)+
  geom_sf(data = platte_geom, color = fg_color, fill = fg_color)+
  geom_sf(data = svc_geom, color = nested_color, fill = nested_color)+
  # geom_sf(data = ssvc_geom, color = "white", fill = bg_color)+
  # scale_color_viridis_c(direction = -1, breaks = my_breaks, labels = my_breaks,
  #                       trans = scales::pseudo_log_trans())+
  scale_color_distiller(name = "Width (m)\n",
                        palette = centerline_palette,
                        direction = -1,
                        breaks = my_breaks,
                        labels = my_breaks,
                        trans = scales::pseudo_log_trans(),
                        limits = c(0,4003))+
  geom_sf(data = centerlines1_platte, aes(color = width_m), linewidth = linewidth)+
  geom_sf(data = centerlines2_platte, aes(color = width_m), linewidth = linewidth)+
  geom_sf(data = centerlines3_platte, aes(color = width_m), linewidth = linewidth)+
  geom_sf(data = centerlines4_platte, aes(color = width_m), linewidth = linewidth)+
  geom_sf(data = centerlines5_platte, aes(color = width_m), linewidth = linewidth)+
  geom_sf(data = centerlines6_platte, aes(color = width_m), linewidth = linewidth)+
  geom_sf(data = centerlines7_platte, aes(color = width_m), linewidth = linewidth)+
  geom_sf(data = centerlines8_platte, aes(color = width_m), linewidth = linewidth)+
  geom_sf(data = centerlines9_platte, aes(color = width_m), linewidth = linewidth)+
  geom_sf(data = centerlines10_platte, aes(color = width_m), linewidth = linewidth)+
  geom_sf(data = centerlines11_platte, aes(color = width_m), linewidth = linewidth)+
  geom_sf(data = centerlines12_platte, aes(color = width_m), linewidth = linewidth)+
  geom_sf(data = centerlines13_platte, aes(color = width_m), linewidth = linewidth)+
  geom_sf(data = platte_lakes_wm, color = "transparent", fill = "#99ccff")+
  labs(color = "Width (m)\n")+
  scale_x_continuous(limits = c(-109.5,-95.5))+
  scale_y_continuous(limits = c(38.7,43.3 #37.51,44.49
  ))+
  # geom_rect(aes(xmin = -105.844, xmax = -104.756, ymin = 39.8, ymax = 40.35), color = "white", fill = NA, linewidth = 2*linewidth)+
  coord_sf(expand = F)+
  annotation_scale(location = "br", width_hint = 0.3, style = "bar",
                   pad_x = unit(0.05, "in"),
                   pad_y = unit(0.05, "in"),
                   text_cex = (5/6),
                   height = unit(0.1, "in"),
                   bar_cols = "black",
                   line_col = "black",
                   text_col = "black")+
  theme_void()+
  theme(text = element_text(family = "sans", size = 12),
        legend.position = "none",
        legend.background = element_rect(color = NA, fill = NA),
        legend.text = element_text(color="white"),
        legend.title = element_text(color = "white"),
        panel.background = element_rect(fill='transparent', color = NA_character_),
        plot.background = element_rect(fill='transparent', color=NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.box.background = element_rect(fill='transparent'))
  # geom_image(
  #   x = -97, y = 39, size = 0.3, aes(image = "C:/Users/Carter/OneDrive - Virginia Tech/Desktop/Projects/2024_01_22_Code_for_Publication/Figure Making/Draft 7 Figures/Fig 1 Maps and Hists/platte_hist.png")
  # )

ggsave(filename = "C:/Users/Carter/OneDrive - Virginia Tech/Desktop/Projects/2024_01_22_Code_for_Publication/Figure Making/Draft 7 Figures/Fig 1 Maps and Hists_v2/platte_width_map_v2.png",
       plot = platte_width_map_v2,
       width = 3.2,
       height = 1.2,
       units = "in",
       dpi = 1200,
       bg = "white")



#st vrain basin (done)
svc_width_map_v2 <- ggplot()+
  # geom_sf(data = coastline, color = bottom_color)+
  # geom_sf(data = countries, color = bottom_color, fill = bottom_color)+
  # geom_sf(data = mississippi_geom, color = "white", fill = bg_color)+
  # geom_sf(data = platte_geom, color = bg_color, fill = bg_color)+
  geom_sf(data = svc_geom_union, color = fg_color, fill = fg_color)+
  geom_sf(data = ssvc_geom, color = nested_color, fill = nested_color)+
  # scale_color_viridis_c(direction = -1, breaks = my_breaks, labels = my_breaks,
  #                       trans = scales::pseudo_log_trans())+
  scale_color_distiller(name = "Width (m)\n",
                        palette = centerline_palette,
                        direction = -1,
                        breaks = my_breaks,
                        labels = my_breaks,
                        trans = scales::pseudo_log_trans(),
                        limits = c(0,4003))+
  geom_sf(data = centerlines1_svc, aes(color = width_m), linewidth = linewidth)+
  geom_sf(data = centerlines2_svc, aes(color = width_m), linewidth = linewidth)+
  geom_sf(data = centerlines3_svc, aes(color = width_m), linewidth = linewidth)+
  geom_sf(data = centerlines4_svc, aes(color = width_m), linewidth = linewidth)+
  geom_sf(data = centerlines5_svc, aes(color = width_m), linewidth = linewidth)+
  geom_sf(data = centerlines6_svc, aes(color = width_m), linewidth = linewidth)+
  geom_sf(data = centerlines7_svc, aes(color = width_m), linewidth = linewidth)+
  geom_sf(data = centerlines8_svc, aes(color = width_m), linewidth = linewidth)+
  geom_sf(data = centerlines9_svc, aes(color = width_m), linewidth = linewidth)+
  geom_sf(data = centerlines10_svc, aes(color = width_m), linewidth = linewidth)+
  geom_sf(data = centerlines11_svc, aes(color = width_m), linewidth = linewidth)+
  geom_sf(data = centerlines12_svc, aes(color = width_m), linewidth = linewidth)+
  geom_sf(data = centerlines13_svc, aes(color = width_m), linewidth = linewidth)+
  geom_sf(data = svc_lakes_wm, color = "transparent", fill = "#7BD9F6")+
  labs(color = "Width (m)\n")+
  scale_x_continuous(limits = c(-105.9,-104.7), breaks = c(-105.7,-105.5,-105.3,-105.1,-104.9))+
  scale_y_continuous(limits = c(39.84,40.31 #51
  ))+
  # geom_rect(aes(xmin = -105.66, xmax = -105.57, ymin = 40.038783, ymax = 40.098217), color = "white", fill = NA, linewidth = 2*linewidth)+
  coord_sf(expand = F)+
  annotation_scale(location = "br", width_hint = 0.25, style = "bar",
                   pad_x = unit(0.05, "in"),
                   pad_y = unit(0.05, "in"),
                   text_cex = (5/6),
                   height = unit(0.1, "in"),
                   bar_cols = "black",
                   line_col = "black",
                   text_col = "black")+
  theme_void()+
  theme(text = element_text(family = "sans", size = 12),
        legend.position = "none",
        legend.background = element_rect(color = NA, fill = NA),
        legend.text = element_text(color="white"),
        legend.title = element_text(color = "white"),
        panel.background = element_rect(fill='transparent', color = NA_character_),
        plot.background = element_rect(fill='transparent', color=NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.box.background = element_rect(fill='transparent'))

ggsave(filename = "C:/Users/Carter/OneDrive - Virginia Tech/Desktop/Projects/2024_01_22_Code_for_Publication/Figure Making/Draft 7 Figures/Fig 1 Maps and Hists_V2/svc_width_map_v2.png",
       plot = svc_width_map_v2,
       width = 3.2,
       height = 1.2,
       units = "in",
       dpi = 1200,
       bg = "white")



#south st vrain basin (done)
ssvc_width_map_v2 <- ggplot()+
  # geom_sf(data = coastline, color = bottom_color)+
  # geom_sf(data = countries, color = bottom_color, fill = bottom_color)+
  # geom_sf(data = mississippi_geom, color = "white", fill = bg_color)+
  # geom_sf(data = platte_geom, color = "white", fill = bg_color)+
  # geom_sf(data = svc_geom_union, color = bg_color, fill = bg_color)+
  geom_sf(data = ssvc_geom, color = fg_color, fill = fg_color)+
  # scale_color_viridis_c(direction = -1, breaks = my_breaks, labels = my_breaks,
  #                       trans = scales::pseudo_log_trans())+
  scale_color_distiller(name = "Width (m)\n",
                       palette = centerline_palette,
                       direction = -1,
                       breaks = my_breaks,
                       labels = my_breaks,
                       trans = scales::pseudo_log_trans(),
                       limits = c(0,4003))+
  geom_sf(data = centerlines1_ssvc, aes(color = width_m), linewidth = linewidth)+
  geom_sf(data = centerlines2_ssvc, aes(color = width_m), linewidth = linewidth)+
  geom_sf(data = centerlines3_ssvc, aes(color = width_m), linewidth = linewidth)+
  geom_sf(data = centerlines4_ssvc, aes(color = width_m), linewidth = linewidth)+
  geom_sf(data = centerlines5_ssvc, aes(color = width_m), linewidth = linewidth)+
  geom_sf(data = centerlines6_ssvc, aes(color = width_m), linewidth = linewidth)+
  geom_sf(data = centerlines7_ssvc, aes(color = width_m), linewidth = linewidth)+
  geom_sf(data = centerlines8_ssvc, aes(color = width_m), linewidth = linewidth)+
  geom_sf(data = centerlines9_ssvc, aes(color = width_m), linewidth = linewidth)+
  geom_sf(data = centerlines10_ssvc, aes(color = width_m), linewidth = linewidth)+
  geom_sf(data = centerlines11_ssvc, aes(color = width_m), linewidth = linewidth)+
  geom_sf(data = centerlines12_ssvc, aes(color = width_m), linewidth = linewidth)+
  geom_sf(data = centerlines13_ssvc, aes(color = width_m), linewidth = linewidth)+
  geom_sf(data = ssvc_lakes_wm, color = "transparent", fill = "#7BD9F6")+
  labs(color = "Width (m)\n")+
  scale_x_continuous(limits = c(-105.66,-105.57)#, breaks = c(-105.7,-105.5,-105.3,-105.1,-104.9)
                     )+
  scale_y_continuous(limits = c(40.048,40.084 #51
  ))+
  # geom_rect(aes(xmin = -105.66, xmax = -105.57, ymin = 40.045, ymax = 40.09), color = "white", fill = NA, linewidth = 1)+
  coord_sf(expand = F)+
  annotation_scale(location = "br", width_hint = 0.3, style = "bar",
                   pad_x = unit(0.06, "in"),
                   pad_y = unit(0.06, "in"),
                   text_cex = (5/6),
                   height = unit(0.1, "in"),
                   bar_cols = "black",
                   line_col = "black",
                   text_col = "black")+
  theme_void()+
  theme(text = element_text(family = "sans", size = 12),
        legend.position = "none",
        legend.background = element_rect(color = NA, fill = NA),
        legend.text = element_text(color="white"),
        legend.title = element_text(color = "white"),
        panel.background = element_rect(fill='transparent', color = NA_character_),
        plot.background = element_rect(fill='transparent', color=NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.box.background = element_rect(fill='transparent'))


ggsave(filename = "C:/Users/Carter/OneDrive - Virginia Tech/Desktop/Projects/2024_01_22_Code_for_Publication/Figure Making/Draft 7 Figures/Fig 1 Maps and Hists_V2/ssvc_width_map_v2.png",
       plot = ssvc_width_map_v2,
       width = 3.2,
       height = 1.2,
       units = "in",
       dpi = 1200,
       bg = "white")
#####


mississippi_width_map
platte_width_map
svc_width_map
ssvc_width_map

# mississippi_width_map + mississippi_hist + platte_width_map + platte_hist + svc_width_map + svc_hist + ssvc_width_map + ssvc_hist +
#   plot_layout(ncol = 2, widths = c(1,1))






