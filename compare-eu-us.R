################################################################################
#
# ikashnitsky.github.io 2018-11-21
# compare colorcoded -- EU, US
# Ilya Kashnitsky, ilya.kashnitsky@gmail.com
# Jonas Schöley, jschoeley@health.sdu.dk 
#                               
################################################################################

# load required packages
library(tidyverse)      # data manipulation and viz
library(janitor)
library(ggtern)         # plot ternary diagrams
library(gridExtra)      # arrange subplots
library(lubridate)      # easy manipulations with dates
library(ggthemes)       # themes for ggplot2
library(extrafont)      # custom font
library(hrbrthemes); import_roboto_condensed()
library(tricolore)


# Additional functions for Tricolore --------------------------------------

# a function to create zooming limits
zoom_limits <- function(
    # 3-columns data frame. ! Oreder is important: L, R, T
    df, 
    # whether to minimize zooming triangle and move the data center
    # or keep the data center at (1/3, 1/3, 1/3)
    keep_center = TRUE, 
    # add 1 percentage point margin to avoid cutting the extreme points
    one_pp_margin = FALSE,
    # the default is to calculate average from the provided data
    # though, I leave a possibility to specify custom center
    # in our case, custom center is the EU pop structure
    center = apply(df, 2, mean, na.rm = T)
) {
    # calculate minimums of the variables
    mins <- apply(df, 2, min)
    # calculate max data span
    span <- max(apply(df, 2, function(x) diff(range(x))))
    # add 1 percentage point margin to avoid cutting the extreme points
    if(one_pp_margin == TRUE & min(mins) > .01){
        mins <- mins - .01
        span <- span + .01
    }
    # keep the center at (1/3, 1/3, 1/3) or not
    if(keep_center == TRUE){
        limits <- rbind(
            center - (1/3)*span/(sqrt(2)/2),
            center + (2/3)*span/(sqrt(2)/2)
        )
    } else {
        limits <- rbind(
            mins,
            c(
                1 - (mins[2] + mins[3]),
                1 - (mins[1] + mins[3]),
                1 - (mins[1] + mins[2])
            )
        ) 
    }
    return(limits)
}


# coordinates and labels for the centered gridlines of a ternary diagram
TernaryCentroidGrid <- function (center) {
    # center percent difference labels
    labels <- seq(-1, 1, 0.1)
    labels <- data.frame(
        L = labels[labels >= -center[1]][1:10],
        T = labels[labels >= -center[2]][1:10],
        R = labels[labels >= -center[3]][1:10]
    )
    
    # breaks of uncentered grid
    breaks = data.frame(
        L = labels$L + center[1],
        T = labels$T + center[2],
        R = labels$R + center[3]
    )
    
    list(labels = labels, breaks = breaks)
}



range_diff <- function(x) x %>% range() %>% diff()

round_and_pad <- function(x) x %>% round(3) %>% paste %>% 
    str_pad(5, side = "right", pad = "0")




# preparation -------------------------------------------------------------
library(eurostat)
library(tigris)
library(tidycensus)
library(sf)
library(rmapshaper)
library(rgdal)


# create subdirectories
ifelse(
    !dir.exists("data"),
    dir.create("data"),
    paste("Directory already exists")
)

ifelse(
    !dir.exists("figures"),
    dir.create("figures"),
    paste("Directory already exists")
)


# get european data -------------------------------------------------------


# Find the needed dataset code 
# http://ec.europa.eu/eurostat/web/regions/data/database

# download the data on broad pop structures at NUTS-3 level
df_eu <- get_eurostat("demo_r_pjanaggr3")

save(df_eu, file = "data/df_eu.rda")

# if the automated download does not work, the data can be grabbed manually at
# http://ec.europa.eu/eurostat/estat-navtree-portlet-prod/BulkDownloadListing

# remote areas to remove (NUTS-2)
remote <- c(paste0('ES',c(63,64,70)),paste('FRA',1:5,sep=''),'PT20','PT30')


f <- tempfile()
download.file("http://ec.europa.eu/eurostat/cache/GISCO/geodatafiles/NUTS_2013_20M_SH.zip", destfile = f)
unzip(f, exdir = "data/.")
NUTS_raw <- readOGR("data/NUTS_2013_20M_SH/data/.", "NUTS_RG_20M_2013")

gd <- get_eurostat_geospatial(resolution = "20", year = "2013", nuts_level = "all")

# geodata
gd_eu3 <-  NUTS_raw %>% 
    st_as_sf() %>% 
    clean_names() %>% 
    filter(stat_levl == 3) %>% 
    st_transform(crs = 3035)

# filter NUTS-3, 2015, both sex, calculate shares
df_eu_str <- df_eu %>% filter(sex=="T", nchar(paste(geo))==5,
                        !str_sub(geo, 1, 4) %in% remote,
                        !str_sub(geo, 1, 2) %in% c("AL", "MK"),
                        year(time)==2015) %>% 
    droplevels() %>% 
    transmute(id = geo, age, value = values) %>% 
    spread(age, value) %>% 
    transmute(
        id,
        sy = Y_LT15 / TOTAL, 
        sw = `Y15-64` / TOTAL, 
        so = 1 - (sy + sw),
        pop_total = TOTAL
    ) %>% 
    left_join(gd_eu3, c("id" = "nuts_id")) %>% 
    st_as_sf()%>% 
    st_transform(crs = 3035)

save(df_eu_str, file = "data/df_eu_str.rda")


# country borders
bord_eu <- NUTS_raw %>% 
    st_as_sf() %>% 
    clean_names() %>% 
    filter(stat_levl == 0) %>% 
    st_transform(crs = 3035) %>% 
    ms_innerlines()

save(bord_eu, file = "data/bord_eu.rda")


# get american data -------------------------------------------------------

df_us <- get_estimates(
    geography = "county",
    product = "characteristics",
    breakdown = "AGEGROUP",
    breakdown_labels = TRUE,
    year = 2015,
    geometry = TRUE,
    shift_geo = TRUE
)

save(df_us, file = "data/df_us.rda")

# state boundaries
bord_us <- tidycensus::state_laea %>% 
    rmapshaper::ms_innerlines()

save(bord_us, file = "data/bord_us.rda")


# calculate ternary composition
df_us_str <- df_us %>% 
    clean_names() %>% 
    rename(id = geoid) %>% 
    group_by(id, name) %>%
    summarise(
        so = value[27] / value[1],
        sw = (value[31]+value[26]) / value[1],
        sy = 1 - (so + sw),
        pop_total = value[1]
    ) %>%
    ungroup()

save(df_us_str, file = "data/df_us_str.rda")



# add ternary colours ---------------------------------------------------------

load("data/df_eu_str.rda")
load("data/df_us_str.rda")

df_hex <- bind_rows(
    eu = df_eu_str %>% select(id, so, sw, sy, pop_total) %>% st_set_geometry(NULL),
    us = df_us_str %>% select(id, so, sw, sy, pop_total) %>% st_set_geometry(NULL),
    .id = "geo"
) 


# Whole data mean 
center <- df_hex %>% 
    select(so, sw, sy) %>% 
    summarise_all(.funs = funs(mean)) %>% 
    gather() %>% pull(value)

# calculate TRUE scaling factor for colors, i.e. the factor of proportionality
# from big tern to zoomed tern
mins <- apply(df_hex %>% select(so, sw, sy), 2, min)
zommed_side <- (1 - (mins[2] + mins[3])) - mins[1]
true_scale <- 1 / zommed_side

tric <- Tricolore(
    df_hex, p1 = 'so', p2 = 'sw', p3 = 'sy',
    center = center, show_data = FALSE, spread = true_scale,
    contrast = .5, lightness = 1, chroma = 1, hue = 2/12,
    breaks = 20,
    crop = TRUE, label_as = "pct_diff"
)

df_hex$hex <- tric$rgb

save(df_hex, file = "data/df_hex.rda")


# add hex to dfs
str_eu <- df_eu_str %>% 
    left_join(
        df_hex %>% filter(geo == "eu") %>% select(id, hex),
        by = "id"
    )

save(str_eu, file = "data/str_eu.rda")

str_us <- df_us_str %>% 
    left_join(
        df_hex %>% filter(geo == "us") %>% select(id, hex),
        by = "id"
    )

save(str_us, file = "data/str_us.rda")



# percent-point difference grid
legend_grid <- TernaryCentroidGrid(center)

# legend limits
legend_limits <- zoom_limits(
    df = df_hex %>% select(so, sw, sy),
    keep_center = FALSE,
    one_pp_margin = TRUE
) # try playing with the params


# -------------------------------------------------------------------------
# READY TO PLAY -----------------------------------------------------------



# map Europe -------------------------------------------------------------

load("data/str_eu.rda")
load("data/bord_eu.rda")
load("data/cities_eu.rda")

str_eu %>% 
    ggplot()+
    geom_sf(aes(fill = hex), color = NA)+
    geom_sf(data = bord_eu, color = "#ffffff", size = .5)+
    geom_sf(data = cities_eu, size = 4.5, shape = 1, stroke = .5, color = "#ffffff")+
    scale_fill_identity()+
    coord_sf(datum = NA)+
    theme_minimal(base_family = font_rc)+
    theme(axis.text = element_blank())

map_eu <- last_plot()

ggsave(plot = map_eu, "figures/map-eu.png", 
       width = 12, height = 10, dpi = 600, type = "cairo-png")


# map America ------------------------------------------------------------

load("data/str_us.rda")
load("data/bord_us.rda")
load("data/cities_us.rda")

str_us %>% 
    ggplot()+
    geom_sf(aes(fill = hex), color = NA)+
    geom_sf(data = bord_us, color = "#ffffff", size = .5)+
    geom_sf(data = cities_us, size = 4.5, shape = 1, stroke = .5, color = "#ffffff")+
    scale_fill_identity()+
    coord_sf(datum = NA)+
    theme_minimal(base_family = font_rc)+
    theme(axis.text = element_blank())

map_us <- last_plot()

ggsave(plot = map_us, "figures/map-us.png", 
       width = 12, height = 10, dpi = 600, type = "cairo-png")

# -------------------------------------------------------------------------



# legend Europe -----------------------------------------------------------

tric$key +
    geom_point(data = df %>% filter(geo == "eu"), aes(so, sw, z = sy), 
               shape = 21, fill = "grey50", size = .5, alpha = .5)+
    geom_point(data = tibble(so = center[1], sw = center[2], sy = center[3]), 
               aes(so, sw, z = sy), 
               shape = 43, color = "white", size = 5)+
    scale_L_continuous(NULL, limits = legend_limits[,1]) +
    scale_T_continuous(NULL, limits = legend_limits[,2]) +
    scale_R_continuous(NULL, limits = legend_limits[,3]) +
    theme_classic() +
    theme(plot.background = element_rect(fill = "white", colour = NA),
          text = element_text(family = font_rc, size = 15, color = "grey20"))

legend_eu <- last_plot()

ggsave(plot = legend_eu, "figures/legend-eu.png", 
       width = 6, height = 5, dpi = 600, type = "cairo-png")



# legend US ---------------------------------------------------------------

tric$key +
    geom_point(data = df %>% filter(geo == "us"), aes(so, sw, z = sy), 
               shape = 21, fill = "grey50", size = .5, alpha = .5)+
    geom_point(data = tibble(so = center[1], sw = center[2], sy = center[3]), 
               aes(so, sw, z = sy), 
               shape = 43, color = "white", size = 5)+
    scale_L_continuous(NULL, limits = legend_limits[,1]) +
    scale_T_continuous(NULL, limits = legend_limits[,2]) +
    scale_R_continuous(NULL, limits = legend_limits[,3]) +
    theme_classic() +
    theme(plot.background = element_rect(fill = "white", colour = NA),
          text = element_text(family = font_rc, size = 15, color = "grey20"))

legend_us <- last_plot()

ggsave(plot = legend_us, "figures/legend-us.png", 
       width = 6, height = 5, dpi = 600, type = "cairo-png")



# ternary ellipses --------------------------------------------------------

df_hex %>% 
    filter(!sy %in% range(sy)) %>% 
    ggtern(aes(so, sw, z = sy, color = geo))+
    geom_point(shape = 21, fill = "grey50", size = .5, alpha = .5)+
    scale_L_continuous("Elderly\n(65+)", limits = legend_limits[,1]) +
    scale_T_continuous("Working age\n(15-64)", limits = legend_limits[,2]) +
    scale_R_continuous("Young\n(0-14)", limits = legend_limits[,3])+
    geom_mean_ellipse(size = 1)+
    scale_color_manual(values = c("grey25", "gold"))+
    labs(x = NULL, y = NULL)+
    Larrowlab("% aged 65+") +
    Tarrowlab("% aged 15-64") +
    Rarrowlab("% aged 0-14") +
    theme(tern.axis.arrow.show = TRUE, 
          tern.axis.ticks.length.major = unit(12, "pt"),
          tern.axis.text = element_text(size = 15, colour = "grey20"),
          tern.axis.title.T = element_text(),
          tern.axis.title.L = element_text(hjust = 0.2, vjust = 0.7, angle = -60),
          tern.axis.title.R = element_text(hjust = 0.8, vjust = 0.6, angle = 60),
          text = element_text(family = font_rc, size = 18, color = "grey20"),
          legend.position = "none")

gg_global <- last_plot()




gg_global_build <- ggplot_build(gg_global)

# table of rahges
gg_global_range <- gg_global_build$data[[2]] %>% 
    as_tibble() %>% 
    group_by(group) %>% 
    summarise_at(.vars = vars(x, y, z), 
                 .funs = funs(range_diff)) %>% 
    #clean up the table a bit
    set_names(c("Region", "Elderly", "Working", "Young")) %>% 
    mutate(Region = c("Europe", "The US") %>% paste)

# table of means
gg_global_avg <- gg_global_build$data[[2]] %>% 
    as_tibble() %>% 
    group_by(group) %>% 
    summarise_at(.vars = vars(x, y, z), 
                 .funs = funs(mean)) %>% 
    #clean up the table a bit
    set_names(c("Region", "Elderly", "Working", "Young")) %>% 
    mutate(Region = c("Europe", "The US") %>% paste)




gg_tab_global_range <- ggplot()+
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = F)+
    annotation_custom(
        tableGrob(
            gg_global_range %>% 
                mutate_if(is.numeric, round_and_pad),
            theme = ttheme_minimal(base_family = font_rc, base_size = 20),
            rows = NULL
            
        ),
        xmin = 0, xmax = 1,
        ymin = 0, ymax = 1
    )

gg_tab_global_avg <- ggplot()+
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = F)+
    annotation_custom(
        tableGrob(
            gg_global_avg %>% 
                mutate_if(is.numeric, round_and_pad),
            theme = ttheme_minimal(base_family = font_rc, base_size = 20),
            rows = NULL
            
        ),
        xmin = 0, xmax = 1,
        ymin = 0, ymax = 1
    )




ggplot()+
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = F)+
    annotation_custom(
        ggplotGrob(gg_global),
        xmin = 0, xmax = 1,
        ymin = 0, ymax = 1
    )+
    annotation_custom(
        ggplotGrob(gg_tab_global_avg),
        xmin = 0, xmax = .4,
        ymin = .7, ymax = .95
    )+
    annotate("text", label = "Average compositions",
             x = .2,  y = .97, hjust = .5, vjust = 1, 
             size = 5, family = font_rc
    )+
    annotation_custom(
        ggplotGrob(gg_tab_global_range),
        xmin = .6, xmax = 1,
        ymin = .7, ymax = .95
    )+
    annotate("text", label = "Ranges of ellipses",
             x = .8,  y = .97, hjust = .5, vjust = 1, 
             size = 5, family = font_rc
    )+
    theme_void(base_family = font_rc) ->
    out_global


ggsave(plot = out_global, "figures/ellipses.png", 
       width = 12, height = 10, dpi = 600, type = "cairo-png")


# compare ecdf ------------------------------------------------------------

options(scipen=10)

df_hex %>% 
    mutate(geo = geo %>% fct_recode(Europe = "eu", `The US` = "us")) %>% 
    ggplot() + 
    geom_hline(yintercept = .5, color = 'grey50') +
    stat_ecdf(aes(pop_total, color = geo), size = 1) + 
    scale_color_manual(NULL, values = c("grey25", "gold"))+
    scale_x_comma(trans = "log10", breaks = 10^(3:7), labels = )+
    labs(x = "Log of total population size of a region",
         y = "Empirical cumulative density")+
    theme_grey(base_size = 20, base_family = font_rc)+
    theme(legend.position = c(.12, .85),
          legend.background = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(hjust = .9, size = 14))

ecdf_pop <- last_plot()

ggsave(plot = ecdf_pop, "figures/ecdf.png", 
       width = 6, height = 6, dpi = 600, type = "cairo-png")



# FINAL POSTER ------------------------------------------------------------

# The poster is designed in inkScape 
