#Code for H5N1 dairy model manscript
#B. Bellotti, M. DeWitt
#2024-09-20

############
#Libraries
############
library(dplyr)
library(sf)
library(ggplot2)
library(here)
library(tidyUSDA)
#library(idstyle)

###############
#Read in USDA data key
###############
#Key from USDA - register for a key and read in ######################INPUT REQUIRED#########################
#readRenviron(".env")

###############
#Functions
###############
#function to read zipped shapefiles from web
read_shape_URL <- function(URL){
  cur_tempfile <- tempfile()
  download.file(url = URL, destfile = cur_tempfile)
  out_directory <- tempfile()
  unzip(cur_tempfile, exdir = out_directory)

  st_read(dsn = out_directory) #read_sf also works here
}

###############
#Data
###############
#H5N1 dairy case data saved from APHIS H5N1 dashboard
df.cases <- read.csv(here::here('data-raw', '2024-09-20_H5N1cases.csv'))

#read in state shapefiles from census:
sf.usa <- read_shape_URL("https://www2.census.gov/geo/tiger/GENZ2018/shp/cb_2018_us_state_500k.zip")%>%
  filter(GEOID < 60 & GEOID != "02" & GEOID != "15") #continental US
usa <-  left_join(sf.usa, df.cases %>% count(State) %>% rename(Cases = n, NAME = State)) #append case counts

#Totals by state
cases_summary <- usa%>%
  group_by(STUSPS)%>%
  summarise(Total = sum(Cases))%>%
  arrange(desc(Total))%>%
  st_drop_geometry()

# Get count of operations with sales in 2017 from USDA quick stats
ops.with.sales <- tidyUSDA::getQuickstat(
  sector = NULL,
  group = NULL,
  commodity = NULL,
  category = NULL,
  domain = NULL,
  county = NULL,
  key = Sys.getenv("KEY"),
  program = "CENSUS",
  data_item = "CATTLE, COWS, MILK - INVENTORY", #dairy operations specific
  geographic_level = "COUNTY",
  year = "2017",
  state = NULL,
  geometry = TRUE,
  lower48 = TRUE)

# Generate centroids for states with outbreaks
z <- left_join(cases_summary, tigris::states()) |>
  #filter(STUSPS %in% c("TX", "SD", "OH", "CO", "NC", "ID", "MI", "NM", "KS", "IA")) |>
  rename("Outbreaks" = "Total") |>
  st_as_sf() |>
  group_by(NAME) |>
  st_centroid()

p_cattle_density <- tigris::counties(state = state.abb[!state.abb %in% c("HI", "AK")]) |>
  select(STATEFP, COUNTYFP, GEOID, ALAND) |>
  mutate(ALAND = ALAND/1e9) |>
  left_join(ops.with.sales |>
              as.data.frame() |>
              select(GEOID, Value, domain_desc) |>
              filter(domain_desc == "TOTAL") , by = "GEOID") |>
  mutate(Value = ifelse(is.na(Value), 0, Value))

###############
#Figures
###############
map1 <- ggplot() +
  theme_void(base_size = 24)+
  geom_sf(data = usa, colour = "grey40", size=0.1, aes(fill = Cases))+
  scale_fill_distiller(direction = 1, palette = "Reds", na.value = "grey60")+
  labs(fill = "Outbreaks \nper state")

map2 <- ggplot(p_cattle_density) +
  geom_sf(aes(fill = Value/ALAND),  lwd = 0) +
  scale_fill_viridis_c(trans = "log",
                       labels= c("0", "1", "10", "100", "1k", "10k", "100k"),
                       breaks =c(0, 1, 10, 100, 1000, 10000, 100000)) +
  theme_void(base_size = 24) +
  #scale_fill_gradient(trans = "log") +
  geom_sf(data = z, aes(size = Outbreaks),
          shape = 19, alpha=0.5, color = "tomato")+
  scale_size_continuous(range = c(3, 20))+
  #theme(legend.position = "top") +
  guides(fill = guide_legend(position = "top", title.position="top"))+
  labs(fill = "Head of dairy cattle per 1,000 sq. km", size = "Outbreaks \nper state")

fig <- ggpubr::ggarrange(map1, map2, ncol=1,
                  widths=c(8, 10),
                  labels = c("A", "B"),
                  vjust=c(25,25))

ggsave("Maps.tiff", plot=fig, path = here::here("figures"), width = 8, height = 8, device='tiff', dpi=300)
