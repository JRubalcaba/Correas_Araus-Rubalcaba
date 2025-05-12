require(NicheMapR)
require(RColorBrewer)
require(ggplot2)

############### Endotherm Model (Kearney et al. 2021 Ecography) ------

####### TRAITS -----
TC_MAX = 45 # max body temperature 
SHAPE_B	= 2 # current ratio between long and short axis (-)
ZFURD = 0.002  # hair layer depth (m)
ZFURV = 0.002
RHOD = 3e+07  # hair density (1/m2)
RHOV = 3e+07
DHAIRD = 3e-05 # fur diameter (m)
DHAIRV = 3e-05  
LHAIRD = 0.0239 # fur length
LHAIRV = 0.0239
KHAIR = 0.209 # hair thermal conductivity (W/m°C) 
REFLD = 0.2 # reflectivity dorsal y ventral (fractional, 0-1) 
REFLV = 0.2
EMISAN = 0.9 # animal emissivity (-)
PCTWET = 0.2 # part of the skin surface that is wet (%)
TC_INC = 0.1 # turns on core temperature elevation, the value being the increment by which TC is increased per iteration
PCTWET_INC = 0.1 # turns on sweating, the value being the increment by which PCTWET is increased per iteration (%)
AK1_INC = 0.1 # turns on thermal conductivity increase (W/mK), the value being the increment by which AK1 is increased per iteration (W/m°C)
AK1_MAX = 2.8 # maximum flesh conductivity (W/mK)
SUBQFAT = 1 # is subcutaneous fat present? (0 is no, 1 is yes)
FATPCT = 10 # % body fat
VEL = 0.1 # wind speed

####### SIMULATED vs OBSERVED  -----

# FMR data
dir <- "C:/Users/juanv/Dropbox/Talento CAM/Field metabolic rates/Assessing the energetic impacts of climate change on biodiversity (2022-T1AMB-23753)/Manuscript/Code and Data"
data <- read.csv(paste0(dir,"/FMR_database.csv"), sep=";")

# Get mean annual temperature
global_climate <- raster::brick("C:/globalclimate/global_climate.nc") # New et al. 1999
TMINN <- global_climate[[38:49]]
TMAXX <- global_climate[[50:61]]
TMINN <- raster::extract(TMINN, data.frame(x=data$Lon_deg, y=data$Lat_deg))/10
TMAXX <- raster::extract(TMAXX, data.frame(x=data$Lon_deg, y=data$Lat_deg))/10
months <- array(1, dim=c(nrow(data),12))
annualT <- (rowMeans(TMINN * months, na.rm=T) + rowMeans(TMAXX * months, na.rm=T)) / 2

# Run endotherm model

FMR_ind_simul <- numeric(nrow(data))
for(i in 1:nrow(data)){ 
  Ta <- data$TALOCmean[i] 
  mass <- data$Mass[i]/1000 
  QSOLR <- data$SOLRAD[i] 
  RH <- data$RHLOCmean[i] 
  mean_annT <- data$TMEAN[i] 
  
  if(is.na(QSOLR)) QSOLR <- 200 
  if(is.na(Ta)) Ta <- 18
  if(is.na(RH)) RH <- 70
  if(is.na(mean_annT)) mean_annT <- 10
  
  # CORE TEMPERATURE
  TC = 35.8 + 0.30 * log10(mass*1000) # White & Seymour 2003 PNAS
  # BASAL METABOLIC RATE
  QBASAL = 4.12 * (mass*1000)^0.69 # mLO2/h (White & Seymour 2003 PNAS)
  QBASAL = QBASAL * 20.9 / 3600 # Transform to Watt
  
  # Run NicheMapR
  endo.out <- endoR_devel(TA=Ta, VEL = VEL, RH = RH, QSOLR = QSOLR,
                          TC = TC, TC_MAX = TC_MAX, QBASAL = QBASAL,
                          AMASS = mass, SHAPE_B = SHAPE_B,
                          ZFURD = ZFURD, ZFURV = ZFURV, DHAIR = DHAIRD, DHAIRV = DHAIRV,
                          LHAIRD = LHAIRD, LHAIRV = LHAIRV,
                          REFLD =  REFLD, REFLV =  REFLD, KHAIR=KHAIR,
                          RHOD = RHOD, RHOV = RHOV,
                          PCTWET_INC = PCTWET_INC, TC_INC = TC_INC, 
                          AK1_INC = AK1_INC, AK1_MAX = AK1_MAX, 
                          SUBQFAT = SUBQFAT, FATPCT = FATPCT,
                          EMISAN = EMISAN)
  
  enbal <- data.frame(endo.out$enbal)
  FMR_ind_simul[i] <- enbal$QGEN + enbal$QEVAP 
}

# Remove cases with failed microclimate estimations
FMR_ind_simul[which(is.na(data$TALOCmean | data$SOLRAD |  data$RHLOCmean | annualT))] <- NA
FMR_M_ind_simul <- residuals(lm(log(FMR_ind_simul) ~ log(data$Mass), na.action = "na.exclude"))
