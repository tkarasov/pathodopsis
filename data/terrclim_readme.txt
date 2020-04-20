Climate data were downloaded from the TERRACLIMATE database 

website: http://www.climatologylab.org/terraclimate.html

cd ~/Dropbox/data/geospatial_datasets/terraclim
for year in {2007..2018}
do
wget -nc -c -nd "https://climate.northwestknowledge.net/TERRACLIMATE-DATA/TerraClimate_aet_${year}.nc" 
wget -nc -c -nd "https://climate.northwestknowledge.net/TERRACLIMATE-DATA/TerraClimate_def_${year}.nc" 
wget -nc -c -nd "https://climate.northwestknowledge.net/TERRACLIMATE-DATA/TerraClimate_pet_${year}.nc" 
wget -nc -c -nd "https://climate.northwestknowledge.net/TERRACLIMATE-DATA/TerraClimate_ppt_${year}.nc" 
wget -nc -c -nd "https://climate.northwestknowledge.net/TERRACLIMATE-DATA/TerraClimate_q_${year}.nc" 
wget -nc -c -nd "https://climate.northwestknowledge.net/TERRACLIMATE-DATA/TerraClimate_soil_${year}.nc" 
wget -nc -c -nd "https://climate.northwestknowledge.net/TERRACLIMATE-DATA/TerraClimate_srad_${year}.nc" 
wget -nc -c -nd "https://climate.northwestknowledge.net/TERRACLIMATE-DATA/TerraClimate_swe_${year}.nc" 
wget -nc -c -nd "https://climate.northwestknowledge.net/TERRACLIMATE-DATA/TerraClimate_tmax_${year}.nc" 
wget -nc -c -nd "https://climate.northwestknowledge.net/TERRACLIMATE-DATA/TerraClimate_tmin_${year}.nc" 
wget -nc -c -nd "https://climate.northwestknowledge.net/TERRACLIMATE-DATA/TerraClimate_vap_${year}.nc" 
wget -nc -c -nd "https://climate.northwestknowledge.net/TERRACLIMATE-DATA/TerraClimate_ws_${year}.nc" 
wget -nc -c -nd "https://climate.northwestknowledge.net/TERRACLIMATE-DATA/TerraClimate_vpd_${year}.nc" 
wget -nc -c -nd "https://climate.northwestknowledge.net/TERRACLIMATE-DATA/TerraClimate_PDSI_${year}.nc" 
done

Citation: Abatzoglou, J.T., S.Z. Dobrowski, S.A. Parks, K.C. Hegewisch, 2018, Terraclimate, a high-resolution global dataset of monthly climate and climatic water balance from 1958-2015, Scientific Data,

Datasets:
tmax:Maximum temperature, 
tmin: minimum temperature, 
vp:vapor pressure, 
ppt:precipitation accumulation, 
srad:downward surface shortwave radiation, 
ws: wind-speed
pet:Reference evapotranspiration (ASCE Penman-Montieth), 
q:Runoff*, 
aet:Actual Evapotranspiration*, 
def:Climate Water Deficit*, 
soil:Soil Moisture*, 
swe:Snow Water Equivalent*, 
PDSI:Palmer Drought Severity Index, 
vpd:Vapor pressure deficit



