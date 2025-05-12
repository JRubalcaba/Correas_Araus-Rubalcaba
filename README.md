# R CODE and DATABASE for the paper: "The energetic impacts of climate change vary along climatic gradients in mammals"

Marta Correas-Araus, Juan G. Rubalcaba
Department of Biodiversity, Ecology and Evolution, Faculty of Biological Sciences, Complutense University of Madrid
Juan G. Rubalcaba, e-mail: jg.rubalcaba@gmail.com 
Marta Correas-Araus, e-mail: marta.coraraus@gmail.com 

Abstract
Extreme temperature fluctuations are becoming more frequent due to climate change, likely increasing the cost of living of endotherms through their impact on the energy and water budgets. 
Yet, we know little about the magnitude of these impacts and whether species from different climatic regions will differ in their response to cold and heat waves. We confronted mechanistic models 
with field metabolic rate data (FMR) of mammals to analyze the response of energy and water balance to temperature anomaly across climatic gradients. FMR increased in response to both positive 
and negative anomalies in colder regions, and this pattern was consistent across both observed and modelled FMR. By contrast, warm-living mammals displayed lower FMR during heat waves, a pattern 
that was not captured by the model, probably due to physiological and behavioral adaptations to prevent water loss in warm, dry environments. The energetic impacts of climate change will differ 
among mammals inhabiting different climates, and therefore, species-specific parameterization is required for future mechanistic models to predict these impacts. 

# R version and required packages
R version 4.4.1 (2024-06-14 ucrt)
Platform: x86_64-w64-mingw32/x64
Running under: Windows 11 x64 (build 22631)

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 MCMCglmm_2.36      coda_0.19-4.1      ggpubr_0.6.0       RColorBrewer_1.1-3 NicheMapR_3.3.2   
 car_3.1-2          carData_3.0-5      nls2_0.3-4         proto_1.0.0        stringr_1.5.1     
 lme4_1.1-35.5      Matrix_1.7-0       ggplot2_3.5.1      nlme_3.1-164      
 phytools_2.4-4     maps_3.4.2.1       ape_5.8            readxl_1.4.3      

 # Code description

1) Run the code "Endotherm_model_simulation.R" to generate simulated FMR
2) The code "PGLS_analysis.R" perform phylogenetic analyses in R. Note that the phylogenetic tree "FritzTree_mammals_consensus.tre" is required for the simulations:
     Fritz, S. A., Bininda‐Emonds, O. R., & Purvis, A. (2009). Geographical variation in predictors of mammalian extinction risk: big is bad, but only in the tropics. Ecology letters, 12(6), 538-549.
     https://onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1111%2Fj.1461-0248.2009.01307.x&file=ELE_1307_sm_SA1.tre

   


 


