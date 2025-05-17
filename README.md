Codes used in the analyses of the manuscript:  
Automated extraction of leaf mass per area from digitized herbarium specimens (submitted) 
  
Authors: Thais Vasconcelos, Will Weaver, Aly Baumgartner, Zoe Bugnaski, James Boyko 

----
Description of folders: 
 
- **data/** 
 This folder contains .csv files with measurements of petiole width obtained through both manual methods and computer vision, as well as leaf area measurements generated using LeafMachine 2. 

- **models/** 
  This folder contains output files from multimodal inference analyses conducted using both standardized and non-standardized variables.

- **plots/** 
  This folder contains plots generated from data analyses, many of which are included as figures in the manuscript.

- **supporting_datasets/** 
  This folder includes .csv files with supporting information such as the list of sampled species, life form classifications, a list of plant families excluded from analyses (non-woody dicots among seed plants), and the finalized list of woody dicot species used in the study.

- **tables/** 
   This folder contains .csv files summarizing results from multimodal inference analyses, along with other tabulated outputs used for reporting and interpretation.

- **trees/** 
  This folder contains phylogenetic trees used in the analyses, provided in formats compatible with downstream comparative methods.  
  

----
Scripts:

> 01_list_of_woody_dicots.R  
Creates a list of woody, non-monocot angiosperm species based on the World Checklist of Vascular Plants (WCVP) dataset. It begins by loading two supporting datasets: one with families to exclude (monocots, gymnosperms, and ferns/lycophytes) and another for life form mapping. The script merges the WCVP distribution and names tables, filters for accepted species-level records, and excludes unwanted plant groups. It then assigns life form categories, identifies woody perennials, and removes hybrids. 
  
> 02_first_pass.R  
Identifies woody dicot species from the list created in *01_list_of_woody_dicots.R* that have at least 10 imaged specimens with valid coordinates on GBIF. It uses the mvh package to search GBIF for specimen metadata, iterating through each species and saving the results in a list. The script then checks how many specimens each species has with valid coordinates, keeping only those with 10 or more records. The final list of target taxa is exported as a .csv file.
  
> 03_pruning_tree.R  
Prunes a large phylogenetic tree (Smith & Brown, 2018) to include only the woody dicot species with at least 10 imaged specimens with valid coordinates on GBIF resulting from script *02_first_pass.R*. It loads the original tree (taxized_GBMB.Rdata) and the list of target species, then uses the ape package to keep only the matching species. The pruned tree is saved in newick format, and a reference table (ref_table.Rsave) is created to map the original taxized tree tip labels to the newick-formatted labels
  
> 04_focal_sample_of_woody_dicots.R  
Reads a phylogenetic tree, applies pre-calculated phylogenetic weights, randomly samples 2000 species in a phylogenetically balanced way, and visualizes the selected taxa on the tree. It then saves the sample as a CSV file.

> 05_specimen_image_download.R  
Downloads herbarium specimen images for a list of target species, prioritizing high-quality and relevant specimens based on geographic, temporal, and institutional preferences. It filters and scores the metadata, downloads up to 5 images per species, checks for image quality (removing low-resolution ones), and saves the images and metadata into organized folders.
  
> 06_post_processing.R  
Processes Leaf Machine output by cleaning and merging it with petiole width data, calculates LMA, removes outliers, verifies geographic accuracy, and enriches the dataset with biome and environmental variables for further analysis.  
  
> 07_pcms_standartized.R  
This script standardizes trait and climate variables, performs phylogenetic linear regressions, and uses model selection and averaging to identify key environmental predictors of LMA variation. 

> 08_r2_comparisons.R
This script compares phylogenetic and non-phylogenetic models to estimate the explanatory power (R²) of climate predictors on leaf mass per area (LMA). It performs model averaging, fits both linear and phylogenetic linear models using the averaged formula, and calculates multiple R² metrics, including likelihood-based and residual-based values, using the rr2 package.

> 09_lma_maps_and_plot_by_biome.R
This script analyzes leaf mass per area (LMA) data by creating global distribution maps and biome-specific violin plots. It aggregates LMA values by filename, then visualizes the global variation in LMA and data point density across 2.5° x 2.5° grid cells. Additionally, it generates violin plots for biomes with at least 100 observations, highlighting the distribution of LMA, sample sizes, and median values within each biome. 
  
> 10_lma_vs_error_SI.R
This script visualizes the relationship between absolute error in petiole width measurement and estimated leaf mass per area (LMA), showing a weak but significant correlation. A linear model fit and color gradient highlight how LMA varies with measurement error.
  
> 11_descriptive_plots_results.R
This script analyzes leaf morphology, focusing on leaf area, petiole width, and LMA. It calculates and visualizes species distribution, generates histograms for key measurements, and compares automated vs. manual petiole width measurements. The results, including histograms and scatter plots, are saved in a PDF for Figure 3.
  
> 12_pcm_plots.R
This script calculates mean environmental variables and fits phylogenetic generalized least squares (PGLS) models to assess the relationship between leaf mass per area (LMA) and each environmental factor. The results, including p-values, are stored and used to annotate scatter plots of LMA against each environmental variable. 
  
> 13_ldg_plot.R
This script processes a dataset to visualize how leaf mass per area (LMA) varies with latitude. It calculates the mean LMA for each 1-degree latitude bin and creates line plots for both mean LMA and the number of observations by latitude. 
  
> 14_plotting_regression_results.R
This script loads a full regression model and a set of candidate models, then performs model averaging using the top models within a 2 AICc unit range. It extracts coefficients and their confidence intervals, filters the results, and then visualizes the effects of different terms in the model. The first plot shows the effect size of each term with confidence intervals, color-coded by the sign of the estimate (positive or negative) and adjusted for statistical significance. The second plot shows the importance values of each term, indicating their relative influence on the model. 
  
> 15_ci_around_estimates_SI.R
This script estimates Leaf Mass per Area (LMA) from petiole width and leaf area using a calibration model and the Royer equation, incorporating error propagation through bootstrap simulations. It generates confidence intervals for each specimen's LMA estimate and quantifies the relative contributions of calibration and model error to overall uncertainty. The script outputs summary statistics and visualizations, including histograms and bar plots of confidence intervals and error source contributions.
  
> 16_comparison_with_other_datasets.R
This script shows geographical areas where species sampled by us, but not in TRY, occur. It identifies regions with missing species-level data and visualizes species and family richness across tropical areas. It highlights data gaps which could guide future sampling efforts. The map generated shows species richness with white areas for regions with zero data.
  
----
  
**Climate layers from:**
  
Karger, et al. (2017). Climatologies at high resolution for the earth’s land surface areas. Scientific data, 4(1), 1-20.  
  
Trabucco, A., & Zomer, R. (2018). Global aridity index and potential evapo- transpiration (ET0) climate database v2. CGIAR Consortium for Spatial Information (CGIAR-CSI). Published online, available from the CGIAR-CSI GeoPortal at https://cgiarcsi.community
  

GBIF POINTS COME FROM:
  GBIF.org (16 May 2025) GBIF Occurrence Download  https://doi.org/10.15468/dl.538h6w

  
on git ignore (files are too large for the repo):
- POWO data
- gbif points from Dorey et al.
- environmental layers

Note: packages rgeos, rgdal, and maptools were unfortunately retired as of October 2023. Much of the code here depends on functions from these three packages and, although the functionalities of these packages have been migrated to other packages, I still have to make changes to the code so that it works in computers that dont already have rgeos, rgdal, and maptools installed. 
