## ysm_htvg: Height Vigour Analysis

Collection of scripts and functions to analyze different htvg estimation methods for TASS projections of YSM samples


### Process

**1. Select samples for analysis.**

- Subset of YSM samples (sample_data.rds) available in [here](https://github.com/bcgov/YSMtechrep/tree/main/shiny_app/data)

Example: Williams Lake TSA

```{r}
library(data.table)
library(dplyr)

ci_wil <- sample_data %>% filter(MGMT_UNIT == "TSA29_WilliamsLake", visit_number_new == 1) %>% pull(CLSTR_ID)
ci_wil_2 <- sample_data %>% filter(MGMT_UNIT == "TSA29_WilliamsLake", visit_number_new == 2) %>% pull(CLSTR_ID)
ci_wil_all <- sample_data %>% filter(MGMT_UNIT == "TSA29_WilliamsLake") %>% pull(CLSTR_ID)
```


**2. Compute height growth for the selected samples.**

- Dataset to compute height growth up to the next visit or the last visit
- Use YSM tree data (tree_data.rds) available in [here](https://github.com/bcgov/YSMtechrep/tree/main/shiny_app/data)

```{r}
train_dat_wi <- tree_data %>% filter(CLSTR_ID %in% ci_wil_all) %>% 
  mutate(SITE_IDENTIFIER = as.numeric(substr(CLSTR_ID, 1, 7))) %>% data.table()

tree1 <- tree_data %>%
  filter(CLSTR_ID %in% ci_wil, DAM_NUM == 1) %>%
  select(SITE_IDENTIFIER, CLSTR_ID, MEAS_YR, TREE_NO, SPECIES, DBH, HEIGHT)

tree2 <- tree_data %>%
  filter(CLSTR_ID %in% ci_wil_2, DAM_NUM == 1) %>%
  select(SITE_IDENTIFIER, CLSTR_ID, MEAS_YR, TREE_NO, SPECIES, DBH, HEIGHT) 

tree_last <- tree_data %>%
  filter(CLSTR_ID %in% ci_wil_last[ci_wil_last$max_meas > 1,]$CLSTR_ID, DAM_NUM == 1) %>%
  select(SITE_IDENTIFIER, CLSTR_ID, MEAS_YR, TREE_NO, SPECIES, DBH, HEIGHT) 

tree_wi_comb <- tree1 %>%
  left_join(tree2, by = c('SITE_IDENTIFIER', 'TREE_NO'), suffix = c("", "_2")) %>%
  left_join(tree_last, by = c('SITE_IDENTIFIER', 'TREE_NO'), suffix = c("", "_3")) 

tree_wi_comb <- tree_wi_comb %>%
  left_join(sample_data %>% select(CLSTR_ID, BEC_ZONE), by = 'CLSTR_ID') %>% data.table()

tree_wi_comb <- tree_wi_comb %>%
  rowwise() %>%
  mutate(DBHG_meas2 = DBH_2 - DBH,
         HTG_meas2 = HEIGHT_2 - HEIGHT,
         DBHG_meas2_annual = DBHG_meas2/(MEAS_YR_2 - MEAS_YR),
         HTG_meas2_annual = HTG_meas2/(MEAS_YR_2 - MEAS_YR),
         DBHG_meas3 = DBH_3 - DBH,
         HTG_meas3 = HEIGHT_3 - HEIGHT,
         DBHG_meas3_annual = DBHG_meas3/(MEAS_YR_3 - MEAS_YR),
         HTG_meas3_annual = HTG_meas3/(MEAS_YR_3 - MEAS_YR)         
         ) %>% data.table()

tree_wi_comb <- tree_wi_comb %>%
  left_join(sample_data %>% select(CLSTR_ID, PROJ_AGE_ADJ), by = 'CLSTR_ID') %>%
  left_join(sample_data %>% select(CLSTR_ID, PROJ_AGE_ADJ_2 = PROJ_AGE_ADJ), by = c('CLSTR_ID_2' = 'CLSTR_ID')) %>%
  left_join(sample_data %>% select(CLSTR_ID, PROJ_AGE_ADJ_3 = PROJ_AGE_ADJ), by = c('CLSTR_ID_3' = 'CLSTR_ID')) %>%
  data.table()

```



**3. Estimate HTVG based on the calculated height growth.**

- Based on height growth up to 1) the second and 2) the last visit
- Surrogate the current height-based HTVG estimation
- Garcia's method from `FAIBCompiler` <https://github.com/bcgov/FAIBCompiler/blob/master/R/prepareTASSInputs.R>


```{r}
# 1. based on second measure
htvgdat_1 <- data.frame()

for (k in unique(tree_wi_comb$CLSTR_ID)){
  
  clstr_id_1 <- k
  clstr_id_2 <- unique(tree_wi_comb[CLSTR_ID == k, ]$CLSTR_ID_2)
  clstr_id_3 <- unique(tree_wi_comb[CLSTR_ID == k, ]$CLSTR_ID_3)
  
  if (!is.na(clstr_id_2)){
    
  htvg0 <- tree_wi_comb %>% filter(CLSTR_ID == clstr_id_1)
  
  ## this method is suggested by Rene to select 6 tallest -> 6 most grown trees in main plot
  ## if can not find 6 trees, adding 1 tallest tree from subplot
  htvg0[DBH >= 9, plot_source := "Main"]
  htvg0[is.na(plot_source), plot_source := "Subplot"]
  
  htvg0 <- htvg0[order(CLSTR_ID, SPECIES, plot_source, -HTG_meas2_annual),]  # order by HTG
  htvg0[, sizeorder := 1:length(TREE_NO),
        by = c("CLSTR_ID", "SPECIES", "plot_source")]
  ## select 6 trees from main plot
  htvg_main <- htvg0[plot_source == "Main" &
                       sizeorder <= 6,]
  htvg_main_smry <- htvg_main[,.(num_site_ht = max(sizeorder)),
                              by = c("CLSTR_ID", "SPECIES")]
  htvg_all_clster <- unique(htvg0[,.(CLSTR_ID, SPECIES)])
  htvg_main_smry <- merge(htvg_all_clster,
                          htvg_main_smry,
                          by = c("CLSTR_ID", "SPECIES"),
                          all.x = TRUE)
  htvg_main_smry[is.na(num_site_ht), num_site_ht := 0]
  clsters_needed <- htvg_main_smry[num_site_ht < 6,]$CLSTR_ID
  # select 1 tree from subplot
  htvg_subplot <- htvg0[plot_source == "Subplot" &
                          sizeorder == 1 &
                          CLSTR_ID %in% clsters_needed,]
  htvg1 <- rbind(htvg_main, htvg_subplot)
  ## round 1: remove btop = Y, cr_cl = I and S, resitual = Y
  ## and calculate site_height
  htvg_smry_garcia <- htvg1[,.(SITE_HT_garcia = mean(HTG_meas2_annual),
                               NO_OF_TREE = length(HTG_meas2_annual)),
                            by = c("CLSTR_ID", "SPECIES")]
  htvg_smry <- merge(htvg_main_smry,
                     htvg_smry_garcia,
                     by = c("CLSTR_ID", "SPECIES"),
                     all.x = TRUE)
  
  ## round 2: ease the conditions and select again for the 0 tree summary
  htvg_smry_garcia_round2 <- htvg1[,.(SITE_HT_garcia_add = mean(HTG_meas2_annual)),
                                   by = c("CLSTR_ID", "SPECIES")]
  htvg_smry <- merge(htvg_smry,
                     htvg_smry_garcia_round2,
                     by = c("CLSTR_ID", "SPECIES"),
                     all.x = TRUE)
  htvg_smry[is.na(SITE_HT_garcia),
            SITE_HT_garcia := SITE_HT_garcia_add]
  htvg_smry <- htvg_smry[,.(CLSTR_ID, SPECIES, SITE_HT = SITE_HT_garcia)]
  
    
  treelist8 <- merge(htvg0, htvg_smry,
                     by = c("CLSTR_ID", "SPECIES"),
                     all.x = TRUE)
  ## cleanup treelist8, name as treelist9
  treelist9 <- treelist8[,HTVG_HTG2 := HTG_meas2_annual/SITE_HT]
  treelist9[HTVG_HTG2 > 1.2,
            HTVG_HTG2 := 1.2]
  
  treelist9[is.na(HTVG_HTG2), HTVG_HTG2 := 0.8]
  
  htvgdat_1 <- rbind(htvgdat_1, treelist9)
} 
}


# 2. based on last measure

tree_wi_comb_filtered2 <- tree_wi_comb %>%
  filter(!is.na(CLSTR_ID_2), !is.na(CLSTR_ID_3)) %>% data.table()


htvgdat_2 <- data.frame()

for (k in unique(tree_wi_comb_filtered2$CLSTR_ID)){
  
  clstr_id_1 <- k
  clstr_id_2 <- unique(tree_wi_comb_filtered2[CLSTR_ID == k, ]$CLSTR_ID_2)
  clstr_id_3 <- unique(tree_wi_comb_filtered2[CLSTR_ID == k, ]$CLSTR_ID_3)
  
  if (clstr_id_2 != clstr_id_3){    
    
    htvg0 <- tree_wi_comb_filtered2 %>% filter(CLSTR_ID == clstr_id_1)
    
    ## this method is suggested by Rene to select 6 tallest -> 6 most grown trees in main plot
    ## if can not find 6 trees, adding 1 tallest tree from subplot
    htvg0[DBH >= 9, plot_source := "Main"]
    htvg0[is.na(plot_source), plot_source := "Subplot"]
    
    htvg0 <- htvg0[order(CLSTR_ID, SPECIES, plot_source, -HTG_meas3_annual),]  # order by HTG
    htvg0[, sizeorder := 1:length(TREE_NO),
          by = c("CLSTR_ID", "SPECIES", "plot_source")]
    ## select 6 trees from main plot
    htvg_main <- htvg0[plot_source == "Main" &
                         sizeorder <= 6,]
    htvg_main_smry <- htvg_main[,.(num_site_ht = max(sizeorder)),
                                by = c("CLSTR_ID", "SPECIES")]
    htvg_all_clster <- unique(htvg0[,.(CLSTR_ID, SPECIES)])
    htvg_main_smry <- merge(htvg_all_clster,
                            htvg_main_smry,
                            by = c("CLSTR_ID", "SPECIES"),
                            all.x = TRUE)
    htvg_main_smry[is.na(num_site_ht), num_site_ht := 0]
    clsters_needed <- htvg_main_smry[num_site_ht < 6,]$CLSTR_ID
    # select 1 tree from subplot
    htvg_subplot <- htvg0[plot_source == "Subplot" &
                            sizeorder == 1 &
                            CLSTR_ID %in% clsters_needed,]
    htvg1 <- rbind(htvg_main, htvg_subplot)
    ## round 1: remove btop = Y, cr_cl = I and S, resitual = Y
    ## and calculate site_height
    htvg_smry_garcia <- htvg1[,.(SITE_HT_garcia = mean(HTG_meas3_annual),
                                 NO_OF_TREE = length(HTG_meas3_annual)),
                              by = c("CLSTR_ID", "SPECIES")]
    htvg_smry <- merge(htvg_main_smry,
                       htvg_smry_garcia,
                       by = c("CLSTR_ID", "SPECIES"),
                       all.x = TRUE)
    
    ## round 2: ease the conditions and select again for the 0 tree summary
    htvg_smry_garcia_round2 <- htvg1[,.(SITE_HT_garcia_add = mean(HTG_meas3_annual)),
                                     by = c("CLSTR_ID", "SPECIES")]
    htvg_smry <- merge(htvg_smry,
                       htvg_smry_garcia_round2,
                       by = c("CLSTR_ID", "SPECIES"),
                       all.x = TRUE)
    htvg_smry[is.na(SITE_HT_garcia),
              SITE_HT_garcia := SITE_HT_garcia_add]
    htvg_smry <- htvg_smry[,.(CLSTR_ID, SPECIES, SITE_HT = SITE_HT_garcia)]
       
    
    treelist8 <- merge(htvg0, htvg_smry,
                       by = c("CLSTR_ID", "SPECIES"),
                       all.x = TRUE)
    ## cleanup treelist8, name as treelist9
    treelist9 <- treelist8[,HTVG_HTG3 := HTG_meas3_annual/SITE_HT]
    treelist9[HTVG_HTG3 > 1.2,
              HTVG_HTG3 := 1.2]
    
    treelist9[is.na(HTVG_HTG3), HTVG_HTG3 := 0.8]
    
    htvgdat_2 <- rbind(htvgdat_2, treelist9)
  } 
}
```



**4. Create TASS input files using HTVG derived from height growth.**

- First, create tree list for TASS input generation using `prepareTASSInputs()` function from `FAIBCompiler` package

```{r}
library(openxlsx)
library(FAIBBase)
library(devtools)
source_url("https://raw.githubusercontent.com/bcgov/FAIBCompiler/refs/heads/master/R/prepareTASSInputs.R")

# YSM site index table
ysm_si_final <- fread("ysm_si_final_20250514.csv")

prepareTASSInputs(inputPath = "Inventory/Compilation/ismc/Archive_nonPSP_20250514",
                  outputPath = "HeightVigour/Williams/output",
                  projectName = "special",
                  clstrIDs = ci_wil_all,
                  siteIndexTable = ysm_si_final,
                  siteIndexTableSource = "Rene",
                  siteIndexMethod = "byvisit",
                  treeVigorMethod = "mainsub",
                  vigorAdjust08 = FALSE,
                  randomSeed = 123)
```

This will create two .csv files `trees_with_xy_all.csv` & `trees_without_xy_all.csv` for generating TASS inputs.

- Second, generate TASS inputs (`.in`, `.xy`) using `tassinputgen.r` function.

```{r}
rust <- 0
default_si <- 15
work_dir <- "HeightVigour/Williams/htvg_htg3"

process_data_xy(filename = trees_with_xy_all, rust_option = "norust", htvg_col = "HTVG_HTG3")
```

**5. Run TASS simulations & summarize the outputs.**

