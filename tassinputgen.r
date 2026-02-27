## Generate TASS inputs (.in, .xy) from the tree list created from prepareTASSInputs() in FAIBCompiler
## Also create batch files to run TASS

process_data_xy <- function(filename, rust_option, htvg_col = "HTVG") {
  
  # Open spatial tree list
  infile <- if ("character" %in% class(filename)) {
    fread(filename)
  } else filename
  
  rust <- if (is.null(rust_option) || is.na(rust_option) || rust_option == "norust") 0 else 1
  
  nci <- length(unique(infile$CLSTR_ID))
  
  gocid <- c()
  residfilelist <- c()
  dedfilelist <- c()
  xyfilelist <- c()
  infilelist <- c()
  
  # Write .xy and .res file by CLSTR_ID
  for (k in 1:nci){
    
    clstr_id <- unique(infile$CLSTR_ID)[k]
    
    cat(sprintf("working for %d out of %d sample measurements.\n", k, nci))
    
    xydat <- infile[infile$CLSTR_ID == clstr_id, ]
    
    tass_exe <- if ("TASS_VERSION" %in% names(xydat) & !is.na(unique(xydat$TASS_VERSION))) {
      unique(xydat$TASS_VERSION) 
    } else {
      ba_dat <- xydat %>%
        group_by(SP) %>%
        summarise(ba = sum(DBH^2* 0.00007854)) %>%
        ungroup() %>%
        mutate(ba_prop = ba/sum(ba)*100)
      
      ifelse(ba_dat %>% 
               filter(SP %in% c("Pli", "Se", "Sw", "Ss")) %>% 
               summarise(ba = sum(ba_prop)) %>%
               pull(ba) >= 80, "TASS3", "TASS2")
    }
    
    tass_exe <- ifelse(is.na(tass_exe) | tass_exe == "", "TASS2", tass_exe)
    
    # Process BEC zone
    bec <- unique(xydat$BEC) 
    
    bec <- fcase(bec == "CMA", "MH",
                 bec == "IMA", "ESSF",
                 bec == "BAFA", "ESSF",
                 bec == "BG", "PP",
                 default = bec)
    
    if (!is.null(bec)){
      
      if (!(bec %in% c("BWBS", "CWH", "ESSF", "ICH", "IDF", "MH", 
                       "MS", "AT", "CDF", "PP", "SWB", "SBPS", "SBS"))){
        
        err_bec <- sprintf("Unknown BEC zone code %s in %s:", bec, clstr_id)
      }
    }
    
    # Determine region based on BEC zone
    region <- ifelse(bec %in% c('CWH','CDF','MH'), 'c', 'i')
    
    
    age <- unique(xydat$AGE) 
    
    
    agent_ht <- 0
    
 
    newid <- xydat$TREE_NO
    
    
    # Process residual trees
    RESID_new <- ifelse(xydat$RESID == "REG", "N", ifelse(xydat$RESID == "R", "Y", NA))
    
          
    err_plot <- ifelse(sum(xydat$PLOT <= 0) > 0, sprintf("plot size 0 in %s:", clstr_id), NA)
    
    err_age <- ifelse(age <= 0, sprintf("tree age 0 in %s:", clstr_id), NA)
    
    next10 <- as.integer(floor((age + 10) / 10) * 10)
    
    
    # Process damage information
    damage_agent <- xydat[,.(DAM_AGNA, DAM_AGNB, DAM_AGNC, DAM_AGND, DAM_AGNE)]
    damage_agent[is.na(damage_agent)] <- ""
    
    damage_ht_class <- xydat[,.(HTC_DAM_A, HTC_DAM_B, HTC_DAM_C, HTC_DAM_D, HTC_DAM_E)]
    damage_ht_class[is.na(damage_ht_class)] <- 0
    
    damage_encircle <- xydat[,.(ECC_DAM_A, ECC_DAM_B, ECC_DAM_C, ECC_DAM_D, ECC_DAM_E)]
    damage_encircle[is.na(damage_encircle)] <- 0
    
    comandra_ht <- rep(0, nrow(xydat))
    gall_ht <- matrix(rep(0, nrow(xydat)*5), ncol = 5, byrow = T)
    gall_encircle <- matrix(rep(0, nrow(xydat)*5), ncol = 5, byrow = T)
    
    if (rust == 1){
      
      map1 <- c("0" = 0.25, "1" = 1.0, "2" = 2.0,
                "3" = 3.0, "4" = 4.0, "5" = 5.0,
                "6" = 6.0, "7" = 7.0, "8" = 8.0, "9" = 9.0)
      
      map2 <- c("0" = 2.5, "1" = 10, "2" = 20,
                "3" = 30, "4" = 40, "5" = 50,
                "6" = 60, "7" = 70, "8" = 80, "9" = 90)
      
      for (kk in 1:5) {
        
        comandra_ht <- ifelse(damage_agent[[kk]] == "DSC",
                              map1[as.character(damage_ht_class[[kk]])],
                              comandra_ht)
        
        gall_ht[, kk] <- ifelse(damage_agent[[kk]] == "DSG",
                               map1[as.character(damage_ht_class[[kk]])],
                               gall_ht[, kk])
        
        gall_encircle[, kk] <- ifelse(damage_agent[[kk]] == "DSG",
                                     map2[as.character(damage_ht_class[[kk]])],
                                     gall_encircle[, kk])
        
      }
    }
    
    
    # Live or dead trees
    xydat$L <- if ("L" %in% names(xydat)) xydat$L else "L"
    
    # Remapping species
    xydat$SP <- ifelse(xydat$SP == "HX", "HW", xydat$SP)
    #spec_okay <- species_remap(xydat$SP, region) 
    spec_okay <- species_remap(xydat$SP, rep(region, length(xydat$SP))) 
    
    
    # Site index
    spsi <- unique(data.table(spec_okay, L = xydat$L, site = xydat$SITE))
    spsi <- spsi[L == "L"]   # only live trees
    spsi[is.na(site), ]$site = default_si
    
    if (sum(duplicated(spsi$spec_okay)) == T){
      
      spsi_check <- unique(xydat[L == "L",c('SP', 'SITE')])
      spsi_check[is.na(SITE), ]$SITE = default_si
      spsi_check$SP <- toupper(spsi_check$SP)
      spsi_check$site <- spsi_check$SITE
      
      spsi$SP <- toupper(spsi$spec_okay)
      
      temp1 <- merge(spsi, spsi_check, by = c('SP', 'site'))
      temp2 <- anti_join(spsi, spsi_check, by = c('SP', 'site'))
      
      for (ii in 1:nrow(temp2)){
        if (temp2$spec_okay[ii] %in% c('Cwc', 'Hwc', 'Fdc')){
          temp2$spec_okay2[ii] <- sub('c$', 'i', temp2$spec_okay[ii])
        } else if (temp2$spec_okay[ii] %in% c('Cwi', 'Hwi', 'Fdi')){
          temp2$spec_okay2[ii] <- sub('i$', 'c', temp2$spec_okay[ii])
        } else temp2$spec_okay2[ii] <- temp2$spec_okay[ii]
      }
      
      spsi <- rbind(temp1[,.(spec_okay, L, site)], temp2[,.(spec_okay = spec_okay2, L, site)])
      
      temp3 <- data.table(spec_okay, xydat$SITE)
      #temp3[spec_okay %in% temp2$spec_okay & V2 %in% temp2$site, spec_okay := temp2$spec_okay2]
      for (iii in 1:nrow(temp2)){
        spec1 <- temp2$spec_okay[iii]
        site1 <- temp2$site[iii]
        temp3[spec_okay == spec1 & V2 == site1, spec_okay := temp2$spec_okay2[iii]]
      }
      
      spec_okay <- temp3$spec_okay
    }
    
    spec_site <- spsi$site
    sp <- spsi$spec_okay
    spec_merch <- ifelse(spsi$spec_okay == "Pli", 12.5, 17.5)
    
    spsi_okay <- species_remap(sp, region) 
    spec_dom_okay <- ifelse(spsi_okay == "Pli", "Pli", spsi_okay)
    spec_dom_okay <- ifelse(spsi_okay == "Sw", "Sw", spsi_okay)
    
    curve_name <- species_site_curve(spsi_okay)
    
    num_spec <- length(unique(spsi_okay))
    
    origin <- "Natural"
    spec_origin <- ifelse(origin == "Natural", 1, 0)
    
    plot <- 1/xydat$PLOT
    
  

    ## Select HTVG column
    htvg_vals <- xydat[[htvg_col]]
    
    
    xy <- data.table(tass_exe, clstr_id, xydat$L, newid, 
                     spec_okay, -xydat$y, xydat$x, xydat$AGE, htvg_vals, 
                     xydat$HT, xydat$HTC, xydat$DBH, comandra_ht, 
                     gall_ht[,1], gall_encircle[,1],
                     gall_ht[,2], gall_encircle[,2],
                     gall_ht[,3], gall_encircle[,3],
                     gall_ht[,4], gall_encircle[,4],
                     gall_ht[,5], gall_encircle[,5],
                     RESID_new)
    
    names(xy) <- c('TASS_VERSION', 'CLSTR_ID', 'L', 'newid', 'SP', 'y', 'x', 'AGE', 
                   'HTVG', 'HT', 'HTC', 'DBH', 'comandra_ht',
                   'gall_ht1', 'gall_encircle1', 'gall_ht2', 'gall_encircle2', 
                   'gall_ht3', 'gall_encircle3', 'gall_ht4', 'gall_encircle4', 
                   'gall_ht5', 'gall_encircle5', 'RESID')
    
 
    
    
    if (any(xy$RESID == "Y")){
      
      resid <- xy[xy$RESID == "Y",]$newid
      resfilename <- sprintf("%s/%s.res", work_dir, clstr_id)
      fwrite(data.table(resid), resfilename, col.names = FALSE, quote = FALSE)
      residfilelist <- c(residfilelist, sprintf("%s.res", clstr_id))
      res_exist <- 1
      
    } else res_exist <- 0
    
    
    if (any(xy$L == "D")){
      
      dead <- paste0("IDCODE ", xy[xy$L == "D",]$newid, " KILL")
      dead <- c("TASS External Process File v1.0", dead)
      dedfilename <- sprintf("%s/%s.ded", work_dir, clstr_id)
      fwrite(data.table(dead), dedfilename, col.names = FALSE, quote = FALSE)
      dedfilelist <- c(dedfilelist, sprintf("%s.ded", clstr_id))
      ded_exist <- 1
      
    } else ded_exist <- 0
    
    
    if (tass_exe == "TASS3"){
      
      xy1 <- xy[,4:23]
      
      xyfilename <- sprintf("%s/%s.xy", work_dir, clstr_id)
      
      write.table(format(xy1, nsmall = 6, trim = T), sep = " ", file = xyfilename, 
                  quote = F, col.names = FALSE, row.names = F)
      
      
    } else {
      
      #num_spec <- length(unique(xy$SP))
      
      for (j in 1:num_spec){
        
        spec_ea <- unique(xy$SP)[j]
        
        xy1 <- xy[xy$SP == spec_ea, c(4,6:23)]
        
        xyfilename <- sprintf("%s/%s%s.xy", work_dir, clstr_id, spec_ea)
        
        write.table(format(xy1, nsmall = 6, trim = T), sep = " ", file = xyfilename, 
                    quote = F, col.names = FALSE, row.names = F)
        
      }
      
    }
    
    xyfilelist <- c(xyfilelist, sprintf("%s*.xy", clstr_id))
    
    
    width <- 100
    
    rustline <- if(rust == 1) "BASTARD_TOADFLAX 0.08" else NULL
    
    if (num_spec <= 0) {
      stop(paste("Problem: num_spec=", num_spec, " clstr_id=", clstr_id, sep=""))
    }
    
    
    
    ## Species lines
    if (tass_exe == "TASS3") {
      
      specline <- c()
      
      for (ii in 1:num_spec){
        
        specline_ea <- c(paste0("SPECIES ", spsi_okay[ii], " USING ", 
                                ifelse(spsi_okay[ii] == "Pli", "Pli", "Sw")),
                         sprintf("  MIN_DBH %f", spec_merch[ii]),
                         sprintf("  SITE_INDEX %f", spec_site[ii]),
                         "")
        
        specline <- c(specline, specline_ea)
      }
      
      specline_all <- c(specline,
                        "",
                        paste0("START_UP ORIGIN ", toupper(origin),
                               " FILE ", clstr_id, 
                               ".xy ID_CODE SPECIES ROW COLUMN AGE HT_VIGOR HEIGHT HT_TO_CROWN DBH COMANDRA_HT GALL_HT GALL_ENCIRCLE GALL_HT GALL_ENCIRCLE GALL_HT GALL_ENCIRCLE GALL_HT GALL_ENCIRCLE GALL_HT GALL_ENCIRCLE END END"),
                        ""
      )
      
    } else {
      
      specline <- c()
      
      for (ii in 1:num_spec){
        
        specline_ea <- c(paste0("SPECIES ", spsi_okay[ii]),
                         paste0("  START_UP ORIGIN ", toupper(origin),
                                " FILE ", clstr_id, spsi_okay[ii],
                                ".xy ID_CODE SPECIES ROW COLUMN AGE HT_VIGOR HEIGHT HT_TO_CROWN DBH COMANDRA_HT GALL_HT GALL_ENCIRCLE GALL_HT GALL_ENCIRCLE GALL_HT GALL_ENCIRCLE GALL_HT GALL_ENCIRCLE GALL_HT GALL_ENCIRCLE END END"),
                         sprintf("  MERCH_SPECS TOP_DIB 10 MIN_DBH %f STUMP 0.3", spec_merch[ii]),
                         paste0("  SITE_CURVE ", curve_name[ii]),
                         sprintf("  SITE_INDEX %f", spec_site[ii]),
                         ""
        )
        
        specline <- c(specline, specline_ea)
        
        specline_all <- specline
      }
      
    }
    
    
    ## Year lines
    #if (tass_exe == "TASS3") {
    
    summline <- c()
    
    for (l in 1:num_spec){
      
      summ_ea <- c(#"  CUSTOM_STAND stand %fressumm.txt",
        paste0("  SELECT SPECIES ", spsi_okay[l], " END"),
        paste0("    CUSTOM_STAND stand %fnewsumm", spsi_okay[l], ".txt"),
        "  END"
      )
      
      summline <- c(summline, summ_ea)
      
    }
    
    summline_all <- c("  CUSTOM_STAND stand %fnewsumm.txt",
                      summline)
    
    if (res_exist == 1){
      
      resline <- c()
      
      for (m in 1:num_spec){
        
        resline_ea <- c(#"  CUSTOM_STAND stand %fressumm.txt",
          paste0("    SELECT SPECIES ", spsi_okay[m], " END"),
          paste0("      CUSTOM_STAND stand %fressumm", spsi_okay[m], ".txt"),
          "    END"
        )
        
        resline <- c(resline, resline_ea)
        
      }
      
      resline_all <- c("  SELECT LIST %f.res",
                       "    CUSTOM_STAND stand %fressumm.txt",
                       resline,
                       "  END",
                       "  SELECT NOT LIST %f.res",
                       summline_all,
                       "END")
      
    } else resline_all <- NULL
    
    if (is.null(resline_all)){
      resline_all <- summline_all
    }
    
    if (ded_exist == 1){
      
      #dedline <- c("  PROCESS %f.ded NON-COMMERCIAL %f.tmp",
      #             resline_all, summline_all)
      dedline <- c("  PROCESS %f.ded NON-COMMERCIAL %f.tmp",
                   resline_all)
      
    } else dedline <- NULL
    
    yearline <- c(paste0("WHEN YEAR = ", age),
                  "  CUSTOM_TREE trees %f.tre",
                  resline_all,
                  #summline_all,
                  dedline,
                  "END",
                  ""
    )
    
    #} else{
    
    #}
    ## Next line
    #if (tass_exe == "TASS3") {
    
    #} else {
    nextline <- c(paste0("WHEN YEAR = ", next10, " REPEAT 10 UNTIL 200"),
                  resline_all,
                  #summline_all, 
                  "END")
    #}
    treeline <- c(
      paste0("WHEN YEAR = ", age, " REPEAT 1 UNTIL 70"),
      "  CUSTOM_TREE trees %f_%4y.tre",
      "END"
    )
    
    
    temp <- c("!TASS Input File v02.10",
              paste0("FILESTUB ", clstr_id),
              "GRAPHICS SIZE 600 600 CROWN_MAP END",
              "UNITS METRIC",
              paste0("DESC clstr_id: ", clstr_id),
              "GRID 0.2",
              sprintf("PLOT %f %f NO_WALLS", width, width),
              "FORCE_AGE 1",
              "ANNUAL_NODES ON",
              "!RANDOM VIGOR TIME",
              "!RANDOM KILL TIME",
              paste0("BEC_ZONE ", bec),
              rustline,
              paste0("START_YEAR ", age),
              "",
              specline_all,
              #"",
              "DEF CUSTOM_STAND stand HEADERS",
              "  YEAR xxxxx",
              "  STAND_AGE xxxx",
              "  NUM_TREES xxxxxxx",
              "  BASAL_AREA xxxx.xx", 
              "  VOLUME xxxxx",
              "  VOLUME xxxxx TOP_DIB 10 MIN_DBH SPECIES STUMP 0.30",
              "  QMEAN_DBH xxxx.xx", 
              "  TOP_HEIGHT98 xxx.xx",
              "END",
              "",
              "DEF CUSTOM_TREE trees HEADERS",
              "  YEAR xxxxx",
              "  STAND_AGE xxxx",
              "  ID_CODE xxxxxxxxxxxx",
              "  SPECIES xxxx", 
              "  DBH xxxx.xx",
              "  HEIGHT xxx.xx",
              "  HT_TO_BASE xxx.xx",
              "  VOLUME xx.xxxxx",
              "  VOLUME xx.xxxxx TOP_DIB 10 MIN_DBH 12.5 STUMP 0.30",
              "  VOLUME xx.xxxxx TOP_DIB 10 MIN_DBH 17.5 STUMP 0.30",
              "END",
              "",
              #paste0("WHEN YEAR = ", age),
              yearline,
              nextline,
              "",
              treeline
    )
    
    inputfile <- as.list(temp)
    
    infilename <- sprintf("%s/%s.in", work_dir, clstr_id)
    
    write.table(inputfile, sep='\n', file = infilename, 
                quote = F, col.names = FALSE, row.names = F)
    
    infilelist <- c(infilelist, sprintf("%s.in", clstr_id))
    
    
    if (sum(xy$L == "L") == 0){
      cat(sprintf("No live species for %s:", clstr_id))
    } else {
      #gocid <- c(gocid, paste0('start "" /b doit.bat ', tass_exe, ' ', clstr_id))
      gocid <- c(gocid, paste0("call doit ", tass_exe, " ", clstr_id))
      #gocid <- c(gocid, paste0('"', tass_exe, ' ', clstr_id,'"'))
    }
    
  }
  
  fwrite(data.table(residfilelist), sprintf("%s/list_At_InputFilesRES.txt", work_dir), 
         col.names = FALSE, quote = FALSE, row.names = F)
  fwrite(data.table(dedfilelist), sprintf("%s/list_At_InputFilesDED.txt", work_dir), 
         col.names = FALSE, quote = FALSE, row.names = F)
  fwrite(data.table(infilelist), sprintf("%s/list_At_InputFilesIN.txt", work_dir), 
         col.names = FALSE, quote = FALSE, row.names = F)
  fwrite(data.table(xyfilelist), sprintf("%s/list_At_InputFilesXY.txt", work_dir), 
         col.names = FALSE, quote = FALSE, row.names = F)
  
  fwrite(data.table(gocid), sprintf("%s/go.bat", work_dir), 
         col.names = FALSE, quote = FALSE, row.names = F)
  
  for (i in seq_along(gofilelist)) {
    a <- data.table(gofilelist[[i]])
    gofilename <- file.path(paste0(work_dir,"/go", i,".bat"))
    write.table(a, file = gofilename, #sep = "", 
                row.names = FALSE, col.names = FALSE,
                quote = FALSE, append = FALSE)
  }
  
  
  doit <- c('@echo off',
            'if exist "%2newsumm.txt" goto end',
            'del "%2." >nul 2>&1',
            'if exist "%2" goto end',
            'echo a > "%2"',
            '"%1" "%2.in"',
            ':end')
  
  write.table(as.list(doit), sep='\n', file = sprintf("%s/doit.bat", work_dir), 
              quote = F, col.names = FALSE, row.names = F)
  
  
}



