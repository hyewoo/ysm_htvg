species_remap <- function(sc1, fiz) {
  sc1 <- toupper(sc1)  # Convert input to uppercase to handle case-insensitive comparisons
  sc3 <- rep(NA, length(sc1))
  
  species_mapping <- list(
    "A" = "At", "ABAL" = "Ba", "ABCO" = "Ba", "AC" = "Acb", "ACB" = "Acb",
    "ACT" = "Act", "AD" = "Act", "AH" = "Act", "AT" = "At", "AX" = "Acb",
    "BA" = "Ba", "BAC" = "Ba", "BAI" = "Ba", "BB" = "Bl", "BG" = "Ba",
    "BI" = "At", "BL" = "Bl", "BM" = "Ba", "BN" = "Bp", "BP" = "Bp",
    "BV" = "At", "COT" = "Act", "CT" = "Act", "CWC" = "Cwc", "CWI" = "Cwi",
    "D" = "Dr", "DG" = "Dr", "DM" = "Dr", "DR" = "Dr", "E" = "Ep",
    "EA" = "Ep", "EB" = "Ep", "EE" = "Ep", "EP" = "Ep", "ES" = "Ep",
    "EW" = "Ep", "EXP" = "Ep", "EY" = "Ep", "FDC" = "Fdc", "FDI" = "Fdi",
    "G" = "Dr", "GP" = "Dr", "GR" = "Dr", "HM" = "Hm", "HWC" = "Hwc",
    "HWI" = "Hwi", "K" = "At", "KC" = "At", "L" = "Lw", "LA" = "Lw",
    "LD" = "Lw", "LE" = "Lw", "LS" = "Lw", "LT" = "Lw", "LW" = "Lw",
    "M" = "Dr", "MB" = "Dr", "ME" = "Dr", "MN" = "Dr", "MR" = "Dr",
    "MS" = "Dr", "MV" = "Dr", "OD" = "At", "OE" = "At", "OF" = "At",
    "OG" = "At", "OH" = "At", "OI" = "At", "OJ" = "At", "OK" = "At",
    "P" = "Pli", "PA" = "Pli", "PF" = "Pli", "PJ" = "Pj", "PL" = "Pli",
    "PLC" = "Pli", "PLI" = "Pli", "PM" = "Pli", "PR" = "Pli", "PS" = "Pli",
    "PV" = "Py", "PW" = "Pw", "PXJ" = "Pli", "PY" = "Py", "Q" = "At",
    "QE" = "At", "QG" = "At", "QW" = "At", "R" = "Dr", "RA" = "Dr",
    "SA" = "Sw", "SB" = "Sb", "SE" = "Se", "SI" = "Sw", "SN" = "Sw",
    "SS" = "Ss", "SW" = "Sw", "SXB" = "Sw", "SXW" = "Sw", 
    "U" = "At", "UA" = "At", "UP" = "At", "V" = "At","VB" = "At",
    "VP" = "At", "VS" = "At", "VV" = "At", "W" = "At", "WA" = "At",
    "WB" = "At", "WD" = "At", "WI" = "At", "WP" = "At", "WS" = "At",
    "WT" = "At", "VB" = "At", "XH" = "At", "ZH" =  "At"
  )
  
  # Check for species that require special handling based on 'fiz'
  special_cases <- list(
    "B" = c("Ba", "Bl"), "BC" = c("Ba", "Bl"), "C" = c("Cwc", "Cwi"),
    "CI" = c("Cwc", "Cwi"), "CP" = c("Cwc", "Cwi"), "CW" = c("Cwc", "Cwi"),
    "CY" = c("Cwc", "Cwi"), "DF" = c("Fdc", "Fdi"), "FD" = c("Fdc", "Fdi"),
    "F" = c("Fdc", "Fdi"), "H" = c("Hwc", "Hwi"), "HW" = c("Hwc", "Hwi"),
    "HXM" = c("Hwc", "Hwi"), "IG" = c("Cwc", "Cwi"), "IS" = c("Cwc", "Cwi"),
    "J" = c("Cwc", "Cwi"), "JR" = c("Cwc", "Cwi"), "JS" = c("Cwc", "Cwi"),
    "OA" = c("Cwc", "Cwi"), "OB" = c("Cwc", "Cwi"), "OC" = c("Cwc", "Cwi"),
    "S" = c("Ss", "Sw"), "SX" = c("Ss", "Sw"), "SXE" = c("Ss", "Sw"),
    "SXL" = c("Ss", "Sw"), "SXS" = c("Ss", "Sw"), "SXX" = c("Ss", "Sw"),
    "T" = c("Hwc", "Hwi"), "TW" = c("Hwc", "Hwi"), "X" = c("Fdc", "Fdi"),
    "XC" = c("Fdc", "Fdi"), "Y" = c("Cwc", "Cwi"), "YC" = c("Cwc", "Cwi"), 
    "YP" = c("Cwc", "Cwi"), "Z" = c("Fdc", "Fdi"), "ZC" = c("Fdc", "Fdi")
  )
  
  # Direct lookup in species_mapping
  #sc3 <- ifelse(sc1 %in% names(species_mapping), unlist(species_mapping[sc1]), NA)
  sc2 <- species_mapping[match(sc1, names(species_mapping))]
  #sc2 <- species_mapping[[sc1]]
  sc2[lengths(sc2) == 0] <- NA
  sc3 <- unlist(sc2, use.names = F)
  
  # Handle special cases based on 'fiz'
  for (i in seq_along(sc1)) {
    if (sc1[i] %in% names(special_cases)) {
      sc3[i] <- ifelse(fiz[i] == "c", special_cases[[sc1[i]]][1], special_cases[[sc1[i]]][2])
    }
  }
  
  return(sc3)
}


species_site_curve <- function(sc1) {
  sc1 <- toupper(sc1)  # Convert input to uppercase to handle case-insensitivity
  
  species_map <- c(
    "SW" = "SW_GOUDIE_PLAAC",
    "PLI" = "PLI_THROWER",
    "FDC" = "FDC_BRUCEAC",
    "HWC" = "HWC_WILEYAC",
    "FDI" = "FDI_THROWERAC",
    "CWC" = "CWC_NIGH",
    "HWI" = "HWI_NIGH",
    "CWI" = "CWI_NIGH",
    "SS" = "SS_NIGH",
    "AT" = "AT_NIGH",
    "DR" = "DR_NIGH",
    "BA" = "BA_NIGH",
    "BL" = "BL_CHENAC",
    "EP" = "EP_NIGH",
    "LW" = "LW_NIGH",
    "SB" = "SB_NIGH",
    "SE" = "SE_NIGH",
    "HM" = "HM_MEANSAC",
    "BP" = "BP_CURTISAC",
    "PJ" = "PJ_HUANGAC",
    "PW" = "PW_CURTISAC",
    "PY" = "PY_NIGH",
    "ACB" = "ACB_HUANGAC",
    "ACT" = "ACT_THROWERAC"
  )
  
  sc2 <- species_map[sc1]
  return(sc2)
}

