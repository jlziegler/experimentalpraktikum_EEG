### EEG Analyse Dorie
###

setwd("~/Documents/UNI/Master/M_Sem 4/experimentalpraktikum/experimentalpraktikum_EEG")
library(eegUtils) 
library(tidyverse)

# list of vhdr files
temp = list.files(path = "./raw", pattern="*.vhdr")

# create empty list
DORlist <- list()

#### BEGIN LOOP 1 ####
# IMPORT DATA
# & save all data in list
for (i in temp){
  filename <- paste0(i)     # filename from temp
  DORlist[[i]] <-           # save the following dfs in a list
    assign(filename,        # assign respective filename to data
           import_raw(i,    # import data set
                      file_path = "./raw",   # path to subdirectory "raw"
                      participant_id = make.names(gsub("*.vhdr$", "", i)))  # set particpant ID
  )
  rm(list = ls(pattern = "vhdr$"))
}

#### END LOOP 1 ####

#### BEGIN LOOP 2 ####
# Preprocessing
DORlist <-      # overwrite list
  lapply(DORlist, function(x){     # apply the following to every df in list
# ELECTRODE LOCATIONS
   x <- electrode_locations(x,
                           overwrite = T,
                           method = "biosemi64")
# FILTER
   x <- eeg_filter(x,
                   low_freq = 0.1,
                   high_freq = 30,
                   method = "fir")
# RE-REF
   x <- eeg_reference(x, 
                      ref_chans = c("A1", "A2"))
# EPOCH
   x <- epoch_data(x, 
                   events = c("S111",
                              "S112",
                              "S212",
                              "S211",
                              "S121",
                              "S123",
                              "S223",
                              "S221"), 
                   epoch_labels = c("oeStan",    # o-e (o-Standard)
                                    "eoStan",    
                                    "oeDev",
                                    "eoDev",     # e-o (o-Deviant)
                                    "ouStan",
                                    "uoStan",
                                    "ouDev",
                                    "uoDev"),
                   time_lim = c(-.2, .8),
                   baseline = c(-.2, 0))
}
)

#### END LOOP 2 ####

#### BEGIN LOOP 3 ####

DORlist2 <- lapply(DORlist, function(x){ica <- run_ICA(x)})
### this currently runs ICA, saves them in a new list, names them like the data, 
### maybe something like this will work:
# paste0("ICA_", gsub("*.vhdr$", "", names(x))) <- run_ICA(x)

# Run ICA
EEG_ICA <- run_ICA(EEG_epo)

# ICA EOG
ICA_comp_eog <- ar_eogcor(decomp = EEG_ICA, 
                          data = EEG_epo, 
                          HEOG = c("HEOGli", "HEOGre"), 
                          VEOG = c("VEOGo", "VEOGu")
)

# ICA Muscle
ICA_comp_mus <- ar_acf(EEG_ICA)

# ICA Channel
ICA_comp_chn <- ar_chanfoc(EEG_ICA)

# ICA Trial
ICA_comp_trl <- ar_trialfoc(EEG_ICA)

# Combine Components
ICA_comp_rem <- c(ICA_comp_eog, ICA_comp_mus, ICA_comp_chn, ICA_comp_trl)

# Remove Components
EEG_cln <- apply_ica(data = dor02_epo, decomp = dor02_ICA, comps = ICA_comp_rem)


##### END LOOP #####


# Combine Data

# Plot Data







