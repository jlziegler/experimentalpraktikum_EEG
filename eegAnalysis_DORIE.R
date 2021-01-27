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
  DORlist[[i]] <-           # save the following elements in a list
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
  lapply(DORlist, function(x){     # apply the following to every element in list
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

ICAlist <- list()

#### BEGIN LOOP 3 ####
# ICA
ICAlist <- 
  lapply(DORlist, function(x){
           run_ICA(x)
    }
    )

#### END LOOP 3 ####

# add prefix "ICA_" to ICA data
names(ICAlist) <- paste("ICA_", names(ICAlist), sep = "")

#### BEGIN LOOP 4 ####
# ICA components -> NOT WORKING PROPERLY!
DORlist <-
  mapply(function(x, y){
# ICA EOG
    ICA_comp_eog <- ar_eogcor(decomp = y, 
                              data = x, 
                              HEOG = c("HEOGli", "HEOGre"), 
                              VEOG = c("VEOGo", "VEOGu"),
                              plot = F)
# ICA MUSCLE
    ICA_comp_mus <- ar_acf(y,
                           plot = F)
# ICA CHANNEL
    ICA_comp_chn <- ar_chanfoc(y,
                               plot = F)
# ICA TRIAL
    ICA_comp_trl <- ar_trialfoc(y,
                                plot = F)
# COMBINE COMPONENTS    
    ICA_comp_rem <- c(ICA_comp_eog, ICA_comp_mus, ICA_comp_chn, ICA_comp_trl)
# REMOVE COMPONENTS
    x <- apply_ica(data = x, 
                   decomp = y, 
                   comps = ICA_comp_rem)
  }, DORlist, ICAlist)

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


# Combine Data

# Plot Data







