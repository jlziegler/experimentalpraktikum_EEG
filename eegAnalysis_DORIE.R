### EEG Analyse Dorie
###

setwd("~/Documents/UNI/Master/M_Sem 4/experimentalpraktikum/experimentalpraktikum_EEG")
library(eegUtils) 
library(tidyverse)

# list of vhdr files
temp = list.files(path = "./raw", pattern="*.vhdr")

# create empty list
DORlist <- list()

# save all data in list
for (i in temp){
# IMPORT DATA
  filename <- paste0(i)     # filename from temp
  DORlist[[i]] <-           # save the following dfs in a list
    assign(filename,        # assign respective filename to data
           import_raw(i,    # import data set
                      file_path = "./raw",   # path to subdirectory "raw"
                      participant_id = make.names(gsub("*.vhdr$", "", i)))  # set particpant ID
  )
}

  
# electrode locations
EEG_chl <- electrode_locations(EEG_raw, 
                                overwrite = T,
                                method = "biosemi64"
                              )

# filter
EEG_flt <- eeg_filter(EEG_chl,
                        low_freq = 0.1,
                        high_freq = 30,
                        method = "fir"
)

# re-ref
EEG_ref <- eeg_reference(EEG_flt, 
                           ref_chans = c("A1", "A2")
)

# epoch
EEG_epo <- epoch_data(EEG_ref,
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







