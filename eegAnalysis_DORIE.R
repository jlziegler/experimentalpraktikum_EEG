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
                                    "ouStan",    # o-u (o-Standard)
                                    "uoStan",    # u-o (u-Standard)
                                    "ouDev",     # o-u (u-Deviant)
                                    "uoDev"),    # u-o (o-Deviant)
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
# ICA components (using map2() from purrr package (tidyverse))
DORlist.final <-
  map2(DORlist, ICAlist, function(x, y){
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
  })

# REMOVE ICAlist
rm(ICAlist)

# CONVERT TO DATA FRAME
DOR_df_FCz <- 
  map(DORlist.final, function(x){
    x %>%
      filter(epoch_labels %in% c("uoDev", "ouStan", "uoStan", "ouDev")) %>%
      as.data.frame(long = T) %>%
      filter(electrode == "FCz")
  })

# COMBINE DATASETS
DOR_gravg_FCz <- bind_rows(DOR_df_FCz)

# Plot o Data
DOR_gravg_FCz %>%
  filter(epoch_labels %in% c("uoDev", "ouStan")) %>%
  ggplot(aes(x = time, y = amplitude)) +
    stat_summary(fun = mean, 
               geom = "line", 
               aes(colour = epoch_labels)) + 
    facet_wrap(~electrode) + # wenn mehrere elektroden
    scale_y_reverse() + 
    theme_light() +
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 0)

#ggsave

# Plot u Data
DOR_gravg_FCz %>%
  filter(epoch_labels %in% c("ouDev", "uoStan")) %>%
  ggplot(aes(x = time, y = amplitude)) +
  stat_summary(fun = mean, 
               geom = "line", 
               aes(colour = epoch_labels)) + 
  facet_wrap(~electrode) + # wenn mehrere elektroden
  scale_y_reverse() + 
  theme_light() +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0)

#ggsave

# PLOT DIFFERENCES
# mean for every time by epoch_label
test2 <-
  DOR_gravg_FCz %>%
  group_by(time, epoch_labels, electrode, participant_id) %>%
  summarise_at("amplitude", mean) %>%
  ungroup()

# diff Dev-Stan
# long to wide
test3 <-
  test2 %>%
  pivot_wider(names_from = epoch_labels,
              values_from = amplitude)
  
# differences
test4 <-
  test3 %>%
  mutate(o_diff = uoDev - ouStan,
         u_diff = ouDev - uoStan,
         .keep = "unused")

# wide to long
test5 <-
  test4 %>%
  pivot_longer(cols = c(o_diff, u_diff), 
               names_to = "epoch_diff", 
               values_to = "amplitude_diff")

# Plot differences
test5 %>%
  ggplot(aes(x = time, y = amplitude_diff)) +
  stat_summary(fun = mean, 
               geom = "line", 
               aes(colour = epoch_diff)) + 
  facet_wrap(~electrode) + # wenn mehrere elektroden
  scale_y_reverse() + 
  theme_light() +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0)

# AVERAGE DATA OVER 120-200 ms
test6 <-
  test5 %>%
    filter(time >= 0.12,  # 120 ms
          time <= 0.2)    # 200 ms

# mean by participant and epoch
test7 <-
  test6 %>% 
  group_by(electrode, participant_id, epoch_diff) %>%  # keep electrode, participant & epoch
  summarise_at("amplitude_diff", mean) %>%             # mean amplitude
  ungroup()

# REPEATED MEASURES ANOVA ON DIFFERENCES





