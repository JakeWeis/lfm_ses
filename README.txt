- to run the script:
once matlab is open, add the LFM_SES_demo folder to the Matlab path, then move to the LFM_SES_demo folder
then follow the script lfm_ses_menu, which calls the scripts in the 001_process folder.

- in the first section of the lfm_ses_menu script ("EDIT SECTION" in "lfm_ses_menu", line 35):
edit the information so that it fits your working environment (ethere should only the "sep" variable to change).
The architecture of the demo folder fits with the code that follows and the paths that are defined in this section but if necessary, the user is free to change them.

- in sesf000_display_tags_names
some information is set by default in the script for usage simplicity (platform type / data folder / deployment metadata)

- in the "STARTERS" ("lfm_ses_menu", line 91):
A short script has been added that removes profiles that are out of the October-January in the specific case of post-reproduction SES profiles. The user is free to comment/change this part to fit with the processed data.

- Note: the scripts that take a bit longer to run in the processing are
sesf032_par_data_preprocess2_dark_signal
sesf036_par_data_preprocess6_functionalFit

- remarks:
1. although there is no direct theoretical obstacle to this code being valid for tags deployed in other regions than Kerguelen, the code has not been tested for data from other regions.
2. 2. If the model changes in the future or is retrained to assimilate post-2020 datasets, the only file to change is the LFMcore.mat file in LFM_SES_demo/00_data/01_seal
