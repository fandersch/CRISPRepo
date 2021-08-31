# CRISPRepo
Shiny app for CRISPR-Cas9 screen repository

global.R: establishes database connection and initialized global variables 
server.R: handles application logic upon user input
ui.R: defines graphical user interface of the application

Neccessary files (indicated with relative path) to deploy the Shiny-app:

*) ./screen.db (SQLite database file mounted with all dropout and drug-modifier screen data)

*) ./screen_facs.db (SQLite database file mounted with all facs-based screen data)

*) ./sgRNAs.db (SQLite database file mounted with all genome-wide sgRNA predictions for mouse and human, including prediction and validation analysis results)

*) dict/dict_joined.txt (Gene dictionary table for mouse/human homologs) 
# CRISPRepo
# CRISPRepo
