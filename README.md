Code and data to support "Dynamics of respired dissolved organic carbon in a stream"

Main code files are (run in this order)
"data process"  Only need to run first few lines. the rest of the processesing is saved to a file
"modeling_rates_stan"  Runs most of the mathematical model of DOC uptake and DIC release
"Creston_metabolism"  Both metabolisma nd CO2 flux is in here.  
"CO2_model"  This file simply models what CO2 would be based soley on metabolism.
"Plots"  plotting routines.

Data to run this code can be this code has been included in this Github repository and is published at: https://doi.org/10.4211/hs.988f0d0aa46249b2b654145cf5fbf895

NOTE: If you are using this Github code with only the data published to Hydroshare, you will need to uncomment lines in "Creston_metabolism.R" to write the "blaine_params.csv" file. Additionally, be sure to generate the "data_doc.csv" intermediate file in the "data_process.R" script to ensure subsequent modeling and plotting scripts run smoothly. These intermediate files are included in the Github repository for ease of use.
