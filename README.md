# scRNA-seq data analysis 🧫
#### Alberti A., Ianniciello S., Magos P.M., Testa I., Tolloso M.

This repository contains the scripts used to perform the data analysis presented in the report [“Comparison between COTAN and state of the art libraries for the analysis of scRNA-seq data”](link) for the “Computational Health Laboratory” course at the University of Pisa (a.y. 2022/2023). 

## How to run

We assume you are working on unix-like systems (Linux, MacOS, Windows Subsystem for Linux).

Colone the repository
```bash
git clone https://github.com/arcotan/single-cell-data-analysis.git
```
Move to the repository folder
```bash
cd single-cell-data-analysis
```
Download the datasets with the following command:
```bash 
wget "https://unipiit-my.sharepoint.com/:u:/g/personal/m_tolloso_studenti_unipi_it/ESX8nFEDsHdMo5_mBAiAgBIB6aDBPwyVSNqgqnDsB2__RQ?e=LhIna1&download=1" -O dataset.zip &&
wget "https://unipiit-my.sharepoint.com/:u:/g/personal/m_tolloso_studenti_unipi_it/ERFzNbWhOq1EkDz5ntlopOIB0VZ8kxKBssaYZTxwVmPqMA?e=mU5RE6&download=1" -O filtered.zip &&
wget "https://unipiit-my.sharepoint.com/:u:/g/personal/m_tolloso_studenti_unipi_it/EUsdmgrR3K1KuNULjP5o6ksBSAU0K0T6-yVvceqpOG5RsA?e=8UdvHK&download=1" -O results.zip && 
unzip dataset.zip &&
unzip filtered.zip &&
unzip results.zip &&
rm dataset.zip &&
rm filtered.zip &&
rm results.zip
```

(All the results are precomputed)
You can still run all the analysis code, executing the following command:
```bash 
chmod +x run.sh &&
./run.sh
```
