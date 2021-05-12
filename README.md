# Yeaseq

💬 Yeaseq is a sequence viewer for four fission yeasts: pombe, japonicus, cryophilus, and octosporus

🐶 The main functions are to retrieve sequence and gene annotation from online database (Ensembl Biomart) or download them from NCBI, and show the basic information and colored sequence of queried gene as user defined

🐍 Yeaseq is developed with python 3.7 and packaged with pyinstaller 3.5 for Windows, MacOS, and Linux (Ubuntu 18.04). wxPython is used to generate graphical interface

## Content

+ [Run](#run)

+ [Sequence viewer for yeast](#viewer)

+ [Search an expected gene](#search)

+ [About this project](#about)

+ [Previous version: Japonicusequence](#previous)

+ [LICENSE](#license)

## <span id="run">Run</span>

🍞 The compiled executable files for Windows, MacOS, and Linux (Ubuntu 18.04) are available [here](https://github.com/Sugiyama-Lab/Yeaseq/releases) (https://github.com/Sugiyama-Lab/Yeaseq/releases)

⛏ To use this package in python, please install 'requests' and 'wxpython', and python >= 3.6 is recommended. Then run `Yeaseq-ui.py` or

```python
import wx
from yeaseq import wx_ui
app = wx.App()
frame_search = wx_ui.SearchFrame(None, title='Yeaseq')
frame_search.Show(True)
app.MainLoop()
```

## <span id="viewer">Sequence viewer for fission yeasts</span>

+ Yeaseq supports four fission yeasts, including japonicus, cryophilus, octosporus, and pombe
+ To make sure the users can get the up-to-date data of expected genes, prior information are avoided to be introduced to compiled application. To do so, two search modes are proposed, online search and offline search:
  + 📶 In online search mode, the sequence and annotaion of searched gene will be queried in Ensembl Biomart database, and the parsed information will be shown in the result window
  + 🖥️ In offline search mode, users can download the genome sequence and annotation in one click, and the parsed information will stored in one folder (yeaseq_ref_seq) of user directory. This ensures users can download the newest datasets NCBI
+ 💡 From version 0.3.6, to take the advantage that over 80% of genes in the fission yeast clade are 1:1:1:1 orthologs, a orthologous gene table was integrated into Yeaseq, and users can search orthologous genes in other three fission yeasts from one. To use this, please cite: Nicholas Rhind et.al. Comparative functional genomics of the fission yeasts. Science 332, 930–6 (2011). https://www.sciencemag.org/lookup/doi/10.1126/science.1203357

## <span id="search">Search an expected gene</span>

#### 1. Enter Ensembl ID or Entrez gene ID
1. Ensembl ID used here is also called stable gene ID, which is usually the systematic name of gene. For example:
    * SPBC11B10.09 for pombe
    * SJAG_01836 for japonicus
    * SPOG_02532 for cryophilus
    * SOCG_02930 for octosporus
2. Entrez gene ID indicates the unique digital gene ID used in NCBI Entrez database. For example:
    * 2539869 is the cdc2 gene in pombe
3. If offline search is open, Enterz gene ID will also support the gene information used in Genbank, including RNA Genbank ID, protein Genbank ID, and the gene name. For example:
    * To search the pombe gene cdc2
    * cdc2, cdk1, swo2, tws1 are of the same effect
    * NM_001356226.1 is the RNA Genbank ID of cdc2
    * NP_001342995.1 is the protein Genbank ID of cdc2
#### 2. Choose expected yeast species in the listbox

#### 3. Fill in suitable number for upstream and downstream bases

+ This design slightly alleviate the network use for text transfer and server use for querying data in online database
+ Here, a rare situation may occur when the expected sequence exceeds the reference sequence. This situation will be detected automatically and the limitation of up- or down-stream number will be changed in the sequence frame

#### 4. Use online or offline mode

+ Press the 'Search offline' button, and all searches are offline now. The search status at the lower left corner will also be changed. Please notice that for the current search status

+ In offline mode, offline data should be downloaded first before searching and this will take some time that depends on the network

+ Offline search is recommended if there is a bad network

+ The application may have no response once it starts downloading data

## <span id="about">About this project</span>

🆕 This project started from January, 2019 at Sugiyama lab

💭 The motivation was there was no suitable sequence viewer for japonicus, and now the supported yeasts have been expanded to four

ℹ️ The scripts have been tested with many genes and no error happened, but some unpredictable issues may still occur. If you meet some bugs, please create an issue. And some error information will be greatly helpful from your generous share

## <span id="previous">Previous version: Japonicusequence</span>

💾 A much simpler sequence viewer developed for japonicus only. Japonicusequence supports online search for Ensembl stable gene id and NCBI entrez gene id, and it has succinct interfaces with (maybe) smaller size

## <span id="license">LICENSE</span>

* This repository contains the source code and compiled application, they are under an **MIT** license with ultra-free use
* This repository also contains a **orthologous data table** from Nicholas Rhind et.al. in their publication in Science 332, 930–6 (2011). To use this (the orthologous gene finder tool in Yeaseq), please cite their excellent work > Comparative functional genomics of the fission yeasts
* When searching a gene
    * the online search will link to biomart in **Ensembl** and enterz in **NCBI**
    * or use the offline data which would be downloaded from the FTP site from **NCBI**
* This package used some packages, they are:
    * The **python** language
    * The python built-in packages: urllib, gzip, os, sys, json, re, xml, platform, shutil, io, setuptools
    * The third-party python packages: **requests**, **wxpython**
* The Yeaseq application is compiled by **pyinstaller** to run on Windows, MacOS, and Linux (Ubuntu)

🏳️‍🌈 Here is all, thanks for your reading

