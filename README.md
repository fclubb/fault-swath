# Analysis of tectonics and topography along the San Andreas Fault

This repository contains code for analysing topographic characteristics and tectonics along the San Andreas Fault. It can be used to reproduce the results of Clubb and Bookhagen (in prep).

## Data 

The data required to run the code can be downloaded from XXXX.

## Installation

We recommend to use a `conda` distribution to install the code, such as `anaconda` or `miniconda`.

First, clone the repository. Then in Linux/Mac open a terminal, or in Windows open an Anaconda Powershell. 
Navigate to the repository:
```
cd fault-swath
```
Then use the `environment.yml` to create a conda environment:
```
conda env create -f environment.yml
```
## Running the code

In order to reproduce the figures from the paper we have provided a series of Jupyter notebooks. To run these, navigate to the `notebooks` folder and type the following command into your terminal or Anaconda Powershell:
```
conda activate fault-swath
jupyter notebook
```
This will open a local server on your machine in your web-browser. You can then click on the corresponding notebook for each Figure and run the code to reproduce it.
