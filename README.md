# PanCanTmmProject

Telomere maintenance mechanism (TMM) dysregulation serves as an early predictor of adjuvant therapy response and a potential therapeutic target in human cancers

## Overview

Purpose of this repository is to share analysis procedure, data and help readers or reviewers to know more detail of this work, reproduce or make use of results they are interested in.

## Repo Contents

* [code](https://github.com/ZhangYajing93/PanCanTmmProject/tree/main/code): tidy R functions and R script.
* [data](https://github.com/ZhangYajing93/PanCanTmmProject/tree/main/data): original and preprocessed data used for analysis and share.
* [result](https://github.com/ZhangYajing93/PanCanTmmProject/tree/main/result): important middle results and final results of manuscript, most of them are in form of .RData, which can be easily loaded and operated by R.

## Instructions for Use
For readers who want to obtain raw/result data, locate data file, then download it with one of following ways:

* In Github, download file by clicking either `Download` button or `Raw` button at corresponding data page

* Use linux command `wget` or `curl`, fo example, you can download TMM gene sets by

  wget `https://github.com/ZhangYajing93/PanCanTmmProject/blob/main/data/TMMgenes.RData`

Or you can download whole respository with one of following ways:

* Clone this repository with `git clone https://github.com/ZhangYajing93/PanCanTmmProject.git`

* Download whole respository by clicking `Download` button at top right of url page https://github.com/ZhangYajing93/PanCanTmmProject

## Reproduce analysis

For readers who want to reproduce analysis shown in manuscript, please [install R](https://cran.r-project.org/) in your computer.

## Test Environment
* System: **Linux**
* Software: **R v4.0.2**
