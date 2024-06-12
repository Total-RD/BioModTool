# BioModTool : Biomass Modeling Tool

### What is BioModTool ?

BioModTool is a Python library allowing easy generation of biomass objective functions for genome scale metabolic models from user data. BioModTool loads biomass composition data in the form of a structured Excel file previously completed by the user, normalizes, converts these data into model compatible unit (mmol.gDW<sup>-1</sup>), and creates a structured biomass objective to update a metabolic model. BioModTool can be run as Python command-lines but it also comes with an interface allowing its use by non-modelers. 

### Installation


Use pip to install BioModTool from [PyPI](https://pypi.org/project/BioModTool/)  


     pip install BioModTool


In case you downloaded or cloned the source code from [GitHub](https://github.com/Total-RD/BioModTool)

     pip install <path-to-BioModTool-repo> 

    

### Documentation

Documentation is available on [GitHub](https://github.com/Total-RD/BioModTool/tree/main/BioModTool_Documentation).



### Application examples

BioModTool was applied to GEM of two bacteria species _Escherichia coli_ (iML1515, [Monk et al. 2017](https://doi.org/10.1038/nbt.3956)) and _Alicyclobacillus acidocaldarius_ (CNA_Alicyclo, [Beck, Hunt et Carlson 2018](https://doi.org/10.3390/pr6050038)), and one microalga _Chlamydomonas reinhardtii_ (iRC1080, [Chang et al. 2011](https://doi.org/10.1038/msb.2011.52)). These examples illustrate the three different contexts in which BioModTool can be used:
-	Add a three level BOF in a model in which metabolic formula are available: iML1515
-	Add a two level BOF in a model in which metabolic formula are not available: CNA_Alicyclo
-	Add a one level BOF in a model in which metabolic formula are not available: iRC1080

All scripts and results of these three application examples are fully available on [GitHub](https://github.com/Total-RD/BioModTool/tree/main/Application_examples).

### Authors/Contributors

Clémence Dupont-Thibert<sup>1,2</sup>, Sylvaine Roy<sup>1</sup>, Gilles Curien<sup>1</sup>, Maxime Durot<sup>2</sup>

1 - Université Grenoble Alpes, CNRS, CEA, INRAE, Interdisciplinary Research Institute of Grenoble, IRIG-Laboratoire de Physiologie Cellulaire et Végétale, 17 Avenue des Martyrs, 38000 Grenoble, France  
2 - TotalEnergies, OneTech, Centre de Recherche de Solaize CRES, Chemin du Canal, 69360 Solaize, France  
