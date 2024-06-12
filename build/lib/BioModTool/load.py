#!/usr/bin/env python
# coding: utf-8
# BioModTool

#------------------------------------------------ Imports --------------------------------------------------------------------
import numpy as np
import pandas as pd
import os
import sys

import cobra
from BioModTool.test import test_Excel_sheets
#------------------------------------------------ Functions ------------------------------------------------------------------

def load_cobra_model(path_to_model):
    """
    Load a cobra model. Try using cobra.io.load_json_model and cobra.io.read_sbml_model
    Parameter: path to a cobra model (str). Model must be given as a json or sbml/xml file.
    Returns the cobra model.
    """
    if os.path.isfile(path_to_model): # Test if file exist
        try:
            return cobra.io.load_json_model(path_to_model) # try to charge json file
        except:
            try:
                return cobra.io.read_sbml_model(path_to_model) # try to charge sbml file
            except:
                raise Exception('Model format not compatible, provide a cobra in smbl/xml or json format.')
    else: 
        exception_directory = "No such file or directory: " + path_to_model
        raise Exception(exception_directory)


def load_user_data_variable_metabolites(path_to_data,sheet_list):
    """
    Load user 'variable' metabolites data from Table 1. Only metabolites from Table 1 are loaded by this function (not Table 2 metabolites).
    Data must be given in a Excel file with a specific format. See data_template.xlsx
            - Coefficient must be a float or int >= 0 (or empty).
            - Unit must be 'g per …' or 'mol per …'.
    
    Parameters:
        - path to the Excel file (str)
        - sheet_list: Excel file sheet names (will be used as reaction and metabolites ID).
        
    Returns:
     - a dictionary:
            - key = string sheet/pool id (expl. DNA, RNA etc) 
            - value = pandas.DataFrame with loaded data
    """
    
    # 1) Load data
    dictdf_raw_data = {}
    
    for sheet in sheet_list:
        df_data = pd.read_excel(path_to_data, sheet, index_col=0, usecols="B:D",skiprows=10) 
        df_data = df_data.dropna() # remove raws with missing metabolite id or missing coefficient
        
        # 2) Test data conformity
        # 2.a Coefficients value must be float (or np.float64 or np.int64)
        for index in df_data.index:
            if type(df_data.loc[index,"Coefficient"]) not in [np.float64,float,np.int64]:
                raise Exception("Error while reading data excel file. In '" + sheet +  "' sheet, coefficient values must be np.float64 or float.")
            if not df_data.loc[index,"Coefficient"] >= 0:
                raise Exception("Error while reading data excel file. In '" + sheet +  "' sheet, coefficient values must be >= 0.")
             
        # 2.b Unit must be 'g per …' or 'mol per …'
        for index in df_data.index:
            if df_data.loc[index,"Unit"] not in ["g per …","mol per …"]:
                raise Exception("Error while reading data excel file. Unit must be 'g per …' or 'mol per …'.")
        
        # Save loaded and tested data in the dictionnary
        dictdf_raw_data[sheet] = df_data
      
    return(dictdf_raw_data)


def load_fixed_metabolites_data(path_to_data,dict_structure): 
    """
    Load data of 'fixed' metabolites from Table 2 (metabolites produced or consumed during macromolecule/pool synthesis having a fixed coefficent (expl. ATP, H2O etc).
    Data must be given in a Excel file with a specific format. See data_template.xlsx
            - Coefficient must be a float or int (or empty).
            - Unit must be mol/mol_pool or mmol/gBIOMASS (NOT CONVERTED BY BioModTool)
    
    Parameters:
        - path to the Excel file (str)
        - dict_structure: keys = Excel file sheet names (will be used as reaction and metabolite ID)
                        value = level of macromolecule/pool ('level_1', 'level_2' of 'level_2_lipids', 'level_3')
        
    Returns:
     - a dictionary:
            - key = string sheet/pool id (expl. DNA, RNA etc) 
            - value = pandas.DataFrame with loaded data
    """
    
    
    # Load data
    dictdf_constant_metabolites = {}
   
    for sheet in dict_structure:
        df_data = pd.read_excel(path_to_data, sheet, index_col=0, usecols="I:K",skiprows=10) 
        df_data = df_data.dropna() # remove raws with missing metabolite id or missing coefficient
        df_data.columns = ["Coefficient", "Unit"]
    
        # 2) Test data conformity
        # 2.a Coefficients value must be float (or np.float64)
        for index in df_data.index:
            try:
                df_data.loc[index,"Coefficient"] = float(df_data.loc[index,"Coefficient"])
            except:
                raise Exception("Error while reading data excel file. In '" + sheet +  "' sheet, for constant metabolites: coefficient values must be np.float64 or float.")
            
        # Save data of the sheet in dict
        dictdf_constant_metabolites[sheet] = df_data
    
    return(dictdf_constant_metabolites)



def load_data(path_to_data, dict_structure):
    """
    Load user data. Function call both load_user_data_variable_metabolites and load_fixed_metabolites_data functions.
    Data must be given in a Excel file with a specific format. See data_template.xlsx
    
    Parameters:
        - path to the Excel file (str)
        - dict_structure: keys = Excel file sheet names (will be used as reaction and metabolites ID)
                        value = level of macromolecule/pool ('level_1', 'level_2' of 'level_2_lipids', 'level_3')
        
    Returns:
        - 2 dictionaries (dictdf_raw_data, dictdf_constant_metabolites):
            - key = string sheet/pool id (expl. DNA, RNA etc) 
            - value = pandas.DataFrame with loadind data     
    """
    
    expected_sheet_list = dict_structure.keys()

    # Test excel structure and dict_structure consistency
    test_Excel_sheets(path_to_data,expected_sheet_list)
    
    # Load biomass composition data
    dictdf_raw_data = load_user_data_variable_metabolites(path_to_data,expected_sheet_list)

    # Load constant metabolite data: metabolites required for pool synthesis (atp, h2o etc)
    dictdf_constant_metabolites = load_fixed_metabolites_data(path_to_data,dict_structure)

    return (dictdf_raw_data,dictdf_constant_metabolites)


def format_data(dictdf_raw_pool_data,suffix,list_all_pool_id):
    """
    Change data format: from a pandas.DataFrame to a dictonary. Additional values are added in the dictionary.
    
    Parameters:
        - data as Dictionary of DataFrame: 
                - key = sheet/pool id (expl. DNA, RNA etc) 
                - value = pandas.DataFrame with loaded data (index = metabolite ID in model, columns = ["Coefficient", "Unit"])
        - suffix (string)
        - list of all pool IDs 
        
    Returns
        - Dictionary of Dictionary of Dictionary:
            - key_1 = sheet/pool id (expl. DNA, RNA etc) 
            - value_1 = composition data (Dictionary of Dicitonnary)
                    - key_2 = metabolite ID in model, 
                    - value_2 = dictionary with keys = ["Initial coefficient","Initial unit"]
    """
    dictdictdict_formated_data = {} 
    
    for pool_id in dictdf_raw_pool_data:
        
        formated_pool = {}  # key = metabolite ID (in model)
                            # value = a dictionary {"Initial coefficient":...,
                            #                       "Initial unit":...}
        for met_id in dictdf_raw_pool_data[pool_id].index:
            met_coeff = dictdf_raw_pool_data[pool_id].loc[met_id,"Coefficient"]
            met_unit = dictdf_raw_pool_data[pool_id].loc[met_id,"Unit"]
            formated_pool[met_id] = {"Initial coefficient":met_coeff,"Initial unit":met_unit} 
            
        dictdictdict_formated_data[pool_id] = formated_pool

    return(dictdictdict_formated_data)
