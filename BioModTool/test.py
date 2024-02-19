#!/usr/bin/env python
# coding: utf-8
# BioModTool

#------------------------------------------------ Imports --------------------------------------------------------------------
import re
import warnings
import pandas as pd

#------------------------------------------------ Functions ------------------------------------------------------------------

#________________ Test suffix
def test_suffix_conformity(suffix):
    """
    Test if the suffix given by the user is conform.
    Expected suffix: string "^[a-zA-Z0-9_]*$" (containing only alphanumeric characters or _)
    Parameters: suffix (str).
    Raises an exception if suffix is nont conform.
    """
    if type(suffix) != str :
        raise Exception('Suffix not conform. Provided suffix must be a string.')  
    if not bool(re.match("^[a-zA-Z0-9_]*$", suffix)):
        raise Exception('Suffix not conform. Provided suffix must contain only alphanumeric characters (underscores are also allowed but not spaces).')  


def test_suffix_in_model(cobra_model,suffix,list_all_pool_id):
    """
    Test if suffix already used in model.
    Parameters: suffix (str).
    Raises an exception if suffix already in model.
    """
    suffix_in_model = False
    list_found_in_model = [] # list of found reactions
    
    for rxn in cobra_model.reactions:
        if rxn.id in list_all_pool_id:
            suffix_in_model = True
            list_found_in_model.append(rxn.id)
    if suffix_in_model:
        raise Exception("Choose another suffix. Suffix already used in model (%s)." %(list_found_in_model))


#________________ Test User Data

def test_Excel_sheets(path_to_data,expected_sheet_list):
    """
    Test Excel structure (sheet names).
    Parameters:
        - path_to_data (str)
        - list of expected sheets in the Excel (= dict_structure keys).
    Raises an exeption if test fails.
    """  
    try:
        xlsx = pd.ExcelFile(path_to_data)
    except Exception as e:
        raise Exception(e)

    excel_sheets = xlsx.sheet_names

    if 'Read_Me' in excel_sheets:
        excel_sheets.remove('Read_Me')
    
    for my_sheet in excel_sheets:
        if "new_pool" in my_sheet:
            excel_sheets.remove(my_sheet)

    if set(excel_sheets) != set(expected_sheet_list):
        raise Exception("Provided dict_structure is not consistent with data structure. Please check that all Excel sheets match a key in dict_structure (and vice-versa).")

def test_metabolite_in_model(cobra_model,met_id):
    """
    Checks if a given metabolite exist in the cobra model.
    Parameters:
        - cobra_model
        - id of the metabolite to test (str)
    Raises an exception if metabolite not in model.
    """  
    try:
        cobra_metabolite = cobra_model.metabolites.get_by_id(met_id)
    except:
        exception_metabolite_not_found = 'Metabolite ' + met_id + ' not found in model.'
        raise Exception(exception_metabolite_not_found)

def are_all_coeff_null_empty(dictdict_pool_data):
    """
    Tests if dict_data is empty or contains only 0.
    Parameter: 
        -dictdict_pool_data : Dictionary of Dictionary  
                - key = metabolite ID in model, 
                - value = dictionary with key_2 = ['Initial coefficient'...]
    Returns all_zero (True/False)
    """
    all_zero = True
    for met in dictdict_pool_data:
        if dictdict_pool_data[met]["Initial coefficient"] != 0 :
            all_zero = False

    all_zero_or_empty = (all_zero or dictdict_pool_data == {})
    return all_zero_or_empty

def is_coeff_null_empty_in_level_2(level_n,dictdict_data_level_n,dictdict_data_level_np1):
    """
    Tests if pseudo-metabolite of level n has a coefficient (not missing and not null in pseudo-reaction of level n+1) 
    Parameters:
        - level_n: (str) each level = a pool/macromolecule, with level n metabolite consumed in level n+1 (np1) reaction
                expl:
                - if level n = DAG, level n+1 = LIPIDS
                - if level n = DNA, level n+1 = BIOMASS
        - 2 dictdict_data_level n & n+1
    Returns coeff_nul_empty_in_level_np1 (True/False)
    """
    if (level_n not in dictdict_data_level_np1):
        coeff_nul_empty_in_level_np1 = True
    else:
        # test if coeff = 0 
        if dictdict_data_level_np1[level_n]['Initial coefficient'] == 0:
            coeff_nul_empty_in_level_np1 = True 
        else:
            coeff_nul_empty_in_level_np1 = False
    return coeff_nul_empty_in_level_np1


def test_continue_add_reaction(level_n,level_np1,dictdict_data_level_n,dictdict_data_level_np1):
    """
    Tests if there are missing values in data given by the user. 
    If values are missing, a message is printed or an exception is raised.
    
    Parameters:
        - level_n and level_np1 : (str) each level = a pool/macromolecule, with level n metabolite consumed in level np1 reaction
                expl:
                - level n = DAG, level n+1 = LIPIDS
                - level n = DNA, level n+1 = BIOMASS
        - 2 dictdict_data_level n & np1 
    Returns test_continue (True/False)
    """    
    
    # 1) Test if all zero or empty
    all_zero_or_empty = are_all_coeff_null_empty(dictdict_data_level_n)
    
   # 2) Test pool coeff in levelnp1 sheet
    coeff_nul_empty_in_level_np1 = is_coeff_null_empty_in_level_2(level_n,dictdict_data_level_n,dictdict_data_level_np1)
    
    # Test if reaction will be added or not:
    if (not all_zero_or_empty) and (not coeff_nul_empty_in_level_np1):
        test_continue = True
    elif (not all_zero_or_empty) and coeff_nul_empty_in_level_np1:
        message = "/!\ Reaction %s added to model but metabolite %s not consummed in Biomass reaction." %(level_n,level_n)
        print(message)
        warnings.warn(message)
        test_continue = True
    elif all_zero_or_empty and (not coeff_nul_empty_in_level_np1):
        raise Exception("Error in data. In sheet %s all metabolite coefficients are missing or null, but in %s sheet, %s coefficient is not null/missing." %(level_n,level_np1,level_n))    
    elif all_zero_or_empty and coeff_nul_empty_in_level_np1:
        test_continue = False
        
    return test_continue

#________________ Test added pseudo-reactions

def identify_BOF_by_suffix(cobra_model,suffix):
    """
    Returns the list IDs of reactions constituting the BOF identifiy by a given suffix : containing "_suffix_" in ID
    Parameters:
        - cobra_model
        - suffix (string)
    Returns: list of strings
    """
    list_BOF_rxn_id = []
    complete_suffix = "_" + suffix + "_"
    for rxn in cobra_model.reactions:
        
        if complete_suffix in rxn.id:
            list_BOF_rxn_id.append(rxn.id)
    return list_BOF_rxn_id

def test_BOF_by_suffix(cobra_model,suffix):
    """
    Tests mass and charge balance of created reactions/pseudoreactions using check_mass_balance from COBRApy librairy.
    Parameters:
        - cobra_model
        - suffix (string)
    Returns a dataframe with check_mass_balance result for each reaction.
    """
    list_BOF_rxn_id = identify_BOF_by_suffix(cobra_model,suffix)
    dict_res_BOF = {}
    
    for rxn_id in list_BOF_rxn_id:
        if rxn_id[0:3] != "EX_":
            rxn = cobra_model.reactions.get_by_id(rxn_id)
            rxn_res = rxn.check_mass_balance()
            dict_res_BOF[rxn_id]=rxn_res
    df_res_BOF = pd.DataFrame.from_dict(dict_res_BOF,orient="index")
    return(df_res_BOF)
    
