#!/usr/bin/env python
# coding: utf-8
# BioModTool

#------------------------------------------------ Imports --------------------------------------------------------------------
from pathlib import Path
import numpy as np
from os.path import join as __join
import pandas as pd

import cobra
import BioModTool.load
import BioModTool.test
import BioModTool.calculate
import BioModTool.update
import BioModTool.save

#------------------------------------------------ Functions ------------------------------------------------------------------

def define_reaction_id(pool_id,suffix,user_compartment):
    """
    Formats complete tag to add to pseudo-metabolites and reactions IDs.
    tag = pool_id + "_" + suffix + "_" + user_compartment
    Parameters: 3 str (pool_id, suffix and user_compartment)
    Returns complete tag (str).
    """
    rxn_id = pool_id + "_" + suffix + "_" + user_compartment
    return rxn_id

def get_ids_by_level(dict_structure,level):
    """
    Identify all pseudo-reactions of a given level.
    levels = "level_1", "level_2", "level_3"
    Parameters: 
    - dict_structure: give the structure of BOF.
        - keys = pseudo-reactions IDs (pool_id) (str)
        - value = level (str)
    Returns a list containing IDs of all pseudo-reactions of the level.
    """
    list_level_ids = []
    for key, value in dict_structure.items(): 
        if value == level:
            list_level_ids.append(key)
    return list_level_ids    

def add_biomass_objective_function(cobra_model, path_to_data, suffix, dict_structure, user_compartment = "c", calculate_charge = True, calculate_formula = True, saving_final_data = True):
    """
    Principal function of BioModTool. Generate and add a biomass objective function to a cobra model based on user data.
    Parameters:
        - cobra_model: cobra.model
        - path_to_data: data must be given in an Excel file with a specific format. See data_template.xlsx
        - suffix: string matching "^[a-zA-Z0-9_]*$" (containing only alphanumeric characters or _). 
                  All added pseudo-reactions and reactions and pseudo-metabolites will contain the following tag: _suffix_compartment.
        - dict_structure: (dictionary) key = pool_id (expl DNA, DAG, BIOMASS etc). Warning: must be the same as Excel sheets. For each pool_id=sheet_id one reaction and one pseudo metabolite will be added.
                                     value = level of the reaction: 'level 3' (for DAG,TAG etc), 'level_2' (for DNA, PROTEINS etc), 'level_1' (for BIOMASS)

        - user_compartment (optional, default value = "c") : cobra model's compartment where biomass reactions and pseudo metabolites will be added.
        - calculate_charge (optional, default value = True) : if True BioModTool will calculate charge of all created pseudo metabolites. Must be set to False if charge is missing for one or some metabolite(s) consumed in biomass.
        - calculate_formula (optional, default value = True) : if True BioModTool will calculate formula of all created pseudo metabolites. Must be set to False if formula is missing for one or some metabolite(s) consumed in biomass.
        - saving_final_data (optional, default value = True) : calculated coefficient are saved as "Bilan_BOF_calculations_suffix.xlsx" in data file folder.
    """
    # Deep copy of cobra model
    cobra_model_copy = cobra_model.copy()
    
    # VARIABLES
    # Dictionary of reaction and pseudo metabolite IDs.
    dict_rxn_met_id = {}
    for pool_id in dict_structure:
        dict_rxn_met_id[pool_id] = define_reaction_id(pool_id, suffix,user_compartment)
        
    # List of reaction IDs by level 
    level_1_id_list = get_ids_by_level(dict_structure, "level_1") 
    
    if len(level_1_id_list) == 1:
        list_all_pool_id = level_1_id_list
        biomass_id = level_1_id_list[0]
        if "level_2" in dict_structure.values():
            list_biomass = get_ids_by_level(dict_structure, "level_2") 
            list_all_pool_id = list_all_pool_id + list_biomass 
            list_all_pool_id_wo_biomass = list_biomass 
            if ("level_3" in dict_structure.values()) and ("level_2_lipid" in dict_structure.values()):
                list_lipids = get_ids_by_level(dict_structure, "level_3")
                list_level_2_lipid_id = get_ids_by_level(dict_structure, "level_2_lipid")
                if len(list_level_2_lipid_id) == 1:
                    level_2_lipid_id = list_level_2_lipid_id[0]
                else:
                    raise Exception("Error in dict_structure, must not contain more than level_2_lipid pseudo-metabolites.")

                list_all_pool_id = list_all_pool_id + list_lipids + list_level_2_lipid_id
                list_all_pool_id_wo_biomass = list_all_pool_id_wo_biomass + list_lipids + list_level_2_lipid_id
            elif ("level_3" not in dict_structure.values()) and ("level_2_lipid" not in dict_structure.values()):
                structure = "OK"
            else:
                raise Exception("Error in dict_structure, if you want to add level 3, you must define level_3 and level_2_lipid pseudo-metabolites.")
        else:
        	list_all_pool_id_wo_biomass = []
    else:
        raise Exception("Error in dict_structure, exactly one pseudo-metabolite must be defined as level_1.")

    biomass_suffix_id = define_reaction_id(biomass_id,suffix,user_compartment)
    # Dictionaries to save all results (to be saved as excel file)
    dictsave_df_calculs = {}
    dictsave_df_constant_metabolites = {}
    dictsave_df_reactions = {}
    dictsave_metabolites = {}
    
    # 1) LOAD  
    # Load and format data (data confirmity is also tested during this step: coeff = int or float, conformed unit)
    dictdf_raw_pool_data, dictdf_raw_constant_metabolites = BioModTool.load.load_data(path_to_data,dict_structure)
    dictdictdict_pool_data = BioModTool.load.format_data(dictdf_raw_pool_data,suffix,list_all_pool_id)
    dictdict_biomass_data = dictdictdict_pool_data[biomass_id] 
    if "level_3" in dict_structure.values():
        dictdict_lipid_data = dictdictdict_pool_data[level_2_lipid_id]
    dictdictdict_constant_metabolites_data = BioModTool.load.format_data(dictdf_raw_constant_metabolites,suffix,list_all_pool_id)

    # 2) TEST USER INPUTS
    # 2.a) Test suffix conformity
    BioModTool.test.test_suffix_conformity(suffix)
    BioModTool.test.test_suffix_in_model(cobra_model_copy,suffix,list_all_pool_id)
    
    # PART 1__________ level 2 and level 3 reactions
    # 2.b) Test Data
    # coefficients in mol/mol
    for pool_id in list_all_pool_id_wo_biomass:
        pool_suffix_id = dict_rxn_met_id[pool_id]
        pool_type = dict_structure[pool_id] 
        dictdict_pool_data = dictdictdict_pool_data[pool_id]
        dictdict_constant_metabolites = dictdictdict_constant_metabolites_data[pool_id]
        
        # 2.b.1) Test if coefficients are missing
        if pool_type == "level_3":
            test_add_rxn = BioModTool.test.test_continue_add_reaction(pool_id,level_2_lipid_id,dictdict_pool_data,dictdict_lipid_data)
        elif pool_type in ["level_2","level_2_lipid"]:
            test_add_rxn = BioModTool.test.test_continue_add_reaction(pool_id,biomass_id,dictdict_pool_data,dictdict_biomass_data)
        elif pool_type == "level_1":
            a = "No test"
        else:
            raise Exception("Error in dict_structure, value must be in ['level_1','level_2','level_2_lipid','level_3']")

        if test_add_rxn :
        # 2.b.2) Test if metabolites are in model
            for dico in [dictdict_pool_data,dictdict_constant_metabolites]:
                for met_id in dico:
                    if met_id in list_all_pool_id:
                        met_id_in_model = define_reaction_id(met_id, suffix,user_compartment)
                    else:
                        met_id_in_model = met_id.strip()

                    BioModTool.test.test_metabolite_in_model(cobra_model_copy,met_id_in_model)
                    cobra_metabolite = cobra_model_copy.metabolites.get_by_id(met_id_in_model)
                    dico[met_id]["cobra_metabolite"] = cobra_metabolite
            
            # 3) CALCULATE STOICHIOMETRIC COEFFICIENTS
            dictdict_pool_data = BioModTool.calculate.calculate_stochoimetric_coefficient_molpermol(dictdict_pool_data,calculate_formula)   

            # 4) CREATE REACTION AND PSEUDO METABOLITE
            dict_pool = BioModTool.update.create_reaction_dict(dictdict_pool_data, dictdict_constant_metabolites)
            (pool_met, pool_rxn) = BioModTool.update.generate_pool_metabolite_and_reaction(cobra_model_copy,dict_pool,pool_suffix_id,user_compartment,calculate_charge,calculate_formula)
            
            # 5) UPDATE MODEL
            cobra_model_copy = BioModTool.update.update_model(cobra_model_copy,pool_met,pool_rxn)
            
            # 6) SAVE RESULTS in dictionaries
            dictsave_df_reactions[pool_id] = pd.DataFrame.from_dict(pool_rxn.metabolites,orient="index",columns=["Stochiometric coefficient in updated model"])
            dictsave_metabolites[pool_met.id] = {"formula":pool_met.formula,"charge":pool_met.charge}
            dictsave_df_constant_metabolites[pool_id] = pd.DataFrame.from_dict(dictdict_constant_metabolites,orient="index")
            dictsave_df_calculs[pool_id] = pd.DataFrame.from_dict(dictdict_pool_data,orient="index")

    # PART 2__________ for level reactions (BIOMASS reaction)
    # Warning: Constant coeff unit in mmol/gDCW  
    dictdict_biomass_constant_metabolites = dictdictdict_constant_metabolites_data[biomass_id]
    # Checked that the metabolite is in cobra_model_copy
        
    # 2.b bis) Look for null/missing coefficients
    all_biomass_coeff_null_empty = BioModTool.test.are_all_coeff_null_empty(dictdict_biomass_data)
    if all_biomass_coeff_null_empty:
        raise Exception("No Biomass reaction added, all coefficients of %s sheet are null or empty." %(biomass_id))
    else:
        for dico in [dictdict_biomass_data,dictdict_biomass_constant_metabolites]:
            keys_to_rm = [] # remove met_id where coeff = 0
            for met_id in dico:
                if dico[met_id]['Initial coefficient'] != 0:
                    if met_id in list_all_pool_id:
                       	met_id_in_model = define_reaction_id(met_id, suffix,user_compartment)
                    else:
                        met_id_in_model = met_id.strip()

                    BioModTool.test.test_metabolite_in_model(cobra_model_copy,met_id_in_model)
                    cobra_metabolite = cobra_model_copy.metabolites.get_by_id(met_id_in_model)
                    dico[met_id]["cobra_metabolite"]=cobra_metabolite
                else:
                    keys_to_rm.append(met_id)
            for key in keys_to_rm:
                del dico[key]
                    

        # 3 bis) CALCULATE STOICHIOMETRIC COEFFICIENTS for BIOMASS

        dictdict_biomass_data = BioModTool.calculate.calculate_stochoimetric_coefficient_mmolpergDCW(dictdict_biomass_data,calculate_formula)

        # 4 bis) CREATE REACTION AND PSEUDO METABOLITE for BIOMASS
        dict_biomass = BioModTool.update.create_reaction_dict(dictdict_biomass_data, dictdict_biomass_constant_metabolites)
        (biomass_met, biomass_rxn) = BioModTool.update.generate_pool_metabolite_and_reaction(cobra_model_copy,dict_biomass,biomass_suffix_id,user_compartment,calculate_charge, calculate_formula)

        # 5 bis) UPDATE MODEL 
        # add Biomass reaction
        cobra_model_copy = BioModTool.update.update_model(cobra_model_copy,biomass_met,biomass_rxn)
        # add EX_biomass reaction
        cobra_model_copy = BioModTool.update.add_EX_biomass_reaction(cobra_model_copy,biomass_met,biomass_suffix_id)

        # 6 bis) SAVE RESULTS as dictionary
        dictsave_df_reactions[biomass_id] = pd.DataFrame.from_dict(biomass_rxn.metabolites,orient="index",columns=["Stochiometric coefficient in updated model"])
        dictsave_metabolites[biomass_met.id] = {"formula":biomass_met.formula,"charge":biomass_met.charge}
        dictsave_df_constant_metabolites[biomass_id] = pd.DataFrame.from_dict(dictdict_biomass_constant_metabolites,orient="index")
        dictsave_df_calculs[biomass_id] = pd.DataFrame.from_dict(dictdict_biomass_data,orient="index")

    # 7) Test mass and charge balance of BOF reaction(s)
    df_res_balance = BioModTool.test.test_BOF_by_suffix(cobra_model_copy,suffix)
    # 8) SAVE DATA TO EXCEL (optional)
    if saving_final_data :
        file_path = __join(Path(path_to_data).parent.absolute(),"Bilan_BOF_calculations_"+suffix+".xlsx")
        BioModTool.save.save_results(file_path,dictsave_df_calculs,dictsave_df_constant_metabolites,dictsave_df_reactions,dictsave_metabolites,df_res_balance)
    return cobra_model_copy 
