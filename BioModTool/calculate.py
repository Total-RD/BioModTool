#!/usr/bin/env python
# coding: utf-8
# BioModTool

#------------------------------------------------ Imports --------------------------------------------------------------------
import numpy as np
import pandas as pd

#------------------------------------------------ Functions ------------------------------------------------------------------

def calculate_coeff_in_mol(dict_metabo,calculate_formula):
    """
    Converts coefficient from initial unit (g/... or mol/...) to mol/...
    
    Parameters: 
    - dictonary (keys = 'Initial unit', 'Initial coefficient' and 'cobra_metabolite')
    - calculate_formula = True or False

    Returns: calculated coefficient in mol/...
    """
    
    # Get information from dictionary
    unit = dict_metabo["Initial unit"]
    initial_value = dict_metabo["Initial coefficient"]
    cobra_metabolite = dict_metabo["cobra_metabolite"]
    
    # 1) Calculate metabolite Molecular Weight (MW)
    met_MW = cobra_metabolite.formula_weight

    # 2) Convert coefficient (several options)
    if unit == "mol per …": # data in mol/something (expl. mol/mol_pool, mol/gDCW, umol/10^6cells etc)
        data_mol_per_smt = initial_value
        
    elif unit == "g per …": # data in g/something (expl. g/g_pool, g/gDCW etc)
        data_mol_per_smt = convert_gram_to_mol(initial_value,met_MW,calculate_formula)

    else:
        raise Exception("Problem with given unit, unit should be 'g per …' or 'mol per …'.")
    
    return (met_MW,data_mol_per_smt)


def calculate_coeff_in_gram(dict_metabo,calculate_formula):
    """
    Converts coefficient from initial unit (g/... or mol/...) to g/...
    
    Parameters: 
    - dictonary (keys = 'Initial unit', 'Initial coefficient' and 'cobra_metabolite')
    - calculate_formula = True or False
    Returns: calculated coefficient in g/...
    """
    # Get information from dict
    unit = dict_metabo["Initial unit"]
    initial_value = dict_metabo["Initial coefficient"]
    cobra_metabolite = dict_metabo["cobra_metabolite"]
    
    # 1) Calculate metabolite MW
    met_MW = cobra_metabolite.formula_weight

    # 2) Convert coefficient (several options)
    if unit == "mol per …": # data in mol/something (expl. mol/mol_pool, mol/gDCW, umol/10^6cells etc)
        data_g_per_smt = convert_mol_to_gram(initial_value,met_MW,calculate_formula)
        
    elif unit == "g per …": # data in g/something (expl. g/g_pool, g/gDCW etc)
        data_g_per_smt = initial_value

    else:
        raise Exception("Problem with given unit, unit should be 'g per …' or 'mol per …'.")
    
    return (met_MW,data_g_per_smt)


def convert_mol_to_gram(value_mol,met_MW,calculate_formula):
    """
    Converts a value from mole to gram using the given molecular weight (mol/... --> g/...).
    Parameters:
    - value_mol: molar value (float)
    - met_MW: molecular weight (float)
    - calculate_formula: True or False
    Returns the converted mass value (float)
    """

    if calculate_formula:
        value_g = value_mol * met_MW

        return value_g
    else:
        raise Exception("calculate_formula = False, no unit conversion (requiring molecular weight) can be performed.")



def convert_gram_to_mol(value_g,MW,calculate_formula):
    """
    Converts a value from gram to mole using the given molecular weight (g/... --> mol/...).
    Parameters:
    - value_g: mass value (float)
    - met_MW: molecular weight (float)
    - calculate_formula: True or False
    Returns the converted mole value (float)
    """
    if calculate_formula:
        value_mol = value_g / MW 
        return(value_mol)
    else:
        raise Exception("calculate_formula = False, no unit conversion (requiring molecular weight) can be performed.")


def normalize_data(coeff_met,sum_coeff): 
    """
    Normalize the coefficient of a metabolite.
    Parameters:
    - coeff_met: coefficient of the metabolite (float)
    - sum_coeff: sum of coefficients of all metabolites of Table 1 (in good unit) 
    Returns: calculated normalized coefficient.
    """
    normalized_value = coeff_met / sum_coeff
    
    return normalized_value


def calculate_pool_mol(dictdict_pool_data):
    """
    Sum all metabolites coefficients in mol/smt.
    Parameter = dictionary (key = metabolite id, value = dictionary (with key = 'Coefficient in mol/smt')).
    Returns the calculated sum.
    """ 
    sum_mol = 0
    
    for met_id in dictdict_pool_data:
        met_coeff = dictdict_pool_data[met_id]["Coefficient in mol/smt"]
        sum_mol = sum_mol + met_coeff
        
    return sum_mol


def calculate_stochoimetric_coefficient_molpermol(dictdict_pool_data,calculate_formula):
    """
    Calculated stoichiometric coefficients in mol/mol_pool (expl/ mol/molDNA, mol/molDAG etc.)
    Parameters:
        - dictdict_pool_data = Dictionary of Dictionary 
                                - key = metabolite ID in model, 
                                - value = dictionary with key_2 = ['Initial coefficient','Initial unit','cobra_model']
    Returns:
        - dictdict_pool_data updated:
            - new key_2 = 'MW','Coefficient in mol/smt','Normalized molar fraction in mol/molpool' 'Final stochiometric coefficient'
    """
     
    # 1) Convert data in mole per something
    for met in dictdict_pool_data:
        (met_MW,coeff_mol_per_smt) = calculate_coeff_in_mol(dictdict_pool_data[met],calculate_formula)
        dictdict_pool_data[met]["MW"] = met_MW
        dictdict_pool_data[met]["Coefficient in mol/smt"] = coeff_mol_per_smt
    
    # 2) Normalize data = molar fraction (mol/molmacromolecule or mol/molpool)
    pool_sum_mol = calculate_pool_mol(dictdict_pool_data)
    for met in dictdict_pool_data:
        final_coeff = normalize_data(dictdict_pool_data[met]["Coefficient in mol/smt"],pool_sum_mol)
        dictdict_pool_data[met]["Normalized molar fraction in mol/molpool"] = final_coeff
        dictdict_pool_data[met]["Final stochiometric coefficient"] = - np.round(final_coeff,5) # consumed metabolites
          
    return dictdict_pool_data


def generate_dictdict_pool_data_skip_calculation(dictdict_pool_data):
    """
    Parameter:
        - dictdict_pool_data = Dictionary of Dictionary 
                            - key = metabolite ID in model, 
                            - value = dictionary with key_2 = ['Initial coefficient','Initial unit','cobra_model']
    Returns:
        - dictdict_pool_data updated:
            -new key : 'Final stochiometric coefficient'
    """

    for met in dictdict_pool_data:
        coeff = dictdict_pool_data[met]
        dictdict_pool_data[met]["Final stochiometric coefficient"] = - coeff # consumed metabolites
          
    return dictdict_pool_data


def calculate_pool_weight(dictdict_pool_data):
    """
    Sum metabolites coefficients in g/smt.
    Parameter = dictionary (key = metabolite id, value = dictionary (with key = 'Coefficient in g/smt')).
    Returns the calculated sum.
    """
    pool_weight = 0
    
    for met_id in dictdict_pool_data:
        met_weight = dictdict_pool_data[met_id]["Coefficient in g/smt"]
        pool_weight = pool_weight + met_weight
        
    return pool_weight


def calculate_stochoimetric_coefficient_mmolpergDCW(dictdict_biomass_data,calculate_formula):
    """
    Calculates the stoichiometric coefficients in mmol/gDCW.
    Parameters:
        - dictdict_pool_data = Dictionary of Dictionary 
                                - key = metabolite ID in model, 
                                - value = dictionary with key_2 = ['Initial coefficient','Initial unit','cobra_model']
    Returns:
        - dictdict_pool_data updated:
            - new key_2 = 'MW','Coefficient in g/smt', 'Normalized massic fraction in g/gpool', 'Normalized massic fraction in g/gpool', 'Coeff in mmol/gpool','Final stochiometric coefficient' 
    """
    # 1) Convert data in g per g something if calculate formula
    if calculate_formula:
        for met in dictdict_biomass_data:

            (met_MW,coeff_g_per_smt) = calculate_coeff_in_gram(dictdict_biomass_data[met],calculate_formula)
            dictdict_biomass_data[met]["MW"] = met_MW
            dictdict_biomass_data[met]["Coefficient in g/smt"] = coeff_g_per_smt
        
        # 2) Normalize data (in g/gpool)
        pool_weight = calculate_pool_weight(dictdict_biomass_data)
        for met in dictdict_biomass_data:
            dictdict_biomass_data[met]["Normalized massic fraction in g/gpool"] = normalize_data(dictdict_biomass_data[met]["Coefficient in g/smt"],pool_weight)
            
         # 3) Convert in mmol/g
            coeff_mol_per_g = convert_gram_to_mol(dictdict_biomass_data[met]["Normalized massic fraction in g/gpool"],dictdict_biomass_data[met]["MW"],calculate_formula)
            final_coeff = coeff_mol_per_g * 1000 
            dictdict_biomass_data[met]["Coeff in mmol/gpool"] = final_coeff
            dictdict_biomass_data[met]["Final stochiometric coefficient"] = - np.round(final_coeff,5) # consumed metabolites
    else:
        for met in dictdict_biomass_data:
            original_coeff = dictdict_biomass_data[met]["Initial coefficient"]
            dictdict_biomass_data[met]["Final stochiometric coefficient"] = - np.round(original_coeff,5)

            
    return dictdict_biomass_data


def calculate_pool_charge(dict_rxn):
    """
	Calculates the charge of the new pseudo metabolite (pool).
        Warning: Requires that all metabolites consumed have a charge in model.
    Parameter: dict_rxn (dictionary : cobra_model.reactions.metabolites format, keys = cobra metabolite, value = its stoichiometric coefficient)
    Returns the calculated charge.
    """
    pool_charge = 0 

    for met in dict_rxn:
        
        met_coeff = dict_rxn[met]
        pool_charge = pool_charge + (- met_coeff * met.charge)

    pool_charge = int(np.round(pool_charge,0))
    if pool_charge == - 0:
        pool_charge = 0

    return(pool_charge)


def calculate_pool_formula(dict_rxn):
    """
    Calculates the forumula of the new pseudo metabolite (pool).
    Warning: Requires that all metabolites consumed have a formula in model.
    Parameter: dict_rxn (dictionary : cobra_model.reactions.metabolites format, keys = cobra metabolite, value = its stoichiometric coefficient)
    Returns the calculated formula.
    """
    dict_pool_formula = {'C': 0,
                         'H': 0,
                         'O': 0,
                         'N': 0,
                         'P': 0,
                         'S': 0,
                         'Mg':0}

    for met in dict_rxn:
        
        met_coeff = dict_rxn[met]
        dict_met_formula = met.elements

        for atom in dict_met_formula:
            dict_pool_formula[atom] = dict_pool_formula[atom] + (- met_coeff*dict_met_formula[atom])
    
    pool_formula = ""
    for atom in dict_pool_formula:
        if dict_pool_formula[atom] != 0 :
            atom_coeff = dict_pool_formula[atom]
            if atom_coeff < 0 :
                if atom_coeff > -1 :
                    print("original NEGATIVE pool formula :", atom_coeff, ". Coefficient rounded to 0.")

                else:
                    err = "Error in pool formula. Coefficient of atom = " + atom + atom_coeff
                    raise Exception(err)
            else :
                pool_formula = pool_formula + atom + str(atom_coeff)

    return(pool_formula)
