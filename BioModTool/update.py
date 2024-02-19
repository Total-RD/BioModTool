#!/usr/bin/env python
# coding: utf-8
# BioModTool

#------------------------------------------------ Imports --------------------------------------------------------------------
import numpy as np
import pandas as pd

import cobra
from BioModTool.calculate import  calculate_pool_formula, calculate_pool_charge

#------------------------------------------------ Functions ------------------------------------------------------------------

def create_reaction_dict(dictdict_pool_data, dictdict_constant_metabolites):
        """
        Create a dictionary with information required to create a new reaction (including constant metabolites)
        Parameters:
            - dictdict_pool_data = Dictionary of Dictionary 
                                    - key = metabolite ID in model, 
                                    - value = dictionary with keys2 = ['cobra_model','Final stochiometric coefficient'...]
            - dictdict_constant_metabolites = Dictionary of Dictionary 
                                    - key = metabolite ID in model, 
                                    - value = dictionary with keys2 = ['cobra_model','Initial coefficient'...]
                                    
        Returns the Dictionary (key = cobra.metabolite, value = stoichiometric coefficient)
        """
        dict_rxn = {}

        # Data from variables metabolites
        for met_id in dictdict_pool_data:
            met = dictdict_pool_data[met_id]["cobra_metabolite"]
            dict_rxn[met] = float(dictdict_pool_data[met_id]["Final stochiometric coefficient"]) 

        # Data for set/constant metabolites
        for met_id in dictdict_constant_metabolites:
            met = dictdict_constant_metabolites[met_id]["cobra_metabolite"]
            
            if met in dict_rxn: # if met already in dict_rxn, addition of constant coefficient to variable coefficient
                dict_rxn[met] = float(dict_rxn[met]) + float(dictdict_constant_metabolites[met_id]["Initial coefficient"])
            else:
                dict_rxn[met] = float(dictdict_constant_metabolites[met_id]["Initial coefficient"])
            
        return dict_rxn


def create_metabolite(cobra_model, met_id, name, compartment, formula=False, charge=False): 
    """
    Create a Cobra metabolite. /!\ Metabolite still needs to be added to the model.
    """
    met = cobra.Metabolite()
    if compartment not in cobra_model.compartments.keys():
        raise Exception('Compartment %s not found' % compartment)
    met.id = str(met_id)
    met.name = str(name)
    met.compartment = str(compartment)
    met.annotation["sbo"] ="SBO:0000247"
    if formula:
        met.formula = str(formula)
    if charge:
        met.charge = int(charge)
    return met

def create_pool_metabolite(cobra_model,dict_pool,pool_suffix_id,user_compartment,calculate_charge, calculate_formula):
    """
    Create pseudo metabolite macromolecule/pool (expl DNA, TAG etc).
    Calculate its charge and formula.
    
    Parameters:
        - dict_pool : dictionary (key = cobra.metabolite, value = stoichiometric coefficient)
        - pool_suffix_id: _suffix_ (str)
        - cobra model
    
    Return the created cobra.metabolite    
    """

    # 1) Create pool metabolite
    pool_met = create_metabolite(cobra_model=cobra_model,
                                 met_id = pool_suffix_id, 
                                 name = pool_suffix_id, 
                                 compartment = user_compartment)
    # 2) Calculate and add metabolite formula and charge
    if calculate_charge:
        pool_charge = calculate_pool_charge(dict_pool)
        pool_met.charge = pool_charge
    if calculate_formula:
        pool_formula = calculate_pool_formula(dict_pool)
        pool_met.formula = pool_formula

    return(pool_met)

def create_reaction(cobra_model, rxn_id, name, metabolites, subsystem='', lower_bound=-1000., upper_bound=1000.,gene_reaction_rule =""):
    """
    Create a Cobra reaction. /!\ Reaction still needs to be added to the model.
    """
    rxn = cobra.Reaction()
    rxn.id = str(rxn_id)  # create_reaction_id(cobra_model, name)
    rxn.name = str(name)
    rxn.subsystem = str(subsystem)
    rxn.annotation["sbo"] = "SBO:0000629"
    if np.isnan(lower_bound) or np.isnan(upper_bound):
        raise Exception('Undefined lower or upper bound when creating '
                        'reaction %s' % rxn_id)
    rxn.lower_bound = lower_bound
    rxn.upper_bound = upper_bound
    rxn.add_metabolites(metabolites)
    rxn.gene_reaction_rule = gene_reaction_rule
    return rxn

def create_pool_synthesis_reaction(cobra_model,pool_met,dict_pool):
    """
    Create the cobra.reaction: synthesis of the pool/macromolecule.
    
    Parameters:
        - cobra.model
        - pool_met: cobra.metabolite, the synthetized metabolite (expl DNA or TAG)
        - dict_pool: dictionary (key = cobra.metabolite, value = stoichiometric coefficient), missing synthetized pool
    Returns: the created cobra.reaction
    """

    # 1) Add pool_met to dict_pool
    dict_pool[pool_met] = 1

    # 2) Create pool reaction
    rxn_id_in_model = pool_met.id
    reaction = create_reaction(cobra_model=cobra_model, 
                               rxn_id = rxn_id_in_model , 
                               name = "Biomass objective function: synthesis of " + rxn_id_in_model, 
                               metabolites = dict_pool, 
                               subsystem='Biomass reaction', 
                               lower_bound=0, 
                               upper_bound=100000.,
                               gene_reaction_rule ="NO_GENE")
    return reaction


def generate_pool_metabolite_and_reaction(cobra_model,dict_pool,pool_suffix_id,user_compartment,calculate_charge, calculate_formula):
    """
    Create cobra.metabolite and cobra.reaction: synthesis of the pool/macromolecule.
    
    Parameters:
        - cobra.model
        - dict_pool: dictionary (key = cobra.metabolite, value = stoichiometric coefficient)
        - suffix id
    Returns: cobra.metabolite, cobra.reaction
    """
    
    # 1) Create Metabolite
    pool_met = create_pool_metabolite(cobra_model,dict_pool,pool_suffix_id,user_compartment,calculate_charge, calculate_formula)
    # 2) Create Reaction
    pool_rxn = create_pool_synthesis_reaction(cobra_model,pool_met,dict_pool)
    
    return(pool_met, pool_rxn)


def add_EX_biomass_reaction(cobra_model,biomass_met,biomass_suffix_id):
    # Create Biomass pseudo metabolite export reaction

    # Reaction informations
    ex_id = "EX_" + biomass_suffix_id[:-1]
    ex_name = "Drain to BIOMASS (" + biomass_suffix_id + ")"
    ex_rxn_dict = {biomass_met: -1}
    
    # Create reaction
    ex_rxn = create_reaction(cobra_model, ex_id, ex_name, ex_rxn_dict, subsystem='Exchange reaction', lower_bound=0, upper_bound=1000.,gene_reaction_rule ="NO_GENE")
    
    # Add reaction
    cobra_model.add_reactions([ex_rxn])
    
    return(cobra_model)


def update_model(cobra_model,met,rxn):
    """
    Add cobra.metabolite and cobra.reaction to cobra.model.
    
    Parameters:
        - cobra.model
        - cobra.metabolites
        - cobra.reactions
        
    Returns: updated cobra.model
    """
    
    # Add metabolite
    cobra_model.add_metabolites(met)
    cobra_model.metabolites.get_by_id(met.id).formula = met.formula
    cobra_model.metabolites.get_by_id(met.id).charge = met.charge
    
    # Add reaction
    cobra_model.add_reactions([rxn])
    
    print("Metabolite %s (formula: %s, charge: %s) added to model. \nReaction %s added to model. \n%s" %(cobra_model.metabolites.get_by_id(met.id).id,cobra_model.metabolites.get_by_id(met.id).formula,cobra_model.metabolites.get_by_id(met.id).charge,cobra_model.reactions.get_by_id(rxn.id).id,cobra_model.reactions.get_by_id(rxn.id).reaction))
    print("---------------------------------------------------------------------------------------------------------------------")
    return(cobra_model)