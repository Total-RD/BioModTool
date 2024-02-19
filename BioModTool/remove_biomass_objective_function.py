#!/usr/bin/env python
# coding: utf-8
# BioModTool

#------------------------------------------------ Imports --------------------------------------------------------------------
import pandas as pd

import BioModTool.main_add_biomass_objective_function

import BioModTool.test

#------------------------------------------------ Functions ------------------------------------------------------------------


def remove_biomass_objective_function(cobra_model,suffix):  
    """
    Remove all reactions/metabolites ID of Biomass Obective Function associated with a given suffix (expl. BIOMASS_suffix_c, DNA_suffix_c, DAG_suffix_c etc.)
    Parameters: 
        - cobra.model
        - suffix (str)
    Returns: updated model without targeted reactions and pseudo metabolites.
    """
    cobra_model_copy = cobra_model.copy()
    list_BOF_ids = BioModTool.test.identify_BOF_by_suffix(cobra_model,suffix)
    met_to_rm = []
    for my_id in list_BOF_ids:
        if my_id[:3] != "EX_":
            met_to_rm.append(cobra_model_copy.metabolites.get_by_id(my_id))
    cobra_model_copy.remove_reactions(list_BOF_ids)
    cobra_model_copy.remove_metabolites(met_to_rm)

    return(cobra_model_copy)