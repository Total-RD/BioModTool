#!/usr/bin/env python
# coding: utf-8
# BioModTool

#------------------------------------------------ Imports --------------------------------------------------------------------
import pandas as pd
from pathlib import Path

import cobra
#------------------------------------------------ Functions ------------------------------------------------------------------

def save_model(cobra_model,path_to_updated_model):
    # save a given cobra model in both json and sbml format.
    file_path = path_to_updated_model.split(".")[0]
    # save as json
    cobra.io.save_json_model(cobra_model,file_path+".json")
    # save as sbml
    cobra.io.write_sbml_model(cobra_model,file_path+".xml")


def save_results(file_path,dict_df_calculs,dict_df_fixed_metabolites,dict_df_rxn,dict_metabolites,df_res_balance):
    """    
    Save all calculations performed by BioModTool and final stoichiometric coefficients in an Excel file.
    The Excel is composed of several sheets:
        - one sheet for each added pseudo-reaction, recapitulated final stoichiometric coefficients and their calculations,
        - one sheet recapitulating all pseudo-metabolites added to model (with their charge and formula),
        - one sheet with results of mass and charge balance test for all created pseudo-reactions.
    """
    with pd.ExcelWriter(file_path) as writer:

        for pool_id in dict_df_calculs:
            df_rxn = dict_df_rxn[pool_id]
            for idx in df_rxn.index:
                df_rxn.loc[idx,"Metabolite_id"] = idx.id
            df_rxn = df_rxn.set_index("Metabolite_id")
            
            # Concatenate df
            df_all_metabolites = pd.concat([dict_df_calculs[pool_id],dict_df_fixed_metabolites[pool_id]])

            df_all_metabolites = df_all_metabolites.set_index("cobra_metabolite")
            for df_all_met_idx in df_all_metabolites.index:
                df_all_metabolites.loc[df_all_met_idx,"Metabolite_id"] = df_all_met_idx.id
            df_all_metabolites = df_all_metabolites.set_index("Metabolite_id")
            df_all_metabolites.index = df_all_metabolites.index.where(~df_all_metabolites.index.duplicated(keep="first"), df_all_metabolites.index + ' (from table 2)')
            

            df_bilan = pd.concat([df_rxn,df_all_metabolites],axis=1)
            df_bilan.index.name="Metabolite_id"
            df_bilan.to_excel(writer, sheet_name=pool_id)
        
        df_metabo = pd.DataFrame.from_dict(dict_metabolites,orient="index")
        df_metabo.to_excel(writer, sheet_name="added_metabolites")
        df_res_balance.to_excel(writer, sheet_name="mass_charge_balance")