#!/usr/bin/env python
# coding: utf-8
# BioModTool : interface

#------------------------------------------------ Imports --------------------------------------------------------------------
from tkinter import *
from tkinter import ttk
from tkinter import filedialog , messagebox
from functools import partial

import cobra
import os
import pandas as pd

from BioModTool.main_add_biomass_objective_function import get_ids_by_level,define_reaction_id,add_biomass_objective_function
from BioModTool.load import load_cobra_model
from BioModTool.save import save_model
from BioModTool.test import test_suffix_conformity, test_suffix_in_model

#------------------------------------------------ Functions ------------------------------------------------------------------
def open_data_file(entry_path_to_data):
    path_to_data = filedialog.askopenfilename(title="Open data (xlsx)",
                                               filetypes=[("Excel file",".xlsx"),("All files",".*")])
    entry_path_to_data.insert(0, path_to_data)

def open_cobra_model_file(entry_path_to_model):
    path_to_model = filedialog.askopenfilename(title="Open Model (sbml or json)",
                                               filetypes=[("SBML model",".xml"),("JSON model",".json"),("All files",".*")])
    entry_path_to_model.insert(0, path_to_model)


def select_formula(variable_formula,label_formula2):
    if variable_formula.get():
        label_formula2.config(text = "Formula will be calculated.")
    else:
        label_formula2.config(text="Formula will not be calculated. Level 1 coefficients must be in mmol/gDW.")


def select_charge(variable_charge,label_charge2):
    if variable_charge.get():
        label_charge2.config(text = "Charge will be calculated.")
    else:
        label_charge2.config(text="Charge will not be calculated.")



def test_load_Model(model_file,variable_formula,variable_charge,window):
    bool_calculate_formula = variable_formula.get()
    bool_calculate_charge = variable_charge.get()
    
    # test to load model
    path_to_model = model_file.get()
    try:
        model = load_cobra_model(path_to_model)
        
    except:
        messagebox.showerror("Error", "Problem to load model. Please check that you selected model in JSON or SBML format.")
        window.destroy()
        interface_step1()
        
     
    window.destroy()
    interface_step2(model,bool_calculate_formula,bool_calculate_charge)


def get_user_compartment(cobra_model,bool_calculate_formula,bool_calculate_charge,user_compartment,window_2):
    
    # Get user compartment
    key_list = [k  for (k, val) in cobra_model.compartments.items() if val == user_compartment.get()]
    my_comp = key_list[0]
    window_2.destroy()
    interface_step3(cobra_model,bool_calculate_formula,bool_calculate_charge,my_comp)


def add_user_BOF(cobra_model,bool_calculate_formula,bool_calculate_charge,my_comp,data_file,dict_structure,suffix_entry,window_5):
    try:
        # Test biomass structure
        # List of reaction IDs by level : to lipid, to biomass
        biomass_id = get_ids_by_level(dict_structure,"level_1") 
        list_all_pool_id = [biomass_id]

        if "level_2" in dict_structure.values():
            list_biomass = get_ids_by_level(dict_structure, "level_2") 
            list_all_pool_id = list_all_pool_id + list_biomass 
            list_all_pool_id_wo_biomass = list_biomass 

            if ("level_3" in dict_structure.values()) and ("level_2_lipid" in dict_structure.values()):
                list_lipids = get_ids_by_level(dict_structure, "level_3")
                level_2_lipid_id = get_ids_by_level(dict_structure, "level_2_lipid")
                list_all_pool_id = list_all_pool_id + list_lipids + [level_2_lipid_id]
                list_all_pool_id_wo_biomass = list_all_pool_id_wo_biomass + list_lipids + [level_2_lipid_id]

            elif ("level_3" not in dict_structure.values()) and ("level_2_lipid" not in dict_structure.values()):
                structure = "OK"
            else:
                messagebox.showerror("Error","Error in choosen, if you want to add level 3, you must define level_3 and level_2_lipid pseudo-metabolites.Please restart BioModTool.")
                window_5.destroy()
                interface_step3(cobra_model,bool_calculate_formula,bool_calculate_charge,my_comp)

        # test suffix
        suffix = suffix_entry.get()
        try:
            test_suffix_conformity(suffix)
        except:
            messagebox.showerror("Error", "Problem with suffix, not conform (must be combinaison of alpha numeric characters).")
            window_5.destroy()
            interface_step3(cobra_model,bool_calculate_formula,bool_calculate_charge,my_comp)
        try:
            test_suffix_in_model(cobra_model,suffix,list_all_pool_id)

        except:
            messagebox.showerror("Error", "Problem with suffix, already in model.")
            window_5.destroy()
            interface_step3(cobra_model,bool_calculate_formula,bool_calculate_charge,my_comp)
        try:
            updated_model = add_biomass_objective_function(cobra_model=cobra_model,
                                                           path_to_data=data_file,
                                                           suffix=suffix,
                                                           dict_structure=dict_structure,
                                                           user_compartment=my_comp,
                                                           calculate_charge=bool_calculate_charge,
                                                           calculate_formula = bool_calculate_formula,
                                                           saving_final_data = True)

            button_save = Button(window_5, text="Save and Quit", command=partial(user_save_model,updated_model,window_5,suffix),bg="darkgrey")
            button_save.grid(row=7,column=0,pady=25) 

        except Exception as e:
            messagebox.showerror("Error",e)
            window_5.destroy()
            
    except Exception as e:
        messagebox.showerror("Error",e)
        window_5.destroy()
        

def user_save_model(cobra_model,window_5,suffix):
    file_name = filedialog.asksaveasfilename(initialfile = "Updated_model_"+suffix+".json",defaultextension=".txt",filetypes=[("All Files","*.*")],title="Save updated Model (.json and .xml)")
    save_model(cobra_model,file_name)
    messagebox.showinfo("Information", "Model was successfully updated and saved.")

    
    window_5.destroy()


def interface_step1():
    i_row = 0

    # Create new window
    window = Tk() 
    window.title("BioModTool ")
    window.geometry("960x450")

    # Label / subtitle (Step 1/5)
    subtitle_label = Label(window, text="STEP 1/5 : Load a Genome Scale Metabolic Model",anchor=CENTER)
    subtitle_label.grid(row=i_row,column=0)
    i_row +=1

    lf_path_to_mod = LabelFrame(window, text="Choose model file")
    lf_path_to_mod.grid(row=i_row,column=0,padx=10,pady=20)
    i_row +=1


    my_label = Label(lf_path_to_mod, text="Path to model (.json or .sbml/xml) :")
    my_label.grid(row=i_row,column=0,padx=10,pady=5)

    model_file = StringVar()

    entry_path_to_model = Entry(lf_path_to_mod, textvariable=model_file, width=45)
    entry_path_to_model.grid(row=i_row,column=1)

    button_browthe_user_path_to_model =Button(lf_path_to_mod, text="Browse...", command=partial(open_cobra_model_file,entry_path_to_model))
    button_browthe_user_path_to_model.grid(row=i_row,column=2,padx=3)
    i_row +=1


    lf_Model_properties = LabelFrame(window, text="Model properties")
    lf_Model_properties.grid(row=i_row,column=0,padx=10,pady=5)
    i_row +=1

    lf_warning_formula = Label(lf_Model_properties, text= "\n /!\ Note that if you select formula = False, the molecular weight cannot be calculated from the metabolite formula preventing unit conversion.")
    lf_warning_formula.config(fg="red")
    lf_warning_formula.grid(row=i_row,column=0,padx=3,stick="w",columnspan=4)
    i_row +=1
    
    lf_warning_formula2 = Label(lf_Model_properties, text= "   - Level 1 data must be given in mmol.gDW-1.")
    lf_warning_formula2.config(fg="red")
    lf_warning_formula2.grid(row=i_row,column=0,padx=3,stick="w",columnspan=4)
    i_row +=1
    
    lf_warning_formula3 = Label(lf_Model_properties, text= "   - Levels 2 and 3 data must be given in mol per ...\n\n ")
    lf_warning_formula3.config(fg="red")
    lf_warning_formula3.grid(row=i_row,column=0,padx=3,stick="w",columnspan=4)
    i_row +=1
    
     
    # Create the new variable
    variable_formula = BooleanVar()
    variable_charge = BooleanVar()

    # Create the radio button
    # Formula
    label_formula1 = Label(lf_Model_properties,text="Are metabolic formula available in model ?")
    label_formula1.grid(row = i_row, column = 0)
    i_row +=1
    label_formula2 = Label(lf_Model_properties,text="                                             ")
    label_formula2.grid(row = i_row, column = 3)

    button_formula_yes = Radiobutton(lf_Model_properties, variable=variable_formula, text="Yes", command=partial(select_formula,variable_formula,label_formula2),value=True)
    button_formula_no = Radiobutton(lf_Model_properties, variable=variable_formula, text="No", command=partial(select_formula,variable_formula,label_formula2), value=False)
    button_formula_yes.grid(row = i_row, column = 1,stick="w",pady=5)
    button_formula_no.grid(row = i_row, column = 2,stick="w")
    i_row +=3

    # Charge
    label_charge1 = Label(lf_Model_properties,text="Are metabolic charge available in model ?")
    label_charge1.grid(row = i_row, column = 0)
    i_row +=1
    label_charge2 = Label(lf_Model_properties)
    label_charge2.grid(row = i_row, column = 3)

    button_charge_yes = Radiobutton(lf_Model_properties, variable=variable_charge, text="Yes", command=partial(select_charge,variable_charge,label_charge2),value=True)
    button_charge_no = Radiobutton(lf_Model_properties, variable=variable_charge, text="No", command=partial(select_charge,variable_charge,label_charge2), value=False)
    button_charge_yes.grid(row = i_row, column = 1,stick="w")
    button_charge_no.grid(row = i_row, column = 2,stick="w")
    i_row +=1

    button_next_to_step2 = Button(window, text="Next Step", command=partial(test_load_Model,model_file,variable_formula,variable_charge,window),bg="darkgrey")
    button_next_to_step2.grid(row=i_row,column=0,pady=25)

    window.mainloop()



def interface_step2(cobra_model,bool_calculate_formula,bool_calculate_charge):
    i_row = 0

    # Create new window
    window_2 = Tk() 
    window_2.title("BioModTool ")
    window_2.geometry("500x200")

    # window_2 Label / subtitle (Step 1/5)
    subtitle_label = Label(window_2, text="STEP 2/5: Define the compartment in which the biomass reaction will be added",anchor=CENTER)
    subtitle_label.grid(row=i_row,column=0,padx=15)
    i_row +=1


    user_compartment = StringVar()
    listeCombo = ttk.Combobox(window_2, textvariable=user_compartment, values=list(cobra_model.compartments.values()))

    # Choose the default item to be displayed
    listeCombo.current(0)

    listeCombo.grid(row=i_row,column=0,pady=5)
    i_row+=1

    button_next_to_step2 = Button(window_2, text="Next Step",bg="darkgrey",command=partial(get_user_compartment,cobra_model,bool_calculate_formula,bool_calculate_charge,user_compartment,window_2))
    button_next_to_step2.grid(row=i_row,column=0,pady=15)

    window_2.mainloop()

def interface_step3(cobra_model,bool_calculate_formula,bool_calculate_charge,my_comp):
    i_row = 0

    # Create new window
    window_3 = Tk()  # nouvelle instance de Tk
    window_3.title("BioModTool ")
    window_3.geometry("800x150")

    # window_3 Label / subtitle (Step 1/5)
    subtitle_label = Label(window_3, text="STEP 3/5: Load Biomass composition data (Excel file)",anchor=CENTER)
    subtitle_label.grid(row=i_row,column=0,padx=10)
    i_row +=1

    my_label = Label(window_3, text="Path to data (.xlsx) :")
    my_label.grid(row=i_row,column=0,padx=10,pady=5)

    data_file = StringVar()

    entry_path_to_data = Entry(window_3, textvariable=data_file, width=60)
    entry_path_to_data.grid(row=i_row,column=1)

    button_browthe_user_path_to_data =Button(window_3, text="Browse...", command=partial(open_data_file,entry_path_to_data))
    button_browthe_user_path_to_data.grid(row=i_row,column=2,padx=3)
    i_row +=1

    button_next_to_step4 = Button(window_3, text="Next Step", command=partial(interface_step4,cobra_model,bool_calculate_formula,bool_calculate_charge,my_comp,data_file,window_3),bg="darkgrey")
    button_next_to_step4.grid(row=i_row,column=1,pady=25)

    window_3.mainloop()


def interface_step4(cobra_model,bool_calculate_formula,bool_calculate_charge,my_comp,var_data_file,window_3):
    data_file = var_data_file.get()

    try:
        xls = pd.ExcelFile(data_file)
    except:
        messagebox.showerror("Error", "Error while loading data. Data must be provided as an Excel file.")
        window_3.destroy()
        interface_step3(cobra_model,bool_calculate_formula,bool_calculate_charge,my_comp)
    
    sheets = xls.sheet_names

    window_3.destroy()
    
    def function_choose_level(user_choice,my_sheet,valeur):
        user_choice[my_sheet] = valeur
        
    i_row = 0

    # Create new window
    window_4 = Tk() # nouvelle instance de Tk
    window_4.title("BioModTool ")
    window_4.geometry("1000x1600")

    # window_4 Label / subtitle (Step 1/5)
    subtitle_label = Label(window_4, text="STEP 4/5: Define structure of the biomass reaction.",anchor=CENTER)
    subtitle_label.grid(row=i_row,column=3)
    i_row +=1

    my_label = Label(window_4, text="Choose a level for each pseudo-metabolite (proposed metabolites correspond to sheet IDs of excel file):")
    my_label.grid(row=i_row,column=3,pady=5)
    i_row +=1
    
    user_choice = {}
    if 'Read_Me' in sheets:
        sheets.remove('Read_Me')
    
    if 'new_pool' in sheets:
        sheets.remove('new_pool')

    for i in range(0,len(sheets)):
        my_sheet = sheets[i]
        # Create the new variable
        variable = IntVar()

        # Create the radio button
        met_label = Label(window_4,text=my_sheet)
        met_label.grid(row = i_row + i, column = 1)
        button_level1 = Radiobutton(window_4, variable=variable, text="level_1", command=partial(function_choose_level,user_choice, my_sheet,"level_1"),value=1)
        button_level2 = Radiobutton(window_4, variable=variable, text="level_2", command=partial(function_choose_level,user_choice, my_sheet,"level_2"),value=2)
        button_level2_lip = Radiobutton(window_4, variable=variable, text="level_2_lipid", command=partial(function_choose_level,user_choice, my_sheet,"level_2_lipid"),value=3)
        button_level3 = Radiobutton(window_4, variable=variable, text="level_3", command=partial(function_choose_level,user_choice, my_sheet,"level_3"),value=4)
        button_level1.grid(row = i_row + i, column = 2, pady=5)
        button_level2.grid(row = i_row + i, column = 3, pady=5)
        button_level2_lip.grid(row = i_row + i, column = 4, pady=5)
        button_level3.grid(row = i_row + i, column = 5, pady=5)

    i_row = i_row + i + 1
    
    button_next_to_step5 = Button(window_4, text="Next Step", command=partial(interface_step5,cobra_model,bool_calculate_formula,bool_calculate_charge,my_comp,data_file,user_choice,window_4),bg="darkgrey")
    button_next_to_step5.grid(row=i_row,column=3,pady=25)
    

    window_4.mainloop() 


def interface_step5(cobra_model,bool_calculate_formula,bool_calculate_charge,my_comp,data_file,user_choice,window_4):
    
    window_4.destroy()
    i_row = 0

    # Create new window
    window_5 = Tk()
    window_5.title("BioModTool ")
    window_5.geometry("500x250")

    # window_5 Label / subtitle (Step 1/5)
    subtitle_label = Label(window_5, text="STEP 5/5: Define the suffix of the biomass reaction",anchor=CENTER)
    subtitle_label.grid(row=i_row,column=0)
    i_row +=1
    
    user_suffix = StringVar()
    label_suffix = Label(window_5, text="Choose a suffix (allowed characters: alphanumeric + _) :")
    label_suffix.grid(row=i_row,column = 0)
    suffix_entry = Entry(window_5,textvariable=user_suffix, width=30)
    suffix_entry.grid(row=i_row,column=1)
    
    i_row +=1 
    
    Add_BOF = Button(window_5, text="Add new BOF", command=partial(add_user_BOF,cobra_model,bool_calculate_formula,bool_calculate_charge,my_comp,data_file,user_choice,suffix_entry,window_5),bg="darkgrey")
    Add_BOF.grid(row=i_row,column=1,pady=25)  

    window_5.mainloop() 


def BioModTool_interface():
    interface_step1()


BioModTool_interface()

