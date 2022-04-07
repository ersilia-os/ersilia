import streamlit as st
from ersilia import ErsiliaModel
import json

from rdkit import Chem
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import Draw

import streamlit as st

def show_output(obj):
    obj = list(obj) 
    d = dict(obj[0])
    d = d['output']
    pair = next(iter((d.items())) )
    st.write(pair[0], " : ", pair[1])

def app():
    st.header("Ersilia Model Hub")
    st.subheader("Calculate and Predict properties for your molecules of interest")

    option = st.selectbox('Which model would you like to run?',('chemprop-antibiotic','molecular-weight'))
    model = ErsiliaModel(option)
    model.serve()

    compound_smiles = st.text_input("Enter the chemical structure", "")
    m = Chem.MolFromSmiles(compound_smiles)
    im=Draw.MolToImage(m,size=(200,200))
    st.image(im)

    apis = tuple(model.get_apis())

    #if only one api available, no need to create radio buttons
    if len(apis) == 1:
        opt = apis[0]
        if opt=='predict':
            if(st.button("Predict")):
                output_obj = model.predict(compound_smiles)
                show_output(output_obj)
        elif opt=='calculate':
            if(st.button("Calculate")):
                output_obj = model.calculate(compound_smiles)
                show_output(output_obj)
    # if model supports more than one api, create 
    else:
        opt = st.radio("Select API", apis)
        if opt=='predict':
            output_obj = model.predict(compound_smiles)
            show_output(output_obj)
        elif opt=='calculate':
            output_obj = model.calculate(compound_smiles)
            show_output(output_obj)
            
    



