import streamlit as st
from ersilia import ErsiliaModel
import json

from rdkit import Chem
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import Draw

import streamlit as st


def app():
    st.header("Ersilia Model Hub")
    st.subheader("Calculate and Predict properties for your molecules of interest")

    option = st.selectbox('Which model would you like to run?',('chemprop-antibiotic','molecular-weight'))
    model = ErsiliaModel(option)
    model.serve()

    compound_smiles = st.text_input("Enter the chemical structure", "")
    m = Chem.MolFromSmiles(compound_smiles)
    im=Draw.MolToImage(m)

    st.image(im)

    opt = st.radio("Select API", ('Predict','Calculate'))
    if opt=='Predict':
        obj = model.predict(compound_smiles)
        obj = list(obj) 
        json_response = obj[0]
        st.write(json_response)
    else:
        obj = model.calculate(compound_smiles)
        obj = list(obj) 
        json_response = obj[0]
        st.write(json_response)



