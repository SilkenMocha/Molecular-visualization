import streamlit as st
import py3Dmol
from stmol import showmol
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

import requests
import json
from io import StringIO

st.title('VISUALIZACIÓN MOLECUALR')

seleccion_molecula = st.selectbox("Seleccione una opción: ", ["SMILES", "Subir un archivo"])

if seleccion_molecula == "Subir un archivo":
    #render xyz
    def render_mol(xyz):
       xyzview = py3Dmol.view(width=400,height=400)
       xyzview.addModel(xyz,'xyz')
       xyzview.addModel(xyz,'sdf')
       xyzview.addModel(xyz,'mol')
       xyzview.setStyle({'stick':{}})
       xyzview.setBackgroundColor('white')#('0xeeeeee')
       xyzview.zoomTo()
       showmol(xyzview, height = 500,width=800)      
      
    uploaded_files = st.sidebar.file_uploader("Choose xyz files", accept_multiple_files=True)
    #file_type = st.sidebar.radio("Tipo de archivo", ("xyz","mol","sdf"))
        
    for uploaded_file in uploaded_files:
        xyz = uploaded_file.getvalue().decode("utf-8")
        render_mol(xyz)

    #if file_type == "xyz": 
        #for uploaded_file in uploaded_files:
            #xyz = uploaded_file.getvalue().decode("utf-8")
            #render_mol(xyz)
            #st.write(xyz)
      
    #render sdf
    #def render_mol(sdf):
        #sdfview = py3Dmol.view(width=400,height=400)
        #sdfview.addModel(sdf,'sdf')
        #sdfview.setStyle({'stick':{}})
        #xyzview.setBackgroundColor('white')#('0xeeeeee')
        #sdfview.zoomTo()
        #sdfview.show()
    
    for uploaded_file in uploaded_files:
        stringio = StringIO(uploaded_file.getvalue().decode("utf-8"))
        string_data = stringio.read()
        


if seleccion_molecula == "SMILES":
    compound_smiles=st.text_input('SMILES please','COc1cccc2cc(C(=O)NCCCCN3CCN(c4cccc5nccnc54)CC3)oc21')
    def makeblock(smi):
        mol = Chem.MolFromSmiles(smi)
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol)
        mblock = Chem.MolToMolBlock(mol)
        return mblock

    def render_mol(xyz):
        xyzview = py3Dmol.view()#(width=1000,height=1000)
        xyzview.addModel(xyz,'mol')
        xyzview.setStyle({'stick':{}})
        xyzview.setBackgroundColor('white')
        xyzview.zoomTo()
        showmol(xyzview,height=500,width=500)

    blk=makeblock(compound_smiles)
    render_mol(blk)
    otros_parametros(compound_smiles)




