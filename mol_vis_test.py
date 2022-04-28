import streamlit as st
import pandas as pd
import math 
#_________________________
import py3Dmol
from stmol import showmol
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
#_________________________
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors

from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D

from rdkit.Chem.Draw import SimilarityMaps
from rdkit.Chem import rdMolDescriptors
#_________________________
import requests
import json
from io import StringIO
#_________________________
#Inicio#

st.title('VISUALIZACIÓN MOLECUALR')

seleccion = st.selectbox("Seleccione una opción: ", ["Visualizacion molecular", "Reactividad"])

if seleccion == "Visualizacion molecular":
    st.title('VISUALIZACIÓN MOLECUALR')
    st.write("Bienvenido. Aquí podrás ver la molecula en su forma tridimensional")
    
    seleccion_molecula = st.selectbox("Seleccione una opción: ", ["SMILES", "Subir un archivo"])
    def otros_parametros(mol):
      
      m = Chem.MolFromSmiles(compound_smiles)
      tpsa = Descriptors.TPSA(m)
      logP = Descriptors.MolLogP(m)
      st.write("TPSA: " + str(tpsa))
      st.write("Log P: "+ str(logP))
      mol = Chem.MolFromSmiles(compound_smiles)
      
      # Gasteiger partial charges
      st.subheader("Gesteiger partial charges")
      AllChem.ComputeGasteigerCharges(mol)
      contribs = [mol.GetAtomWithIdx(i).GetDoubleProp('_GasteigerCharge') for i in range(mol.GetNumAtoms())]
      fig = SimilarityMaps.GetSimilarityMapFromWeights(mol, contribs, colorMap='jet', contourLines=10)
      st.pyplot(fig)

      # Crippen contributions to logP
      st.subheader("Crippen contributions to log P")
      contribs = rdMolDescriptors._CalcCrippenContribs(mol)
      fig2 = SimilarityMaps.GetSimilarityMapFromWeights(mol,[x for x,y in contribs], colorMap='jet', contourLines=10)
      st.pyplot(fig2)

    if seleccion_molecula == "Subir un archivo":
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
      file_type = st.sidebar.radio("Tipo de archivo", ("xyz","mol","sdf"))

      for uploaded_file in uploaded_files:
        xyz = uploaded_file.getvalue().decode("utf-8")
        render_mol(xyz)
        #st.write(xyz)
      
      #xyz to SMILES
      if file_type == "xyz":
        def xyz_to_smi(str_input):
          webserver_url = "https://www.cheminfo.org/webservices/babel"
          options='{"ph":"","hydrogens":"No change","coordinates":"None","inputFormat":"xyz -- XYZ cartesian coordinates format","outputFormat":"smi -- SMILES format"}'
          req = requests.post(webserver_url,data = {'options':options,'input':str_input})
      
          if req.status_code == 200:
            return json.loads(req.text)['result'].split()[0]
          else:
            return None
    
      if file_type == "mol":
        def xyz_to_smi(str_input):
          webserver_url = "https://www.cheminfo.org/webservices/babel"
          options='{"ph":"","hydrogens":"No change","coordinates":"None","inputFormat":"mol -- MDL MOL format","outputFormat":"smi -- SMILES format"}'
          req = requests.post(webserver_url,data = {'options':options,'input':str_input})
      
          if req.status_code == 200:
            return json.loads(req.text)['result'].split()[0]
          else:
            return None

      if file_type== "sdf":
        def xyz_to_smi(str_input):
          webserver_url = "https://www.cheminfo.org/webservices/babel"
          options='{"ph":"","hydrogens":"No change","coordinates":"None","inputFormat":"sdf -- MDL MOL format","outputFormat":"smi -- SMILES format"}'
          req = requests.post(webserver_url,data = {'options':options,'input':str_input})
      
          if req.status_code == 200:
            return json.loads(req.text)['result'].split()[0]
          else:
            return None

      for uploaded_file in uploaded_files:
        stringio = StringIO(uploaded_file.getvalue().decode("utf-8"))
        string_data = stringio.read()

        compound_smiles = xyz_to_smi(string_data)
        st.subheader("SMILES: " + xyz_to_smi(string_data))

        otros_parametros(compound_smiles)


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







