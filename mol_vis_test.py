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
import pubchempy

#Inicio#
st.title ("FENÓMENOS CUÁNTICOS")
st.subheader("Erick López Saldviar 348916")





seleccion = st.selectbox("Seleccione una opción: ", ["Visualizacion molecular", "Reactividad"])

#__________________________________________________________________________________________
#Calculo reactividad#
if seleccion == "Reactividad":
  with st.form(key='calc_react'):
    st.write("Bienvenido. Este programa te calculara parámetros de reactividad")
    st.subheader('Hartress')
    ht = st.number_input("Hartress: ", format="%.4f", step=1e-4)
    ht0 = st.number_input("Hartress 0: ", format="%.4f", step=1e-4)
    ht1p = st.number_input("Hartress +1: ", format="%.4f", step=1e-4)
    ht1m = st.number_input("Hartress -1: ", format="%.4f", step=1e-4)
    st.subheader('Nucleofilicidad')
    homo = st.number_input("HOMO:", format="%.4f", step=1e-4)
    lumo = st.number_input("LUMO: ", format="%.4f", step=1e-4)
    homo_o = st.number_input("HOMO (O):", format="%.4f", step=1e-4)
    lumo_o = st.number_input("LUMO (V): ", format="%.4f", step=1e-4)

    calcular = st.form_submit_button('Calcular')
 
 
  if calcular:
    eV0 = ht0 * 27.2116
    eV1p = ht1p * 27.2116
    eV1m = ht1m * 27.2116
  
    kcal = ht*627.5

    A = eV0-eV1m
    I = eV1p - eV0

    n = (I-A)/2
    u = (I+A)/2
    w = (pow(u,2))/(2*n)


    col1, col2 = st.columns(2)
    col1.metric(label="H 0", value=round(ht0,5))
    col2.metric(label="kCal/mol", value=str(round(kcal,5)))

    col1, col2 = st.columns(2)
    col1.metric(label="H -1", value=round(ht1m,5))
    col2.metric(label="H +1", value=round(ht1p,5))


    col1, col2, col3 = st.columns(3)
    col1.metric(label="eV 0", value=str(round(eV0,5)))
    col2.metric(label="eV -1", value=str(round(eV1p,5)))
    col3.metric(label="eV +1", value=str(round(eV1m,5)))
  
    st.subheader("Aproximación de energias")

    col1, col2 = st.columns(2)
    col1.metric(label="Afinidad electrónica", value=str(round(A,5)))
    col2.metric(label="Potencial de ionización", value=str(round(I,5)))

    col1, col2, col3 = st.columns(3)
    col1.metric(label="Dureza", value=str(round(n,5)))
    col2.metric(label="Electronegatividad", value=str(round(u,5)))
    col3.metric(label="Electrofilicidad", value=str(round(w,5)))

    st.subheader("Aproximación Orbital")

    A_orb = -27*lumo
    I_orb = -27*homo

    n_orb = (I_orb - A_orb) / 2
    u_orb = (I_orb + A_orb) / 2
    w_orb = pow(u_orb, 2) / (2 * n_orb)
  
    col1, col2 = st.columns(2)
    col1.metric(label="Afinidad electrónica", value=str(round(A_orb,5)))
    col2.metric(label="Potencial de ionización", value=str(round(I_orb,5)))

    col1, col2, col3 = st.columns(3)
    col1.metric(label="Dureza", value=str(round(n_orb,5)))
    col2.metric(label="Electronegatividad", value=str(round(u_orb,5)))
    col3.metric(label="Electrofilicidad", value=str(round(w_orb,5)))

    gap = lumo_o - homo_o
    gap_eV = gap*27.2116
  
    col1, col2 = st.columns(2)
    col1.metric(label="GAP", value=str(round(gap,5)))
    col2.metric(label="GAP (eV", value=str(round(gap_eV,5)))
#__________________________________________________________________________________________
#Visualización molecular#
if seleccion == "Visualizacion molecular":
    st.title('VISUALIZACIÓN MOLECULAR')
    st.write("Bienvenido. Aquí podrás ver la molecula en su forma tridimensional")
    
    seleccion_molecula = st.selectbox("Seleccione una opción: ", ["SMILES", "Subir un archivo"])
    def otros_parametros(mol):
      
      m = Chem.MolFromSmiles(compound_smiles)
      #tpsa = Descriptors.TPSA(m)
      logP = Descriptors.MolLogP(m)
      #st.write("TPSA: " + str(tpsa))
      st.metric("Log P", str(logP))
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
        def sdf_data(content):
          lines = content.splitlines()
          res = {}
          aux = ""
          aux2 = ""
          app_line = False
          for l in lines:
              if app_line:
                  if l == "":
                      app_line = False
                      res[aux] = aux2[1:]
                      aux2 = ""
                  else:
                      aux2 = aux2 + "\n" + l
              elif len(l) > 0 and l[0] == ">":
                  aux = l[3:-1]
                  app_line = True
          return res 

        for uploaded_file in uploaded_files:
          stringio = StringIO(uploaded_file.getvalue().decode("utf-8"))
          string_data = stringio.read()          
          propiedades = sdf_data(string_data)

          #st.write(propiedades)

          if "ROTATIONAL.CONSTANTS" in propiedades:
            st.metric("Constantes de rotacion: ", propiedades['ROTATIONAL.CONSTANTS'])
        
          if "ELECTRONIC.ENERGY" in propiedades: 
            st.metric("Energía electronica: ", propiedades['ELECTRONIC.ENERGY'])
        
          if "DIPOLE.MOMENT" in propiedades:
            st.metric("Momento dipolar: ", propiedades['DIPOLE.MOMENT'])        
        
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

        compounds = pubchempy.get_compounds(compound_smiles, namespace='smiles')
        match = compounds[0]
        st.subheader(match.iupac_name)

        MolecularFormula = pubchempy.get_properties('MolecularFormula',compound_smiles, namespace='smiles')
        mol_weight = pubchempy.get_properties('MolecularWeight',compound_smiles, namespace='smiles')
        ExactMass = pubchempy.get_properties('ExactMass',compound_smiles, namespace='smiles')
        MonoisotopicMass = pubchempy.get_properties('MonoisotopicMass',compound_smiles, namespace='smiles')
        tpsa1 = pubchempy.get_properties('TPSA',compound_smiles, namespace='smiles')
        RotatableBondCount = pubchempy.get_properties('RotatableBondCount',compound_smiles, namespace='smiles')
        HBondDonorCount = pubchempy.get_properties('HBondDonorCount',compound_smiles, namespace='smiles')
        HBondAcceptorCount = pubchempy.get_properties('HBondAcceptorCount',compound_smiles, namespace='smiles')
        FeatureHydrophobeCount3D = pubchempy.get_properties('FeatureHydrophobeCount3D',compound_smiles, namespace='smiles')
        Complexity = pubchempy.get_properties('Complexity',compound_smiles, namespace='smiles')
        Charge = pubchempy.get_properties('Charge',compound_smiles, namespace='smiles')
        CovalentUnitCount = pubchempy.get_properties('CovalentUnitCount',compound_smiles, namespace='smiles')
        FeatureAcceptorCount3D = pubchempy.get_properties('FeatureAcceptorCount3D',compound_smiles, namespace='smiles')
        FeatureDonorCount3D = pubchempy.get_properties('FeatureDonorCount3D',compound_smiles, namespace='smiles')
        FeatureRingCount3D = pubchempy.get_properties('FeatureRingCount3D',compound_smiles, namespace='smiles')

        col1, col2, col3 = st.columns(3)
        col1.metric("Molecular Formula", MolecularFormula[0]['MolecularFormula'])
        col2.metric("Exact Mass", ExactMass[0]['ExactMass'])
        col3.metric("Monoisotopic Mass", MonoisotopicMass[0]['MonoisotopicMass'])

        col1, col2, col3 = st.columns(3)
        col1.metric("Molecular weight", mol_weight[0]['MolecularWeight'])
        col2.metric("TPSA", str(round(tpsa1[0]['TPSA'], 4)))
        col3.metric("Complexity", str(Complexity[0]['Complexity']))

        col1, col2, col3 = st.columns(3)
        col1.metric("Rotable bond", str(RotatableBondCount[0]['RotatableBondCount']))
        col2.metric("H bond donor", str(HBondDonorCount[0]['HBondDonorCount']))
        col3.metric("H bond acceptor", str(HBondAcceptorCount[0]['HBondAcceptorCount']))

        col1, col2, col3 = st.columns(3)
        col1.metric("Ring count", str(FeatureRingCount3D[0]['FeatureRingCount3D']))
        col2.metric("Acceptor count", str(FeatureAcceptorCount3D[0]['FeatureAcceptorCount3D']))
        col3.metric("Donor count", str(FeatureDonorCount3D[0]['FeatureDonorCount3D']))

        col1, col2, col3= st.columns(3)
        col1.metric("Covalent unit", str(CovalentUnitCount[0]['CovalentUnitCount']))
        col2.metric("Hydrophobe count", str(FeatureHydrophobeCount3D[0]['FeatureHydrophobeCount3D']))
        col3.metric("Charge", Charge[0]['Charge'])

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

      compounds = pubchempy.get_compounds(compound_smiles, namespace='smiles')
      match = compounds[0]
      st.subheader(match.iupac_name)

      MolecularFormula = pubchempy.get_properties('MolecularFormula',compound_smiles, namespace='smiles')
      mol_weight = pubchempy.get_properties('MolecularWeight',compound_smiles, namespace='smiles')
      ExactMass = pubchempy.get_properties('ExactMass',compound_smiles, namespace='smiles')
      MonoisotopicMass = pubchempy.get_properties('MonoisotopicMass',compound_smiles, namespace='smiles')
      tpsa1 = pubchempy.get_properties('TPSA',compound_smiles, namespace='smiles')
      RotatableBondCount = pubchempy.get_properties('RotatableBondCount',compound_smiles, namespace='smiles')
      HBondDonorCount = pubchempy.get_properties('HBondDonorCount',compound_smiles, namespace='smiles')
      HBondAcceptorCount = pubchempy.get_properties('HBondAcceptorCount',compound_smiles, namespace='smiles')
      FeatureHydrophobeCount3D = pubchempy.get_properties('FeatureHydrophobeCount3D',compound_smiles, namespace='smiles')
      Complexity = pubchempy.get_properties('Complexity',compound_smiles, namespace='smiles')
      Charge = pubchempy.get_properties('Charge',compound_smiles, namespace='smiles')
      CovalentUnitCount = pubchempy.get_properties('CovalentUnitCount',compound_smiles, namespace='smiles')
      FeatureAcceptorCount3D = pubchempy.get_properties('FeatureAcceptorCount3D',compound_smiles, namespace='smiles')
      FeatureDonorCount3D = pubchempy.get_properties('FeatureDonorCount3D',compound_smiles, namespace='smiles')
      FeatureRingCount3D = pubchempy.get_properties('FeatureRingCount3D',compound_smiles, namespace='smiles')

      blk=makeblock(compound_smiles)
      render_mol(blk)
      
      col1, col2, col3 = st.columns(3)
      col1.metric("Molecular Formula", MolecularFormula[0]['MolecularFormula'])
      col2.metric("Exact Mass", ExactMass[0]['ExactMass'])
      col3.metric("Monoisotopic Mass", MonoisotopicMass[0]['MonoisotopicMass'])

      col1, col2, col3 = st.columns(3)
      col1.metric("Molecular weight", mol_weight[0]['MolecularWeight'])
      col2.metric("TPSA", str(round(tpsa1[0]['TPSA'], 4)))
      col3.metric("Complexity", str(Complexity[0]['Complexity']))

      col1, col2, col3 = st.columns(3)
      col1.metric("Rotable bond", str(RotatableBondCount[0]['RotatableBondCount']))
      col2.metric("H bond donor", str(HBondDonorCount[0]['HBondDonorCount']))
      col3.metric("H bond acceptor", str(HBondAcceptorCount[0]['HBondAcceptorCount']))

      col1, col2, col3 = st.columns(3)
      col1.metric("Ring count", str(FeatureRingCount3D[0]['FeatureRingCount3D']))
      col2.metric("Acceptor count", str(FeatureAcceptorCount3D[0]['FeatureAcceptorCount3D']))
      col3.metric("Donor count", str(FeatureDonorCount3D[0]['FeatureDonorCount3D']))

      col1, col2, col3= st.columns(3)
      col1.metric("Covalent unit", str(CovalentUnitCount[0]['CovalentUnitCount']))
      col2.metric("Hydrophobe count", str(FeatureHydrophobeCount3D[0]['FeatureHydrophobeCount3D']))
      col3.metric("Charge", Charge[0]['Charge'])

      otros_parametros(compound_smiles)

