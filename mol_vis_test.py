import streamlit as st
import py3Dmol
from stmol import showmol
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

st.title('VISUALIZACIÓN MOLECUALR')

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

compound_smiles=st.text_input('SMILES please','CC')
blk=makeblock(compound_smiles)
render_mol(blk)

