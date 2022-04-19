#$ conda activate myenv
#$ source activate myenv

import streamlit as st
#import py3Dmol
#from stmol import showmol
#from psikit import Psikit
#import numpy as np
#import psi4

#from rdkit import Chem
#from rdkit.Chem import AllChem

#pk = Psikit()
#pk.read_from_smiles('O')
#pk.energy()
#pk.create_cube_files(gridspace=0.5)

#homo_voldata = open("Psi_a_5_1-A\"_HOMO.cube", "r").read()
#lumo_voldata = open("Psi_a_6_5-A\'_LUMO.cube", "r").read()

#v = py3Dmol.view()
#v.addVolumetricData(homo_voldata, "cube", {'isoval': -0.03, 'color': "red", 'opacity': 0.75})
#v.addVolumetricData(homo_voldata, "cube", {'isoval': 0.03, 'color': "blue", 'opacity': 0.75})
#v.addModel(Chem.MolToMolBlock(pk.mol), 'mol')
#v.setStyle({'stick':{}})
#v.zoomTo()
#v.show()

st.title('RDKit + Py3DMOL 😀')
#st.write('Hello world')


#def makeblock(smi):
    #mol = Chem.MolFromSmiles(smi)
    #mol = Chem.AddHs(mol)
    #AllChem.EmbedMolecule(mol)
    #mblock = Chem.MolToMolBlock(mol)
    #return mblock

#def render_mol(xyz):
    #xyzview = py3Dmol.view()(width=1000,height=1000)
    #xyzview.addModel(xyz,'mol')
    #xyzview.setStyle({'stick':{}})
    #xyzview.setBackgroundColor('white')
    #xyzview.zoomTo()
    #showmol(xyzview,height=500,width=500)

#compound_smiles=st.text_input('SMILES please','CC')
#blk=makeblock(compound_smiles)
#render_mol(blk)
