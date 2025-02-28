import streamlit as st
from rdkit import Chem
from rdkit.Chem import AllChem
import py3Dmol
st.set_page_config(layout="wide")

def visualize_molecule(smiles, size=(400, 400)):
    """
    Generates a 3D molecular visualization from a SMILES string.

    Parameters:
      smiles (str): The SMILES representation of the molecule.
      size (tuple): Width and height for the visualization.

    Returns:
      view (py3Dmol.view): The interactive 3D view of the molecule.
    """
    # Create an RDKit molecule object from SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        st.error("Invalid SMILES string!")
        return None

    # Generate 3D coordinates
    mol = Chem.AddHs(mol)  # Add hydrogen atoms
    AllChem.EmbedMolecule(mol, randomSeed=42)  # Generate 3D coordinates
    AllChem.MMFFOptimizeMolecule(mol)  # Energy minimize the structure

    # Create a Py3DMol view
    view = py3Dmol.view(width=size[0], height=size[1])

    # Convert the molecule to PDB format
    pdb = Chem.MolToPDBBlock(mol)

    # Add the molecule to the view
    view.addModel(pdb, "pdb")

    # Set the style of the visualization and zoom
    view.setStyle({'stick': {}})
    view.zoomTo()

    return view

# Streamlit UI
st.title("Interactive Molecule Visualizer")
st.write("Enter a SMILES string to generate an interactive 3D molecular visualization.")

# Use a default SMILES string for aspirin
default_smiles = "CC(=O)OC1=CC=CC=C1C(=O)O"
smiles_input = st.text_input("SMILES Input", default_smiles)

if st.button("Visualize Molecule"):
    view = visualize_molecule(smiles_input)
    if view:
        # Get the HTML representation of the Py3DMol view.
        # Note: _make_html() is used here; depending on your version, you may need to adjust this.
        html = view._make_html()
        st.components.v1.html(html, width=500, height=500)

