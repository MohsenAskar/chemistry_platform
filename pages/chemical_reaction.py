import streamlit as st
from rdkit import Chem
from rdkit.Chem import AllChem, rdChemReactions
import streamlit.components.v1 as components

st.set_page_config(layout="wide")

def mol_to_html(mol, size=(400, 400)):
    """
    Generates an interactive 3D HTML view for an RDKit molecule.
    """
    try:
        # Create a copy of the molecule to avoid modifying the original
        mol_copy = Chem.Mol(mol)
        
        # Sanitize the molecule
        Chem.SanitizeMol(mol_copy)
        
        # Calculate implicit valence before adding hydrogens
        for atom in mol_copy.GetAtoms():
            atom.UpdatePropertyCache(strict=False)
        
        # Add hydrogens with coordinates
        mol_copy = Chem.AddHs(mol_copy, addCoords=True)
        
        # Generate 3D coordinates with error handling
        try:
            # Try ETKDG first
            AllChem.EmbedMolecule(mol_copy, randomSeed=42, useRandomCoords=True)
        except:
            # If ETKDG fails, try with different parameters
            AllChem.EmbedMolecule(mol_copy, randomSeed=42, useRandomCoords=True, maxAttempts=5000)
        
        # Energy minimization
        try:
            AllChem.MMFFOptimizeMolecule(mol_copy)
        except:
            # If MMFF fails, try UFF
            AllChem.UFFOptimizeMolecule(mol_copy)
        
        # Create py3Dmol view
        import py3Dmol
        view = py3Dmol.view(width=size[0], height=size[1])
        pdb = Chem.MolToPDBBlock(mol_copy)
        view.addModel(pdb, "pdb")
        view.setStyle({'stick':{}, 'atoms':{'showLabels':True}})
        view.zoomTo()
        return view._make_html()
    except Exception as e:
        st.error(f"Error generating 3D structure: {str(e)}")
        return None

def run_reaction(r1_smiles, r2_smiles, reaction_smarts):
    """
    Runs a chemical reaction with improved error handling.
    """
    try:
        # Create reaction from SMARTS
        reaction = rdChemReactions.ReactionFromSmarts(reaction_smarts)
        
        # Create molecules from SMILES with sanitization
        r1 = Chem.MolFromSmiles(r1_smiles, sanitize=True)
        r2 = Chem.MolFromSmiles(r2_smiles, sanitize=True)
        
        if not r1 or not r2:
            return None, "One or both reactant SMILES are invalid."
            
        # Update property cache for both reactants
        for mol in [r1, r2]:
            for atom in mol.GetAtoms():
                atom.UpdatePropertyCache(strict=False)
        
        # Run the reaction
        product_sets = reaction.RunReactants((r1, r2))
        
        if not product_sets:
            return None, "Reaction did not yield any products."
            
        # Sanitize and optimize products
        sanitized_products = []
        for product in product_sets[0]:
            try:
                Chem.SanitizeMol(product)
                sanitized_products.append(product)
            except Exception as e:
                st.warning(f"Warning: Could not sanitize a product: {str(e)}")
                
        return sanitized_products if sanitized_products else None, None
        
    except Exception as e:
        return None, f"Error running reaction: {str(e)}"

# Streamlit App Layout
st.title("Interactive Chemical Reaction Visualizer")
st.write("This app demonstrates a chemical reaction by visualizing the reactants and products in 3D.")


with st.expander("About SMILES Notation"):
    st.write("""
    SMILES (Simplified Molecular Input Line Entry System) is a way to represent chemical structures using text.
    Some basic examples:
    - CC = ethane
    - CCO = ethanol
    - CC(=O)O = acetic acid
    - c1ccccc1 = benzene
    
    If you're entering custom molecules, make sure to use correct SMILES notation.
    """)
# Default values for an esterification reaction
default_r1 = "CC(=O)O"  # Acetic acid
default_r2 = "CCO"      # Ethanol
default_reaction_smarts = "[C:1](=[O:2])[O:3].[O:4][C:5]>>[C:1](=[O:2])[O:4][C:5]"

r1_smiles = st.text_input("Reactant 1 SMILES", default_r1)
r2_smiles = st.text_input("Reactant 2 SMILES", default_r2)
reaction_smarts = st.text_input("Reaction SMARTS", default_reaction_smarts)

if st.button("Run Reaction"):
    st.subheader("Reactants")
    col1, col2 = st.columns(2)
    
    # Display reactants
    with col1:
        st.markdown("**Reactant 1:**")
        mol1 = Chem.MolFromSmiles(r1_smiles, sanitize=True)
        if mol1:
            html1 = mol_to_html(mol1)
            if html1:
                st.components.v1.html(html1, width=500, height=500)
        else:
            st.error("Invalid Reactant 1 SMILES")
            
    with col2:
        st.markdown("**Reactant 2:**")
        mol2 = Chem.MolFromSmiles(r2_smiles, sanitize=True)
        if mol2:
            html2 = mol_to_html(mol2)
            if html2:
                st.components.v1.html(html2, width=500, height=500)
        else:
            st.error("Invalid Reactant 2 SMILES")
    
    # Run and display reaction products
    st.subheader("Products")
    product_mols, error = run_reaction(r1_smiles, r2_smiles, reaction_smarts)
    
    if error:
        st.error(error)
    elif product_mols:
        num_products = len(product_mols)
        cols = st.columns(num_products)
        for i, prod in enumerate(product_mols):
            with cols[i]:
                st.markdown(f"**Product {i+1}:**")
                html_prod = mol_to_html(prod)
                if html_prod:
                    st.components.v1.html(html_prod, width=500, height=500)
    else:
        st.error("No products generated.")


st.title("Try these reactions ⌬")
st.markdown("""
**Some tips for using these:**
- Copy and paste the SMILES strings carefully
- Make sure to include all brackets and special characters
- The 3D visualization might take a few seconds to load
- You can rotate the molecules by clicking and dragging
- Use the scroll wheel to zoom in/out

**Esterification (Formation of Ethyl Acetate)**
- Reactant 1 (Acetic acid): CC(=O)O
- Reactant 2 (Ethanol): CCO
- Reaction SMARTS: [C:1](=[O:2])[O:3].[O:4][C:5]>>[C:1](=[O:2])[O:4][C:5]


**Amide Formation**
- Reactant 1 (Acetic acid): CC(=O)O
- Reactant 2 (Methylamine): CN
- Reaction SMARTS: [C:1](=[O:2])[O:3].[N:4][C:5]>>[C:1](=[O:2])[N:4][C:5]


**Aldol Condensation**
- Reactant 1 (Acetaldehyde): CC=O
- Reactant 2 (Acetaldehyde): CC=O
- Reaction SMARTS: [C:1][C:2](=[O:3]).[C:4][C:5](=[O:6])>>[C:1][C:2](=[O:3])[C:4][C:5](O)[C:6]


**Formation of Aspirin**
- Reactant 1 (Salicylic acid): O=C(O)c1ccccc1O
- Reactant 2 (Acetic anhydride): CC(=O)OC(=O)C
- eaction SMARTS: [O:1][c:2]1[c:3][c:4][c:5][c:6][c:7]1[C:8](=[O:9])[O:10].[C:11](=[O:12])[O:13][C:14](=[O:15])[C:16]>>[C:11](=[O:12])[O:1][c:2]1[c:3][c:4][c:5][c:6][c:7]1[C:8](=[O:9])[O:10]


**Simple Alcohol Oxidation**
- Reactant 1 (Ethanol): CCO
- Reactant 2 (Oxygen): O=O
- Reaction SMARTS: [C:1][C:2][O:3].[O:4]=[O:5]>>[C:1][C:2]=[O:3]


**Formation of Ibuprofen Salt**
- Reactant 1 (Ibuprofen): CC(C)Cc1ccc(cc1)[C@H](C)C(=O)O
- Reactant 2 (Sodium hydroxide): [Na+][OH-]
- Reaction SMARTS: [C:1](=[O:2])[O:3].[Na+][OH-]>>[C:1](=[O:2])[O-].[Na+]


**Simple Alkene Formation (Dehydration)**
- Reactant 1 (2-butanol): CCC(C)O
- Reactant 2 (Acid catalyst): [H+]
- Reaction SMARTS: [C:1][C:2][C:3]([C:4])[O:5].[H+]>>[C:1][C:2]=[C:3][C:4]


**Formation of Methyl Benzoate**
- Reactant 1 (Benzoic acid): O=C(O)c1ccccc1
- Reactant 2 (Methanol): CO
- Reaction SMARTS: [C:1](=[O:2])[O:3].[O:4][C:5]>>[C:1](=[O:2])[O:4][C:5]
""")