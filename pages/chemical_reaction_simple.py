import streamlit as st
from rdkit import Chem
from rdkit.Chem import AllChem, rdChemReactions
import py3Dmol
st.set_page_config(layout="wide")

# Dictionary of common reactions with their SMARTS patterns and descriptions
REACTIONS = {
    "Esterification": {
        "smarts": "[C:1](=[O:2])[O:3].[O:4][C:5]>>[C:1](=[O:2])[O:4][C:5]",
        "description": "Formation of an ester from a carboxylic acid and an alcohol",
        "reactant1_type": "Carboxylic Acid",
        "reactant2_type": "Alcohol"
    },
    "Acid-Base": {
        "smarts": "[C:1](=[O:2])[O:3].[Na+][OH-]>>[C:1](=[O:2])[O-].[Na+]",
        "description": "Formation of a salt from an acid and a base",
        "reactant1_type": "Acid",
        "reactant2_type": "Base"
    },
    "Amide Formation": {
        "smarts": "[C:1](=[O:2])[O:3].[N:4][C:5]>>[C:1](=[O:2])[N:4][C:5]",
        "description": "Formation of an amide from a carboxylic acid and an amine",
        "reactant1_type": "Carboxylic Acid",
        "reactant2_type": "Amine"
    }
}

# Common molecules dictionary with names, SMILES, and types
COMMON_MOLECULES = {
    "Acetic acid": {"smiles": "CC(=O)O", "type": "Carboxylic Acid"},
    "Ethanol": {"smiles": "CCO", "type": "Alcohol"},
    "Methanol": {"smiles": "CO", "type": "Alcohol"},
    "Propanoic acid": {"smiles": "CCC(=O)O", "type": "Carboxylic Acid"},
    "Methylamine": {"smiles": "CN", "type": "Amine"},
    "Ethylamine": {"smiles": "CCN", "type": "Amine"},
    "Sodium hydroxide": {"smiles": "[Na+][OH-]", "type": "Base"},
    "Benzoic acid": {"smiles": "O=C(O)c1ccccc1", "type": "Carboxylic Acid"},
    "Isopropanol": {"smiles": "CC(C)O", "type": "Alcohol"}
}

def mol_to_html(mol, size=(400, 400)):
    """Generate 3D visualization of molecule with atom labels."""
    try:
        mol_copy = Chem.Mol(mol)
        Chem.SanitizeMol(mol_copy)
        
        for atom in mol_copy.GetAtoms():
            atom.UpdatePropertyCache(strict=False)
        
        mol_copy = Chem.AddHs(mol_copy, addCoords=True)
        
        try:
            AllChem.EmbedMolecule(mol_copy, randomSeed=42, useRandomCoords=True)
        except:
            AllChem.EmbedMolecule(mol_copy, randomSeed=42, useRandomCoords=True, maxAttempts=5000)
        
        try:
            AllChem.MMFFOptimizeMolecule(mol_copy)
        except:
            AllChem.UFFOptimizeMolecule(mol_copy)
        
        # Create the py3Dmol view
        view = py3Dmol.view(width=size[0], height=size[1])
                
        # Convert to PDB format and add model
        pdb = Chem.MolToPDBBlock(mol_copy)
        view.addModel(pdb, "pdb")

        # Basic style for entire molecule
        view.setStyle({'model': -1}, {'stick':{}, 'sphere':{'scale':0.3}})

        # Add atom labels only for non-carbon atoms
        for atom in mol_copy.GetAtoms():
            if atom.GetSymbol() != "C":  # Skip carbons to reduce clutter
                pos = mol_copy.GetConformer().GetAtomPosition(atom.GetIdx())
                view.addLabel(atom.GetSymbol(), 
                            {'position': {'x':pos.x,'y':pos.y,'z':pos.z},
                            'fontSize': 14,
                            'fontColor':'black',
                            'backgroundOpacity': 0.0})

        view.zoomTo()
        return view._make_html()
    except Exception as e:
        st.error(f"Error generating 3D structure: {str(e)}")
        return None

def run_reaction(r1_smiles, r2_smiles, reaction_smarts):
    """Run chemical reaction with error handling."""
    try:
        reaction = rdChemReactions.ReactionFromSmarts(reaction_smarts)
        r1 = Chem.MolFromSmiles(r1_smiles, sanitize=True)
        r2 = Chem.MolFromSmiles(r2_smiles, sanitize=True)
        
        if not r1 or not r2:
            return None, "One or both reactant SMILES are invalid."
            
        for mol in [r1, r2]:
            for atom in mol.GetAtoms():
                atom.UpdatePropertyCache(strict=False)
        
        product_sets = reaction.RunReactants((r1, r2))
        
        if not product_sets:
            return None, "Reaction did not yield any products."
            
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

# Streamlit UI
st.title("Chemistry Student's Reaction Visualizer")
st.write("Select a reaction type and reactants to see the 3D visualization of the reaction.")

# Reaction type selection
reaction_type = st.selectbox(
    "Select Reaction Type",
    list(REACTIONS.keys()),
    help="Choose the type of reaction you want to study"
)

# Show reaction description
st.info(REACTIONS[reaction_type]["description"])

# Create two columns for reactant selection
col1, col2 = st.columns(2)

with col1:
    st.subheader(f"Select {REACTIONS[reaction_type]['reactant1_type']}")
    # Filter molecules by type
    r1_molecules = {name: info for name, info in COMMON_MOLECULES.items() 
                   if info["type"] == REACTIONS[reaction_type]["reactant1_type"]}
    r1_name = st.selectbox("Reactant 1", list(r1_molecules.keys()))
    r1_smiles = r1_molecules[r1_name]["smiles"]
    # Option for custom SMILES
    use_custom_r1 = st.checkbox("Use custom molecule for Reactant 1")
    if use_custom_r1:
        r1_smiles = st.text_input("Enter SMILES for Reactant 1", r1_smiles)

with col2:
    st.subheader(f"Select {REACTIONS[reaction_type]['reactant2_type']}")
    r2_molecules = {name: info for name, info in COMMON_MOLECULES.items() 
                   if info["type"] == REACTIONS[reaction_type]["reactant2_type"]}
    r2_name = st.selectbox("Reactant 2", list(r2_molecules.keys()))
    r2_smiles = r2_molecules[r2_name]["smiles"]
    # Option for custom SMILES
    use_custom_r2 = st.checkbox("Use custom molecule for Reactant 2")
    if use_custom_r2:
        r2_smiles = st.text_input("Enter SMILES for Reactant 2", r2_smiles)

if st.button("Run Reaction"):
    st.subheader("Reactants")
    col1, col2 = st.columns(2)
    
    with col1:
        st.markdown(f"**{REACTIONS[reaction_type]['reactant1_type']}:**")
        mol1 = Chem.MolFromSmiles(r1_smiles, sanitize=True)
        if mol1:
            html1 = mol_to_html(mol1)
            if html1:
                st.components.v1.html(html1, width=500, height=500)
        else:
            st.error("Invalid Reactant 1 SMILES")
            
    with col2:
        st.markdown(f"**{REACTIONS[reaction_type]['reactant2_type']}:**")
        mol2 = Chem.MolFromSmiles(r2_smiles, sanitize=True)
        if mol2:
            html2 = mol_to_html(mol2)
            if html2:
                st.components.v1.html(html2, width=500, height=500)
        else:
            st.error("Invalid Reactant 2 SMILES")
    
    st.subheader("Products")
    product_mols, error = run_reaction(r1_smiles, r2_smiles, REACTIONS[reaction_type]["smarts"])
    
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

# Add helpful information
with st.expander("How to Use This Tool"):
    st.write("""
    1. Select the type of reaction you want to study from the dropdown menu
    2. Choose your reactants from the provided lists of common molecules
    3. Optionally, you can enter custom molecules using SMILES notation
    4. Click 'Run Reaction' to see the 3D visualization of reactants and products
    5. You can rotate the molecules by clicking and dragging
    6. Use the scroll wheel to zoom in/out
    """)

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