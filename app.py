import streamlit as st
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D
import pandas as pd
import py3Dmol

# -------------------------------
# PAGE CONFIG
# -------------------------------
st.set_page_config(page_title="Epothilone B Analyzer", layout="wide")

# -------------------------------
# THEME TOGGLE
# -------------------------------
theme = st.toggle("🌗 Dark / Light Mode", value=True)

bg = "#0e1117" if theme else "#ffffff"
text = "white" if theme else "black"

st.markdown(f"""
<style>
body {{ background-color: {bg}; color: {text}; }}

.card {{
    padding: 20px;
    border-radius: 12px;
    background-color: {"#1e1e1e" if theme else "#f5f5f5"};
    margin-bottom: 20px;
    transition: transform 0.3s ease;
}}

.card:hover {{
    transform: scale(1.05);
}}

.title {{
    font-size: 30px;
    font-weight: bold;
    animation: fadeIn 1s ease-in;
}}

@keyframes fadeIn {{
    from {{opacity: 0; transform: translateY(-10px);}}
    to {{opacity: 1; transform: translateY(0);}}
}}
</style>
""", unsafe_allow_html=True)

# -------------------------------
# TITLE
# -------------------------------
st.markdown('<div class="title">🧬 Epothilone B Stereochemistry Analyzer</div>', unsafe_allow_html=True)

# -------------------------------
# INFO
# -------------------------------
col1, col2 = st.columns(2)

with col1:
    st.markdown('<div class="card">', unsafe_allow_html=True)
    st.subheader("🔬 Stereochemistry")
    st.write("3D arrangement of atoms affecting drug behavior.")
    st.markdown('</div>', unsafe_allow_html=True)

with col2:
    st.markdown('<div class="card">', unsafe_allow_html=True)
    st.subheader("💊 Epothilone B")
    st.write("Anticancer drug with multiple chiral centers.")
    st.markdown('</div>', unsafe_allow_html=True)

# -------------------------------
# SMILES
# -------------------------------
smiles = "CC1=C[C@@H]2[C@@H](O)[C@H](OC(=O)C[C@H](C)[C@H](O)C(=O)N[C@@H](C)C(=O)O)[C@@H](O)[C@H](OC)[C@H](C)[C@@H](O)[C@H](C)C(=O)O[C@H]2O1"

# -------------------------------
# BUTTON
# -------------------------------
if st.button("🔍 Analyze Molecule"):

    mol = Chem.MolFromSmiles(smiles)

    if mol is None:
        st.error("Invalid SMILES")
    else:
        mol = Chem.AddHs(mol)
        Chem.AssignStereochemistry(mol, force=True)

        # 3D coordinates
        AllChem.EmbedMolecule(mol)
        AllChem.UFFOptimizeMolecule(mol)

        centers = Chem.FindMolChiralCenters(mol, includeUnassigned=False)
        highlight_atoms = [idx for idx, _ in centers]

        col1, col2 = st.columns([1, 1])

        # -------------------------------
        # LEFT: DATA
        # -------------------------------
        with col1:
            st.subheader("🧪 Chiral Centers")

            data = []
            for idx, config in centers:
                atom = mol.GetAtomWithIdx(idx)
                st.write(f"Atom {idx} ({atom.GetSymbol()}): **{config}**")

                data.append({
                    "Atom Index": idx,
                    "Element": atom.GetSymbol(),
                    "Configuration": config
                })

            st.success(f"Total chiral centers: {len(centers)}")

            df = pd.DataFrame(data)
            st.dataframe(df)

        # -------------------------------
        # RIGHT: 2D + 3D
        # -------------------------------
        with col2:
            st.subheader("🧬 2D Structure")

            drawer = rdMolDraw2D.MolDraw2DSVG(350, 350)
            drawer.DrawMolecule(mol, highlightAtoms=highlight_atoms)
            drawer.FinishDrawing()

            svg = drawer.GetDrawingText()
            st.image(svg)

            # -------------------------------
            # 3D VIEW WITH ANIMATION
            # -------------------------------
            st.subheader("🌐 3D Animated View")

            mol_block = Chem.MolToMolBlock(mol)

            view = py3Dmol.view(width=400, height=400)
            view.addModel(mol_block, "mol")

            # base style
            view.setStyle({"stick": {}})

            # highlight chiral atoms (animated spheres)
            for atom_idx in highlight_atoms:
                view.addStyle(
                    {"serial": atom_idx + 1},
                    {"sphere": {"color": "red", "radius": 0.6}}
                )

            # rotate animation
            view.spin(True)

            view.zoomTo()

            st.components.v1.html(view._make_html(), height=400)

        # -------------------------------
        # EXPLANATION
        # -------------------------------
        st.markdown('<div class="card">', unsafe_allow_html=True)
        st.subheader("📘 R/S Configuration")
        st.write("R = Clockwise, S = Counterclockwise.")
        st.markdown('</div>', unsafe_allow_html=True)

        st.markdown('<div class="card">', unsafe_allow_html=True)
        st.subheader("📘 Chirality")
        st.write("Chiral molecules are non-superimposable mirror images.")
        st.markdown('</div>', unsafe_allow_html=True)

# -------------------------------
# FOOTER
# -------------------------------
st.write("---")
st.markdown("""
### 👨‍🎓 Student Details

**Name:** I MOHAMMED ABIDH  
**Register Number:** RA2511026050042  
**Class:** AIML - A  

---
Developed using RDKit + Streamlit
""")