# app/streamlit_app.py

import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import streamlit as st
from core.functional_groups import detect_functional_groups
from core.mcq_generator import generate_mcqs
from core.jee_parser import parse_reaction_string
from core.molecule_renderer import get_rdkit_mol, draw_2d, draw_3d
from core.properties import get_molecular_properties
from core.explainer import classify_mechanism, get_explanation
from core.reaction_classifier import classify_reaction
from core.mechanism_visuals import get_mechanism_image_path

# --- Initialize session state ---
if 'reaction_history' not in st.session_state:
    st.session_state.reaction_history = []

# --- Page config ---
st.set_page_config(page_title="JeeChemSage.AI", layout="wide")
st.title(" ⌬ JeeChemSage.AI")
st.subheader("Visual Chemistry Tutor for IIT-JEE")

# --- Sidebar: Reaction History ---
st.sidebar.markdown("## 🕘 Reaction History")
for past in reversed(st.session_state.reaction_history):
    if st.sidebar.button(past, key=past):
        st.session_state.reaction_input = past

# --- Chapter-wise samples ---
chapters = {
    "Haloalkanes & Haloarenes": {
        "SN2 Reaction": "ethyl bromide + KOH -> ethene + HBr",
        "SN1 Reaction": "tert-butyl bromide + H2O -> tert-butyl alcohol + HBr",
        "Finkelstein Reaction": "CH3Cl + NaI -> CH3I + NaCl"
    },
    "Alcohols, Phenols & Ethers": {
        "Dehydration (E1)": "ethanol + H2SO4 -> ethene + H2O",
        "Williamson Ether Synthesis": "ethyl bromide + sodium ethoxide -> ethyl ethyl ether + NaBr"
    },
    "Amines": {
        "Acylation": "aniline + acetyl chloride -> acetanilide + HCl",
        "Hoffmann Rearrangement": "benzamide + Br2 + KOH -> aniline + KBr + CO2"
    }
}

st.markdown("### Choose a Chapter + Sample Reaction")
selected_chapter = st.selectbox("Select Chapter", list(chapters.keys()))
sample_rxn = st.selectbox("Sample Reaction", list(chapters[selected_chapter].values()))

# Input
reaction_input = st.text_input("Or enter your own reaction", st.session_state.get('reaction_input', sample_rxn)).strip()

if reaction_input:
    try:
        # Save history
        if reaction_input not in st.session_state.reaction_history:
            st.session_state.reaction_history.append(reaction_input)

        # Parse
        reactants, products = parse_reaction_string(reaction_input)

        if not reactants or not products:
            st.warning("Could not parse reactants/products. Please check the format.")
            st.stop()

        # Display reactants
        st.markdown("### Reactants")
        for smi in reactants:
            if not smi: continue
            c1, c2 = st.columns(2)
            mol = get_rdkit_mol(smi)
            with c1:
                st.image(draw_2d(mol), caption=smi)
                gr = detect_functional_groups(smi)
                if gr: st.markdown(f"**FG:** {', '.join(gr)}")
            with c2:
                st.components.v1.html(draw_3d(smi)._repr_html_(), height=300)

        # Display products
        st.markdown("### Products")
        for smi in products:
            if not smi: continue
            c1, c2 = st.columns(2)
            mol = get_rdkit_mol(smi)
            with c1:
                st.image(draw_2d(mol), caption=smi)
                gr = detect_functional_groups(smi)
                if gr: st.markdown(f"**FG:** {', '.join(gr)}")
            with c2:
                st.components.v1.html(draw_3d(smi)._repr_html_(), height=300)

        # Properties
        st.markdown("### Properties: Reactants")
        for smi in reactants:
            props = get_molecular_properties(smi)
            if props:
                st.markdown(f"**{smi}**")
                st.json(props)

        st.markdown("### Properties: Products")
        for smi in products:
            props = get_molecular_properties(smi)
            if props:
                st.markdown(f"**{smi}**")
                st.json(props)

        # Mechanism classification
        st.markdown("### Mechanism Prediction")
        mechanism = classify_mechanism(reactants, products)
        st.success(f"**Mechanism:** {mechanism}")

        # Explanation toggle
        simple = st.checkbox("Explain like I'm 15")
        st.markdown("### 💡 Explanation")
        explanation = get_explanation(mechanism, simple)
        st.info(explanation if explanation else "No simple analogy available yet.")

        # Reaction classifier ML + fallback
        st.markdown("### ML-Predicted Reaction Type")
        reaction_type = classify_reaction(reactants, products)
        st.success(f"**Reaction Class:** {reaction_type}")

        # Visual mechanism
        st.markdown("### 📷 Visual Mechanism Diagram")
        img = get_mechanism_image_path(mechanism)
        st.text(f"[DEBUG] Mechanism image path: {img}")
        if img and os.path.exists(img):
            st.image(img, caption=mechanism)
        else:
            st.warning("No diagram available for this mechanism.")

        # MCQs within same block
        st.markdown("### Test Yourself!")
        mcqs = generate_mcqs(reactants, products, mechanism, reaction_type)
        for idx, q in enumerate(mcqs):
            ans = st.radio(q['question'], q['options'], key=f"mcq_{idx}")
            if ans:
                if ans == q['answer']:
                    st.success("Correct ☑")
                else:
                    st.error(f"Incorrect ✖ — Correct: {q['answer']}")

    except Exception as e:
        st.error(f"Error: {e}")
