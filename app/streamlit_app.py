# app/streamlit_app.py
import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
import streamlit as st
from core.jee_parser import parse_reaction_string
from core.molecule_renderer import get_rdkit_mol, draw_2d, draw_3d
from core.properties import get_molecular_properties
from core.explainer import classify_mechanism, get_explanation

st.set_page_config(page_title="JeeChemSage.AI", layout="wide")
st.title("JeeChemSage.AI")
st.subheader("Visual Chemistry Tutor for IIT-JEE")

reaction_input = st.text_input("Enter a reaction (SMILES)", "CCBr + [OH-] -> C=C + Br")

if reaction_input:
    try:
        reactants, products = parse_reaction_string(reaction_input)

        st.markdown("### Reactants")
        for smi in reactants:
            col1, col2 = st.columns(2)
            mol = get_rdkit_mol(smi)
            with col1:
                st.image(draw_2d(mol), caption=smi)
            with col2:
                st.components.v1.html(draw_3d(smi)._repr_html_(), height=300)

        st.markdown("### Products")
        for smi in products:
            col1, col2 = st.columns(2)
            mol = get_rdkit_mol(smi)
            with col1:
                st.image(draw_2d(mol), caption=smi)
            with col2:
                st.components.v1.html(draw_3d(smi)._repr_html_(), height=300)

        st.markdown("### Properties (Reactants Only)")
        for smi in reactants:
            props = get_molecular_properties(smi)
            if props:
                st.markdown(f"**{smi}**")
                st.json(props)

        st.markdown("### Reaction Mechanism")
        mechanism = classify_mechanism(reactants, products)
        st.success(f"**Detected:** {mechanism}")

        st.markdown("### Explanation")
        st.info(get_explanation(mechanism))

    except Exception as e:
        st.error(f"Error: {e}")
