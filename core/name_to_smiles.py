# core/name_to_smiles.py

from difflib import get_close_matches
from pubchempy import get_compounds
import functools

# Extended Mappings (✅ Add more as needed)
NAME_TO_SMILES = {
    # Acids/Bases/Reagents
    "hcl": "[H+].[Cl-]",
    "hbr": "[H+].[Br-]",
    "h2so4": "OS(=O)(=O)O",
    "koh": "[K+].[OH-]",
    "naoh": "[Na+].[OH-]",
    "hno3": "O=N(=O)O",
    "nh3": "N",
    "br2": "BrBr",

    # Alkenes/Alkanes
    "ethene": "C=C",
    "ethylene": "C=C",
    "ethane": "CC",

    # Alcohols
    "ethanol": "CCO",
    "methanol": "CO",
    "tert-butyl alcohol": "CC(C)(C)O",

    # Alkyl Halides
    "ethyl bromide": "CCBr",
    "methyl bromide": "CBr",
    "bromoethane": "CCBr",
    "chloroethane": "CCCl",
    "ch3cl": "CCl",
    "ch3i": "CI",
    "tert-butyl bromide": "CC(C)(C)CBr",

    # Aromatics
    "benzene": "c1ccccc1",
    "phenol": "c1ccc(cc1)O",
    "aniline": "c1ccc(cc1)N",
    "acetanilide": "CC(=O)Nc1ccccc1",
    "toluene": "Cc1ccccc1",
    "benzamide": "NC(=O)c1ccccc1",

    # Ethers
    "ethyl ethyl ether": "CCOCC",
    "diethyl ether": "CCOCC",
    "ether": "COC",
    "sodium ethoxide": "CC[O-].[Na+]",

    # Others
    "acetyl chloride": "CC(=O)Cl",
    "water": "O",
    "h2o": "O",
    "co2": "O=C=O",
}

# Local match (first-pass)
def local_lookup(name: str):
    name = name.lower().strip()
    if name in NAME_TO_SMILES:
        print(f"[INFO] Translated name '{name}' to SMILES (local): {NAME_TO_SMILES[name]}")
        return NAME_TO_SMILES[name]
    matches = get_close_matches(name, NAME_TO_SMILES.keys(), n=1, cutoff=0.8)
    if matches:
        print(f"[INFO] Approx match '{name}' → '{matches[0]}'")
        return NAME_TO_SMILES[matches[0]]
    return None

# Cached PubChem fallback
@functools.lru_cache(maxsize=128)
def pubchem_lookup(name: str):
    try:
        compounds = get_compounds(name, 'name')
        if compounds:
            return compounds[0].canonical_smiles
    except Exception as e:
        print(f"[PubChem Error] {e}")
    return None

# Combined interface
def name_to_smiles(name: str):
    name = name.strip().lower()
    return local_lookup(name) or pubchem_lookup(name)
