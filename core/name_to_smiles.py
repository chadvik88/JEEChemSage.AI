from rdkit import Chem
from difflib import get_close_matches

# Expandable dictionary
NAME_TO_SMILES = {
    "ethyl bromide": "CCBr",
    "bromoethane": "CCBr",
    "ethene": "C=C",
    "ethylene": "C=C",
    "hydroxide": "[OH-]",
    "koh": "[K+].[OH-]",
    "hbr": "Br",
    "bromide": "Br",
    "phenol": "c1ccc(cc1)O",
    "benzene": "c1ccccc1"
}

def name_to_smiles(name):
    name = name.strip().lower()
    if name in NAME_TO_SMILES:
        return NAME_TO_SMILES[name]
    close = get_close_matches(name, NAME_TO_SMILES.keys(), n=1, cutoff=0.75)
    if close:
        return NAME_TO_SMILES[close[0]]
    return None
