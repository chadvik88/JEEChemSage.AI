# core/functional_groups.py

from rdkit import Chem

# Optional: dictionary if you want to show descriptions/UI labels
FUNCTIONAL_GROUPS = {
    "alcohol": "[CX4][OH]",
    "ether": "COC",
    "carboxylic acid": "C(=O)[OH]",
    "amine": "[NX3;H2,H1;!$(NC=O)]",
    "amide": "C(=O)N",
    "ketone": "C(=O)C",
    "aldehyde": "[CX3H1](=O)[#6]",
    "halide": "[F,Cl,Br,I]",
    "nitro": "[$([NX3](=O)=O)]",
    "benzene ring": "c1ccccc1",
    "alkyl halide": "[CX4][Br,Cl,I,F]",
    "alkene": "C=C",
}

def detect_functional_groups(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return []

    groups = []
    for group, smarts in FUNCTIONAL_GROUPS.items():
        if mol.HasSubstructMatch(Chem.MolFromSmarts(smarts)):
            groups.append(group)

    return groups
