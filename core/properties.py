# core/properties.py

from rdkit import Chem
from rdkit.Chem import Descriptors

def get_molecular_properties(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    return {
        "Molecular Weight": round(Descriptors.MolWt(mol), 2),
        "LogP": round(Descriptors.MolLogP(mol), 2),
        "Num H-Donors": Descriptors.NumHDonors(mol),
        "Num H-Acceptors": Descriptors.NumHAcceptors(mol),
        "TPSA": round(Descriptors.TPSA(mol), 2),
        "Rotatable Bonds": Descriptors.NumRotatableBonds(mol)
    }
