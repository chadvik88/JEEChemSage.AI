# core/molecule_renderer.py

from rdkit import Chem
from rdkit.Chem import Draw
import py3Dmol

def get_rdkit_mol(smiles):
    return Chem.MolFromSmiles(smiles)

def draw_2d(mol):
    """Returns a PIL image of the 2D structure."""
    if mol is None:
        raise ValueError("Invalid SMILES or null molecule passed to draw_2d()")
    return Draw.MolToImage(mol, size=(300, 300))

def draw_3d(smiles):
    """Returns a Py3Dmol view for interactive 3D display."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    mblock = Chem.MolToMolBlock(mol)
    view = py3Dmol.view(width=400, height=300)
    view.addModel(mblock, 'mol')
    view.setStyle({'stick': {}})
    view.zoomTo()
    return view
