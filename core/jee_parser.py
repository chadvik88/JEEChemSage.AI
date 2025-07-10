# core/jee_parser.py

from core.name_to_smiles import name_to_smiles
from rdkit import Chem

def _translate(term):
    """
    Attempts to convert a chemical name or SMILES into valid SMILES.
    """
    term = term.strip()
    smiles = name_to_smiles(term)

    if smiles:
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                print(f"[INFO] Translated name '{term}' to SMILES: {smiles}")
                return smiles
            else:
                print(f"[❌ FAIL] RDKit could not parse SMILES: {smiles}")
        except:
            print(f"[❌ EXCEPTION] RDKit failed on: {smiles}")
    else:
        print(f"[❌ FAIL] Could not resolve: '{term}'")

    return None


def parse_reaction_string(reaction_str):
    """
    Parses a reaction string like 'ethyl bromide + KOH -> ethene + HBr'
    into lists of valid SMILES strings for reactants and products.
    Returns: ([reactant_smiles], [product_smiles])
    """
    try:
        print(f"[DEBUG] Parsing reaction: {reaction_str}")
        if '→' in reaction_str:
            reactants_str, products_str = reaction_str.split('→')
        elif '->' in reaction_str:
            reactants_str, products_str = reaction_str.split('->')
        else:
            raise ValueError("Reaction must contain '->' or '→'")

        print(f"[DEBUG] Reactants: {reactants_str.strip()}, Products: {products_str.strip()}")

        reactants = [_translate(r) for r in reactants_str.split('+')]
        products  = [_translate(p) for p in products_str.split('+')]

        # Filter out anything that failed translation
        reactants = [r for r in reactants if r]
        products  = [p for p in products if p]

        print(f"[DEBUG] Parsed Reactants SMILES: {reactants}")
        print(f"[DEBUG] Parsed Products SMILES: {products}")

        return reactants, products

    except Exception as e:
        print(f"[Parser Error] {e}")
        return [], []
