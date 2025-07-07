# core/jee_parser.py

from core.name_to_smiles import name_to_smiles

def _translate(term):
    return name_to_smiles(term.strip()) or term.strip()

def parse_reaction_string(reaction_str):
    """
    Parses a reaction string into reactants and products.
    Accepts natural language or SMILES.
    Example: 'ethyl bromide + KOH -> ethene + HBr'
    Returns: (reactants_list, products_list)
    """
    if '→' in reaction_str:
        reactants, products = reaction_str.split('→')
    elif '->' in reaction_str:
        reactants, products = reaction_str.split('->')
    else:
        raise ValueError("Reaction must contain '->' or '→'")

    reactants = [_translate(r) for r in reactants.split('+')]
    products  = [_translate(p) for p in products.split('+')]

    return reactants, products

