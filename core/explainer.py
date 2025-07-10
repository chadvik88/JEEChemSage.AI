from core.functional_groups import detect_functional_groups

# Classifies reaction into known mechanism names

def classify_mechanism(reactants, products):
    fg_reactants = []
    for smi in reactants:
        fg_reactants += detect_functional_groups(smi)

    fg_products = []
    for smi in products:
        fg_products += detect_functional_groups(smi)

    fg_reactants = set(fg_reactants)
    fg_products = set(fg_products)

    # === Rule-Based Mechanism Classification ===
    if "alkyl halide" in fg_reactants and "alkene" in fg_products:
        return "E2 Elimination Reaction"

    if "alkyl halide" in fg_reactants and "alcohol" in fg_products:
        return "SN1 Substitution Reaction" if "H2O" in reactants or "water" in reactants else "SN2 Substitution Reaction"

    if "amide" in fg_reactants and "amine" in fg_products:
        return "Hoffmann Rearrangement"

    if "amine" in fg_reactants and "amide" in fg_products:
        return "Acylation"

    if "alcohol" in fg_reactants and "alkene" in fg_products:
        return "E1 Dehydration"

    if "alkyl halide" in fg_reactants and "ether" in fg_products:
        return "Williamson Ether Synthesis"

    return "Unknown Organic Reaction (rule-based fallback)"

# Gives child-level explanation
def get_explanation(mechanism, simple=False):
    explanations = {
        "E2 Elimination Reaction": {
            "normal": (
                "E2 is a single-step elimination reaction where a strong base removes a β-hydrogen as the leaving group exits. "
                "This leads to the formation of a double bond without forming a carbocation."
            ),
            "simple": (
                "Imagine a Lego model where you pull off one brick (hydrogen) and another brick (Br) pops off at the same time, "
                "and the rest snaps into a double bond that too, all in one smooth move."
            )
        },
        "SN2 Substitution Reaction": {
            "normal": (
                "SN2 is a one-step reaction where the nucleophile attacks from the opposite side of the leaving group, "
                "causing a direct substitution and inversion of stereochemistry."
            ),
            "simple": (
                "It’s like someone enters from the back door (nucleophile) while someone else exits the front (leaving group). "
                "Quick, one motion with the house flipping direction!"
            )
        },
        "SN1 Substitution Reaction": {
            "normal": (
                "SN1 is a two-step reaction where the leaving group first departs, forming a carbocation. "
                "Then the nucleophile attacks, possibly from either side."
            ),
            "simple": (
                "Like a chair gets empty (Br leaves), and someone else (H2O) walks in and sits down with a short pause in between!"
            )
        },
        "E1 Dehydration": {
            "normal": (
                "E1 involves two steps: the loss of water to form a carbocation, followed by elimination of a β-hydrogen "
                "to form a double bond. It’s typically acid-catalyzed."
            ),
            "simple": (
                "Think of water leaving first like a spilled drink, then the molecule ‘dries up’ into a double bond."
            )
        }
    }

    # Exact match needed for your reaction labels
    entry = explanations.get(mechanism.strip())
    if not entry:
        return None

    return entry["simple"] if simple else entry["normal"]
