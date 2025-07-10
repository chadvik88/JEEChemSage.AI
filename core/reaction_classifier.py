from core.functional_groups import detect_functional_groups

def classify_reaction(reactants, products):
    try:
        all_smiles = ' + '.join(reactants) + ' -> ' + ' + '.join(products)
        print(f"[DEBUG] Classifying Reaction: {all_smiles}")

        fg_reactants = []
        for smi in reactants:
            fg_reactants += detect_functional_groups(smi)

        fg_products = []
        for smi in products:
            fg_products += detect_functional_groups(smi)

        fg_reactants = set(fg_reactants)
        fg_products = set(fg_products)

        print(f"[DEBUG] FG Reactants: {fg_reactants}")
        print(f"[DEBUG] FG Products: {fg_products}")

        # --- Rule-Based Classification ---
        if 'alkyl halide' in fg_reactants and 'alkene' in fg_products:
            return "E2 Elimination Reaction"

        if 'alkyl halide' in fg_reactants and 'alcohol' in fg_products:
            if 'H2O' in reactants or 'water' in reactants:
                return "SN1 Substitution Reaction"
            else:
                return "SN2 Substitution Reaction"

        if 'amide' in fg_reactants and 'amine' in fg_products:
            return "Hoffmann Rearrangement"

        if 'amine' in fg_reactants and 'amide' in fg_products:
            return "Acylation Reaction"

        if 'alcohol' in fg_reactants and 'alkene' in fg_products:
            return "E1 Dehydration"

        if 'alkyl halide' in fg_reactants and 'ether' in fg_products:
            return "Williamson Ether Synthesis"

        return "Unknown Organic Reaction (rule-based fallback)"

    except Exception as e:
        print(f"[ERROR] classify_reaction failed: {e}")
        return "Unknown (error)"
