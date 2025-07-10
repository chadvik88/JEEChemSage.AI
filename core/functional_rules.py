from core.functional_groups import detect_functional_groups

def rule_based_classify(reactants, products):
    r_groups = sum([detect_functional_groups(r) for r in reactants], [])
    p_groups = sum([detect_functional_groups(p) for p in products], [])

    r = set(r_groups)
    p = set(p_groups)

    def contains(group): return group in r
    def forms(group): return group in p

    # --- Substitution ---
    if contains("Halide") and forms("Alcohol"):
        return "Nucleophilic Substitution (SN1/SN2)"
    if contains("Halide") and forms("Amine"):
        return "Ammonolysis"
    if contains("Alcohol") and forms("Ether"):
        return "Williamson Ether Synthesis"
    if contains("Aromatic Ring") and forms("Substituted Aromatic Ring"):
        return "Electrophilic Aromatic Substitution"

    # --- Elimination ---
    if contains("Halide") and forms("Alkene"):
        return "Dehydrohalogenation (E1/E2)"
    if contains("Alcohol") and forms("Alkene"):
        return "Dehydration (E1 Mechanism)"

    # --- Addition ---
    if contains("Alkene") and forms("Alcohol"):
        return "Electrophilic Addition (Hydration)"
    if contains("Alkene") and forms("Halide"):
        return "Electrophilic Addition (Halogenation)"
    if contains("Alkyne") and forms("Alkene"):
        return "Partial Reduction of Alkyne"

    # --- Oxidation / Reduction ---
    if contains("Alcohol") and forms("Ketone"):
        return "Oxidation of Secondary Alcohol"
    if contains("Alcohol") and forms("Aldehyde"):
        return "Oxidation of Primary Alcohol"
    if contains("Ketone") and forms("Alcohol"):
        return "Reduction of Ketone"
    if contains("Aldehyde") and forms("Alcohol"):
        return "Reduction of Aldehyde"
    if contains("Carboxylic Acid") and forms("Alcohol"):
        return "Reduction of Carboxylic Acid"

    # --- Rearrangements ---
    if contains("Amide") and forms("Amine") and "CO2" in p:
        return "Hoffmann Rearrangement"
    if contains("Alcohol") and forms("Alkene"):
        return "Pinacol Rearrangement"

    # --- Named Reactions ---
    if contains("Aldehyde") and forms("Alcohol") and "Carboxylic Acid" in p:
        return "Cannizzaro Reaction"
    if contains("Carbonyl") and forms("β-Hydroxy Ketone"):
        return "Aldol Condensation"
    if contains("Amide") and forms("Amine") and not "CO2" in p:
        return "Gabriel Synthesis"
    if contains("Aromatic Ring") and forms("Acylated Aromatic Ring"):
        return "Friedel-Crafts Acylation"
    if contains("Aromatic Ring") and "Diazonium Salt" in p:
        return "Diazotization"
    if contains("Diazonium Salt") and forms("Phenol"):
        return "Hydrolysis of Diazonium"
    if contains("Phenol") and forms("Salicylaldehyde"):
        return "Reimer-Tiemann Reaction"
    if contains("Methyl Ketone") and "CHI3" in p:
        return "Haloform Reaction"
    if contains("Alkyl Halide") and forms("Alkane"):
        return "Wurtz Reaction"
    if contains("Aromatic Amine") and forms("Diazonium Salt"):
        return "Sandmeyer Reaction"

    # --- Functional Group Interconversions ---
    if contains("Carboxylic Acid") and forms("Ester"):
        return "Esterification (Fischer)"
    if contains("Carboxylic Acid") and forms("Amide"):
        return "Amidation"
    if contains("Acid Chloride") and forms("Amide"):
        return "Acylation of Amine"
    if contains("Aldehyde") and forms("Carboxylic Acid"):
        return "Oxidation of Aldehyde"
    if contains("Nitrile") and forms("Carboxylic Acid"):
        return "Nitrile Hydrolysis"

    return "Unknown Organic Reaction (rule-based fallback)"
