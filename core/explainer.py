# core/explainer.py

def classify_mechanism(reactants, products):
    """
    Basic logic for classifying mechanisms.
    """
    if "CCBr" in reactants and "[OH-]" in reactants and "C=C" in products:
        return "E2 Elimination"
    if "CCBr" in reactants and "[OH-]" in reactants and "CCOH" in products:
        return "SN2 Substitution"
    return "Unknown"

def get_explanation(mechanism):
    if mechanism == "SN2 Substitution":
        return (
            "This is an SN2 reaction. Hydroxide ion attacks the electrophilic carbon, "
            "causing the bromide to leave. It proceeds in one step with inversion of configuration."
        )
    if mechanism == "E2 Elimination":
        return (
            "This is an E2 reaction. The base (OH⁻) abstracts a proton while the leaving group (Br⁻) departs, "
            "resulting in the formation of a double bond."
        )
    return "Mechanism unknown. Try a supported reaction."
