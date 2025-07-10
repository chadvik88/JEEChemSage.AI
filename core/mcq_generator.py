def generate_mcqs(reactants, products, mechanism, reaction_type):
    mcqs = []
    mechanism = mechanism.strip().lower()

    if "e2" in mechanism:
        mcqs.append({
            "question": "What type of base is typically involved in E2 reactions?",
            "options": ["Weak base", "Strong base", "Acid catalyst", "Radical initiator"],
            "answer": "Strong base"
        })
        mcqs.append({
            "question": "How many steps are involved in an E2 reaction?",
            "options": ["One", "Two", "Three", "Depends on substrate"],
            "answer": "One"
        })

    elif "sn2" in mechanism:
        mcqs.append({
            "question": "What is the stereochemical outcome of an SN2 reaction?",
            "options": ["Racemization", "Retention", "Inversion", "No change"],
            "answer": "Inversion"
        })

    elif "sn1" in mechanism:
        mcqs.append({
            "question": "What intermediate is formed in an SN1 reaction?",
            "options": ["Carbocation", "Carbanion", "Radical", "None"],
            "answer": "Carbocation"
        })

    elif "e1" in mechanism:
        mcqs.append({
            "question": "What is the first step in an E1 reaction?",
            "options": ["Loss of leaving group", "Base attack", "Radical formation", "None"],
            "answer": "Loss of leaving group"
        })

    return mcqs
