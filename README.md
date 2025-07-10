# ⌬ JeeChemSage.AI  
**A Visual Chemistry Tutor for IIT-JEE Aspirants**

JeeChemSage.AI is an interactive chemistry learning tool built to help students master **organic reaction mechanisms** for competitive exams like the **IIT-JEE**. Designed with high school learners in mind, the app visually breaks down complex reactions, identifies functional groups, predicts reaction classes, and even tests your understanding with MCQs — all within a clean Streamlit interface.

---

## 📸 Live Demo

> 🧬 [Try the app (Streamlit link)](https://huggingface.co/spaces/charvik/JEEChemAI)

---

## Features

### 1. Smart Reaction Parsing
- Accepts custom or sample organic reactions like `ethyl bromide + KOH -> ethene + HBr`
- Automatically translates IUPAC/common names to SMILES using RDKit and PubChem fallback
- Parses reactants/products using custom `jee_parser.py`

### 2. Reaction Classification
- Identifies **reaction mechanism** (e.g., E2, SN1, SN2, etc.)
- Detects key **functional groups** (alcohols, alkyl halides, amines, ethers, etc.)
- Displays both **ML-predicted reaction class** and rule-based interpretation

### 3. Molecular Properties (Reactants & Products)
- Molecular weight  
- LogP  
- Hydrogen donors/acceptors  
- Topological Polar Surface Area (TPSA)  
- Rotatable bonds  

### 4. Mechanism Explanation (ELI15 Mode)
- Option to **"Explain like I’m 15"** for each predicted mechanism
- Helps beginners understand the *why* behind each transformation
- Example:
  > “Imagine a Lego model where you pull off one brick (hydrogen) and another (Br) pops off at the same time...”

### 5. Visual Mechanism Diagrams
- Loads mechanism-specific PNG diagrams from `assets/mechanisms/` (like `e2.png`, `sn2.png`)
- Visual walkthrough of how bonds break and form

### 6. "Test Yourself" MCQs
- Auto-generates multiple-choice questions based on the input reaction
- Provides **immediate feedback** on correctness
- Covers mechanism identification, byproducts, nucleophile/electrophile roles, etc.

### 7. Reaction History Sidebar
- Easily re-run or revisit previous inputs via clickable buttons in the sidebar

---

## 🛠 Tech Stack

| Tool            | Use                                       |
|-----------------|--------------------------------------------|
| **Python 3.10+** | Core application logic                   |
| **Streamlit**    | Frontend + UI rendering                  |
| **RDKit**        | Chemistry parsing, SMILES, and drawings  |
| **PubChemPy**    | Backup molecular name resolution         |
| **Pillow/3Dmol.js** | 2D and 3D Molecule Visualization    |

---

## Project Structure

```bash
JeeChemSage.AI/
│
├── app/
│   └── streamlit_app.py         # Main app interface
│
├── core/                        # Core logic modules
│   ├── jee_parser.py            # Reaction parser
│   ├── name_to_smiles.py        # Name to SMILES resolver (RDKit + PubChem fallback)
│   ├── properties.py            # Molecular property calculator
│   ├── reaction_classifier.py   # Rule-based reaction classification
│   ├── functional_groups.py     # SMARTS pattern detection
│   ├── mechanism_visuals.py     # Mechanism image resolver
│   └── mcq_generator.py         # MCQ generator
│
├── assets/
│   └── mechanisms/              # PNG diagrams (e.g., e2.png, sn1.png)
│
└── README.md
```

## Sample Reactions:
  - reaction: "ethyl bromide + KOH -> ethene + HBr"
    type: "E2 Elimination"
  - reaction: "tert-butyl bromide + H2O -> tert-butyl alcohol + HBr"
    type: "SN1"
  - reaction: "ethanol + H2SO4 -> ethene + H2O"
    type: "E1"
  - reaction: "benzamide + Br2 + KOH -> aniline + KBr + CO2"
    type: "Hoffmann Rearrangement"

## Future Ideas:
  - "ML-powered reaction classification (currently rule-based fallback)"
  - "Intelligent MCQ difficulty levels"
  - "Integration with NCERT & Allen chemistry syllabus"
  - "Voice input (e.g., via Whisper API)"
  - "Online hosting via Streamlit Cloud or Hugging Face Spaces"

## Installation:
  - step: "Clone the repository"
    command: "git clone https://github.com/yourusername/JeeChemSage.AI.git"
  - step: "Navigate into the project directory"
    command: "cd JeeChemSage.AI"
  - step: "Install dependencies (recommended inside a virtual environment)"
    command: "pip install -r requirements.txt"
  - step: "Run the Streamlit app"
    command: "streamlit run app/streamlit_app.py"
  - tip: "If you see RDKit or SMILES parsing errors, make sure you’ve installed RDKit correctly"
    command: "conda install -c rdkit rdkit"

## Tagline: "For students, by a student."
