from typing import List, Dict, Optional

from rdkit import Chem
from rdkit.Chem import DataStructs

from .safety_db import BANNED_PESTICIDES


def smiles_to_fp(smiles: str):
    """
    Convert SMILES to an RDKit fingerprint.
    Returns None if parsing fails.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return Chem.RDKFingerprint(mol)


def similarity(sm1: str, sm2: str) -> float:
    """
    Tanimoto similarity between two SMILES strings.
    Returns 0.0 if either can't be parsed.
    """
    fp1 = smiles_to_fp(sm1)
    fp2 = smiles_to_fp(sm2)
    if fp1 is None or fp2 is None:
        return 0.0
    return DataStructs.FingerprintSimilarity(fp1, fp2)


def safety_scan(user_smiles: str, threshold: float = 0.7) -> List[Dict]:
    """
    Compare user SMILES to banned/restricted pesticides list.
    Returns a list of dicts: {id, name, similarity}, sorted by similarity desc.
    """
    hits: List[Dict] = []

    for key, info in BANNED_PESTICIDES.items():
        sim = similarity(user_smiles, info["smiles"])
        if sim >= threshold:
            hits.append(
                {
                    "id": key,
                    "name": info["name"],
                    "similarity": float(sim),
                }
            )

    hits.sort(key=lambda x: x["similarity"], reverse=True)
    return hits


# ============================
# Structural toxicophore alerts
# ============================

_TOX_ALERT_DEFS: List[Dict] = [
    {
        "id": "organophosphate",
        "name": "Organophosphate / phosphate ester",
        "smarts": "P(=O)(O)(O)",
        "description": (
            "Phosphate ester motif often seen in compounds with "
            "acetylcholinesterase inhibition and acute neurotoxicity."
        ),
    },
    {
        "id": "organothiophosphate",
        "name": "Organothiophosphate",
        "smarts": "P(=S)(O)(O)",
        "description": (
            "Thio-phosphate motif found in many legacy insecticides with "
            "nervous-system toxicity and metabolic activation risks."
        ),
    },
    {
        "id": "aromatic_nitro",
        "name": "Aromatic nitro group",
        "smarts": "[N+](=O)[O-]",
        "description": (
            "Nitroaromatic groups can be associated with mutagenicity and "
            "long-term toxicity in some chemical classes."
        ),
    },
    {
        "id": "halogenated_aromatic",
        "name": "Halogenated aromatic ring",
        "smarts": "[Cl,Br,F,I]c1ccccc1",
        "description": (
            "Halogenated aromatics can be persistent and bioaccumulative in "
            "the environment, raising eco-toxicity concerns."
        ),
    },
]

# Pre-compile SMARTS patterns once at import time
for alert in _TOX_ALERT_DEFS:
    pattern = Chem.MolFromSmarts(alert["smarts"])
    alert["pattern"] = pattern


def analyze_structural_alerts(user_smiles: str) -> List[Dict]:
    """
    Runs a simple toxicophore screen on the input SMILES.
    Returns a list of dicts: {id, name, description} for each alert that hits.
    If the SMILES cannot be parsed, returns [].
    """
    mol = Chem.MolFromSmiles(user_smiles)
    if mol is None:
        return []

    hits: List[Dict] = []
    for alert in _TOX_ALERT_DEFS:
        patt = alert.get("pattern")
        if patt is None:
            continue
        if mol.HasSubstructMatch(patt):
            hits.append(
                {
                    "id": alert["id"],
                    "name": alert["name"],
                    "description": alert["description"],
                }
            )

    return hits


def summarize_safety(
    user_smiles: str,
    banned_hits: List[Dict],
    alert_hits: List[Dict],
) -> Dict:
    """
    Combine banned-similarity hits + structural alerts into a simple, robust
    safety summary:
      - flag: short human-readable string
      - level: High / Medium / Low / Minimal / Unknown
      - score: 0–1 heuristic (None if SMILES invalid)
      - max_similarity: max similarity to banned list
    """
    mol = Chem.MolFromSmiles(user_smiles)
    if mol is None:
        # Don't pretend we're confident if we can't even parse the SMILES
        return {
            "flag": "⚠ SMILES could not be parsed – safety screen not reliable",
            "level": "Unknown",
            "score": None,
            "max_similarity": 0.0,
        }

    max_sim = max((h["similarity"] for h in banned_hits), default=0.0)
    num_alerts = len(alert_hits)

    # Simple heuristic: similarity + alert count => risk score [0,1]
    base = max_sim  # structural similarity to banned list
    score = min(1.0, base + 0.2 * num_alerts)  # each alert adds 0.2, capped at 1.0

    if max_sim >= 0.9 or num_alerts >= 2:
        level = "High"
        flag = "⚠ High structural toxicity concern"
    elif max_sim >= 0.8 or num_alerts == 1:
        level = "Medium"
        flag = "⚠ Moderate structural toxicity concern"
    elif max_sim >= 0.7:
        level = "Low"
        flag = "⚠ Weak similarity to restricted chemistries"
    else:
        level = "Minimal"
        flag = "No strong structural toxicity alerts found"

    return {
        "flag": flag,
        "level": level,
        "score": float(score),
        "max_similarity": float(max_sim),
    }
