from rdkit import Chem
from rdkit.Chem import DataStructs
from .safety_db import BANNED_PESTICIDES

def smiles_to_fp(smiles: str):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return Chem.RDKFingerprint(mol)

def similarity(sm1: str, sm2: str) -> float:
    fp1 = smiles_to_fp(sm1)
    fp2 = smiles_to_fp(sm2)
    if fp1 is None or fp2 is None:
        return 0.0
    return DataStructs.FingerprintSimilarity(fp1, fp2)

def safety_scan(user_smiles: str, threshold: float = 0.7):
    hits = []
    for key, info in BANNED_PESTICIDES.items():
        sim = similarity(user_smiles, info["smiles"])
        if sim >= threshold:
            hits.append({
                "id": key,
                "name": info["name"],
                "similarity": sim,
            })
    # sort strongest matches first
    hits.sort(key=lambda x: x["similarity"], reverse=True)
    return hits
