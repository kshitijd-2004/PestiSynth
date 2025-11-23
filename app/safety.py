from rdkit import Chem
from rdkit.Chem import DataStructs

def similarity(sm1: str, sm2: str) -> float:
    m1 = Chem.MolFromSmiles(sm1)
    m2 = Chem.MolFromSmiles(sm2)
    if m1 is None or m2 is None:
        return 0.0
    fp1 = Chem.RDKFingerprint(m1)
    fp2 = Chem.RDKFingerprint(m2)
    return float(DataStructs.FingerprintSimilarity(fp1, fp2))
