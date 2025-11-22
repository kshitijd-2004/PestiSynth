# main.py
import os
import sys 

from functools import lru_cache
from typing import List, Optional

from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PLAPT_DIR = os.path.join(BASE_DIR, "WELP-PLAPT")

if PLAPT_DIR not in sys.path:
    sys.path.append(PLAPT_DIR)

from plapt import Plapt


# =========================
# Domain config (in-memory)
# =========================

PESTS = {
    "mosquito": {
        "name": "Mosquito",
        "sequence": "MNSAPRNAGAARSVDRRGFIAACGLLVLLVVRMLGAADATRIFDIENIPVVKYVAVAGRNVTLNCPGVTEHSLVDTLVWRTSQTIAEYVDGGLPLVSSLRITLLPDNFSLHFNPAFASDTDEYSCLVNDRHSPDALVDLLVQDVPDPPGRPLVLSFTSRTVNLSWAYTQDPRNAPINNFVIETRVGENGEWDQVDPLYTNSNEAFYQVTGLVPFTVYSFRVIAVNELGHSPPSKESYYFVTLREAPTGKPVTTIAHNTSATSVYISWKPPPAETILGEFLGYRITYRARDKGTDDDVKEIYIRDSTVESHEIHHLETYTQYLASIQVFNPEGLGPPTTVLVMTDEGVPSKPLNLSVLEVTSTTIKITWREPEKLNGAIHGYRVYYVHQNQTLLHLPILKADAAVNSVYTYTLSNLKPYTDYKIIVAAFTKKFDGEPSEVSQRTDISGPSAPKVVNLTCHSQHELLFGWHIPQTYYNTIDYYIISYRNLEHSEYKDIRLTANSSIVETSMIIPNLTTNMVYEVKVRAASASVINPKQIILGSYSEPKKISLQLHCEKIPQPSQRQVYDDYNLAVLGGIVFSCFGLLLIVLSFLLWKKCFHAAYYYLDDPPACQGANTAGLIDWEAPCEVAGEVRGPIPVTEFAKHVASLHVDGDIGFSKEYEAIQGEALNDEYPSEHSQHPDNKGKNRYLNVIAYDHSRVHLRQVPGQKKHLDYINANYIDGYQRPRSFIGTQGPLPGTFDCFWRMIWEQRVAIIVMITNLVERGRRKCDMYWPKDGAEMYGVIQVRLIREDVMATYTVRTLQIKHTKLKKKKASQSEKLVYQYHYTNWPDHGTPDHPLPVINFVKKSTAANPSDAGPIVVHCSAGVGRTGTYIVLDAMLKQIESKGMLNVFGFLRYIRAQRNYLVQTEEQYIFIHDALVEAIDSGETNIKMDAIGGLVNNIDFIDNQYKLITSYQPKEINLTSALKSVNAIKNRSSLVPLEGSRVHLTPKPGVEGSDYINASWLHGFRRLRDFVVTQHPLIETFKDFWQMVWDHNAQTVVLLSSADNMSFLQFWPNESEPMESDYYRIRMVSETSENNYIVRNFVIQSIQDDYELSVRMFENPMWPDMANPRSIYDFAVRVHERCAQYRNGPIVVVDRYGGFQACQFCCISSLAMQLEYDQTANVYTYAKLYHNKRPGIWSSYEDIRQIYRILSYMPKDLGLLKCTELRTEFDDAAIMTATPDLYSKICSNGSINTHLNSGDGGGNGNDGVPTGNGTNGGLPMSGGGTTTAATIQNGGTVIVKMNGEDNDELSVVVATSNHLNLDHNQS"
    },
    "locust": {
        "name": "Locust",
        "sequence": "MSEQUENCEFORLOCUSTACHE..."
    }
}

PESTICIDES = {
    "imidacloprid": {
        "name": "Imidacloprid",
        "smiles": "CC1=NC(N)=NC(=N1)N2C=NC3=CC(Cl)=C(Cl)C=C23",
    },
    "fipronil": {
        "name": "Fipronil",
        "smiles": "CN=C(NCC1CCOC1)N[N+](=O)[O-]",
    },
    "malathion": {
        "name": "Malathion",
        "smiles": "CCOC(=O)CC(OP(=S)(OC(C)C)OC(C)C)C(=O)OCC",
    },
}


# =========================
# Pydantic models
# =========================

class PestInfo(BaseModel):
    id: str
    name: str


class PesticideInfo(BaseModel):
    id: str
    name: str
    smiles: str


class MoleculeInput(BaseModel):
    name: Optional[str] = None
    pesticide_id: Optional[str] = None
    smiles: Optional[str] = None


class BatchScoreRequest(BaseModel):
    pest_id: str
    molecules: List[MoleculeInput]


class MoleculeResult(BaseModel):
    name: str
    smiles: str
    score: float
    label: str
    interpretation: str


class BatchScoreResponse(BaseModel):
    pest: str
    results: List[MoleculeResult]


# =========================
# PLAPT wrapper
# =========================

@lru_cache(maxsize=1)
def get_plapt_model() -> Plapt:
    model_path = os.path.join(PLAPT_DIR, "models", "affinity_predictor.onnx")
    return Plapt(prediction_module_path=model_path)


def score_single_molecule(protein_seq: str, smiles: str) -> tuple[float, str, float]:
    plapt_model = get_plapt_model()

    results = plapt_model.score_candidates(protein_seq, [smiles])
    result = results[0]

    affinity_uM = float(result["affinity_uM"])
    # neg_log10 = float(result["neg_log10_affinity_M"])

    if affinity_uM <= 0.1:
        label = "High"
    elif affinity_uM <= 1.0:
        label = "Medium"
    else:
        label = "Low"

    return affinity_uM, label



# =========================
# FastAPI app
# =========================

app = FastAPI(
    title="Pest Pesticide Affinity API",
    description="Predicts binding affinity between pest proteins and pesticide molecules using PLAPT.",
    version="0.1.0",
)

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)


@app.get("/health")
def health():
    return {"status": "ok"}


@app.get("/pests", response_model=List[PestInfo])
def list_pests():
    return [
        PestInfo(id=pest_id, name=pest["name"])
        for pest_id, pest in PESTS.items()
    ]


@app.get("/pesticides", response_model=List[PesticideInfo])
def list_pesticides():
    return [
        PesticideInfo(id=pid, name=p["name"], smiles=p["smiles"])
        for pid, p in PESTICIDES.items()
    ]


@app.post("/score/batch", response_model=BatchScoreResponse)
def score_batch(req: BatchScoreRequest):
    pest = PESTS.get(req.pest_id)
    if pest is None:
        raise HTTPException(status_code=400, detail="Unknown pest_id")

    protein_seq = pest["sequence"]

    if not req.molecules:
        raise HTTPException(status_code=400, detail="At least one molecule is required.")

    if len(req.molecules) > 5:
        raise HTTPException(status_code=400, detail="Maximum of 5 molecules allowed.")

    results: List[MoleculeResult] = []

    for mol in req.molecules:
        smiles = None
        display_name = mol.name

        if mol.pesticide_id:
            pesticide = PESTICIDES.get(mol.pesticide_id)
            if pesticide is None:
                raise HTTPException(status_code=400, detail=f"Unknown pesticide_id: {mol.pesticide_id}")
            smiles = pesticide["smiles"]
            if display_name is None:
                display_name = pesticide["name"]

        if mol.smiles:
            smiles = mol.smiles.strip()

        if not smiles:
            raise HTTPException(status_code=400, detail="Each molecule needs pesticide_id or smiles.")

        if not display_name:
            display_name = "Unnamed molecule"

        score_value, label = score_single_molecule(protein_seq, smiles)

        interpretation = (
            f"Predicted {label.lower()} binding to the protein sequence associated with {pest['name']}."
        )

        results.append(
            MoleculeResult(
                name=display_name,
                smiles=smiles,
                score=score_value,
                label=label,
                interpretation=interpretation,
            )
        )

    return BatchScoreResponse(
        pest=pest["name"],
        results=results,
    )
