import os
import sys

from functools import lru_cache
from typing import List, Optional

from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware

from pydantic import BaseModel

# Local modules
from .pests import PESTS
from .pesticides import PESTICIDES
from .safety_db import BANNED_PESTICIDES
from .safety import safety_scan      

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PLAPT_DIR = os.path.join(BASE_DIR, "WELP-PLAPT")

if PLAPT_DIR not in sys.path:
    sys.path.append(PLAPT_DIR)

from plapt import Plapt


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


class SafetyMatch(BaseModel):
    name: str
    similarity: float


class MoleculeResult(BaseModel):
    name: str
    smiles: str
    score: float
    label: str
    interpretation: str
    safety_matches: List[SafetyMatch]     
    safety_flag: str                      


class BatchScoreRequest(BaseModel):
    pest_id: str
    molecules: List[MoleculeInput]


class BatchScoreResponse(BaseModel):
    pest: str
    results: List[MoleculeResult]


# =========================
# PLAPT Model Loader
# =========================

@lru_cache(maxsize=1)
def get_plapt_model() -> Plapt:
    model_path = os.path.join(PLAPT_DIR, "models", "affinity_predictor.onnx")
    return Plapt(prediction_module_path=model_path)


def score_single_molecule(protein_seq: str, smiles: str) -> tuple[float, str]:
    plapt_model = get_plapt_model()

    results = plapt_model.score_candidates(protein_seq, [smiles])
    result = results[0]

    affinity_uM = float(result["affinity_uM"])

    # Simple classification
    if affinity_uM <= 0.1:
        label = "Low"
    elif affinity_uM <= 1.0:
        label = "Medium"
    else:
        label = "High"

    return affinity_uM, label


# =========================
# FastAPI app
# =========================

app = FastAPI(
    title="PestiSynth API",
    description="Predicts pesticide–protein affinity + safety scan using PLAPT and chemical similarity.",
    version="1.1.0",
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


# =========================
# Main Scoring Endpoint
# =========================

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

        # Predefined pesticide
        if mol.pesticide_id:
            pesticide = PESTICIDES.get(mol.pesticide_id)
            if pesticide is None:
                raise HTTPException(
                    status_code=400,
                    detail=f"Unknown pesticide_id: {mol.pesticide_id}",
                )
            smiles = pesticide["smiles"]
            if display_name is None:
                display_name = pesticide["name"]

        # Custom SMILES
        if mol.smiles:
            smiles = mol.smiles.strip()

        if not smiles:
            raise HTTPException(
                status_code=400,
                detail="Each molecule requires pesticide_id or smiles input.",
            )

        if not display_name:
            display_name = "Unnamed molecule"

        # 1) Affinity Prediction
        affinity_value, label = score_single_molecule(protein_seq, smiles)

        # 2) SafetyScan – RDKit similarity to banned list
        safety_hits = safety_scan(smiles)       # list of dicts

        if safety_hits:
            safety_flag = "⚠ Similar to hazardous or restricted chemicals"
            safety_models = [SafetyMatch(name=h["name"], similarity=h["similarity"])
                             for h in safety_hits]
        else:
            safety_flag = "No known toxicity matches detected"
            safety_models = []

        # 3) Assemble result
        interpretation = (
            f"Predicted {label.lower()} binding to {pest['name']} protein target. "
            "SafetyScan performed for banned/restricted pesticide similarity."
        )

        results.append(
            MoleculeResult(
                name=display_name,
                smiles=smiles,
                score=affinity_value,
                label=label,
                interpretation=interpretation,
                safety_matches=safety_models,
                safety_flag=safety_flag,
            )
        )

    return BatchScoreResponse(
        pest=pest["name"],
        results=results,
    )
