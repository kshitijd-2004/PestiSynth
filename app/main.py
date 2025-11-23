import os
import sys
from functools import lru_cache
from typing import List, Optional

from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel, Field

# Local modules
from .pests import PESTS, PEST_TARGET_PATHWAY
from .pesticides import PESTICIDES
from .safety_db import BANNED_PESTICIDES
from .safety import (
    safety_scan,
    analyze_structural_alerts,
    summarize_safety,
)

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PLAPT_DIR = os.path.join(BASE_DIR, "WELP-PLAPT")

if PLAPT_DIR not in sys.path:
    sys.path.append(PLAPT_DIR)

from plapt import Plapt  # type: ignore


# =========================
# Pydantic models
# =========================

class PestInfo(BaseModel):
    id: str
    name: str
    pathway: str


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


class SafetyAlert(BaseModel):
    id: str
    name: str
    description: str


class MoleculeResult(BaseModel):
    name: str
    smiles: str
    score: float
    label: str
    interpretation: str

    # Old safety fields (still used by frontend)
    safety_matches: List[SafetyMatch]
    safety_flag: str

    # New, more detailed safety fields
    safety_level: Optional[str] = None       # High / Medium / Low / Minimal / Unknown
    safety_score: Optional[float] = None     # 0–1 heuristic score
    safety_alerts: List[SafetyAlert] = Field(default_factory=list)

    # Selectivity metrics
    selectivity_index: Optional[float] = None
    best_off_target_pest: Optional[str] = None
    best_off_target_affinity: Optional[float] = None


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
    title="PestiSynth API",
    description="Predicts pesticide–protein affinity + safety scan using PLAPT and chemical similarity.",
    version="1.2.0",
)

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],  # relax for local dev / demo
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
        PestInfo(
            id=pest_id,
            name=pest["name"],
            pathway=PEST_TARGET_PATHWAY.get(pest_id, "Unknown"),
        )
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
        raise HTTPException(
            status_code=400,
            detail="At least one molecule is required.",
        )

    if len(req.molecules) > 5:
        raise HTTPException(
            status_code=400,
            detail="Maximum of 5 molecules allowed.",
        )

    results: List[MoleculeResult] = []

    for mol in req.molecules:
        smiles: Optional[str] = None
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

        # Custom SMILES overrides (or provides) the structure
        if mol.smiles:
            smiles = mol.smiles.strip()

        if not smiles:
            raise HTTPException(
                status_code=400,
                detail="Each molecule requires pesticide_id or smiles input.",
            )

        if not display_name:
            display_name = "Unnamed molecule"

        # 1) Affinity Prediction for the selected target pest
        affinity_value, label = score_single_molecule(protein_seq, smiles)

        # 1b) Selectivity analysis: compare to other pests
        best_off_target_affinity: Optional[float] = None
        best_off_target_name: Optional[str] = None

        for other_pest_id, other_pest in PESTS.items():
            if other_pest_id == req.pest_id:
                continue  # skip chosen target pest

            other_protein_seq = other_pest["sequence"]
            other_affinity, _ = score_single_molecule(other_protein_seq, smiles)

            if best_off_target_affinity is None or other_affinity < best_off_target_affinity:
                best_off_target_affinity = other_affinity
                best_off_target_name = other_pest["name"]

        # Compute selectivity index, if we evaluated any off-targets
        if best_off_target_affinity is not None and affinity_value > 0:
            selectivity_index = best_off_target_affinity / affinity_value
        else:
            selectivity_index = None

        # 2) Safety – similarity to banned list + structural alerts
        safety_hits = safety_scan(smiles)  # [{name, similarity}, ...]
        alert_hits = analyze_structural_alerts(smiles)  # [{id, name, description}, ...]

        summary = summarize_safety(smiles, safety_hits, alert_hits)

        safety_flag = summary.get("flag", "No safety assessment available")
        safety_level = summary.get("level")
        safety_score = summary.get("score")

        safety_models = [
            SafetyMatch(name=h["name"], similarity=h["similarity"])
            for h in safety_hits
        ]

        safety_alert_models = [
            SafetyAlert(
                id=a["id"],
                name=a["name"],
                description=a["description"],
            )
            for a in alert_hits
        ]

        # 3) Assemble result
        interpretation = (
            f"Predicted {label.lower()} binding to {pest['name']} protein target. "
            "SafetyScan includes similarity to banned/restricted pesticides plus "
            "structural toxicophore alerts."
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
                safety_level=safety_level,
                safety_score=safety_score,
                safety_alerts=safety_alert_models,
                selectivity_index=selectivity_index,
                best_off_target_pest=best_off_target_name,
                best_off_target_affinity=best_off_target_affinity,
            )
        )

    return BatchScoreResponse(
        pest=pest["name"],
        results=results,
    )
