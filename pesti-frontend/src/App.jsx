// src/App.jsx
import React, { useEffect, useState } from "react";
import imgFallArmyworm from "./assets/pest_img/Fall_Armyworm.png";
import imgAphids from "./assets/pest_img/Aphids.png";
import imgColoradoPotatoBeetle from "./assets/pest_img/Colorado_Potato_Beetle.png";
import imgBrownPlanthopper from "./assets/pest_img/Brown_Planthopper.png";
import imgMosquitoLarvae from "./assets/pest_img/Mosquito_Larvae.png";
import imgCottonBollworm from "./assets/pest_img/Cotton_Ballworm.png";
import imgDiamondbackMoth from "./assets/pest_img/Diamondback_Moth.png";
import imgRiceStemBorer from "./assets/pest_img/RiceStem_Borer.png";

import imgAcetamiprid from "./assets/pesticide_img/Acetamiprid.png";
import imgBifenthrin from "./assets/pesticide_img/Bifenthrin.png";
import imgChlorpyrifos from "./assets/pesticide_img/Chlorpyrifos.png";
import imgClothianidin from "./assets/pesticide_img/Clothianidin.png";
import imgCypermethrin from "./assets/pesticide_img/Cypermethrin.png";
import imgDeltamethrin from "./assets/pesticide_img/Deltamethrin.png";
import imgDiazinon from "./assets/pesticide_img/Diazinon.png";
import imgDichlorvos from "./assets/pesticide_img/Dichlorvos.png";
import imgDinotefuran from "./assets/pesticide_img/Dinotefuran.png";
import imgFenvalerate from "./assets/pesticide_img/Fenvalerate.png";
import imgImidacloprid from "./assets/pesticide_img/Imidacloprid.png";
import imgLambdaCyhalothrin from "./assets/pesticide_img/Lambda-cyhalothrin.png";
import imgNitenpyram from "./assets/pesticide_img/Nitenpyram.png";
import imgParathion from "./assets/pesticide_img/Parathion.png";
import imgPermethrin from "./assets/pesticide_img/Permethrin.png";
import imgThiacloprid from "./assets/pesticide_img/Thiacloprid.png";
import imgThiamethoxam from "./assets/pesticide_img/Thiamethoxam.png";
import imgMalathion from "./assets/pesticide_img/Malathion.png";

const API_BASE_URL = "http://localhost:8000"; // change if needed

// Optional: map pest/pesticide IDs to image URLs
const PEST_IMAGE_MAP = {
  fall_armyworm: imgFallArmyworm,
  aphid: imgAphids,
  colorado_potato_beetle: imgColoradoPotatoBeetle,
  brown_planthopper: imgBrownPlanthopper,
  mosquito_larvae: imgMosquitoLarvae,
  cotton_bollworm: imgCottonBollworm,
  diamondback_moth: imgDiamondbackMoth,
  rice_stem_borer: imgRiceStemBorer,
};

const PESTICIDE_IMAGE_MAP = {
  imidacloprid: imgImidacloprid,
  thiamethoxam: imgThiamethoxam,
  clothianidin: imgClothianidin,
  acetamiprid_metabolite: imgAcetamiprid,
  nitenpyram: imgNitenpyram,
  dinotefuran: imgDinotefuran,
  thiacloprid_amide: imgThiacloprid,
  permethrin: imgPermethrin,
  cypermethrin: imgCypermethrin,
  deltamethrin: imgDeltamethrin,
  lambda_cyhalothrin: imgLambdaCyhalothrin,
  fenvalerate: imgFenvalerate,
  bifenthrin: imgBifenthrin,
  chlorpyrifos: imgChlorpyrifos,
  diazinon: imgDiazinon,
  dichlorvos: imgDichlorvos,
  parathion: imgParathion,
  phosmet: imgMalathion,
};

const VIEW_SELECTION = "selection";
const VIEW_REPORT = "report";

// max TOTAL molecules per run (library + custom)
const MAX_MOLECULES = 5;

function App() {
  const [pests, setPests] = useState([]);
  const [pesticides, setPesticides] = useState([]);

  const [pestsLoading, setPestsLoading] = useState(true);
  const [pesticidesLoading, setPesticidesLoading] = useState(true);
  const [error, setError] = useState(null);

  const [selectedPestId, setSelectedPestId] = useState(null);
  const [selectedPesticideIds, setSelectedPesticideIds] = useState([]);

  // custom SMILES
  const [customMolecules, setCustomMolecules] = useState([]); // { id, name, smiles }
  const [customInput, setCustomInput] = useState("");
  const [customCounter, setCustomCounter] = useState(1);
  const [fileImportMessage, setFileImportMessage] = useState("");

  const [scoreLoading, setScoreLoading] = useState(false);
  const [scoreError, setScoreError] = useState(null);
  const [scoreResponse, setScoreResponse] = useState(null);

  const [view, setView] = useState(VIEW_SELECTION);

  // search inputs
  const [pestSearch, setPestSearch] = useState("");
  const [pesticideSearch, setPesticideSearch] = useState("");

  useEffect(() => {
    const fetchPests = async () => {
      try {
        setPestsLoading(true);
        const res = await fetch(`${API_BASE_URL}/pests`);
        if (!res.ok) throw new Error(`Failed to fetch pests: ${res.status}`);
        const data = await res.json();
        setPests(data);
      } catch (err) {
        console.error(err);
        setError(err.message || "Failed to load pests.");
      } finally {
        setPestsLoading(false);
      }
    };

    const fetchPesticides = async () => {
      try {
        setPesticidesLoading(true);
        const res = await fetch(`${API_BASE_URL}/pesticides`);
        if (!res.ok)
          throw new Error(`Failed to fetch pesticides: ${res.status}`);
        const data = await res.json();
        setPesticides(data);
      } catch (err) {
        console.error(err);
        setError(err.message || "Failed to load pesticides.");
      } finally {
        setPesticidesLoading(false);
      }
    };

    fetchPests();
    fetchPesticides();
  }, []);

  const totalSelectedMolecules =
    selectedPesticideIds.length + customMolecules.length;

  const handlePestSelect = (pestId) => {
    setSelectedPestId(pestId);
    setScoreResponse(null);
    setScoreError(null);
  };

  const handlePesticideToggle = (pesticideId) => {
    setScoreResponse(null);
    setScoreError(null);

    if (selectedPesticideIds.includes(pesticideId)) {
      setSelectedPesticideIds((prev) =>
        prev.filter((id) => id !== pesticideId)
      );
    } else {
      const totalIfAdded = totalSelectedMolecules + 1;
      if (totalIfAdded > MAX_MOLECULES) return;
      setSelectedPesticideIds((prev) => [...prev, pesticideId]);
    }
  };

  const handleAddCustomSmiles = () => {
    const trimmed = customInput.trim();
    if (!trimmed) return;

    if (totalSelectedMolecules >= MAX_MOLECULES) return;

    const name = `Custom molecule ${customCounter}`;
    const id = `custom-${customCounter}`;

    setCustomMolecules((prev) => [...prev, { id, name, smiles: trimmed }]);
    setCustomCounter((c) => c + 1);
    setCustomInput("");
    setScoreResponse(null);
    setScoreError(null);
  };

  const handleRemoveCustom = (id) => {
    setCustomMolecules((prev) => prev.filter((m) => m.id !== id));
    setScoreResponse(null);
    setScoreError(null);
  };

  const handleFileUpload = (event) => {
    const file = event.target.files && event.target.files[0];
    if (!file) return;

    const reader = new FileReader();
    reader.onload = (e) => {
      const text = String(e.target.result || "");
      const rawTokens = text
        .split(",")
        .map((s) => s.trim())
        .filter((s) => s.length > 0);

      if (rawTokens.length === 0) {
        setFileImportMessage("No valid SMILES found in file.");
        return;
      }

      setCustomMolecules((prev) => {
        let next = [...prev];
        let added = 0;
        let skipped = 0;
        let counter = customCounter;

        for (const smiles of rawTokens) {
          const already =
            next.some((m) => m.smiles === smiles) ||
            selectedPesticideIds.length + next.length >= MAX_MOLECULES;

          if (already) {
            skipped += 1;
            continue;
          }

          const name = `Custom molecule ${counter}`;
          next.push({ id: `custom-${counter}`, name, smiles });
          counter += 1;
          added += 1;

          if (selectedPesticideIds.length + next.length >= MAX_MOLECULES) {
            break;
          }
        }

        setCustomCounter(counter);

        if (added === 0) {
          setFileImportMessage(
            "No additional SMILES added (duplicates or limit reached)."
          );
        } else {
          let msg = `Added ${added} custom molecule${
            added > 1 ? "s" : ""
          } from file.`;
          if (skipped > 0) {
            msg += ` Skipped ${skipped} (duplicates or beyond limit).`;
          }
          setFileImportMessage(msg);
        }

        return next;
      });
    };

    reader.readAsText(file);
  };

  const handleRunScoring = async () => {
    if (!selectedPestId || totalSelectedMolecules === 0) return;

    setScoreLoading(true);
    setScoreError(null);

    try {
      const pesticideMolecules = selectedPesticideIds.map((id) => ({
        pesticide_id: id,
      }));

      const customInputs = customMolecules.map((m) => ({
        name: m.name,
        smiles: m.smiles,
      }));

      const payload = {
        pest_id: selectedPestId,
        molecules: [...pesticideMolecules, ...customInputs],
      };

      const res = await fetch(`${API_BASE_URL}/score/batch`, {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify(payload),
      });

      if (!res.ok) {
        let message = `Scoring request failed (${res.status})`;
        try {
          const errData = await res.json();
          if (errData?.detail) message = errData.detail;
        } catch (_) {}
        throw new Error(message);
      }

      const data = await res.json();
      setScoreResponse(data);
      setView(VIEW_REPORT);
    } catch (err) {
      console.error(err);
      setScoreError(err.message || "Failed to run scoring.");
    } finally {
      setScoreLoading(false);
    }
  };

  const handleBackToSelection = () => {
    setView(VIEW_SELECTION);
  };

  const selectedPest = pests.find((p) => p.id === selectedPestId) || null;

  // link pesticide metadata (for report)
  const selectedPesticideMetas = selectedPesticideIds.map(
    (id) => pesticides.find((p) => p.id === id) || null
  );

  // search/filtering
  const pestQuery = pestSearch.trim().toLowerCase();
  const pesticideQuery = pesticideSearch.trim().toLowerCase();

  const filteredPests = pests.filter((p) =>
    p.name.toLowerCase().includes(pestQuery)
  );

  const filteredPesticides = pesticides.filter((p) =>
    p.name.toLowerCase().includes(pesticideQuery)
  );

  return (
    <div className="app">
      <header className="app-header">
        <div className="app-header-inner">
          <h1 className="app-title">PestiSynth</h1>
          <p className="app-subtitle">
            Predict pesticide–protein affinity and safety with PLAPT.
          </p>
        </div>
      </header>

      {view === VIEW_SELECTION && (
        <>
          <main className="app-main">
            {error && (
              <div className="error-banner">
                <span>⚠</span>
                <span>{error}</span>
              </div>
            )}

            {/* Summary / Selection Status */}
            <section className="summary-section">
              <div className="summary-card">
                <div>
                  <h2>Selection Overview</h2>
                  <p>
                    Choose a target pest and up to {MAX_MOLECULES} molecules
                    (library or custom SMILES) to evaluate.
                  </p>
                </div>

                <div className="summary-pill-row">
                  <div className="summary-pill">
                    <span className="summary-label">Target Pest</span>
                    <span className="summary-value">
                      {selectedPest ? selectedPest.name : "None selected"}
                    </span>
                  </div>

                  <div className="summary-pill">
                    <span className="summary-label">Selected Molecules</span>
                    <span className="summary-value">
                      {totalSelectedMolecules} / {MAX_MOLECULES}
                    </span>
                  </div>
                </div>
              </div>
            </section>

            {/* Pest Selection */}
            <section className="section">
              <div className="section-header">
                <h2>1. Choose a Target Pest</h2>
                <p>
                  Select the pest whose protein you want to evaluate. Only one
                  pest can be active at a time.
                </p>
              </div>

              <div className="search-row">
                <input
                  type="text"
                  className="search-input"
                  placeholder="Search pests by name…"
                  value={pestSearch}
                  onChange={(e) => setPestSearch(e.target.value)}
                />
              </div>

              {pestsLoading ? (
                <div className="loading">Loading pests…</div>
              ) : filteredPests.length === 0 ? (
                <div className="no-results">No pests match your search.</div>
              ) : (
                <div className="card-grid">
                  {filteredPests.map((pest) => (
                    <PestCard
                      key={pest.id}
                      pest={pest}
                      isSelected={pest.id === selectedPestId}
                      onSelect={() => handlePestSelect(pest.id)}
                    />
                  ))}
                </div>
              )}
            </section>

            {/* Pesticide Selection */}
            <section className="section">
              <div className="section-header">
                <h2>2. Choose Molecules from Library</h2>
                <p>
                  Pick molecules from your predefined pesticide library. You can
                  combine these with custom SMILES below.
                </p>
              </div>

              {!selectedPestId && (
                <div className="section-hint">
                  Please select a pest first to continue with molecule selection.
                </div>
              )}

              <div className="search-row">
                <input
                  type="text"
                  className="search-input"
                  placeholder="Search pesticides by name…"
                  value={pesticideSearch}
                  onChange={(e) => setPesticideSearch(e.target.value)}
                />
              </div>

              {pesticidesLoading ? (
                <div className="loading">Loading pesticides…</div>
              ) : filteredPesticides.length === 0 ? (
                <div className="no-results">
                  No pesticides match your search.
                </div>
              ) : (
                <div
                  className={`card-grid ${
                    !selectedPestId ? "card-grid-disabled" : ""
                  }`}
                >
                  {filteredPesticides.map((p) => (
                    <PesticideCard
                      key={p.id}
                      pesticide={p}
                      isSelected={selectedPesticideIds.includes(p.id)}
                      disabled={!selectedPestId}
                      onToggle={() => handlePesticideToggle(p.id)}
                      maxReached={
                        totalSelectedMolecules >= MAX_MOLECULES &&
                        !selectedPesticideIds.includes(p.id)
                      }
                    />
                  ))}
                </div>
              )}
            </section>

            {/* Custom SMILES input */}
            <section className="section">
              <div className="section-header">
                <h2>3. Add Custom SMILES (optional)</h2>
                <p>
                  Include novel molecules by typing a SMILES string or importing
                  a comma-separated list from a .txt file.
                </p>
              </div>

              <div className="custom-card">
                <div className="custom-input-row">
                  <input
                    type="text"
                    className="custom-input"
                    placeholder="e.g. CCOP(=S)(OCC)OC1=CC=C(C=C1)[N+](=O)[O-]"
                    value={customInput}
                    onChange={(e) => setCustomInput(e.target.value)}
                    onKeyDown={(e) => {
                      if (e.key === "Enter") {
                        e.preventDefault();
                        handleAddCustomSmiles();
                      }
                    }}
                  />
                  <button
                    type="button"
                    className="secondary-button custom-add-button"
                    onClick={handleAddCustomSmiles}
                    disabled={
                      !customInput.trim() ||
                      totalSelectedMolecules >= MAX_MOLECULES
                    }
                  >
                    Add
                  </button>
                </div>

                <div className="custom-file-row">
                  <label className="custom-file-label">
                    <span className="custom-file-title">
                      Upload SMILES from .txt
                    </span>
                    <span className="custom-file-subtitle">
                      Use comma-separated SMILES. We’ll import up to{" "}
                      {MAX_MOLECULES - selectedPesticideIds.length} additional
                      molecules.
                    </span>
                    <input
                      type="file"
                      accept=".txt"
                      onChange={handleFileUpload}
                    />
                  </label>
                </div>

                {fileImportMessage && (
                  <p className="custom-file-message">{fileImportMessage}</p>
                )}

                {customMolecules.length > 0 && (
                  <div className="custom-list">
                    {customMolecules.map((m) => (
                      <div key={m.id} className="custom-chip">
                        <div className="custom-chip-main">
                          <span className="custom-chip-name">{m.name}</span>
                          <code className="custom-chip-smiles">
                            {m.smiles}
                          </code>
                        </div>
                        <button
                          type="button"
                          className="custom-chip-remove"
                          onClick={() => handleRemoveCustom(m.id)}
                        >
                          ×
                        </button>
                      </div>
                    ))}
                  </div>
                )}

                <p className="custom-limit-note">
                  Total molecules this run:{" "}
                  <strong>{totalSelectedMolecules}</strong> / {MAX_MOLECULES}
                </p>
              </div>

              {totalSelectedMolecules >= MAX_MOLECULES && (
                <p className="limit-hint" style={{ marginTop: "10px" }}>
                  You’ve reached the maximum of {MAX_MOLECULES} molecules per
                  run (library + custom).
                </p>
              )}

              {scoreError && (
                <div className="error-banner" style={{ marginTop: "12px" }}>
                  <span>⚠</span>
                  <span>{scoreError}</span>
                </div>
              )}
            </section>
          </main>

          {/* Bottom bar for selection view */}
          <section className="bottom-bar">
            <div className="bottom-bar-inner">
              <div>
                <div className="bottom-bar-title">
                  Ready for affinity analysis
                </div>
                <div className="bottom-bar-subtitle">
                  Pest:{" "}
                  <strong>
                    {selectedPest ? selectedPest.name : "None selected"}
                  </strong>{" "}
                  · Molecules:{" "}
                  <strong>{totalSelectedMolecules}</strong>
                </div>
              </div>
              <button
                className="primary-button"
                disabled={
                  !selectedPestId ||
                  totalSelectedMolecules === 0 ||
                  scoreLoading
                }
                onClick={handleRunScoring}
              >
                {scoreLoading ? "Running analysis…" : "Open Affinity Report"}
              </button>
            </div>
          </section>
        </>
      )}

      {view === VIEW_REPORT && scoreResponse && (
        <main className="app-main report-layout">
          <ReportView
            scoreResponse={scoreResponse}
            onBack={handleBackToSelection}
            pest={selectedPest}
            pesticideMetas={selectedPesticideMetas}
            customMolecules={customMolecules}
          />
        </main>
      )}
    </div>
  );
}

/* ---------- Selection view cards ---------- */

function PestCard({ pest, isSelected, onSelect }) {
  const imgSrc = PEST_IMAGE_MAP[pest.id];

  return (
    <button
      type="button"
      className={`card pest-card ${isSelected ? "card-selected" : ""}`}
      onClick={onSelect}
    >
      <div className="card-image-wrapper">
        {imgSrc ? (
          <img src={imgSrc} alt={pest.name} className="card-image" />
        ) : (
          <div className="card-image-placeholder">
            <span>{pest.name.charAt(0).toUpperCase()}</span>
          </div>
        )}
      </div>
      <div className="card-content">
        <h3 className="card-title">{pest.name}</h3>
        <p className="card-subtitle">Click to set as target organism.</p>
      </div>
    </button>
  );
}

function PesticideCard({
  pesticide,
  isSelected,
  disabled,
  onToggle,
  maxReached,
}) {
  const imgSrc = PESTICIDE_IMAGE_MAP[pesticide.id];

  const handleClick = () => {
    if (disabled || maxReached) return;
    onToggle();
  };

  return (
    <button
      type="button"
      className={`card pesticide-card ${
        isSelected ? "card-selected" : ""
      } ${disabled ? "card-disabled" : ""} ${
        maxReached ? "card-max-reached" : ""
      }`}
      onClick={handleClick}
    >
      <div className="card-image-wrapper">
        {imgSrc ? (
          <img src={imgSrc} alt={pesticide.name} className="card-image" />
        ) : (
          <div className="card-image-placeholder">
            <span>{pesticide.name.charAt(0).toUpperCase()}</span>
          </div>
        )}
      </div>
      <div className="card-content">
        <h3 className="card-title">{pesticide.name}</h3>
        <p className="card-subtitle">
          SMILES: <code>{pesticide.smiles}</code>
        </p>
        <div className="card-tag-row">
          <span className="card-tag">
            {isSelected ? "Selected" : maxReached ? "Limit reached" : "Available"}
          </span>
        </div>
      </div>
    </button>
  );
}

/* ---------- Report view (new page) ---------- */

function ReportView({
  scoreResponse,
  onBack,
  pest,
  pesticideMetas,
  customMolecules,
}) {
  const pestName = pest?.name || scoreResponse.pest;
  const pestImgSrc = pest ? PEST_IMAGE_MAP[pest.id] : null;

  // Build metadata array in the same order as payload:
  // first library pesticides, then custom molecules
  const allMeta = [
    ...pesticideMetas.map((meta) => ({ type: "library", meta })),
    ...customMolecules.map((cm) => ({ type: "custom", meta: cm })),
  ];

  const paired = scoreResponse.results.map((res, idx) => ({
    result: res,
    context: allMeta[idx] || { type: "unknown", meta: null },
  }));

  const sorted = [...paired].sort(
    (a, b) => Number(a.result.score) - Number(b.result.score)
  );

  const maxScore =
    sorted.length > 0
      ? Math.max(...sorted.map((p) => Number(p.result.score)))
      : 0;

  return (
    <div className="report-page">
      <div className="report-header-row">
        <button className="secondary-button" onClick={onBack}>
          ← Back to selection
        </button>
        <div className="report-header-text">
          <h2>Affinity & Safety Report</h2>
          <p>
            Target pest: <strong>{pestName}</strong> ·{" "}
            {sorted.length} molecule{sorted.length !== 1 ? "s" : ""} evaluated
          </p>
        </div>
      </div>

      {/* Hero pest profile */}
      <section className="report-section">
        <div className="report-hero">
          <div className="report-hero-image-wrapper">
            {pestImgSrc ? (
              <img
                src={pestImgSrc}
                alt={pestName}
                className="report-hero-image"
              />
            ) : (
              <div className="report-hero-image-placeholder">
                <span>{pestName.charAt(0).toUpperCase()}</span>
              </div>
            )}
          </div>

          <div className="report-hero-content">
            <h3 className="report-hero-title">{pestName}</h3>
            <p className="report-hero-subtitle">
              Protein affinity profile generated using PLAPT, based on the
              selected candidate molecules.
            </p>

            <div className="report-hero-meta-row">
              <div className="report-hero-meta">
                <span className="report-hero-meta-label">
                  Molecules evaluated
                </span>
                <span className="report-hero-meta-value">
                  {sorted.length}
                </span>
              </div>
              <div className="report-hero-meta">
                <span className="report-hero-meta-label">
                  Best affinity (µM)
                </span>
                <span className="report-hero-meta-value">
                  {sorted.length > 0
                    ? Number(sorted[0].result.score).toFixed(4)
                    : "—"}
                </span>
              </div>
              <div className="report-hero-meta">
                <span className="report-hero-meta-label">
                  Worst affinity (µM)
                </span>
                <span className="report-hero-meta-value">
                  {sorted.length > 0
                    ? Number(sorted[sorted.length - 1].result.score).toFixed(4)
                    : "—"}
                </span>
              </div>
            </div>
          </div>
        </div>
      </section>

      {/* Two-column layout: molecules + chart */}
      <section className="report-section">
        <div className="report-two-column">
          <div className="report-column">
            <h3 className="report-section-title">
              Molecule profiles (ranked best to worst)
            </h3>
            <p className="report-section-subtitle">
              Lower affinity values (µM) indicate stronger binding for this model.
            </p>

            <div className="report-list">
              {sorted.map((pair, idx) => (
                <ReportCard key={pair.result.name + idx} pair={pair} rank={idx + 1} />
              ))}
            </div>
          </div>

          {sorted.length > 1 && (
            <div className="report-column report-column-chart">
              <h3 className="report-section-title">
                Affinity comparison chart
              </h3>
              <p className="report-section-subtitle">
                Bar height represents predicted affinity in µM (0 to max across
                selected molecules).
              </p>

              <BarChartAffinity pairs={sorted} maxScore={maxScore} />
            </div>
          )}
        </div>
      </section>
    </div>
  );
}

function ReportCard({ pair, rank }) {
  const { result, context } = pair;
  const meta = context?.meta;
  const isCustom = context?.type === "custom";

  const displayName =
    result.name || meta?.name || (isCustom ? "Custom molecule" : "Unnamed");

  const displaySmiles = result.smiles || meta?.smiles || "—";

  const imgSrc =
    !isCustom && meta && meta.id ? PESTICIDE_IMAGE_MAP[meta.id] : null;

  const { score, label, interpretation } = result;

  const labelColorClass =
    label === "High"
      ? "chip-high"
      : label === "Medium"
      ? "chip-medium"
      : "chip-low";

  return (
    <div className="report-card">
      <div className="report-card-main">
        <div className="report-card-image-wrapper">
          {imgSrc ? (
            <img src={imgSrc} alt={displayName} className="report-card-image" />
          ) : (
            <div className="report-card-image-placeholder">
              <span>{displayName.charAt(0).toUpperCase()}</span>
            </div>
          )}
        </div>

        <div className="report-card-body">
          <div className="report-card-header">
            <div className="report-card-title-row">
              <span className="report-rank">#{rank}</span>
              <h4 className="report-card-title">{displayName}</h4>
              <span className={`chip ${labelColorClass}`}>{label} affinity</span>
            </div>
            <div className="report-card-smiles">
              <span>SMILES: </span>
              <code>{displaySmiles}</code>
            </div>
          </div>

          <div className="report-metrics">
            <div className="metric">
              <div className="metric-label">Predicted affinity (µM)</div>
              <div className="metric-value">
                {typeof score === "number" ? score.toFixed(4) : score}
              </div>
            </div>
          </div>

          <p className="report-interpretation">{interpretation}</p>
        </div>
      </div>
    </div>
  );
}

/* ---------- Bar chart ---------- */

function BarChartAffinity({ pairs, maxScore }) {
  const safeMax = maxScore || 1;
  const axisMaxLabel = safeMax.toFixed(2);

  return (
    <div className="affinity-chart">
      <div className="affinity-chart-inner">
        {/* Y axis numeric labels */}
        <div className="affinity-chart-y">
          <span>{axisMaxLabel}</span>
          <span>0</span>
        </div>

        {/* Right side: plot (bars + labels) */}
        <div className="affinity-chart-plot">
          {/* Bars with axes */}
          <div className="affinity-chart-bars">
            {pairs.map(({ result }) => {
              const scoreNum = Number(result.score);

              // Keep very small affinities visible with a minimum height
              const rawPct = (scoreNum / safeMax) * 100;
              const heightPct = Math.max((safeMax - scoreNum) / safeMax * 100, 3);

              const label = result.name || result.smiles;

              return (
                <div
                  className="affinity-chart-bar-wrapper"
                  key={result.name + result.smiles}
                >
                  <div className="affinity-chart-bar">
                    <div
                      className="affinity-chart-bar-fill"
                      style={{ height: `${heightPct}%` }}
                      title={`${label} – ${scoreNum.toFixed(4)} µM`}
                    />
                  </div>
                </div>
              );
            })}
          </div>

          {/* X axis category labels (compound names) under the axis */}
          <div className="affinity-chart-x-categories">
            {pairs.map(({ result }) => {
              const label = result.name || result.smiles;
              return (
                <div
                  className="affinity-chart-x-category"
                  key={`label-${result.name + result.smiles}`}
                >
                  <span title={label}>{label}</span>
                </div>
              );
            })}
          </div>
        </div>
      </div>

      {/* Axis titles */}
      <div className="affinity-chart-axis-labels">
        <span className="axis-y-label">Affinity (µM)</span>
        <span className="axis-x-label">Compound</span>
      </div>
    </div>
  );
}

export default App;
