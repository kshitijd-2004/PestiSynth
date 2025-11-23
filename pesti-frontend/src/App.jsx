// src/App.jsx
import React, { useEffect, useState } from "react";

const API_BASE_URL = "http://localhost:8000"; // change if needed

// Optional: map pest/pesticide IDs to image URLs
const PEST_IMAGE_MAP = {
  // fall_armyworm: "/images/pests/fall_armyworm.jpg",
  // mosquito_larvae: "/images/pests/mosquito_larvae.jpg",
};

const PESTICIDE_IMAGE_MAP = {
  // "2_4_D": "/images/pesticides/2_4_D.png",
  // "glyphosate": "/images/pesticides/glyphosate.png",
};

const VIEW_SELECTION = "selection";
const VIEW_REPORT = "report";

function App() {
  const [pests, setPests] = useState([]);
  const [pesticides, setPesticides] = useState([]);

  const [pestsLoading, setPestsLoading] = useState(true);
  const [pesticidesLoading, setPesticidesLoading] = useState(true);
  const [error, setError] = useState(null);

  const [selectedPestId, setSelectedPestId] = useState(null);
  const [selectedPesticideIds, setSelectedPesticideIds] = useState([]);

  const MAX_PESTICIDES = 5;

  const [scoreLoading, setScoreLoading] = useState(false);
  const [scoreError, setScoreError] = useState(null);
  const [scoreResponse, setScoreResponse] = useState(null);

  const [view, setView] = useState(VIEW_SELECTION);

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
        if (!res.ok) throw new Error(`Failed to fetch pesticides: ${res.status}`);
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

  const handlePestSelect = (pestId) => {
    setSelectedPestId(pestId);
    setScoreResponse(null);
    setScoreError(null);
  };

  const handlePesticideToggle = (pesticideId) => {
    setScoreResponse(null);
    setScoreError(null);

    if (selectedPesticideIds.includes(pesticideId)) {
      setSelectedPesticideIds((prev) => prev.filter((id) => id !== pesticideId));
    } else {
      if (selectedPesticideIds.length >= MAX_PESTICIDES) return;
      setSelectedPesticideIds((prev) => [...prev, pesticideId]);
    }
  };

  const handleRunScoring = async () => {
    if (!selectedPestId || selectedPesticideIds.length === 0) return;

    setScoreLoading(true);
    setScoreError(null);

    try {
      const payload = {
        pest_id: selectedPestId,
        molecules: selectedPesticideIds.map((id) => ({
          pesticide_id: id,
        })),
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

  const selectedPest = pests.find((p) => p.id === selectedPestId);

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
                    Choose a target pest and up to {MAX_PESTICIDES} library
                    molecules to evaluate.
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
                      {selectedPesticideIds.length} / {MAX_PESTICIDES}
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
                  Select the pest whose protein you want to evaluate. Only one pest
                  can be active at a time.
                </p>
              </div>

              {pestsLoading ? (
                <div className="loading">Loading pests…</div>
              ) : (
                <div className="card-grid">
                  {pests.map((pest) => (
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
                  Pick up to {MAX_PESTICIDES} pesticides from your predefined
                  library.
                </p>
              </div>

              {!selectedPestId && (
                <div className="section-hint">
                  Please select a pest first to continue with molecule selection.
                </div>
              )}

              {pesticidesLoading ? (
                <div className="loading">Loading pesticides…</div>
              ) : (
                <div
                  className={`card-grid ${
                    !selectedPestId ? "card-grid-disabled" : ""
                  }`}
                >
                  {pesticides.map((p) => (
                    <PesticideCard
                      key={p.id}
                      pesticide={p}
                      isSelected={selectedPesticideIds.includes(p.id)}
                      disabled={!selectedPestId}
                      onToggle={() => handlePesticideToggle(p.id)}
                      maxReached={
                        selectedPesticideIds.length >= MAX_PESTICIDES &&
                        !selectedPesticideIds.includes(p.id)
                      }
                    />
                  ))}
                </div>
              )}

              {selectedPesticideIds.length >= MAX_PESTICIDES && (
                <p className="limit-hint">
                  You’ve reached the maximum of {MAX_PESTICIDES} molecules per run.
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
                  <strong>{selectedPesticideIds.length}</strong>
                </div>
              </div>
              <button
                className="primary-button"
                disabled={
                  !selectedPestId ||
                  selectedPesticideIds.length === 0 ||
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
          />
        </main>
      )}
    </div>
  );
}

/* ---------- Cards ---------- */

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

/* ---------- Report View (new page) ---------- */

function ReportView({ scoreResponse, onBack }) {
  const { pest, results } = scoreResponse;

  // Sort: best (lowest affinity µM) to worst
  const sorted = [...results].sort((a, b) => a.score - b.score);
  const maxScore =
    sorted.length > 0 ? Math.max(...sorted.map((r) => Number(r.score))) : 0;

  return (
    <div className="report-page">
      <div className="report-header-row">
        <button className="secondary-button" onClick={onBack}>
          ← Back to selection
        </button>
        <div className="report-header-text">
          <h2>Affinity & Safety Report</h2>
          <p>
            Target pest: <strong>{pest}</strong> ·{" "}
            {sorted.length} molecule{sorted.length !== 1 ? "s" : ""} evaluated
          </p>
        </div>
      </div>

      {/* 1) Pest profile */}
      <section className="section report-section">
        <h3 className="report-section-title">1. Target pest profile</h3>
        <div className="report-pest-card">
          <div className="report-pest-icon">
            <span>{pest.charAt(0).toUpperCase()}</span>
          </div>
          <div>
            <div className="report-pest-name">{pest}</div>
            <div className="report-pest-subtitle">
              Protein affinity predictions generated using PLAPT.
            </div>
          </div>
        </div>
      </section>

      {/* 2) SMILES profiles sorted best→worst */}
      <section className="section report-section">
        <h3 className="report-section-title">
          2. Molecules ranked by predicted affinity
        </h3>
        <p className="report-section-subtitle">
          Lower affinity values (µM) indicate stronger binding for this model.
        </p>

        <div className="report-list">
          {sorted.map((res, idx) => (
            <ReportCard key={res.name + res.smiles + idx} result={res} rank={idx + 1} />
          ))}
        </div>
      </section>

      {/* 3) Bar graph: affinities vs pest */}
      {sorted.length > 1 && (
        <section className="section report-section">
          <h3 className="report-section-title">
            3. Affinity comparison chart
          </h3>
          <p className="report-section-subtitle">
            Bar height represents predicted affinity in µM (0 to max across
            selected molecules).
          </p>

          <BarChartAffinity results={sorted} maxScore={maxScore} />
        </section>
      )}
    </div>
  );
}

function ReportCard({ result, rank }) {
  const { name, smiles, score, label, interpretation } = result;

  const labelColorClass =
    label === "High"
      ? "chip-high"
      : label === "Medium"
      ? "chip-medium"
      : "chip-low";

  return (
    <div className="report-card">
      <div className="report-card-header">
        <div className="report-card-title-row">
          <span className="report-rank">#{rank}</span>
          <h4 className="report-card-title">{name}</h4>
          <span className={`chip ${labelColorClass}`}>{label} affinity</span>
        </div>
        <div className="report-card-smiles">
          <span>SMILES: </span>
          <code>{smiles}</code>
        </div>
      </div>

      <div className="report-card-body">
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
  );
}

/* ---------- Simple bar chart (CSS-only) ---------- */

function BarChartAffinity({ results, maxScore }) {
  const safeMax = maxScore || 1;

  return (
    <div className="affinity-chart">
      <div className="affinity-chart-inner">
        <div className="affinity-chart-y">
          <span>{safeMax.toFixed(2)}</span>
          <span>0</span>
        </div>

        <div className="affinity-chart-bars">
          {results.map((res) => {
            const heightPct = (Number(res.score) / safeMax) * 100;
            const shortSmiles =
              res.smiles.length > 16
                ? res.smiles.slice(0, 15) + "…"
                : res.smiles;

            return (
              <div className="affinity-chart-bar-wrapper" key={res.name + res.smiles}>
                <div className="affinity-chart-bar">
                  <div
                    className="affinity-chart-bar-fill"
                    style={{ height: `${heightPct}%` }}
                    title={`${res.name} – ${res.score.toFixed(4)} µM`}
                  />
                </div>
                <div className="affinity-chart-x-label">
                  <span title={res.smiles}>{shortSmiles}</span>
                </div>
              </div>
            );
          })}
        </div>
      </div>

      <div className="affinity-chart-axis-labels">
        <span className="axis-y-label">Affinity (µM)</span>
        <span className="axis-x-label">SMILES</span>
      </div>
    </div>
  );
}

export default App;
