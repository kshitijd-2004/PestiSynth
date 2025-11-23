// src/App.jsx
import React, { useEffect, useState } from "react";

const API_BASE_URL = "http://localhost:8000"; // <-- change if needed

// Optional: map pest/pesticide IDs to image URLs
const PEST_IMAGE_MAP = {
  // fall_armyworm: "/images/pests/fall_armyworm.jpg",
  // mosquito_larvae: "/images/pests/mosquito_larvae.jpg",
};

const PESTICIDE_IMAGE_MAP = {
  // "2_4_D": "/images/pesticides/2_4_D.png",
  // "glyphosate": "/images/pesticides/glyphosate.png",
};

function App() {
  const [pests, setPests] = useState([]);
  const [pesticides, setPesticides] = useState([]);

  const [pestsLoading, setPestsLoading] = useState(true);
  const [pesticidesLoading, setPesticidesLoading] = useState(true);
  const [error, setError] = useState(null);

  const [selectedPestId, setSelectedPestId] = useState(null);
  const [selectedPesticideIds, setSelectedPesticideIds] = useState([]);

  const MAX_PESTICIDES = 5;

  // Score /report state
  const [scoreLoading, setScoreLoading] = useState(false);
  const [scoreError, setScoreError] = useState(null);
  const [scoreResponse, setScoreResponse] = useState(null);

  useEffect(() => {
    // Fetch pests
    const fetchPests = async () => {
      try {
        setPestsLoading(true);
        const res = await fetch(`${API_BASE_URL}/pests`);
        if (!res.ok) {
          throw new Error(`Failed to fetch pests: ${res.status}`);
        }
        const data = await res.json();
        setPests(data);
      } catch (err) {
        console.error(err);
        setError(err.message || "Failed to load pests.");
      } finally {
        setPestsLoading(false);
      }
    };

    // Fetch pesticides
    const fetchPesticides = async () => {
      try {
        setPesticidesLoading(true);
        const res = await fetch(`${API_BASE_URL}/pesticides`);
        if (!res.ok) {
          throw new Error(`Failed to fetch pesticides: ${res.status}`);
        }
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
    // Clear previous results when pest changes
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
      if (selectedPesticideIds.length >= MAX_PESTICIDES) return;
      setSelectedPesticideIds((prev) => [...prev, pesticideId]);
    }
  };

  const handleRunScoring = async () => {
    if (!selectedPestId || selectedPesticideIds.length === 0) return;

    setScoreLoading(true);
    setScoreError(null);
    setScoreResponse(null);

    try {
      const payload = {
        pest_id: selectedPestId,
        molecules: selectedPesticideIds.map((id) => ({
          pesticide_id: id,
        })),
      };

      const res = await fetch(`${API_BASE_URL}/score/batch`, {
        method: "POST",
        headers: {
          "Content-Type": "application/json",
        },
        body: JSON.stringify(payload),
      });

      if (!res.ok) {
        let message = `Scoring request failed (${res.status})`;
        try {
          const errData = await res.json();
          if (errData?.detail) message = errData.detail;
        } catch (_) {
          // ignore JSON parse error
        }
        throw new Error(message);
      }

      const data = await res.json();
      setScoreResponse(data);
    } catch (err) {
      console.error(err);
      setScoreError(err.message || "Failed to run scoring.");
    } finally {
      setScoreLoading(false);
    }
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
                Start by choosing a target pest, then pick up to {MAX_PESTICIDES}{" "}
                molecules from your pesticide library.
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

        {/* Pest Selection Section */}
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

        {/* Pesticide Selection Section */}
        <section className="section">
          <div className="section-header">
            <h2>2. Choose Molecules from Library</h2>
            <p>
              Pick up to {MAX_PESTICIDES} pesticides from your predefined
              library. You can also add custom SMILES later in the workflow.
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
        </section>

        {/* Affinity & Safety Report Section */}
        <section className="section report-section">
          <div className="section-header">
            <h2>3. Affinity & Safety Report</h2>
            <p>
              Run the PLAPT model to generate affinity scores, qualitative labels,
              and a SafetyScan summary for each selected molecule.
            </p>
          </div>

          {scoreLoading && (
            <div className="report-status">Running analysis…</div>
          )}

          {scoreError && (
            <div className="error-banner" style={{ marginTop: "12px" }}>
              <span>⚠</span>
              <span>{scoreError}</span>
            </div>
          )}

          {scoreResponse && (
            <div className="report-results">
              <div className="report-header">
                <div>
                  <div className="report-title">
                    Target: {scoreResponse.pest}
                  </div>
                  <div className="report-subtitle">
                    {scoreResponse.results.length} molecule
                    {scoreResponse.results.length !== 1 ? "s" : ""} evaluated.
                  </div>
                </div>
              </div>

              <div className="report-grid">
                {scoreResponse.results.map((res) => (
                  <ReportCard key={res.name + res.smiles} result={res} />
                ))}
              </div>
            </div>
          )}
        </section>
      </main>

      {/* Bottom bar action */}
      <section className="bottom-bar">
        <div className="bottom-bar-inner">
          <div>
            <div className="bottom-bar-title">Ready for affinity analysis</div>
            <div className="bottom-bar-subtitle">
              Pest:{" "}
              <strong>{selectedPest ? selectedPest.name : "None selected"}</strong>{" "}
              · Molecules: <strong>{selectedPesticideIds.length}</strong>
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
            {scoreLoading ? "Running analysis…" : "Run Affinity + Safety Scan"}
          </button>
        </div>
      </section>
    </div>
  );
}

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

function ReportCard({ result }) {
  const { name, smiles, score, label, interpretation, safety_flag, safety_matches } =
    result;

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
          <h3 className="report-card-title">{name}</h3>
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
          <div className="metric">
            <div className="metric-label">Safety flag</div>
            <div className="metric-value metric-value-secondary">
              {safety_flag}
            </div>
          </div>
        </div>

        <p className="report-interpretation">{interpretation}</p>

        {safety_matches && safety_matches.length > 0 && (
          <div className="report-safety-matches">
            <div className="safety-matches-title">
              Similar to hazardous/restricted:
            </div>
            <div className="safety-matches-list">
              {safety_matches.map((m) => (
                <span key={m.name} className="chip chip-safety">
                  {m.name} · sim {m.similarity.toFixed(2)}
                </span>
              ))}
            </div>
          </div>
        )}
      </div>
    </div>
  );
}

export default App;
