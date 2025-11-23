// src/App.jsx
import React, { useEffect, useState } from "react";

const API_BASE_URL = "http://localhost:8000"; // <-- change if needed

// Optional: map pest/pesticide IDs to image URLs
// Replace these with your real images when you have them.
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
  };

  const handlePesticideToggle = (pesticideId) => {
    if (selectedPesticideIds.includes(pesticideId)) {
      setSelectedPesticideIds((prev) =>
        prev.filter((id) => id !== pesticideId)
      );
    } else {
      if (selectedPesticideIds.length >= MAX_PESTICIDES) return;
      setSelectedPesticideIds((prev) => [...prev, pesticideId]);
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

        {/* Placeholder footer action (future: call /score/batch) */}
        <section className="bottom-bar">
          <div className="bottom-bar-inner">
            <div>
              <div className="bottom-bar-title">Ready for affinity analysis</div>
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
              disabled={!selectedPestId || selectedPesticideIds.length === 0}
              onClick={() => {
                // This is where you’ll later hook in /score/batch
                // For now, we just log to verify the selections.
                console.log("Selected pest:", selectedPestId);
                console.log("Selected pesticide IDs:", selectedPesticideIds);
                alert(
                  "Selections captured. Next step: hook this button up to /score/batch."
                );
              }}
            >
              Continue to Affinity Report
            </button>
          </div>
        </section>
      </main>
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

export default App;
