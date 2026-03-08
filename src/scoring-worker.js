// Web Worker for parallel compound processing (feature extraction + CPU scoring)

const FTYPE = { hba: 0, hbd: 1, hydrophobic: 2, positive: 3, negative: 4 };

function extractFeaturesFromSMILES(smiles, props) {
  const features = [];
  const hash = (s, i) => {
    let h = 0;
    for (let c = 0; c < s.length; c++) h = ((h << 5) - h + s.charCodeAt(c) + i * 31) | 0;
    return (h & 0x7fffffff) / 0x7fffffff;
  };
  const patterns = [
    { re: /\[NH2?\]|N\(/g, type: "hbd", weight: 1.3 },
    { re: /O\(|OH|O\)/g, type: "hba", weight: 1.4 },
    { re: /C\(=O\)O/g, type: "negative", weight: 1.5 },
    { re: /C\(=O\)N/g, type: "hba", weight: 1.2 },
    { re: /c1.*c.*c.*c/g, type: "hydrophobic", weight: 1.0 },
    { re: /N\+|NH3/g, type: "positive", weight: 1.1 },
    { re: /S\(|s1/g, type: "hydrophobic", weight: 0.7 },
    { re: /F|Cl|Br/g, type: "hydrophobic", weight: 0.6 },
  ];
  let fIdx = 0;
  for (const p of patterns) {
    const matches = smiles.match(p.re);
    if (matches) {
      for (let m = 0; m < Math.min(matches.length, 3); m++) {
        const radius = Math.sqrt(parseFloat(props.full_mwt) || 300) * 0.3;
        features.push({
          type: p.type,
          x: (hash(smiles, fIdx * 3) - 0.5) * radius,
          y: (hash(smiles, fIdx * 3 + 1) - 0.5) * radius,
          z: (hash(smiles, fIdx * 3 + 2) - 0.5) * radius,
          weight: p.weight,
        });
        fIdx++;
      }
    }
  }
  if (!features.some((f) => f.type === "hba") && (parseInt(props.hba) || 0) > 0)
    features.push({ type: "hba", x: (hash(smiles, 100) - 0.5) * 5, y: (hash(smiles, 101) - 0.5) * 5, z: (hash(smiles, 102) - 0.5) * 5, weight: 1.0 });
  if (!features.some((f) => f.type === "hbd") && (parseInt(props.hbd) || 0) > 0)
    features.push({ type: "hbd", x: (hash(smiles, 200) - 0.5) * 5, y: (hash(smiles, 201) - 0.5) * 5, z: (hash(smiles, 202) - 0.5) * 5, weight: 1.0 });
  return features;
}

function scoreCompound(compound, pharmacophore) {
  if (!pharmacophore?.features) return 0;
  let score = 0;
  for (const pf of pharmacophore.features) {
    let bestMatch = 0;
    for (const cf of compound.features) {
      let ts = -0.1;
      if (pf.type === cf.type) ts = 1.0;
      else if ((pf.type === "hba" && cf.type === "hbd") || (pf.type === "hbd" && cf.type === "hba")) ts = 0.8;
      else if ((pf.type === "positive" && cf.type === "negative") || (pf.type === "negative" && cf.type === "positive")) ts = 1.2;
      const dx = pf.x - cf.x, dy = pf.y - cf.y, dz = pf.z - cf.z;
      bestMatch = Math.max(bestMatch, ts * Math.exp(-(dx * dx + dy * dy + dz * dz) / 4.0) * pf.weight);
    }
    score += bestMatch;
  }
  if (!compound.lipinski) score *= 0.7;
  score *= Math.max(0.3, 1.0 - Math.abs(compound.mw - 350) / 400);
  score -= compound.tpsa * 0.003;
  score -= compound.rotBonds * 0.08;
  return +score.toFixed(4);
}

function processMolecule(m) {
  const p = m.molecule_properties;
  if (!p) return null;
  const smiles = m.molecule_structures?.canonical_smiles;
  if (!smiles) return null;
  const features = extractFeaturesFromSMILES(smiles, p);
  if (features.length < 2) return null;
  return {
    id: m.molecule_chembl_id,
    name: m.pref_name || m.molecule_chembl_id,
    smiles,
    mw: +parseFloat(p.full_mwt).toFixed(1),
    hbd: parseInt(p.hbd) || 0,
    hba: parseInt(p.hba) || 0,
    logP: +parseFloat(p.alogp).toFixed(2),
    tpsa: +parseFloat(p.psa).toFixed(1),
    rotBonds: parseInt(p.rtb) || 0,
    aromRings: parseInt(p.aromatic_rings) || 0,
    heavyAtoms: parseInt(p.heavy_atoms) || 0,
    features,
    lipinski: (parseFloat(p.full_mwt) || 0) <= 500 && (parseInt(p.hbd) || 0) <= 5 && (parseInt(p.hba) || 0) <= 10 && (parseFloat(p.alogp) || 0) <= 5,
    source: "ChEMBL",
  };
}

function generateSyntheticBatch(startIdx, count) {
  const mols = [];
  for (let i = 0; i < count; i++) {
    const idx = startIdx + i;
    const r = (n) => { let x = Math.sin(idx * 9301 + n * 49297) * 49297; return x - Math.floor(x); };
    const mw = 200 + r(1) * 330, hbd = Math.floor(r(2) * 5), hba = Math.floor(r(3) * 8) + 1;
    const logP = -0.5 + r(4) * 5, tpsa = 30 + r(5) * 100, rotBonds = Math.floor(r(6) * 8), aromRings = Math.floor(r(7) * 4);
    const nF = 3 + Math.floor(r(10) * 5), features = [];
    for (let f = 0; f < nF; f++) {
      const types = ["hba", "hbd", "hydrophobic", "positive", "negative"];
      features.push({ type: types[Math.floor(r(20 + f) * types.length)], x: (r(30 + f) - 0.5) * 8, y: (r(40 + f) - 0.5) * 8, z: (r(50 + f) - 0.5) * 8, weight: 0.8 + r(60 + f) * 0.8 });
    }
    mols.push({
      id: `SYN${String(idx).padStart(8, "0")}`, name: `Synthetic-${idx}`, smiles: null,
      mw: +mw.toFixed(1), hbd, hba, logP: +logP.toFixed(2), tpsa: +tpsa.toFixed(1), rotBonds, aromRings,
      heavyAtoms: Math.floor(mw / 14), features, lipinski: mw <= 500 && hbd <= 5 && hba <= 10 && logP <= 5, source: "synthetic",
    });
  }
  return mols;
}

self.onmessage = (e) => {
  const { type, id } = e.data;

  if (type === "process_molecules") {
    // Extract features from raw ChEMBL molecules
    const { molecules } = e.data;
    const compounds = molecules.map(processMolecule).filter(Boolean);
    self.postMessage({ type: "processed", id, compounds });
  }

  if (type === "score_cpu") {
    // CPU-parallel scoring fallback
    const { compounds, pharmacophore } = e.data;
    const scores = compounds.map((c) => scoreCompound(c, pharmacophore));
    self.postMessage({ type: "scored", id, scores });
  }

  if (type === "generate_synthetic") {
    // Generate synthetic compounds
    const { startIdx, count } = e.data;
    const compounds = generateSyntheticBatch(startIdx, count);
    self.postMessage({ type: "processed", id, compounds });
  }
};
