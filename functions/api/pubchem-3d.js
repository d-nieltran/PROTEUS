export async function onRequest({ request }) {
  if (request.method !== "POST") {
    return Response.json({ error: "POST required" }, { status: 405 });
  }

  const { smiles } = await request.json();
  if (!Array.isArray(smiles) || smiles.length === 0) {
    return Response.json({ error: "smiles array required" }, { status: 400 });
  }

  // Limit batch size to prevent abuse
  const batch = smiles.slice(0, 20);
  const results = [];

  for (const smi of batch) {
    try {
      // Step 1: SMILES → CID
      const cidUrl = `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/${encodeURIComponent(smi)}/cids/JSON`;
      const cidRes = await fetch(cidUrl);
      if (!cidRes.ok) { results.push({ smiles: smi, atoms: null }); continue; }
      const cidData = await cidRes.json();
      const cid = cidData?.IdentifierList?.CID?.[0];
      if (!cid) { results.push({ smiles: smi, atoms: null }); continue; }

      // Step 2: CID → 3D conformer (JSON format)
      const confUrl = `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/${cid}/record/JSON?record_type=3d`;
      const confRes = await fetch(confUrl);
      if (!confRes.ok) { results.push({ smiles: smi, atoms: null }); continue; }
      const confData = await confRes.json();

      // Extract 3D coordinates from PubChem compound record
      const atoms = extractAtomsFromRecord(confData);
      results.push({ smiles: smi, cid, atoms });
    } catch {
      results.push({ smiles: smi, atoms: null });
    }
  }

  return Response.json(results, {
    headers: { "Cache-Control": "public, max-age=86400" },
  });
}

function extractAtomsFromRecord(data) {
  try {
    const compound = data?.PC_Compounds?.[0];
    if (!compound) return null;

    const atoms = compound.atoms;
    const coords = compound.coords?.[0];
    if (!atoms || !coords) return null;

    const elements = atoms.element || [];
    const conformer = coords.conformers?.[0];
    if (!conformer) return null;

    const xs = conformer.x || [];
    const ys = conformer.y || [];
    const zs = conformer.z || [];

    // PubChem uses atomic numbers for elements
    const ELEMENT_MAP = {
      1: "H", 6: "C", 7: "N", 8: "O", 9: "F", 15: "P", 16: "S",
      17: "Cl", 35: "Br", 53: "I", 5: "B", 14: "Si", 34: "Se",
    };

    const result = [];
    for (let i = 0; i < elements.length; i++) {
      if (i >= xs.length || i >= ys.length || i >= zs.length) break;
      const el = ELEMENT_MAP[elements[i]] || "C";
      // Skip hydrogens for pharmacophore features (too many, not relevant)
      if (el === "H") continue;
      result.push({ x: xs[i], y: ys[i], z: zs[i], element: el });
    }
    return result.length > 0 ? result : null;
  } catch {
    return null;
  }
}
