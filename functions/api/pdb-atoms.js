export async function onRequest({ request }) {
  const url = new URL(request.url);
  const pdbId = url.searchParams.get("pdbId");
  const ligandId = url.searchParams.get("ligandId");

  if (!pdbId || !ligandId) {
    return Response.json({ error: "pdbId and ligandId required" }, { status: 400 });
  }

  const upstream = `https://models.rcsb.org/v1/${encodeURIComponent(pdbId)}/atoms?label_comp_id=${encodeURIComponent(ligandId)}&encoding=cif`;

  for (let attempt = 0; attempt < 2; attempt++) {
    try {
      const res = await fetch(upstream);
      if (!res.ok) continue;
      const cif = await res.text();
      const atoms = parseCifAtoms(cif);
      return Response.json(atoms, {
        headers: { "Cache-Control": "public, max-age=86400" },
      });
    } catch {
      if (attempt === 0) await new Promise((r) => setTimeout(r, 500));
    }
  }

  return Response.json({ error: "PDB upstream failed" }, { status: 502 });
}

function parseCifAtoms(cif) {
  const atoms = [];
  const lines = cif.split("\n");
  let inLoop = false;
  let columns = [];
  let colMap = {};

  for (let i = 0; i < lines.length; i++) {
    const line = lines[i].trim();

    if (line === "loop_") {
      inLoop = true;
      columns = [];
      colMap = {};
      continue;
    }

    if (inLoop && line.startsWith("_atom_site.")) {
      const colName = line.replace("_atom_site.", "");
      colMap[colName] = columns.length;
      columns.push(colName);
      continue;
    }

    // If we were collecting _atom_site columns and hit a non-column line, parse data
    if (inLoop && columns.length > 0 && !line.startsWith("_")) {
      if (line === "" || line.startsWith("#") || line.startsWith("loop_") || line.startsWith("data_")) {
        // End of this loop's data
        if (atoms.length > 0) break;
        inLoop = false;
        columns = [];
        continue;
      }

      // Parse data row — split by whitespace
      const vals = splitCifRow(line);
      if (vals.length < columns.length) continue;

      const xIdx = colMap["Cartn_x"];
      const yIdx = colMap["Cartn_y"];
      const zIdx = colMap["Cartn_z"];
      const elIdx = colMap["type_symbol"];

      if (xIdx === undefined || yIdx === undefined || zIdx === undefined) continue;

      const x = parseFloat(vals[xIdx]);
      const y = parseFloat(vals[yIdx]);
      const z = parseFloat(vals[zIdx]);
      if (isNaN(x) || isNaN(y) || isNaN(z)) continue;

      atoms.push({
        x,
        y,
        z,
        type_symbol: elIdx !== undefined ? vals[elIdx] : "C",
      });
    }
  }

  return atoms;
}

function splitCifRow(line) {
  // Handle quoted strings in CIF format
  const parts = [];
  let i = 0;
  while (i < line.length) {
    if (line[i] === "'" || line[i] === '"') {
      const quote = line[i];
      let j = i + 1;
      while (j < line.length && line[j] !== quote) j++;
      parts.push(line.slice(i + 1, j));
      i = j + 1;
    } else if (line[i] === " " || line[i] === "\t") {
      i++;
    } else {
      let j = i;
      while (j < line.length && line[j] !== " " && line[j] !== "\t") j++;
      parts.push(line.slice(i, j));
      i = j;
    }
  }
  return parts;
}
