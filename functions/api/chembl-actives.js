// ChEMBL target IDs for PROTEUS NTD targets
const TARGET_MAP = {
  PfDHFR: "CHEMBL1939",   // P. falciparum bifunctional DHFR-TS
  LmPTR1: "CHEMBL6194",   // L. major pteridine reductase 1
  TcTR: "CHEMBL5131",     // T. cruzi trypanothione reductase
};

export async function onRequest({ request }) {
  const url = new URL(request.url);
  const targetId = url.searchParams.get("target");

  if (!targetId || !TARGET_MAP[targetId]) {
    return Response.json(
      { error: "target required (PfDHFR, LmPTR1, or TcTR)", targets: Object.keys(TARGET_MAP) },
      { status: 400 }
    );
  }

  const chemblTarget = TARGET_MAP[targetId];
  // Fetch known actives with pChEMBL >= 5 (IC50 <= 10 µM)
  const upstream =
    `https://www.ebi.ac.uk/chembl/api/data/activity.json` +
    `?target_chembl_id=${chemblTarget}` +
    `&pchembl_value__gte=5` +
    `&limit=1000`;

  for (let attempt = 0; attempt < 2; attempt++) {
    try {
      const res = await fetch(upstream, {
        headers: { Accept: "application/json" },
      });
      if (!res.ok) continue;
      const data = await res.json();
      const activities = data?.activities || [];

      // Deduplicate by molecule ChEMBL ID, keep best pChEMBL value
      const byMol = new Map();
      for (const a of activities) {
        const id = a.molecule_chembl_id;
        if (!id) continue;
        const pchembl = parseFloat(a.pchembl_value) || 0;
        const existing = byMol.get(id);
        if (!existing || pchembl > existing.pchembl) {
          byMol.set(id, {
            chemblId: id,
            smiles: a.canonical_smiles || null,
            pchembl,
            assayType: a.assay_type || "",
            activityType: a.standard_type || "",
          });
        }
      }

      const actives = Array.from(byMol.values())
        .sort((a, b) => b.pchembl - a.pchembl);

      return Response.json(
        {
          target: targetId,
          chemblTarget,
          count: actives.length,
          actives,
        },
        {
          headers: { "Cache-Control": "public, max-age=86400" },
        }
      );
    } catch {
      if (attempt === 0) await new Promise((r) => setTimeout(r, 500));
    }
  }

  return Response.json({ error: "ChEMBL upstream failed" }, { status: 502 });
}
