export async function onRequest({ request }) {
  const url = new URL(request.url);
  const offset = url.searchParams.get("offset") || "0";
  const limit = url.searchParams.get("limit") || "100";

  const upstream =
    `https://www.ebi.ac.uk/chembl/api/data/molecule.json` +
    `?limit=${limit}&offset=${offset}` +
    `&molecule_properties__full_mwt__lte=550` +
    `&molecule_properties__full_mwt__gte=200` +
    `&molecule_properties__alogp__gte=-1` +
    `&molecule_properties__alogp__lte=6` +
    `&molecule_properties__hbd__lte=5` +
    `&molecule_properties__hba__lte=10` +
    `&molecule_structures__canonical_smiles__isnull=false`;

  for (let attempt = 0; attempt < 2; attempt++) {
    try {
      const res = await fetch(upstream, {
        headers: { Accept: "application/json" },
      });
      if (!res.ok) continue;
      const data = await res.text();
      return new Response(data, {
        headers: {
          "Content-Type": "application/json",
          "Cache-Control": "public, max-age=3600",
        },
      });
    } catch {
      if (attempt === 0) await new Promise((r) => setTimeout(r, 500));
    }
  }

  return Response.json({ error: "ChEMBL upstream failed" }, { status: 502 });
}
