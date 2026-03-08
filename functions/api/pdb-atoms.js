export async function onRequest({ request }) {
  const url = new URL(request.url);
  const pdbId = url.searchParams.get("pdbId");
  const ligandId = url.searchParams.get("ligandId");

  if (!pdbId || !ligandId) {
    return Response.json({ error: "pdbId and ligandId required" }, { status: 400 });
  }

  const upstream = `https://models.rcsb.org/v1/${encodeURIComponent(pdbId)}/atoms?label_comp_id=${encodeURIComponent(ligandId)}&encoding=json`;

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
          "Cache-Control": "public, max-age=86400",
        },
      });
    } catch {
      if (attempt === 0) await new Promise((r) => setTimeout(r, 500));
    }
  }

  return Response.json({ error: "PDB upstream failed" }, { status: 502 });
}
