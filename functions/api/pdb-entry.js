export async function onRequest({ request }) {
  const url = new URL(request.url);
  const pdbId = url.searchParams.get("pdbId");

  if (!pdbId) {
    return Response.json({ error: "pdbId required" }, { status: 400 });
  }

  const upstream = `https://data.rcsb.org/rest/v1/core/entry/${encodeURIComponent(pdbId)}`;

  try {
    const res = await fetch(upstream, {
      headers: { Accept: "application/json" },
    });
    const data = await res.text();
    return new Response(data, {
      status: res.status,
      headers: {
        "Content-Type": "application/json",
        "Cache-Control": "public, max-age=86400",
      },
    });
  } catch {
    return Response.json({ error: "PDB entry upstream failed" }, { status: 502 });
  }
}
