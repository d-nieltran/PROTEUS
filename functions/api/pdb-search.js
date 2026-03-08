export async function onRequest({ request }) {
  if (request.method !== "POST") {
    return Response.json({ error: "POST required" }, { status: 405 });
  }

  const body = await request.text();
  const upstream = "https://search.rcsb.org/rcsbsearch/v2/query";

  try {
    const res = await fetch(upstream, {
      method: "POST",
      headers: {
        "Content-Type": "application/json",
        Accept: "application/json",
      },
      body,
    });
    const data = await res.text();
    return new Response(data, {
      status: res.status,
      headers: {
        "Content-Type": "application/json",
        "Cache-Control": "public, max-age=3600",
      },
    });
  } catch {
    return Response.json({ error: "PDB search upstream failed" }, { status: 502 });
  }
}
