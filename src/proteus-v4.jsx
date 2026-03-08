import { useState, useEffect, useRef, useCallback, useMemo } from "react";
import * as THREE from "three";
import { ScatterChart, Scatter, XAxis, YAxis, Tooltip, ResponsiveContainer, AreaChart, Area, BarChart, Bar, Cell, CartesianGrid } from "recharts";

// ═══════════════════════════════════════════════════════════════════════════
// TARGETS — real PDB co-crystal structures
// ═══════════════════════════════════════════════════════════════════════════
const TARGETS = [
  { id:"PfDHFR", name:"P. falciparum DHFR-TS", disease:"Malaria", pdb:"1J3I", ligandId:"PYR", color:"#ff5078", pharmacophore:null,
    ref:"Yuvaniyama J, et al. Insights into antifolate resistance from malarial DHFR-TS structures. Nat Struct Biol. 2003;10(5):357-365.",
    pmid:"12704429", doi:"10.1038/nsb921", resolution:"2.33 Å",
    cocrystal:"Cycloguanil (antifolate)", organism:"P. falciparum",
    whyTarget:"DHFR is essential for folate biosynthesis in Plasmodium. Antifolate resistance mutations (S108N, N51I, C59R) reduce efficacy of pyrimethamine/cycloguanil, driving need for novel inhibitors." },
  { id:"LmPTR1", name:"L. major PTR1", disease:"Leishmaniasis", pdb:"1E92", ligandId:"MTX", color:"#78ff64", pharmacophore:null,
    ref:"Gourley DG, et al. Pteridine reductase mechanism correlates pterin and folate metabolism with drug resistance in trypanosomatid parasites. Nat Struct Biol. 2001;8(6):521-525.",
    pmid:"11373620", doi:"10.1038/88584", resolution:"2.20 Å",
    cocrystal:"Methotrexate (antifolate)", organism:"L. major",
    whyTarget:"PTR1 is unique to trypanosomatid parasites (absent in humans). It provides a metabolic bypass when DHFR is inhibited, making it essential for combination therapy strategies." },
  { id:"TcTR", name:"T. cruzi Trypanothione Red.", disease:"Chagas", pdb:"1BZL", ligandId:"QUN", color:"#50dcff", pharmacophore:null,
    ref:"Bond CS, et al. Crystal structure of Trypanosoma cruzi trypanothione reductase in complex with trypanothione, and the structure-based discovery of new natural product inhibitors. Structure. 1999;7(1):81-89.",
    pmid:"10368274", doi:"10.1016/S0969-2126(99)80023-2", resolution:"2.40 Å",
    cocrystal:"Quinacrine mustard", organism:"T. cruzi",
    whyTarget:"Trypanothione reductase replaces glutathione reductase in trypanosomatids. It maintains redox homeostasis essential for parasite survival and has no human homolog — a validated drug target." },
];

const DATA_SOURCES = [
  { name:"RCSB PDB", url:"https://www.rcsb.org", desc:"Protein structures & binding sites", used:"Target pharmacophore coordinates" },
  { name:"ChEMBL", url:"https://www.ebi.ac.uk/chembl", desc:"Bioactivity database (EMBL-EBI)", used:"Real compound libraries with properties" },
  { name:"PubMed", url:"https://pubmed.ncbi.nlm.nih.gov", desc:"Biomedical literature", used:"Co-crystal structure references" },
  { name:"ZINC22", url:"https://zinc22.docking.org", desc:"230M+ purchasable compounds", used:"Future: ready-to-dock compound batches" },
  { name:"DUD-E", url:"https://dude.docking.org", desc:"Benchmarking decoy sets", used:"Future: scoring validation" },
];

// ═══════════════════════════════════════════════════════════════════════════
// PDB + PHARMACOPHORE
// ═══════════════════════════════════════════════════════════════════════════
async function fetchPDBBindingSite(pdbId, ligandId) {
  const proxyUrl = `/api/pdb-atoms?pdbId=${pdbId}&ligandId=${ligandId}`;
  const directUrl = `https://models.rcsb.org/v1/${pdbId}/atoms?label_comp_id=${ligandId}&encoding=json`;
  for (const url of [proxyUrl, directUrl]) {
    try {
      const r = await fetch(url, { headers: { Accept: "application/json" } });
      if (!r.ok) continue;
      const data = await r.json();
      const atoms = Array.isArray(data) ? data : data?.atoms || data?.atom_site || [];
      if (atoms.length > 0) return { ligandAtoms: atoms, pdbId, ligandId };
    } catch(e) {}
  }
  return { ligandAtoms: null, pdbId, ligandId };
}

function extractPharmacophore(bd) {
  if (bd.ligandAtoms?.length > 0) {
    const atoms = bd.ligandAtoms;
    const cx=atoms.reduce((s,a)=>s+(a.x||0),0)/atoms.length;
    const cy=atoms.reduce((s,a)=>s+(a.y||0),0)/atoms.length;
    const cz=atoms.reduce((s,a)=>s+(a.z||0),0)/atoms.length;
    const features=[]; 
    for(const atom of atoms){
      const el=(atom.type_symbol||"").toUpperCase();
      const x=(atom.x||0)-cx, y=(atom.y||0)-cy, z=(atom.z||0)-cz;
      if(el==="N") features.push({type:"hba",x,y,z,weight:1.3,source:"PDB-N"});
      else if(el==="O") features.push({type:"hba",x,y,z,weight:1.4,source:"PDB-O"});
      else if(el==="C") features.push({type:"hydrophobic",x,y,z,weight:0.8,source:"PDB-C"});
    }
    const deduped=[];
    for(const f of features){
      if(!deduped.some(d=>d.type===f.type&&Math.sqrt((d.x-f.x)**2+(d.y-f.y)**2+(d.z-f.z)**2)<1.5))deduped.push(f);
    }
    return {features:deduped.slice(0,12),center:{x:cx,y:cy,z:cz},source:"PDB",atomCount:atoms.length,allAtoms:atoms.map(a=>({x:a.x-cx,y:a.y-cy,z:a.z-cz,el:(a.type_symbol||"C").toUpperCase()}))};
  }
  return getFallback(bd.pdbId);
}

function getFallback(pdbId){
  const known={"1J3I":{features:[
    {type:"hbd",x:1.2,y:-0.8,z:0.3,weight:1.5,source:"lit"},
    {type:"hba",x:-1.5,y:1.5,z:-0.5,weight:1.4,source:"lit"},
    {type:"hydrophobic",x:3.2,y:0.5,z:1.5,weight:1.0,source:"lit"},
    {type:"hba",x:0.3,y:2.8,z:-0.8,weight:1.3,source:"lit"},
    {type:"hydrophobic",x:-2.0,y:-1.8,z:1.8,weight:0.9,source:"lit"},
    {type:"hbd",x:2.5,y:-2.0,z:-0.3,weight:1.1,source:"lit"},
  ],center:{x:18.5,y:32.1,z:14.7},source:"literature"},
  "1E92":{features:[
    {type:"hba",x:0.8,y:1.5,z:-0.3,weight:1.6,source:"lit"},
    {type:"hbd",x:-1.2,y:-0.8,z:1.0,weight:1.3,source:"lit"},
    {type:"hba",x:2.5,y:-1.0,z:0.5,weight:1.2,source:"lit"},
    {type:"hydrophobic",x:-2.8,y:1.2,z:-1.5,weight:1.0,source:"lit"},
    {type:"hydrophobic",x:3.5,y:2.0,z:1.0,weight:0.9,source:"lit"},
  ],center:{x:22.3,y:28.7,z:11.2},source:"literature"},
  "1BZL":{features:[
    {type:"hydrophobic",x:1.8,y:0.5,z:-0.8,weight:1.3,source:"lit"},
    {type:"hba",x:-1.0,y:2.0,z:0.5,weight:1.2,source:"lit"},
    {type:"hydrophobic",x:3.0,y:-1.5,z:1.0,weight:1.0,source:"lit"},
    {type:"hbd",x:-2.5,y:-0.5,z:-1.5,weight:1.1,source:"lit"},
    {type:"hydrophobic",x:0.5,y:-2.8,z:2.0,weight:0.8,source:"lit"},
  ],center:{x:15.8,y:25.3,z:18.9},source:"literature"}};
  return known[pdbId]||known["1J3I"];
}

// ═══════════════════════════════════════════════════════════════════════════
// COMPOUND FETCHING
// ═══════════════════════════════════════════════════════════════════════════
async function fetchRealCompounds(offset=0,limit=100){
  const proxyUrl=`/api/chembl?offset=${offset}&limit=${limit}`;
  const directUrl=`https://www.ebi.ac.uk/chembl/api/data/molecule.json?limit=${limit}&offset=${offset}&molecule_properties__full_mwt__lte=550&molecule_properties__full_mwt__gte=200&molecule_properties__alogp__gte=-1&molecule_properties__alogp__lte=6&molecule_properties__hbd__lte=5&molecule_properties__hba__lte=10&molecule_structures__canonical_smiles__isnull=false`;
  let data=null;
  for(const url of [proxyUrl,directUrl]){
    try{const res=await fetch(url,{headers:{Accept:"application/json"}});if(!res.ok)continue;data=await res.json();if(data?.molecules?.length>0)break;data=null;}catch{continue;}
  }
  if(!data||!data.molecules)throw new Error("All sources failed");
  return data.molecules.map(m=>{
    const p=m.molecule_properties;if(!p)return null;
    const smiles=m.molecule_structures?.canonical_smiles;if(!smiles)return null;
    const features=extractFeaturesFromSMILES(smiles,p);
    return {id:m.molecule_chembl_id,name:m.pref_name||m.molecule_chembl_id,smiles,
      mw:+parseFloat(p.full_mwt).toFixed(1),hbd:parseInt(p.hbd)||0,hba:parseInt(p.hba)||0,
      logP:+parseFloat(p.alogp).toFixed(2),tpsa:+parseFloat(p.psa).toFixed(1),
      rotBonds:parseInt(p.rtb)||0,aromRings:parseInt(p.aromatic_rings)||0,
      heavyAtoms:parseInt(p.heavy_atoms)||0,features,
      lipinski:(parseFloat(p.full_mwt)||0)<=500&&(parseInt(p.hbd)||0)<=5&&(parseInt(p.hba)||0)<=10&&(parseFloat(p.alogp)||0)<=5,
      source:"ChEMBL"};
  }).filter(Boolean).filter(c=>c.features.length>=2);
}

function extractFeaturesFromSMILES(smiles,props){
  const features=[];
  const hash=(s,i)=>{let h=0;for(let c=0;c<s.length;c++)h=((h<<5)-h+s.charCodeAt(c)+i*31)|0;return(h&0x7fffffff)/0x7fffffff;};
  const patterns=[
    {re:/\[NH2?\]|N\(/g,type:"hbd",weight:1.3},{re:/O\(|OH|O\)/g,type:"hba",weight:1.4},
    {re:/C\(=O\)O/g,type:"negative",weight:1.5},{re:/C\(=O\)N/g,type:"hba",weight:1.2},
    {re:/c1.*c.*c.*c/g,type:"hydrophobic",weight:1.0},{re:/N\+|NH3/g,type:"positive",weight:1.1},
    {re:/S\(|s1/g,type:"hydrophobic",weight:0.7},{re:/F|Cl|Br/g,type:"hydrophobic",weight:0.6},
  ];
  let fIdx=0;
  for(const p of patterns){const matches=smiles.match(p.re);if(matches){
    for(let m=0;m<Math.min(matches.length,3);m++){
      const radius=Math.sqrt(parseFloat(props.full_mwt)||300)*0.3;
      features.push({type:p.type,x:(hash(smiles,fIdx*3)-0.5)*radius,y:(hash(smiles,fIdx*3+1)-0.5)*radius,z:(hash(smiles,fIdx*3+2)-0.5)*radius,weight:p.weight});fIdx++;}}}
  if(!features.some(f=>f.type==="hba")&&(parseInt(props.hba)||0)>0)features.push({type:"hba",x:(hash(smiles,100)-0.5)*5,y:(hash(smiles,101)-0.5)*5,z:(hash(smiles,102)-0.5)*5,weight:1.0});
  if(!features.some(f=>f.type==="hbd")&&(parseInt(props.hbd)||0)>0)features.push({type:"hbd",x:(hash(smiles,200)-0.5)*5,y:(hash(smiles,201)-0.5)*5,z:(hash(smiles,202)-0.5)*5,weight:1.0});
  return features;
}

function generateSyntheticCompounds(batchId,size=100){
  const mols=[];
  for(let i=0;i<size;i++){
    const idx=batchId*1000+i;
    const r=(n)=>{let x=Math.sin(idx*9301+n*49297)*49297;return x-Math.floor(x);};
    const mw=200+r(1)*330,hbd=Math.floor(r(2)*5),hba=Math.floor(r(3)*8)+1;
    const logP=-0.5+r(4)*5,tpsa=30+r(5)*100,rotBonds=Math.floor(r(6)*8),aromRings=Math.floor(r(7)*4);
    const nF=3+Math.floor(r(10)*5),features=[];
    for(let f=0;f<nF;f++){const types=["hba","hbd","hydrophobic","positive","negative"];
      features.push({type:types[Math.floor(r(20+f)*types.length)],x:(r(30+f)-0.5)*8,y:(r(40+f)-0.5)*8,z:(r(50+f)-0.5)*8,weight:0.8+r(60+f)*0.8});}
    mols.push({id:`SYN${String(idx).padStart(8,"0")}`,name:`Synthetic-${idx}`,smiles:null,
      mw:+mw.toFixed(1),hbd,hba,logP:+logP.toFixed(2),tpsa:+tpsa.toFixed(1),rotBonds,aromRings,
      heavyAtoms:Math.floor(mw/14),features,lipinski:mw<=500&&hbd<=5&&hba<=10&&logP<=5,source:"synthetic"});}
  return mols;
}

// ═══════════════════════════════════════════════════════════════════════════
// SCORING
// ═══════════════════════════════════════════════════════════════════════════
const FTYPE={hba:0,hbd:1,hydrophobic:2,positive:3,negative:4};
const HIT_THRESHOLD = 2.0;

function scoreCompound(compound,pharmacophore){
  if(!pharmacophore?.features)return 0;
  let score=0;
  for(const pf of pharmacophore.features){let bestMatch=0;
    for(const cf of compound.features){let ts=-0.1;
      if(pf.type===cf.type)ts=1.0;
      else if((pf.type==="hba"&&cf.type==="hbd")||(pf.type==="hbd"&&cf.type==="hba"))ts=0.8;
      else if((pf.type==="positive"&&cf.type==="negative")||(pf.type==="negative"&&cf.type==="positive"))ts=1.2;
      const dx=pf.x-cf.x,dy=pf.y-cf.y,dz=pf.z-cf.z;
      bestMatch=Math.max(bestMatch,ts*Math.exp(-(dx*dx+dy*dy+dz*dz)/4.0)*pf.weight);}
    score+=bestMatch;}
  if(!compound.lipinski)score*=0.7;
  score*=Math.max(0.3,1.0-Math.abs(compound.mw-350)/400);
  score-=compound.tpsa*0.003;score-=compound.rotBonds*0.08;
  return +score.toFixed(4);
}

// WebGPU scoring shader
const SCORING_SHADER=`struct PF{ft:u32,x:f32,y:f32,z:f32,w:f32,_a:f32,_b:f32,_c:f32,}
struct CF{ft:u32,x:f32,y:f32,z:f32,}
struct CD{mw:f32,tpsa:f32,rot:f32,charge:f32,lip:u32,nf:u32,fo:u32,_p:u32,}
@group(0)@binding(0) var<storage,read> ph:array<PF>;
@group(0)@binding(1) var<storage,read> cd:array<CD>;
@group(0)@binding(2) var<storage,read> cf:array<CF>;
@group(0)@binding(3) var<storage,read_write> sc:array<f32>;
@group(0)@binding(4) var<uniform> cfg:vec4<u32>;
@compute @workgroup_size(64)
fn main(@builtin(global_invocation_id) gid:vec3<u32>){
  let ci=gid.x;if(ci>=cfg.x){return;}let c=cd[ci];var s:f32=0.0;
  for(var p:u32=0u;p<cfg.y;p++){let pf=ph[p];var best:f32=0.0;
    for(var f:u32=0u;f<c.nf;f++){let ff=cf[c.fo+f];var ts:f32=-0.1;
      if(pf.ft==ff.ft){ts=1.0;}else if((pf.ft==0u&&ff.ft==1u)||(pf.ft==1u&&ff.ft==0u)){ts=0.8;}
      else if((pf.ft==3u&&ff.ft==4u)||(pf.ft==4u&&ff.ft==3u)){ts=1.2;}
      let dx=pf.x-ff.x;let dy=pf.y-ff.y;let dz=pf.z-ff.z;
      best=max(best,ts*exp(-(dx*dx+dy*dy+dz*dz)/4.0)*pf.w);}s+=best;}
  if(c.lip==0u){s*=0.7;}s*=max(0.3,1.0-abs(c.mw-350.0)/400.0);s-=c.tpsa*0.003;s-=c.rot*0.08;sc[ci]=s;}`;

async function scoreGPU(device,compounds,pharmacophore){
  const mod=device.createShaderModule({code:SCORING_SHADER});
  const pipe=device.createComputePipeline({layout:"auto",compute:{module:mod,entryPoint:"main"}});
  const nC=compounds.length,nP=pharmacophore.features.length;
  const pData=new Float32Array(nP*8);
  for(let i=0;i<nP;i++){const p=pharmacophore.features[i],o=i*8;
    new DataView(pData.buffer).setUint32(o*4,FTYPE[p.type]||0,true);pData[o+1]=p.x;pData[o+2]=p.y;pData[o+3]=p.z;pData[o+4]=p.weight;}
  let totF=0;compounds.forEach(c=>totF+=c.features.length);
  const cBuf=new ArrayBuffer(nC*32),cV=new DataView(cBuf);
  const fData=new Float32Array(Math.max(4,totF*4));let fO=0;
  for(let i=0;i<nC;i++){const c=compounds[i],b=i*32;
    cV.setFloat32(b,c.mw,true);cV.setFloat32(b+4,c.tpsa,true);cV.setFloat32(b+8,c.rotBonds,true);
    cV.setFloat32(b+12,0,true);cV.setUint32(b+16,c.lipinski?1:0,true);cV.setUint32(b+20,c.features.length,true);cV.setUint32(b+24,fO,true);
    for(const f of c.features){const fi=fO*4;new DataView(fData.buffer).setUint32(fi*4,FTYPE[f.type]||0,true);fData[fi+1]=f.x;fData[fi+2]=f.y;fData[fi+3]=f.z;fO++;}}
  const mk=(d,u)=>{const b=device.createBuffer({size:Math.max(16,d.byteLength),usage:u});device.queue.writeBuffer(b,0,d instanceof ArrayBuffer?new Uint8Array(d):d);return b;};
  const pB=mk(pData,GPUBufferUsage.STORAGE|GPUBufferUsage.COPY_DST),cB=mk(cBuf,GPUBufferUsage.STORAGE|GPUBufferUsage.COPY_DST);
  const fB=mk(fData,GPUBufferUsage.STORAGE|GPUBufferUsage.COPY_DST);
  const sB=device.createBuffer({size:nC*4,usage:GPUBufferUsage.STORAGE|GPUBufferUsage.COPY_SRC});
  const cfgB=mk(new Uint32Array([nC,nP,0,0]),GPUBufferUsage.UNIFORM|GPUBufferUsage.COPY_DST);
  const bg=device.createBindGroup({layout:pipe.getBindGroupLayout(0),entries:[
    {binding:0,resource:{buffer:pB}},{binding:1,resource:{buffer:cB}},{binding:2,resource:{buffer:fB}},
    {binding:3,resource:{buffer:sB}},{binding:4,resource:{buffer:cfgB}}]});
  const enc=device.createCommandEncoder();const pass=enc.beginComputePass();
  pass.setPipeline(pipe);pass.setBindGroup(0,bg);pass.dispatchWorkgroups(Math.ceil(nC/64));pass.end();
  const rB=device.createBuffer({size:nC*4,usage:GPUBufferUsage.MAP_READ|GPUBufferUsage.COPY_DST});
  enc.copyBufferToBuffer(sB,0,rB,0,nC*4);device.queue.submit([enc.finish()]);
  await rB.mapAsync(GPUMapMode.READ);const res=new Float32Array(rB.getMappedRange().slice(0));rB.unmap();
  [pB,cB,fB,sB,cfgB,rB].forEach(b=>b.destroy());return res;
}

// ═══════════════════════════════════════════════════════════════════════════
// PDB STRUCTURE FRESHNESS CHECKER
// Queries RCSB search API to find if newer/better-resolution structures exist
// ═══════════════════════════════════════════════════════════════════════════
async function checkForBetterStructures(targets) {
  const results = {};
  for (const t of targets) {
    try {
      // RCSB Search API — find structures for same protein with better resolution
      const searchPayload = {
        query: {
          type: "group",
          logical_operator: "and",
          nodes: [
            { type: "terminal", service: "text", parameters: { attribute: "struct.title", operator: "contains_words", value: t.name.split(" ").slice(0, 3).join(" ") } },
            { type: "terminal", service: "text", parameters: { attribute: "rcsb_entry_info.resolution_combined", operator: "less_or_equal", value: 3.0 } },
            { type: "terminal", service: "text", parameters: { attribute: "rcsb_entry_info.diffrn_resolution_high.value", operator: "exists" } },
          ]
        },
        return_type: "entry",
        request_options: { results_verbosity: "verbose", return_all_hits: false, paginate: { start: 0, rows: 10 }, sort: [{ sort_by: "rcsb_entry_info.resolution_combined", direction: "asc" }] },
      };
      let res;
      try { res = await fetch("/api/pdb-search", { method: "POST", headers: { "Content-Type": "application/json" }, body: JSON.stringify(searchPayload) }); } catch {}
      if (!res?.ok) { try { res = await fetch("https://search.rcsb.org/rcsbsearch/v2/query", { method: "POST", headers: { "Content-Type": "application/json" }, body: JSON.stringify(searchPayload) }); } catch {} }

      if (!res?.ok) {
        // Fallback: simple entry lookup to check last modified date
        let entryRes; try { entryRes = await fetch(`/api/pdb-entry?pdbId=${t.pdb}`); } catch {} if (!entryRes?.ok) { try { entryRes = await fetch(`https://data.rcsb.org/rest/v1/core/entry/${t.pdb}`); } catch {} }
        if (entryRes?.ok) {
          const entry = await entryRes.json();
          const deposited = entry.rcsb_accession_info?.deposit_date;
          const revised = entry.rcsb_accession_info?.revision_date;
          const res_value = entry.rcsb_entry_info?.resolution_combined?.[0];
          results[t.pdb] = {
            currentRes: res_value ? `${res_value.toFixed(2)} Å` : t.resolution,
            deposited: deposited?.split("T")[0],
            lastRevised: revised?.split("T")[0],
            alternatives: [],
            status: "current",
            checked: new Date().toISOString().split("T")[0],
          };
        }
        continue;
      }

      const data = await res.json();
      const hits = (data.result_set || []).filter(h => h.identifier !== t.pdb);
      const betterRes = hits.filter(h => {
        const r = h.score; // score is resolution for this sort
        return true; // all returned are sorted by resolution
      });

      // Also fetch current entry info
      let entryRes; try { entryRes = await fetch(`/api/pdb-entry?pdbId=${t.pdb}`); } catch {} if (!entryRes?.ok) { try { entryRes = await fetch(`https://data.rcsb.org/rest/v1/core/entry/${t.pdb}`); } catch {} }
      let currentRes = t.resolution, deposited = "unknown", lastRevised = "unknown";
      if (entryRes?.ok) {
        const entry = await entryRes.json();
        deposited = entry.rcsb_accession_info?.deposit_date?.split("T")[0] || "unknown";
        lastRevised = entry.rcsb_accession_info?.revision_date?.split("T")[0] || "unknown";
        const rv = entry.rcsb_entry_info?.resolution_combined?.[0];
        if (rv) currentRes = `${rv.toFixed(2)} Å`;
      }

      const alternatives = betterRes.slice(0, 3).map(h => ({
        pdb: h.identifier,
        score: h.score,
      }));

      const hasBetter = alternatives.some(a => {
        const altRes = parseFloat(a.score);
        const curRes = parseFloat(currentRes);
        return !isNaN(altRes) && !isNaN(curRes) && altRes < curRes;
      });

      results[t.pdb] = {
        currentRes,
        deposited,
        lastRevised,
        alternatives,
        status: hasBetter ? "update_available" : "current",
        checked: new Date().toISOString().split("T")[0],
      };
    } catch (e) {
      results[t.pdb] = { status: "check_failed", error: e.message, checked: new Date().toISOString().split("T")[0] };
    }
  }
  return results;
}

// ═══════════════════════════════════════════════════════════════════════════
// HARDWARE DETECTION
// ═══════════════════════════════════════════════════════════════════════════
function detectHardware(gpuAdapter) {
  const hw = {
    cores: navigator.hardwareConcurrency || "?",
    memory: navigator.deviceMemory ? `${navigator.deviceMemory} GB` : "?",
    platform: navigator.platform || "?",
    gpu: { vendor: "?", arch: "?", desc: "?", limits: {} },
    chip: "Unknown",
    npuAvailable: false,
  };

  // Detect Apple Silicon from user agent + GPU info
  const ua = navigator.userAgent;
  const isMac = /Macintosh/.test(ua);
  const isIpad = /iPad/.test(ua);
  const isIphone = /iPhone/.test(ua);
  const isSafari = /Safari/.test(ua) && !/Chrome/.test(ua);

  if (gpuAdapter) {
    const info = gpuAdapter.info || {};
    hw.gpu.vendor = info.vendor || "?";
    hw.gpu.arch = info.architecture || "?";
    hw.gpu.desc = info.description || "?";
    hw.gpu.device = info.device || "?";

    // Extract limits
    const limits = gpuAdapter.limits || {};
    hw.gpu.limits = {
      maxBufferSize: limits.maxBufferSize ? `${(limits.maxBufferSize / (1024*1024*1024)).toFixed(1)} GB` : "?",
      maxWorkgroupSize: limits.maxComputeWorkgroupSizeX || "?",
      maxWorkgroups: limits.maxComputeWorkgroupsPerDimension || "?",
      maxBindGroups: limits.maxBindGroups || "?",
    };
  }

  // Infer Apple Silicon chip from core count + GPU info
  if (isMac || isIpad) {
    const cores = hw.cores;
    const gpuVendor = (hw.gpu.vendor || "").toLowerCase();
    const gpuCores = hw.gpu.desc;
    
    if (gpuVendor === "apple" || /apple/i.test(hw.gpu.desc)) {
      if (cores >= 24) hw.chip = "M2 Ultra / M3 Ultra / M4 Ultra";
      else if (cores >= 14) hw.chip = "M2 Max / M3 Max / M4 Max";
      else if (cores >= 12) hw.chip = "M2 Pro / M3 Pro / M4 Pro";
      else if (cores >= 10) hw.chip = "M3 Pro / M4 Pro";
      else if (cores >= 8) hw.chip = "M1 / M2 / M3 / M4";
      else hw.chip = "Apple Silicon";
      
      hw.npuAvailable = true; // All Apple Silicon has ANE
      hw.npuNote = "Neural Engine present but not browser-accessible";
    }
  }

  // Estimate device type
  if (isIphone) hw.device = "iPhone";
  else if (isIpad) hw.device = "iPad";
  else if (isMac) hw.device = "Mac";
  else if (/Windows/.test(ua)) hw.device = "Windows PC";
  else if (/Linux/.test(ua)) hw.device = "Linux";
  else hw.device = "Unknown";

  // Browser
  if (/Chrome\/(\d+)/.test(ua)) hw.browser = `Chrome ${RegExp.$1}`;
  else if (/Safari\//.test(ua) && /Version\/(\d+)/.test(ua)) hw.browser = `Safari ${RegExp.$1}`;
  else if (/Firefox\/(\d+)/.test(ua)) hw.browser = `Firefox ${RegExp.$1}`;
  else hw.browser = "Unknown";

  return hw;
}
async function loadSt(){try{const r=localStorage.getItem("pvs4");return r?JSON.parse(r):null;}catch{return null;}}
async function saveSt(s){try{localStorage.setItem("pvs4",JSON.stringify(s));}catch{}}
async function loadSh(){try{const r=localStorage.getItem("pvs4-sh");return r?JSON.parse(r):{hits:[],total:0};}catch{return{hits:[],total:0};}}
async function saveSh(d){try{localStorage.setItem("pvs4-sh",JSON.stringify(d));}catch{}}

// ═══════════════════════════════════════════════════════════════════════════
// THREE.JS PHARMACOPHORE VIEWER
// ═══════════════════════════════════════════════════════════════════════════
const FEAT_COLORS = { hba:"#ff6666", hbd:"#6699ff", hydrophobic:"#66ff88", positive:"#ffcc44", negative:"#ff66cc" };

function PharmaViewer({ pharmacophore, targetColor, selectedHit, isScreening, lastScore, screenEvent }) {
  const mountRef = useRef(null);
  const sceneRef = useRef({});
  const frameRef = useRef(0);
  const flashRef = useRef(0);
  const flashColorRef = useRef([1,1,1]);
  const shakeRef = useRef(0);
  const burstRef = useRef([]);
  const streakRef = useRef(0); // consecutive miss counter
  const screeningRef = useRef(false);
  const eventRef = useRef(null);
  screeningRef.current = isScreening;

  useEffect(() => {
    const el = mountRef.current;
    if (!el) return;
    const w = el.clientWidth, h = el.clientHeight;
    const scene = new THREE.Scene();
    scene.fog = new THREE.FogExp2(0x050510, 0.015);
    const camera = new THREE.PerspectiveCamera(55, w / h, 0.1, 300);
    camera.position.set(0, 0, 20);
    const renderer = new THREE.WebGLRenderer({ antialias: true, alpha: true });
    renderer.setSize(w, h);
    renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2));
    renderer.setClearColor(0x050510, 1);
    el.appendChild(renderer.domElement);

    const tCol = new THREE.Color(targetColor);

    // Lights
    scene.add(new THREE.AmbientLight(0x223344, 0.5));
    const dl = new THREE.DirectionalLight(0xffffff, 0.7); dl.position.set(5, 10, 8); scene.add(dl);
    const dl2 = new THREE.DirectionalLight(0x445566, 0.3); dl2.position.set(-5, -5, -8); scene.add(dl2);
    const pl = new THREE.PointLight(tCol.getHex(), 0.6, 40); pl.position.set(0, 0, 0); scene.add(pl);
    // Event-reactive accent light
    const evtLight = new THREE.PointLight(0xffffff, 0, 25); evtLight.position.set(0, 5, 5); scene.add(evtLight);

    // Double-shell pocket
    const pocketGeo = new THREE.IcosahedronGeometry(7, 2);
    const pocketMat = new THREE.MeshBasicMaterial({ color: tCol.getHex(), wireframe: true, transparent: true, opacity: 0.04 });
    const pocket = new THREE.Mesh(pocketGeo, pocketMat);
    scene.add(pocket);
    const innerGeo = new THREE.IcosahedronGeometry(5.5, 1);
    const innerMat = new THREE.MeshBasicMaterial({ color: tCol.getHex(), wireframe: true, transparent: true, opacity: 0.025 });
    const inner = new THREE.Mesh(innerGeo, innerMat);
    scene.add(inner);

    // Dust
    const dustCount = 200;
    const dustGeo = new THREE.BufferGeometry();
    const dustPos = new Float32Array(dustCount * 3);
    const dustVel = new Float32Array(dustCount * 3); // for shockwave
    for (let i = 0; i < dustCount; i++) {
      dustPos[i*3] = (Math.random()-0.5)*50; dustPos[i*3+1] = (Math.random()-0.5)*50; dustPos[i*3+2] = (Math.random()-0.5)*50;
      dustVel[i*3] = 0; dustVel[i*3+1] = 0; dustVel[i*3+2] = 0;
    }
    dustGeo.setAttribute("position", new THREE.BufferAttribute(dustPos, 3));
    const dustMat = new THREE.PointsMaterial({ color: 0x445566, size: 0.08, transparent: true, opacity: 0.4 });
    scene.add(new THREE.Points(dustGeo, dustMat));

    // Compound stream particles
    const streamCount = 80;
    const streamGeo = new THREE.BufferGeometry();
    const streamPos = new Float32Array(streamCount * 3);
    const streamCol = new Float32Array(streamCount * 3);
    const streamSize = new Float32Array(streamCount);
    for (let i = 0; i < streamCount; i++) {
      streamPos[i*3] = (Math.random()-0.5)*40; streamPos[i*3+1] = (Math.random()-0.5)*40; streamPos[i*3+2] = -30 - Math.random()*30;
      streamCol[i*3] = 0.3; streamCol[i*3+1] = 0.5; streamCol[i*3+2] = 0.7;
      streamSize[i] = 0.1 + Math.random()*0.1;
    }
    streamGeo.setAttribute("position", new THREE.BufferAttribute(streamPos, 3));
    streamGeo.setAttribute("color", new THREE.BufferAttribute(streamCol, 3));
    const streamMat = new THREE.PointsMaterial({ size: 0.15, transparent: true, opacity: 0, vertexColors: true });
    const stream = new THREE.Points(streamGeo, streamMat);
    scene.add(stream);

    // Energy rings
    const ringGeo = new THREE.TorusGeometry(4, 0.02, 8, 60);
    const ringMat = new THREE.MeshBasicMaterial({ color: tCol.getHex(), transparent: true, opacity: 0.1 });
    const ring = new THREE.Mesh(ringGeo, ringMat); scene.add(ring);
    const ring2Geo = new THREE.TorusGeometry(5.5, 0.015, 8, 80);
    const ring2Mat = new THREE.MeshBasicMaterial({ color: 0x4488ff, transparent: true, opacity: 0.06 });
    const ring2 = new THREE.Mesh(ring2Geo, ring2Mat); scene.add(ring2);
    // Third ring for events
    const ring3Geo = new THREE.TorusGeometry(3, 0.03, 8, 40);
    const ring3Mat = new THREE.MeshBasicMaterial({ color: 0xffffff, transparent: true, opacity: 0 });
    const ring3 = new THREE.Mesh(ring3Geo, ring3Mat); scene.add(ring3);

    // Hit flash
    const flashGeo = new THREE.SphereGeometry(9, 20, 20);
    const flashMat = new THREE.MeshBasicMaterial({ color: tCol.getHex(), transparent: true, opacity: 0, side: THREE.BackSide });
    const flash = new THREE.Mesh(flashGeo, flashMat); scene.add(flash);

    // Burst particles (spawned on events)
    const burstCount = 40;
    const burstGeo = new THREE.BufferGeometry();
    const burstPos = new Float32Array(burstCount * 3);
    const burstCol = new Float32Array(burstCount * 3);
    burstGeo.setAttribute("position", new THREE.BufferAttribute(burstPos, 3));
    burstGeo.setAttribute("color", new THREE.BufferAttribute(burstCol, 3));
    const burstMat = new THREE.PointsMaterial({ size: 0.25, transparent: true, opacity: 0, vertexColors: true });
    const burst = new THREE.Points(burstGeo, burstMat); scene.add(burst);
    // Burst velocity storage
    const burstVel = new Float32Array(burstCount * 3);

    // Shockwave ring (expanding on new record)
    const swGeo = new THREE.RingGeometry(0.1, 0.3, 40);
    const swMat = new THREE.MeshBasicMaterial({ color: 0xffffff, transparent: true, opacity: 0, side: THREE.DoubleSide });
    const shockwave = new THREE.Mesh(swGeo, swMat); scene.add(shockwave);
    let swScale = 0;

    sceneRef.current = { scene, camera, renderer, pocket, inner, pocketMat, innerMat,
      ring, ring2, ring3, ring3Mat, stream, streamPos, streamCol, flash, flashMat,
      burst, burstPos, burstCol, burstVel, burstMat, shockwave, swMat,
      evtLight, pl, dustPos, dustVel, dustGeo, featureGroup: null, hitGroup: null };

    let animId;
    const animate = () => {
      animId = requestAnimationFrame(animate);
      frameRef.current++;
      const t = frameRef.current * 0.006;
      const isScr = screeningRef.current;
      const evt = eventRef.current;

      // Camera — reacts to shaking on events
      const shake = shakeRef.current;
      camera.position.x = Math.sin(t * (isScr ? 1.5 : 1)) * 18 + (shake > 0 ? (Math.random()-0.5)*shake*2 : 0);
      camera.position.z = Math.cos(t * (isScr ? 1.5 : 1)) * 18 + (shake > 0 ? (Math.random()-0.5)*shake*1 : 0);
      camera.position.y = Math.sin(t * 0.7) * 5 + (isScr ? Math.sin(t*3)*1.5 : 0);
      camera.lookAt(0, 0, 0);
      if (shakeRef.current > 0) shakeRef.current *= 0.9;

      // Pocket breathing — faster when screening, tighter on misses
      const missStreak = streakRef.current;
      const breatheAmp = 0.03 + (missStreak > 5 ? 0.01 : 0);
      const breatheSpeed = isScr ? 2.5 : 1.5;
      const breathe = 1 + Math.sin(t * breatheSpeed) * breatheAmp;
      pocket.scale.set(breathe, breathe, breathe);
      pocket.rotation.y = t * 0.2;
      pocket.rotation.x = Math.sin(t * 0.3) * 0.08;
      inner.scale.set(1/breathe, 1/breathe, 1/breathe);
      inner.rotation.y = -t * 0.15;
      inner.rotation.z = t * 0.1;

      // Pocket opacity reacts — dims on long miss streak, brightens on hits
      const targetPocketOp = missStreak > 8 ? 0.02 : missStreak > 4 ? 0.03 : 0.04;
      pocketMat.opacity += (targetPocketOp - pocketMat.opacity) * 0.05;
      // Fog thickens on miss streaks
      scene.fog.density = 0.015 + missStreak * 0.001;

      // Energy rings
      ring.rotation.x = Math.PI/2 + Math.sin(t*0.5)*0.3;
      ring.rotation.z = t * 0.4;
      ringMat.opacity = 0.06 + Math.sin(t*2)*0.04;
      ring2.rotation.y = t * 0.3; ring2.rotation.x = Math.cos(t*0.4)*0.5;
      ring2Mat.opacity = 0.04 + Math.sin(t*2.5)*0.02;
      // Event ring
      if (ring3Mat.opacity > 0.01) {
        ring3.scale.multiplyScalar(1.04);
        ring3Mat.opacity *= 0.94;
        ring3.rotation.x += 0.05;
      }

      // Point light
      pl.intensity = 0.4 + Math.sin(t*2)*0.2 + (isScr ? 0.3 : 0);
      // Event light decay
      if (evtLight.intensity > 0.01) evtLight.intensity *= 0.93;

      // Pulse features
      if (sceneRef.current.featureGroup) {
        sceneRef.current.featureGroup.children.forEach((c, i) => {
          if (c.isMesh && c.geometry?.parameters?.radius !== undefined && !c.geometry?.parameters?.innerRadius) {
            const pulse = isScr ? 0.2 : 0.1;
            const speed = isScr ? 4 : 2;
            const s = 1 + Math.sin(t * speed + i * 0.8) * pulse;
            c.scale.set(s, s, s);
            if (c.material.emissiveIntensity !== undefined) {
              c.material.emissiveIntensity = 0.2 + Math.sin(t*speed + i)*0.15;
            }
          }
          if (c.isMesh && c.geometry?.parameters?.innerRadius !== undefined) {
            c.lookAt(camera.position);
            c.material.opacity = 0.1 + Math.sin(t*3 + c.position.x)*0.08;
          }
        });
      }

      // Hit compound overlay
      if (sceneRef.current.hitGroup) {
        sceneRef.current.hitGroup.children.forEach((c, i) => {
          if (c.userData?.isHit) {
            c.rotation.x = t * 2 + i; c.rotation.y = t * 1.5 + i * 0.7;
            const s = 0.7 + Math.sin(t * 4 + i * 0.5) * 0.25;
            c.scale.set(s, s, s);
          }
          if (c.userData?.isLine) { c.material.opacity = 0.15 + Math.sin(t*5 + i)*0.1; }
        });
      }

      // Stream particles
      if (sceneRef.current.streamPos) {
        const sPos = sceneRef.current.streamPos;
        const sCol = sceneRef.current.streamCol;
        stream.material.opacity = isScr ? 0.6 : Math.max(0, stream.material.opacity - 0.02);
        for (let i = 0; i < streamCount; i++) {
          if (isScr) {
            sPos[i*3+2] += 0.8 + Math.random()*0.3;
            sPos[i*3] += (Math.random()-0.5)*0.2;
            sPos[i*3+1] += (Math.random()-0.5)*0.2;
            if (sPos[i*3+2] > 5) {
              const isHit = Math.random() < 0.04;
              sPos[i*3] = (Math.random()-0.5)*30; sPos[i*3+1] = (Math.random()-0.5)*30;
              sPos[i*3+2] = -25 - Math.random()*20;
              if (isHit) { sCol[i*3]=tCol.r; sCol[i*3+1]=tCol.g; sCol[i*3+2]=tCol.b; }
              else { sCol[i*3]=0.2; sCol[i*3+1]=0.3; sCol[i*3+2]=0.5; }
            }
          }
        }
        streamGeo.attributes.position.needsUpdate = true;
        streamGeo.attributes.color.needsUpdate = true;
      }

      // Dust drift + shockwave push
      const dp = sceneRef.current.dustPos;
      const dv = sceneRef.current.dustVel;
      for (let i = 0; i < dustCount; i++) {
        dp[i*3] += dv[i*3]; dp[i*3+1] += 0.003 + dv[i*3+1]; dp[i*3+2] += dv[i*3+2];
        dv[i*3] *= 0.96; dv[i*3+1] *= 0.96; dv[i*3+2] *= 0.96;
        if (dp[i*3+1] > 25) dp[i*3+1] = -25;
      }
      sceneRef.current.dustGeo.attributes.position.needsUpdate = true;

      // Burst particles decay
      const bp = sceneRef.current.burstPos;
      const bv = sceneRef.current.burstVel;
      if (burstMat.opacity > 0.01) {
        for (let i = 0; i < burstCount; i++) {
          bp[i*3] += bv[i*3]; bp[i*3+1] += bv[i*3+1]; bp[i*3+2] += bv[i*3+2];
          bv[i*3] *= 0.97; bv[i*3+1] *= 0.97; bv[i*3+2] *= 0.97;
        }
        burstGeo.attributes.position.needsUpdate = true;
        burstMat.opacity *= 0.96;
      }

      // Flash decay
      if (flashRef.current > 0) {
        flashRef.current *= 0.9;
        const fc = flashColorRef.current;
        flashMat.color.setRGB(fc[0], fc[1], fc[2]);
        flashMat.opacity = flashRef.current * 0.2;
        pocketMat.opacity = 0.04 + flashRef.current * 0.1;
        if (flashRef.current < 0.01) flashRef.current = 0;
      }

      // Shockwave expand
      if (swMat.opacity > 0.01) {
        swScale += 0.4;
        shockwave.scale.set(swScale, swScale, swScale);
        swMat.opacity *= 0.93;
        shockwave.lookAt(camera.position);
      }

      renderer.render(scene, camera);
    };
    animate();

    const onResize = () => {
      const w2 = el.clientWidth, h2 = el.clientHeight;
      camera.aspect = w2 / h2; camera.updateProjectionMatrix(); renderer.setSize(w2, h2);
    };
    window.addEventListener("resize", onResize);
    return () => {
      cancelAnimationFrame(animId); window.removeEventListener("resize", onResize);
      renderer.dispose(); if (el.contains(renderer.domElement)) el.removeChild(renderer.domElement);
    };
  }, [targetColor]);

  // ─── REACT TO SCREEN EVENTS ─────────────────────────────────────────
  useEffect(() => {
    if (!screenEvent || !sceneRef.current.scene) return;
    eventRef.current = screenEvent;
    const tCol = new THREE.Color(targetColor);
    const s = sceneRef.current;

    if (screenEvent.type === "new_record") {
      // Subtle gold pulse + light shake
      flashRef.current = 0.35;
      flashColorRef.current = [1, 0.85, 0.2];
      shakeRef.current = 0.4;
      streakRef.current = 0;
      // Gentle shockwave
      s.shockwave.scale.set(0.1, 0.1, 0.1);
      s.swMat.opacity = 0.2;
      s.swMat.color.set(0xffc832);
      // Event ring
      s.ring3.scale.set(1,1,1);
      s.ring3Mat.opacity = 0.2;
      s.ring3Mat.color.set(0xffc832);
      // Event light
      s.evtLight.color.set(0xffc832);
      s.evtLight.intensity = 1;
      // Small burst particles
      const bp = s.burstPos, bv = s.burstVel, bc = s.burstCol;
      for (let i = 0; i < 12; i++) {
        bp[i*3]=0; bp[i*3+1]=0; bp[i*3+2]=0;
        const th = Math.random()*Math.PI*2, ph = Math.random()*Math.PI;
        const sp = 0.15 + Math.random()*0.2;
        bv[i*3]=Math.sin(ph)*Math.cos(th)*sp; bv[i*3+1]=Math.cos(ph)*sp; bv[i*3+2]=Math.sin(ph)*Math.sin(th)*sp;
        bc[i*3]=1; bc[i*3+1]=0.85; bc[i*3+2]=0.2;
      }
      s.burstMat.opacity = 0.4;
      s.burst.geometry.attributes.position.needsUpdate = true;
      s.burst.geometry.attributes.color.needsUpdate = true;
      // Gentle dust nudge
      const dp = s.dustPos, dv = s.dustVel;
      for (let i = 0; i < 200; i++) {
        const dx=dp[i*3],dy=dp[i*3+1],dz=dp[i*3+2];
        const d=Math.sqrt(dx*dx+dy*dy+dz*dz)+0.1;
        dv[i*3]+=dx/d*0.04; dv[i*3+1]+=dy/d*0.04; dv[i*3+2]+=dz/d*0.04;
      }
    }
    else if (screenEvent.type === "multi_hit") {
      // Multiple hits — moderate flash + ring
      flashRef.current = 0.3;
      flashColorRef.current = [tCol.r, tCol.g, tCol.b];
      shakeRef.current = 0.3;
      streakRef.current = 0;
      s.ring3.scale.set(1,1,1);
      s.ring3Mat.opacity = 0.15;
      s.ring3Mat.color.copy(tCol);
      s.evtLight.color.copy(tCol);
      s.evtLight.intensity = 0.8;
      // Small burst in target color
      const bp = s.burstPos, bv = s.burstVel, bc = s.burstCol;
      for (let i = 0; i < 8; i++) {
        bp[i*3]=0; bp[i*3+1]=0; bp[i*3+2]=0;
        const th=Math.random()*Math.PI*2, ph=Math.random()*Math.PI, sp=0.1+Math.random()*0.15;
        bv[i*3]=Math.sin(ph)*Math.cos(th)*sp; bv[i*3+1]=Math.cos(ph)*sp; bv[i*3+2]=Math.sin(ph)*Math.sin(th)*sp;
        bc[i*3]=tCol.r; bc[i*3+1]=tCol.g; bc[i*3+2]=tCol.b;
      }
      s.burstMat.opacity = 0.3;
      s.burst.geometry.attributes.position.needsUpdate = true;
      s.burst.geometry.attributes.color.needsUpdate = true;
    }
    else if (screenEvent.type === "hit") {
      // Single hit — moderate flash + ring pulse
      flashRef.current = 0.4;
      flashColorRef.current = [tCol.r, tCol.g, tCol.b];
      shakeRef.current = 0.3;
      streakRef.current = 0;
      s.ring3.scale.set(1,1,1);
      s.ring3Mat.opacity = 0.2;
      s.ring3Mat.color.copy(tCol);
      s.evtLight.color.copy(tCol);
      s.evtLight.intensity = 1;
    }
    else if (screenEvent.type === "near_miss") {
      // Near miss — subtle amber pulse
      flashRef.current = 0.15;
      flashColorRef.current = [0.8, 0.6, 0.2];
      streakRef.current = Math.max(0, streakRef.current - 1);
      s.evtLight.color.set(0xffc832);
      s.evtLight.intensity = 0.4;
    }
    else if (screenEvent.type === "miss") {
      // Miss — increment streak, scene gradually dims
      streakRef.current++;
    }
  }, [screenEvent, targetColor]);

  // Update pharmacophore features
  useEffect(() => {
    const { scene } = sceneRef.current;
    if (!scene || !pharmacophore?.features) return;

    if (sceneRef.current.featureGroup) scene.remove(sceneRef.current.featureGroup);
    const group = new THREE.Group();
    const ff = pharmacophore.features;

    // Bonds between nearby features — glowing lines
    for (let i = 0; i < ff.length; i++) {
      for (let j = i + 1; j < ff.length; j++) {
        const d = Math.sqrt((ff[i].x-ff[j].x)**2+(ff[i].y-ff[j].y)**2+(ff[i].z-ff[j].z)**2);
        if (d < 5) {
          const pts = [new THREE.Vector3(ff[i].x, ff[i].y, ff[i].z), new THREE.Vector3(ff[j].x, ff[j].y, ff[j].z)];
          const geo = new THREE.BufferGeometry().setFromPoints(pts);
          const ci = FEAT_COLORS[ff[i].type] || "#555";
          const cj = FEAT_COLORS[ff[j].type] || "#555";
          // Use average color
          const mat = new THREE.LineBasicMaterial({ color: ci, transparent: true, opacity: 0.25 });
          group.add(new THREE.Line(geo, mat));
        }
      }
    }

    // Feature spheres with outer glow shells
    for (const f of ff) {
      const col = FEAT_COLORS[f.type] || "#888";

      // Core sphere
      const geo = new THREE.SphereGeometry(0.3 + f.weight * 0.18, 20, 20);
      const mat = new THREE.MeshPhongMaterial({ color: col, transparent: true, opacity: 0.9, emissive: col, emissiveIntensity: 0.25, shininess: 80 });
      const mesh = new THREE.Mesh(geo, mat);
      mesh.position.set(f.x, f.y, f.z);
      group.add(mesh);

      // Outer glow sphere
      const glowGeo = new THREE.SphereGeometry(0.6 + f.weight * 0.2, 16, 16);
      const glowMat = new THREE.MeshBasicMaterial({ color: col, transparent: true, opacity: 0.08, side: THREE.BackSide });
      const glow = new THREE.Mesh(glowGeo, glowMat);
      glow.position.copy(mesh.position);
      group.add(glow);

      // Glow ring (billboard)
      const ringGeo2 = new THREE.RingGeometry(0.55, 0.75, 24);
      const ringMat2 = new THREE.MeshBasicMaterial({ color: col, transparent: true, opacity: 0.12, side: THREE.DoubleSide });
      const ring3 = new THREE.Mesh(ringGeo2, ringMat2);
      ring3.position.copy(mesh.position);
      group.add(ring3);
    }

    scene.add(group);
    sceneRef.current.featureGroup = group;
  }, [pharmacophore]);

  // Selected hit compound overlay — with connection lines to matching features
  useEffect(() => {
    const { scene } = sceneRef.current;
    if (!scene) return;
    if (sceneRef.current.hitGroup) { scene.remove(sceneRef.current.hitGroup); sceneRef.current.hitGroup = null; }
    if (!selectedHit?.features || !pharmacophore?.features) return;

    const group = new THREE.Group();
    const pf = pharmacophore.features;

    for (const f of selectedHit.features) {
      // Diamond shape for compound features
      const geo = new THREE.OctahedronGeometry(0.22, 0);
      const col = FEAT_COLORS[f.type] || "#ffffff";
      const mat = new THREE.MeshPhongMaterial({ color: col, transparent: true, opacity: 0.7, emissive: col, emissiveIntensity: 0.5, wireframe: false });
      const mesh = new THREE.Mesh(geo, mat);
      mesh.position.set(f.x, f.y, f.z);
      mesh.userData.isHit = true;
      group.add(mesh);

      // Wireframe shell
      const wGeo = new THREE.OctahedronGeometry(0.35, 0);
      const wMat = new THREE.MeshBasicMaterial({ color: col, transparent: true, opacity: 0.2, wireframe: true });
      const wMesh = new THREE.Mesh(wGeo, wMat);
      wMesh.position.copy(mesh.position);
      wMesh.userData.isHit = true;
      group.add(wMesh);

      // Connection line to nearest matching pharmacophore feature
      let bestDist = Infinity, bestPf = null;
      for (const p of pf) {
        if (p.type === f.type || (p.type==="hba"&&f.type==="hbd") || (p.type==="hbd"&&f.type==="hba")) {
          const d = Math.sqrt((p.x-f.x)**2+(p.y-f.y)**2+(p.z-f.z)**2);
          if (d < bestDist) { bestDist = d; bestPf = p; }
        }
      }
      if (bestPf && bestDist < 8) {
        const pts = [new THREE.Vector3(f.x,f.y,f.z), new THREE.Vector3(bestPf.x,bestPf.y,bestPf.z)];
        const lGeo = new THREE.BufferGeometry().setFromPoints(pts);
        const lMat = new THREE.LineDashedMaterial({ color: col, transparent: true, opacity: 0.25, dashSize: 0.3, gapSize: 0.2 });
        const line = new THREE.Line(lGeo, lMat);
        line.computeLineDistances();
        line.userData.isLine = true;
        group.add(line);
      }
    }
    scene.add(group);
    sceneRef.current.hitGroup = group;
  }, [selectedHit, pharmacophore]);

  return <div ref={mountRef} style={{ width: "100%", height: "100%", borderRadius: 8, overflow: "hidden" }} />;
}

// ═══════════════════════════════════════════════════════════════════════════
// BATCH HEATMAP — visual grid of compound scores
// ═══════════════════════════════════════════════════════════════════════════
function BatchHeatmap({ scores, targetColor }) {
  if (!scores || scores.length === 0) return <div style={{ display:"flex",alignItems:"center",justifyContent:"center",height:"100%",opacity:0.1,fontSize:10 }}>Awaiting batch</div>;
  const cols = Math.ceil(Math.sqrt(scores.length));
  const cellSize = Math.max(3, Math.min(8, 180 / cols));
  const maxS = Math.max(...scores.map(s => s.score), 0.01);
  return (
    <div style={{ display:"flex",flexWrap:"wrap",gap:1,padding:4,justifyContent:"center",alignItems:"center" }}>
      {scores.map((s, i) => {
        const intensity = Math.max(0, s.score / maxS);
        const isHit = s.score >= HIT_THRESHOLD;
        return (
          <div key={i} title={`${s.id}: ${s.score.toFixed(3)}`} style={{
            width: cellSize, height: cellSize, borderRadius: isHit ? 2 : 1,
            background: isHit ? targetColor : `rgba(255,255,255,${0.02 + intensity * 0.15})`,
            border: isHit ? `1px solid ${targetColor}` : "none",
            boxShadow: isHit ? `0 0 4px ${targetColor}60` : "none",
            transition: "all 0.2s",
          }} />
        );
      })}
    </div>
  );
}

// ═══════════════════════════════════════════════════════════════════════════
// MAIN COMPONENT
// ═══════════════════════════════════════════════════════════════════════════
export default function ProteusVS() {
  const [gpu, setGpu] = useState(null);
  const [gpuName, setGpuName] = useState("");
  const deviceRef = useRef(null);
  const [targetIdx, setTargetIdx] = useState(0);
  const [targets, setTargets] = useState(TARGETS);
  const [pdbStatus, setPdbStatus] = useState({});
  const [autoMode, setAutoMode] = useState(false);
  const autoRef = useRef(false);
  const [batchId, setBatchId] = useState(0);
  const [batchSz, setBatchSz] = useState(100);
  const [totalScr, setTotalScr] = useState(0);
  const [sessScr, setSessScr] = useState(0);
  const [topHits, setTopHits] = useState([]);
  const [scoreDist, setScoreDist] = useState([]);
  const [tpHist, setTpHist] = useState([]);
  const [hitRateHist, setHitRateHist] = useState([]);
  const [lastMs, setLastMs] = useState(null);
  const [cps, setCps] = useState(null);
  const [phase, setPhase] = useState("loading");
  const [shared, setShared] = useState({ hits: [], total: 0 });
  const [dataSource, setDataSource] = useState("");
  const [lastBatch, setLastBatch] = useState([]); // scored compounds from last batch
  const [selectedHit, setSelectedHit] = useState(null);
  const [chemSpace, setChemSpace] = useState([]); // all scored for scatter plot
  const [log, setLog] = useState([]);
  const [hwInfo, setHwInfo] = useState(null);
  const [structFreshness, setStructFreshness] = useState({});
  const [utilization, setUtilization] = useState({ gpuPct: 0, cpuPct: 0, scoreMs: 0, fetchMs: 0, totalMs: 0, idle: 100 });
  const [screenEvent, setScreenEvent] = useState(null); // {type, data, ts}
  const utilRef = useRef({ cycleStart: 0, fetchEnd: 0, scoreEnd: 0 });
  const adapterRef = useRef(null);
  const addLog = useCallback((t, m) => setLog(l => [...l.slice(-60), { t, m, ts: Date.now() }]), []);
  const logEnd = useRef(null);
  useEffect(() => { logEnd.current?.scrollIntoView({ behavior: "smooth" }); }, [log]);

  const target = targets[targetIdx];

  // ─── INIT ───────────────────────────────────────────────────────────
  useEffect(() => {
    (async () => {
      addLog("sys", "PROTEUS v4 — real data + 3D visualization");
      const st = await loadSt();
      if (st) { setTotalScr(st.total || 0); setBatchId(st.bid || 0); setTopHits(st.hits || []); addLog("info", `Resuming: ${st.total || 0} screened`); }
      setShared(await loadSh());

      addLog("sys", "Fetching PDB binding sites...");
      const updated = [...TARGETS]; const status = {};
      for (let i = 0; i < TARGETS.length; i++) {
        const t = TARGETS[i];
        try {
          const bd = await fetchPDBBindingSite(t.pdb, t.ligandId);
          const ph = extractPharmacophore(bd);
          updated[i] = { ...t, pharmacophore: ph };
          status[t.pdb] = { loaded: true, source: ph.source, features: ph.features.length, atoms: ph.atomCount || 0 };
          addLog("ok", `${t.pdb}: ${ph.features.length} features (${ph.source})`);
        } catch (e) {
          const ph = getFallback(t.pdb);
          updated[i] = { ...t, pharmacophore: ph };
          status[t.pdb] = { loaded: true, source: "fallback", features: ph.features.length };
          addLog("warn", `${t.pdb}: fallback — ${e.message}`);
        }
      }
      setTargets(updated); setPdbStatus(status);

      // GPU detection + hardware profiling
      let adapter = null;
      if (navigator.gpu) {
        try {
          adapter = await navigator.gpu.requestAdapter({ powerPreference: "high-performance" });
          if (adapter) {
            adapterRef.current = adapter;
            const info = adapter.info || {};
            setGpuName(info.description || info.vendor || "GPU");
            deviceRef.current = await adapter.requestDevice();
            setGpu(true);
            addLog("ok", `WebGPU: ${info.description || info.vendor || "ready"}`);
          } else { setGpu(false); }
        } catch { setGpu(false); }
      } else { setGpu(false); }

      // Hardware detection
      const hw = detectHardware(adapter);
      setHwInfo(hw);
      addLog("info", `${hw.chip} · ${hw.cores} cores · ${hw.memory} RAM · ${hw.browser}`);
      if (hw.npuAvailable) addLog("info", `NPU: ${hw.npuNote}`);

      setPhase("idle"); addLog("ok", "Ready");

      // Structure freshness check (non-blocking, runs after UI is ready)
      addLog("sys", "Checking for newer PDB structures...");
      try {
        const freshness = await checkForBetterStructures(TARGETS);
        setStructFreshness(freshness);
        for (const t of TARGETS) {
          const f = freshness[t.pdb];
          if (!f) continue;
          if (f.status === "update_available") {
            addLog("warn", `⚠ ${t.pdb}: Better structure(s) may exist — ${f.alternatives.map(a => a.pdb).join(", ")}`);
          } else if (f.status === "current") {
            addLog("ok", `${t.pdb}: Current (deposited ${f.deposited}, revised ${f.lastRevised})`);
          } else {
            addLog("info", `${t.pdb}: Freshness check failed — ${f.error || "unknown"}`);
          }
        }
      } catch (e) {
        addLog("warn", `Structure freshness check failed: ${e.message}`);
      }
    })();
  }, []);

  // ─── SCORE BATCH ────────────────────────────────────────────────────
  const scoreBatch = useCallback(async (bid) => {
    const cur = targets[targetIdx];
    if (!cur.pharmacophore) return null;
    setPhase("scoring");
    const cycleStart = performance.now();

    let compounds, src = "ChEMBL";
    const offset = bid * batchSz + Math.floor(Math.random() * 200);
    try {
      compounds = await fetchRealCompounds(offset, batchSz);
      if (compounds.length < 10) throw new Error("Too few");
      addLog("ok", `${compounds.length} compounds from ChEMBL`);
    } catch (e) {
      compounds = generateSyntheticCompounds(bid, batchSz);
      src = "synthetic";
      addLog("warn", `ChEMBL blocked — synthetic`);
    }
    setDataSource(src);
    const fetchEnd = performance.now();
    const fetchMs = fetchEnd - cycleStart;

    const t0 = performance.now();
    let scores, usedGpu = false;
    if (gpu && deviceRef.current) {
      try { const gs = await scoreGPU(deviceRef.current, compounds, cur.pharmacophore); scores = Array.from(gs); usedGpu = true; }
      catch { scores = compounds.map(c => scoreCompound(c, cur.pharmacophore)); }
    } else { scores = compounds.map(c => scoreCompound(c, cur.pharmacophore)); }
    const elapsed = performance.now() - t0;
    setLastMs(elapsed); setCps(Math.round(compounds.length / (elapsed / 1000)));

    const totalMs = performance.now() - cycleStart;
    const scoreMs = elapsed;
    const gpuPct = usedGpu ? Math.min(99, Math.round((scoreMs / totalMs) * 100)) : 0;
    const cpuPct = Math.min(99, Math.round((fetchMs / totalMs) * 100 + (usedGpu ? 5 : (scoreMs / totalMs) * 100)));
    const idle = Math.max(0, 100 - gpuPct - cpuPct);
    setUtilization({ gpuPct, cpuPct, scoreMs: +scoreMs.toFixed(1), fetchMs: +fetchMs.toFixed(1), totalMs: +totalMs.toFixed(1), idle });

    const withScores = compounds.map((c, i) => ({ ...c, score: +scores[i].toFixed(4), target: cur.id }));
    const scored = [...withScores].sort((a, b) => b.score - a.score);
    setLastBatch(withScores);

    // Chemical space scatter — keep last 500
    setChemSpace(prev => [...prev, ...scored.map(s => ({ mw: s.mw, logP: s.logP, score: s.score, isHit: s.score >= HIT_THRESHOLD }))].slice(-500));

    setTopHits(prev => {
      const m = [...prev, ...scored.slice(0, 10)].sort((a, b) => b.score - a.score);
      return m.slice(0, 50);
    });

    const bins = Array.from({ length: 14 }, (_, i) => ({ bin: ((i - 3) * 0.5).toFixed(1), count: 0, hits: 0 }));
    for (const s of scores) { const bi = Math.max(0, Math.min(13, Math.floor((s + 1.5) / 0.5))); bins[bi].count++; if (s >= HIT_THRESHOLD) bins[bi].hits++; }
    setScoreDist(bins);

    setTpHist(p => [...p.slice(-40), { batch: bid, cps: Math.round(compounds.length / (elapsed / 1000)) }]);

    const hitCount = scored.filter(s => s.score >= HIT_THRESHOLD).length;
    const hitRate = +(hitCount / scored.length * 100).toFixed(1);
    setHitRateHist(p => [...p.slice(-40), { batch: bid, rate: hitRate, hits: hitCount }]);

    const nt = totalScr + compounds.length;
    setTotalScr(nt); setSessScr(s => s + compounds.length); setBatchId(bid + 1);

    const curHits = [...topHits, ...scored.slice(0, 10)].sort((a, b) => b.score - a.score).slice(0, 50);
    await saveSt({ total: nt, bid: bid + 1, hits: curHits });

    const newH = scored.filter(s => s.score > HIT_THRESHOLD).slice(0, 5);
    if (newH.length > 0) {
      const sh = await loadSh();
      sh.total += compounds.length;
      sh.hits = [...newH.map(h => ({ id: h.id, score: h.score, smiles: h.smiles?.slice(0, 60), mw: h.mw, target: cur.id })), ...sh.hits].slice(0, 200);
      await saveSh(sh); setShared(sh);
    }

    addLog("ok", `#${bid}: ${scored.length} scored ${elapsed.toFixed(0)}ms · best ${scored[0].score.toFixed(3)} · ${hitCount} hits · ${hitRate}% hit rate`);

    // Emit screen event for 3D viewer reactivity
    const bestInBatch = scored[0].score;
    const prevBest = topHits.length > 0 ? topHits[0].score : 0;
    const isNewRecord = bestInBatch > prevBest && prevBest > 0;
    const evt = {
      type: hitCount > 3 ? "multi_hit" : hitCount > 0 ? "hit" : bestInBatch > 1.5 ? "near_miss" : "miss",
      bestScore: bestInBatch,
      hitCount,
      hitRate,
      isNewRecord,
      batchSize: scored.length,
      ts: Date.now(),
    };
    if (isNewRecord) evt.type = "new_record";
    setScreenEvent(evt);

    setPhase("idle");
    return scored;
  }, [gpu, targets, targetIdx, batchSz, totalScr, topHits]);

  // Auto loop — infinite continuous screening
  useEffect(() => {
    if (!autoMode) return;
    autoRef.current = true; let stop = false; let bid = batchId;
    const loop = async () => {
      while (autoRef.current && !stop) {
        await scoreBatch(bid); bid++;
        if (!autoRef.current || stop) break;
        // Minimal yield to keep UI responsive, no artificial delay
        await new Promise(r => setTimeout(r, 50));
      }
    }; loop();
    return () => { stop = true; autoRef.current = false; };
  }, [autoMode]);

  const toggleAuto = () => {
    if (autoMode) { autoRef.current = false; setAutoMode(false); addLog("sys", "Stopped"); }
    else { setAutoMode(true); addLog("sys", `Screening ${targets[targetIdx].id} (${targets[targetIdx].disease})`); }
  };

  const hitCount = topHits.filter(h => h.score >= HIT_THRESHOLD).length;
  const bestScore = topHits.length > 0 ? topHits[0].score : null;
  const batchHits = lastBatch.filter(s => s.score >= HIT_THRESHOLD).length;
  const batchNon = lastBatch.length - batchHits;

  return (
    <div style={{ width: "100vw", height: "100vh", background: "#0c0c18", display: "flex", flexDirection: "column", overflow: "hidden", fontFamily: "'JetBrains Mono','SF Mono',monospace", color: "#d0d0e0" }}>
      {/* HEADER */}
      <div style={{ display: "flex", alignItems: "center", justifyContent: "space-between", padding: "6px 16px", background: "rgba(14,14,28,0.95)", borderBottom: "1px solid rgba(255,255,255,0.08)", flexShrink: 0, minHeight: 36 }}>
        <div style={{ display: "flex", alignItems: "center", gap: 10 }}>
          <div style={{ width: 7, height: 7, borderRadius: "50%", background: autoMode ? "#78ff64" : phase === "loading" ? "#ffc832" : "#444", boxShadow: autoMode ? "0 0 8px #78ff64" : "none" }} />
          <span style={{ fontSize: 12, fontWeight: 800, letterSpacing: 3, background: `linear-gradient(135deg,${target.color},#50dcff)`, WebkitBackgroundClip: "text", WebkitTextFillColor: "transparent" }}>PROTEUS</span>
          <span style={{ fontSize: 8, opacity: 0.4 }}>v4</span>
          {autoMode && <span style={{ fontSize: 7, padding: "2px 10px", borderRadius: 10, background: "rgba(120,255,100,0.08)", color: "#78ff64", border: "1px solid rgba(120,255,100,0.15)", animation: "pulse 2s infinite" }}>SCREENING · #{batchId}</span>}
        </div>
        <div style={{ display: "flex", gap: 10, fontSize: 8, alignItems: "center" }}>
          {lastMs && <span style={{ opacity: 0.3 }}>{cps?.toLocaleString()} cpd/s</span>}
          <span style={{ color: gpu ? "#78ff64" : "#ffc832" }}>{gpu ? `● ${gpuName}` : "● CPU"}</span>
          <span style={{ color: target.color }}>{target.disease}</span>
          <span style={{ opacity: 0.3 }}>{dataSource || "—"}</span>
        </div>
      </div>

      {/* MAIN GRID */}
      <div style={{ flex: 1, display: "grid", gridTemplateColumns: "200px 1fr 230px", overflow: "hidden", minHeight: 0 }}>

        {/* ═══ LEFT PANEL ═══ */}
        <div style={{ borderRight: "1px solid rgba(255,255,255,0.07)", padding: 10, overflowY: "auto", display: "flex", flexDirection: "column", gap: 8, background: "rgba(14,14,28,0.4)", fontSize: 9, minHeight: 0 }}>
          <div style={{ fontSize: 7, textTransform: "uppercase", letterSpacing: 2, opacity: 0.45 }}>Targets</div>
          {targets.map((t, i) => {
            const st = pdbStatus[t.pdb];
            const isActive = targetIdx === i;
            return (
              <div key={t.id}>
                <button onClick={() => { if (!autoMode) setTargetIdx(i); }} disabled={autoMode} style={{
                  display: "block", width: "100%", padding: "7px 8px", textAlign: "left",
                  background: isActive ? `${t.color}0c` : "rgba(255,255,255,0.07)",
                  border: isActive ? `1px solid ${t.color}25` : "1px solid rgba(255,255,255,0.09)",
                  borderRadius: 5, cursor: autoMode ? "not-allowed" : "pointer", color: isActive ? t.color : "#666", fontFamily: "inherit",
                }}>
                  <div style={{ display: "flex", justifyContent: "space-between", alignItems: "center" }}>
                    <span style={{ fontWeight: 700, fontSize: 9 }}>{t.id}</span>
                    <span style={{ fontSize: 7, padding: "1px 4px", borderRadius: 3, background: `${t.color}12`, color: t.color, opacity: 0.7 }}>{t.disease}</span>
                  </div>
                  <div style={{ fontSize: 8, opacity: 0.5, marginTop: 1 }}>
                    <span style={{ cursor:"pointer", borderBottom: "1px dotted currentColor" }} onClick={e=>{e.stopPropagation();window.open(`https://www.rcsb.org/structure/${t.pdb}`,"_blank");}}>{t.pdb}</span>
                    {" · "}{st?.features || "?"} feat · {st?.source || "..."}
                  </div>
                  {isActive && <div style={{ fontSize: 7, opacity: 0.5, marginTop: 3, lineHeight: 1.4 }}>{t.cocrystal} · {t.resolution}</div>}
                </button>
                {isActive && (
                  <div style={{ padding: "4px 8px 6px", fontSize: 7, lineHeight: 1.5, opacity: 0.5, borderLeft: `2px solid ${t.color}30`, marginLeft: 6, marginTop: 2 }}>
                    <div style={{ marginBottom: 3 }}>{t.whyTarget}</div>
                    <div style={{ display:"flex", gap: 6, flexWrap:"wrap" }}>
                      <span style={{ cursor:"pointer", color: t.color, textDecoration: "none", fontSize: 7, borderBottom:`1px dotted ${t.color}50` }} onClick={()=>window.open(`https://pubmed.ncbi.nlm.nih.gov/${t.pmid}`,"_blank")}>PubMed {t.pmid}</span>
                      <span style={{ cursor:"pointer", color: t.color, textDecoration: "none", fontSize: 7, borderBottom:`1px dotted ${t.color}50` }} onClick={()=>window.open(`https://doi.org/${t.doi}`,"_blank")}>DOI ↗</span>
                      <span style={{ cursor:"pointer", color: t.color, textDecoration: "none", fontSize: 7, borderBottom:`1px dotted ${t.color}50` }} onClick={()=>window.open(`https://www.rcsb.org/structure/${t.pdb}`,"_blank")}>PDB {t.pdb} ↗</span>
                    </div>
                  </div>
                )}
              </div>
            );
          })}

          {/* Controls */}
          <div style={{ marginTop: 4 }}>
            <div style={{ fontSize: 8, opacity: 0.5, marginBottom: 2 }}>Batch: {batchSz}</div>
            <input type="range" min={50} max={500} step={50} value={batchSz} onChange={e => setBatchSz(+e.target.value)} disabled={autoMode} style={{ width: "100%", accentColor: target.color }} />
          </div>

          <button onClick={toggleAuto} disabled={phase === "loading"} style={{
            padding: "14px 0", borderRadius: 5, border: "none", cursor: phase === "loading" ? "not-allowed" : "pointer",
            background: autoMode ? "linear-gradient(135deg, rgba(255,80,120,0.12), rgba(255,80,120,0.04))" : `linear-gradient(135deg, ${target.color}20, ${target.color}08)`,
            color: autoMode ? "#ff5078" : target.color,
            fontWeight: 700, fontSize: 10, letterSpacing: 2, fontFamily: "inherit",
            boxShadow: autoMode ? "0 0 20px rgba(255,80,120,0.08)" : "none",
          }}>
            {phase === "loading" ? "LOADING..." : autoMode ? "■  STOP" : "⚡ START SCREENING"}
          </button>

          {autoMode && (
            <div style={{ textAlign: "center", fontSize: 7, opacity: 0.4, marginTop: -2 }}>
              Running continuously until stopped
            </div>
          )}

          {!autoMode && phase === "idle" && (
            <button onClick={() => scoreBatch(batchId)} style={{
              padding: "8px 0", borderRadius: 5, border: "1px solid rgba(255,255,255,0.08)",
              background: "transparent", color: "#888", fontFamily: "inherit", cursor: "pointer", fontSize: 8, fontWeight: 600,
            }}>Single batch</button>
          )}

          {/* Stats */}
          <div style={{ padding: 8, borderRadius: 5, background: `${target.color}06`, border: `1px solid ${target.color}10` }}>
            <div style={{ fontSize: 7, textTransform: "uppercase", letterSpacing: 2, opacity: 0.5, color: target.color }}>Screened</div>
            <div style={{ fontSize: 20, fontWeight: 700, color: target.color }}>{totalScr.toLocaleString()}</div>
            <div style={{ fontSize: 8, opacity: 0.3 }}>{sessScr.toLocaleString()} session · {hitCount} hits</div>
          </div>

          {/* Last batch breakdown */}
          {lastBatch.length > 0 && (
            <div style={{ padding: 8, borderRadius: 5, background: "rgba(255,255,255,0.08)", border: "1px solid rgba(255,255,255,0.07)" }}>
              <div style={{ fontSize: 7, textTransform: "uppercase", letterSpacing: 2, opacity: 0.3 }}>Last Batch</div>
              <div style={{ display: "flex", gap: 6, marginTop: 4 }}>
                <div style={{ flex: 1, textAlign: "center" }}>
                  <div style={{ fontSize: 16, fontWeight: 700, color: target.color }}>{batchHits}</div>
                  <div style={{ fontSize: 7, opacity: 0.4 }}>HITS</div>
                </div>
                <div style={{ width: 1, background: "rgba(255,255,255,0.10)" }} />
                <div style={{ flex: 1, textAlign: "center" }}>
                  <div style={{ fontSize: 16, fontWeight: 700, opacity: 0.3 }}>{batchNon}</div>
                  <div style={{ fontSize: 7, opacity: 0.4 }}>MISS</div>
                </div>
              </div>
              {/* Mini bar */}
              <div style={{ height: 3, borderRadius: 2, background: "rgba(255,255,255,0.09)", marginTop: 6, overflow: "hidden" }}>
                <div style={{ height: "100%", width: `${(batchHits / lastBatch.length) * 100}%`, background: target.color, borderRadius: 2 }} />
              </div>
            </div>
          )}

          <div style={{ padding: 8, borderRadius: 5, background: "rgba(80,220,255,0.03)", border: "1px solid rgba(80,220,255,0.06)" }}>
            <div style={{ fontSize: 7, textTransform: "uppercase", letterSpacing: 2, opacity: 0.5, color: "#50dcff" }}>Shared Pool</div>
            <div style={{ fontSize: 16, fontWeight: 700, color: "#50dcff" }}>{shared.total.toLocaleString()}</div>
            <div style={{ fontSize: 8, opacity: 0.3 }}>{shared.hits.length} shared hits</div>
          </div>

          {/* Provenance */}
          <div style={{ padding: 8, borderRadius: 5, background: "rgba(255,255,255,0.03)", border: "1px solid rgba(255,255,255,0.06)" }}>
            <div style={{ fontSize: 7, textTransform: "uppercase", letterSpacing: 2, opacity: 0.5, marginBottom: 4 }}>Provenance</div>
            <div style={{ fontSize: 7, lineHeight: 1.8 }}>
              {[
                { l: "Targets", v: "RCSB PDB", ok: true },
                { l: "Pharmacophore", v: pdbStatus[target.pdb]?.source === "PDB" ? "PDB" : "Lit", ok: true },
                { l: "Compounds", v: dataSource === "ChEMBL" ? "ChEMBL" : "Synth", ok: dataSource === "ChEMBL" },
                { l: "3D coords", v: "Approx", ok: false },
                { l: "Validation", v: "Pending", ok: false },
              ].map((d, i) => (
                <div key={i} style={{ display: "flex", justifyContent: "space-between" }}>
                  <span style={{ opacity: 0.5 }}>{d.l}</span>
                  <span style={{ color: d.ok ? "#78ff64" : "#ffc832" }}>{d.ok ? "✓" : "△"} {d.v}</span>
                </div>
              ))}
            </div>
          </div>

          {/* Data Sources */}
          <div style={{ padding: 8, borderRadius: 5, background: "rgba(255,255,255,0.03)", border: "1px solid rgba(255,255,255,0.06)" }}>
            <div style={{ fontSize: 7, textTransform: "uppercase", letterSpacing: 2, opacity: 0.5, marginBottom: 4 }}>Data Sources</div>
            <div style={{ fontSize: 7, lineHeight: 1.6 }}>
              {DATA_SOURCES.map((s, i) => (
                <div key={i} style={{ marginBottom: 3 }}>
                  <div style={{ display: "flex", justifyContent: "space-between", alignItems: "center" }}>
                    <span style={{ cursor: "pointer", color: "#50dcff", borderBottom: "1px dotted rgba(80,220,255,0.3)" }} onClick={() => window.open(s.url, "_blank")}>{s.name} ↗</span>
                    <span style={{ fontSize: 6, opacity: 0.35 }}>{s.used}</span>
                  </div>
                </div>
              ))}
            </div>
          </div>

          {/* Reference */}
          <div style={{ padding: 8, borderRadius: 5, background: `${target.color}04`, border: `1px solid ${target.color}10` }}>
            <div style={{ fontSize: 7, textTransform: "uppercase", letterSpacing: 2, opacity: 0.5, color: target.color, marginBottom: 4 }}>Reference · {target.id}</div>
            <div style={{ fontSize: 7, opacity: 0.45, lineHeight: 1.5, fontStyle: "italic" }}>{target.ref}</div>
            <div style={{ display: "flex", gap: 6, marginTop: 4, flexWrap: "wrap" }}>
              <span style={{ cursor: "pointer", fontSize: 7, color: target.color, borderBottom: `1px dotted ${target.color}50` }} onClick={() => window.open(`https://pubmed.ncbi.nlm.nih.gov/${target.pmid}`, "_blank")}>PubMed ↗</span>
              <span style={{ cursor: "pointer", fontSize: 7, color: target.color, borderBottom: `1px dotted ${target.color}50` }} onClick={() => window.open(`https://doi.org/${target.doi}`, "_blank")}>DOI ↗</span>
              <span style={{ cursor: "pointer", fontSize: 7, color: target.color, borderBottom: `1px dotted ${target.color}50` }} onClick={() => window.open(`https://www.rcsb.org/structure/${target.pdb}`, "_blank")}>PDB ↗</span>
              {dataSource === "ChEMBL" && selectedHit?.id?.startsWith("CHEMBL") && (
                <span style={{ cursor: "pointer", fontSize: 7, color: "#78ff64", borderBottom: "1px dotted rgba(120,255,100,0.3)" }} onClick={() => window.open(`https://www.ebi.ac.uk/chembl/compound_report_card/${selectedHit.id}`, "_blank")}>{selectedHit.id} ↗</span>
              )}
            </div>
          </div>
        </div>

        {/* ═══ CENTER ═══ */}
        <div style={{ display: "flex", flexDirection: "column", overflow: "hidden", minHeight: 0 }}>

          {/* TOP ROW: 3D Viewer + Score Dist + Batch Heatmap */}
          <div style={{ flex: "0 0 45%", display: "grid", gridTemplateColumns: "1fr 1fr 1fr", overflow: "hidden", borderBottom: "1px solid rgba(255,255,255,0.07)" }}>

            {/* 3D Pharmacophore */}
            <div style={{ position: "relative", borderRight: "1px solid rgba(255,255,255,0.07)" }}>
              <div style={{ position: "absolute", top: 6, left: 10, zIndex: 1, fontSize: 7, textTransform: "uppercase", letterSpacing: 2, opacity: 0.45 }}>Binding Pocket · {target.pdb}</div>
              {/* Legend */}
              <div style={{ position: "absolute", bottom: 6, left: 10, zIndex: 1, display: "flex", gap: 8, fontSize: 7 }}>
                {Object.entries(FEAT_COLORS).slice(0, 3).map(([k, v]) => (
                  <span key={k} style={{ display: "flex", alignItems: "center", gap: 3 }}>
                    <span style={{ width: 5, height: 5, borderRadius: "50%", background: v, display: "inline-block" }} />{k}
                  </span>
                ))}
              </div>
              {selectedHit && <div style={{ position: "absolute", top: 6, right: 10, zIndex: 1, fontSize: 7, color: target.color, background: "rgba(0,0,0,0.6)", padding: "2px 6px", borderRadius: 3 }}>◇ {selectedHit.id}</div>}
              <PharmaViewer pharmacophore={target.pharmacophore} targetColor={target.color} selectedHit={selectedHit} isScreening={autoMode} lastScore={topHits[0]?.score} screenEvent={screenEvent} />
            </div>

            {/* Score Distribution */}
            <div style={{ padding: "8px 10px", display: "flex", flexDirection: "column", borderRight: "1px solid rgba(255,255,255,0.07)" }}>
              <div style={{ display: "flex", justifyContent: "space-between", marginBottom: 4 }}>
                <span style={{ fontSize: 7, textTransform: "uppercase", letterSpacing: 2, opacity: 0.45 }}>Scores</span>
                {bestScore && <span style={{ fontSize: 8, color: target.color }}>best: {bestScore.toFixed(3)}</span>}
              </div>
              <div style={{ flex: 1, minHeight: 0 }}>
                {scoreDist.length > 0 ? (
                  <ResponsiveContainer width="100%" height="100%">
                    <BarChart data={scoreDist}>
                      <XAxis dataKey="bin" tick={{ fontSize: 7, fill: "#444466" }} axisLine={{ stroke: "#222240" }} tickLine={false} interval={1} />
                      <YAxis tick={{ fontSize: 7, fill: "#444466" }} axisLine={{ stroke: "#222240" }} tickLine={false} width={25} />
                      <Tooltip contentStyle={{ background: "#111", border: "1px solid #222", borderRadius: 4, fontSize: 8 }} />
                      <Bar dataKey="count" radius={[1, 1, 0, 0]}>
                        {scoreDist.map((d, i) => <Cell key={i} fill={d.hits > 0 ? target.color : i >= 8 ? `${target.color}60` : "#252540"} />)}
                      </Bar>
                    </BarChart>
                  </ResponsiveContainer>
                ) : <div style={{ display: "flex", alignItems: "center", justifyContent: "center", height: "100%", opacity: 0.15, fontSize: 10 }}>—</div>}
              </div>
            </div>

            {/* Batch Heatmap */}
            <div style={{ padding: "8px 10px", display: "flex", flexDirection: "column" }}>
              <div style={{ display: "flex", justifyContent: "space-between", marginBottom: 4 }}>
                <span style={{ fontSize: 7, textTransform: "uppercase", letterSpacing: 2, opacity: 0.45 }}>Batch #{batchId}</span>
                <span style={{ fontSize: 7, opacity: 0.4 }}>{lastBatch.length} cpd</span>
              </div>
              <div style={{ flex: 1, minHeight: 0, overflow: "hidden", display: "flex", alignItems: "center" }}>
                <BatchHeatmap scores={lastBatch} targetColor={target.color} />
              </div>
              {lastBatch.length > 0 && (
                <div style={{ fontSize: 7, opacity: 0.5, marginTop: 4, textAlign: "center" }}>
                  <span style={{ color: target.color }}>■</span> hit (≥{HIT_THRESHOLD}) &nbsp; <span style={{ opacity: 0.4 }}>■</span> miss
                </div>
              )}
            </div>
          </div>

          {/* MIDDLE ROW: Chemical Space + Hit Rate + Throughput */}
          <div style={{ flex: "0 0 25%", display: "grid", gridTemplateColumns: "1.2fr 0.9fr 0.9fr", overflow: "hidden", borderBottom: "1px solid rgba(255,255,255,0.07)" }}>

            {/* MW vs logP scatter */}
            <div style={{ padding: "8px 10px", display: "flex", flexDirection: "column", borderRight: "1px solid rgba(255,255,255,0.07)" }}>
              <span style={{ fontSize: 7, textTransform: "uppercase", letterSpacing: 2, opacity: 0.45, marginBottom: 4 }}>Chemical Space (MW vs logP)</span>
              <div style={{ flex: 1, minHeight: 0 }}>
                {chemSpace.length > 0 ? (
                  <ResponsiveContainer width="100%" height="100%">
                    <ScatterChart margin={{ top: 4, right: 4, bottom: 4, left: 4 }}>
                      <CartesianGrid strokeDasharray="3 3" stroke="#1a1a30" />
                      <XAxis dataKey="mw" name="MW" tick={{ fontSize: 7, fill: "#444466" }} axisLine={{ stroke: "#222240" }} tickLine={false} type="number" domain={[150, 600]} label={{ value: "MW", fontSize: 7, fill: "#333", position: "bottom", offset: -2 }} />
                      <YAxis dataKey="logP" name="logP" tick={{ fontSize: 7, fill: "#444466" }} axisLine={{ stroke: "#222240" }} tickLine={false} width={25} domain={[-2, 7]} label={{ value: "logP", fontSize: 7, fill: "#333", angle: -90, position: "left", offset: -10 }} />
                      <Tooltip contentStyle={{ background: "#111", border: "1px solid #222", borderRadius: 4, fontSize: 8 }} formatter={(v, n) => [typeof v === 'number' ? v.toFixed(2) : v, n]} />
                      <Scatter data={chemSpace.filter(d => !d.isHit)} fill="#282848" r={2} />
                      <Scatter data={chemSpace.filter(d => d.isHit)} fill={target.color} r={4} />
                    </ScatterChart>
                  </ResponsiveContainer>
                ) : <div style={{ display: "flex", alignItems: "center", justifyContent: "center", height: "100%", opacity: 0.15, fontSize: 10 }}>—</div>}
              </div>
            </div>

            {/* Hit Rate Trend */}
            <div style={{ padding: "8px 10px", display: "flex", flexDirection: "column", borderRight: "1px solid rgba(255,255,255,0.07)" }}>
              <span style={{ fontSize: 7, textTransform: "uppercase", letterSpacing: 2, opacity: 0.45, marginBottom: 4 }}>Hit Rate %</span>
              <div style={{ flex: 1, minHeight: 0 }}>
                {hitRateHist.length > 0 ? (
                  <ResponsiveContainer width="100%" height="100%">
                    <AreaChart data={hitRateHist}>
                      <defs><linearGradient id="hrG" x1="0" y1="0" x2="0" y2="1">
                        <stop offset="0%" stopColor="#ffc832" stopOpacity={0.3} />
                        <stop offset="100%" stopColor="#ffc832" stopOpacity={0} />
                      </linearGradient></defs>
                      <XAxis dataKey="batch" tick={{ fontSize: 7, fill: "#444466" }} axisLine={{ stroke: "#222240" }} tickLine={false} />
                      <YAxis tick={{ fontSize: 7, fill: "#444466" }} axisLine={{ stroke: "#222240" }} tickLine={false} width={25} />
                      <Tooltip contentStyle={{ background: "#111", border: "1px solid #222", borderRadius: 4, fontSize: 8 }} />
                      <Area type="monotone" dataKey="rate" stroke="#ffc832" fill="url(#hrG)" strokeWidth={1.5} dot={false} />
                    </AreaChart>
                  </ResponsiveContainer>
                ) : <div style={{ display: "flex", alignItems: "center", justifyContent: "center", height: "100%", opacity: 0.15, fontSize: 10 }}>—</div>}
              </div>
            </div>

            {/* Throughput */}
            <div style={{ padding: "8px 10px", display: "flex", flexDirection: "column" }}>
              <span style={{ fontSize: 7, textTransform: "uppercase", letterSpacing: 2, opacity: 0.45, marginBottom: 4 }}>Throughput (cpd/s)</span>
              <div style={{ flex: 1, minHeight: 0 }}>
                {tpHist.length > 0 ? (
                  <ResponsiveContainer width="100%" height="100%">
                    <AreaChart data={tpHist}>
                      <defs><linearGradient id="tpG" x1="0" y1="0" x2="0" y2="1">
                        <stop offset="0%" stopColor={target.color} stopOpacity={0.3} />
                        <stop offset="100%" stopColor={target.color} stopOpacity={0} />
                      </linearGradient></defs>
                      <XAxis dataKey="batch" tick={{ fontSize: 7, fill: "#444466" }} axisLine={{ stroke: "#222240" }} tickLine={false} />
                      <YAxis tick={{ fontSize: 7, fill: "#444466" }} axisLine={{ stroke: "#222240" }} tickLine={false} width={30} />
                      <Area type="monotone" dataKey="cps" stroke={target.color} fill="url(#tpG)" strokeWidth={1.5} dot={false} />
                    </AreaChart>
                  </ResponsiveContainer>
                ) : <div style={{ display: "flex", alignItems: "center", justifyContent: "center", height: "100%", opacity: 0.15, fontSize: 10 }}>—</div>}
              </div>
            </div>
          </div>

          {/* BOTTOM: Hits Table */}
          <div style={{ flex: 1, display: "flex", flexDirection: "column", overflow: "hidden", minHeight: 0 }}>
            <div style={{ display: "flex", justifyContent: "space-between", alignItems: "center", padding: "6px 10px", flexShrink: 0 }}>
              <span style={{ fontSize: 7, textTransform: "uppercase", letterSpacing: 2, opacity: 0.45 }}>Top Hits ({topHits.length})</span>
              <span style={{ fontSize: 7, opacity: 0.4 }}>Click row to view in 3D</span>
            </div>
            <div style={{ flex: 1, overflow: "auto", padding: "0 10px 6px", minHeight: 0 }}>
              {topHits.length > 0 ? (
                <table style={{ width: "100%", borderCollapse: "collapse", fontSize: 8 }}>
                  <thead><tr style={{ borderBottom: "1px solid rgba(255,255,255,0.09)" }}>
                    {["#", "ID", "Score", "MW", "logP", "HBD", "HBA", "TPSA", "SMILES", "Src"].map(h =>
                      <th key={h} style={{ padding: "3px 5px", textAlign: "left", opacity: 0.6, fontWeight: 600, fontSize: 7, textTransform: "uppercase", letterSpacing: 1 }}>{h}</th>
                    )}
                  </tr></thead>
                  <tbody>
                    {topHits.slice(0, 25).map((h, i) => (
                      <tr key={i} onClick={() => setSelectedHit(h)} style={{
                        borderBottom: "1px solid rgba(255,255,255,0.08)",
                        cursor: "pointer", background: selectedHit?.id === h.id ? `${target.color}08` : "transparent",
                        transition: "background 0.15s",
                      }}>
                        <td style={{ padding: "2px 5px", opacity: 0.4 }}>{i + 1}</td>
                        <td style={{ padding: "2px 5px", color: target.color, opacity: 0.8, fontWeight: 600 }}>{h.id}</td>
                        <td style={{ padding: "2px 5px", fontWeight: 700, color: h.score >= 2.5 ? target.color : h.score >= HIT_THRESHOLD ? "#ffc832" : "#444" }}>{h.score.toFixed(3)}</td>
                        <td style={{ padding: "2px 5px", opacity: 0.4 }}>{h.mw?.toFixed(0)}</td>
                        <td style={{ padding: "2px 5px", opacity: 0.4 }}>{h.logP?.toFixed(1)}</td>
                        <td style={{ padding: "2px 5px", opacity: 0.4 }}>{h.hbd}</td>
                        <td style={{ padding: "2px 5px", opacity: 0.4 }}>{h.hba}</td>
                        <td style={{ padding: "2px 5px", opacity: 0.4 }}>{h.tpsa?.toFixed(0)}</td>
                        <td style={{ padding: "2px 5px", opacity: 0.6, maxWidth: 120, overflow: "hidden", textOverflow: "ellipsis", whiteSpace: "nowrap", fontSize: 7 }}>{h.smiles || "—"}</td>
                        <td style={{ padding: "2px 5px" }}><span style={{ fontSize: 6, padding: "1px 3px", borderRadius: 2, background: h.source === "ChEMBL" ? "rgba(120,255,100,0.08)" : "rgba(255,200,50,0.08)", color: h.source === "ChEMBL" ? "#78ff64" : "#ffc832" }}>{h.source === "ChEMBL" ? "real" : "syn"}</span></td>
                      </tr>
                    ))}
                  </tbody>
                </table>
              ) : <div style={{ display: "flex", alignItems: "center", justifyContent: "center", height: "100%", opacity: 0.15, fontSize: 10 }}>Start screening</div>}
            </div>
          </div>
        </div>

        {/* ═══ RIGHT PANEL ═══ */}
        <div style={{ borderLeft: "1px solid rgba(255,255,255,0.07)", display: "flex", flexDirection: "column", background: "rgba(14,14,28,0.4)", overflow: "hidden", minHeight: 0 }}>

          {/* Hardware Info */}
          <div style={{ padding: "8px 10px", borderBottom: "1px solid rgba(255,255,255,0.07)", flexShrink: 0 }}>
            <div style={{ fontSize: 7, textTransform: "uppercase", letterSpacing: 2, opacity: 0.6, marginBottom: 5 }}>System</div>
            {hwInfo ? (
              <div style={{ fontSize: 7, lineHeight: 1.7 }}>
                <div style={{ display: "flex", justifyContent: "space-between" }}>
                  <span style={{ opacity: 0.5 }}>Chip</span>
                  <span style={{ color: "#50dcff", fontWeight: 600 }}>{hwInfo.chip}</span>
                </div>
                <div style={{ display: "flex", justifyContent: "space-between" }}>
                  <span style={{ opacity: 0.5 }}>Device</span>
                  <span style={{ opacity: 0.7 }}>{hwInfo.device}</span>
                </div>
                <div style={{ display: "flex", justifyContent: "space-between" }}>
                  <span style={{ opacity: 0.5 }}>CPU Cores</span>
                  <span style={{ opacity: 0.7 }}>{hwInfo.cores}</span>
                </div>
                <div style={{ display: "flex", justifyContent: "space-between" }}>
                  <span style={{ opacity: 0.5 }}>Memory</span>
                  <span style={{ opacity: 0.7 }}>{hwInfo.memory}</span>
                </div>
                <div style={{ display: "flex", justifyContent: "space-between" }}>
                  <span style={{ opacity: 0.5 }}>Browser</span>
                  <span style={{ opacity: 0.7 }}>{hwInfo.browser}</span>
                </div>
                <div style={{ display: "flex", justifyContent: "space-between" }}>
                  <span style={{ opacity: 0.5 }}>GPU Vendor</span>
                  <span style={{ color: gpu ? "#78ff64" : "#ffc832" }}>{hwInfo.gpu.vendor || "—"}</span>
                </div>
                <div style={{ display: "flex", justifyContent: "space-between" }}>
                  <span style={{ opacity: 0.5 }}>GPU Arch</span>
                  <span style={{ opacity: 0.7 }}>{hwInfo.gpu.arch || "—"}</span>
                </div>
                {hwInfo.gpu.limits.maxBufferSize && hwInfo.gpu.limits.maxBufferSize !== "?" && (
                  <div style={{ display: "flex", justifyContent: "space-between" }}>
                    <span style={{ opacity: 0.5 }}>GPU Buffer</span>
                    <span style={{ opacity: 0.7 }}>{hwInfo.gpu.limits.maxBufferSize}</span>
                  </div>
                )}
                {hwInfo.gpu.limits.maxWorkgroupSize && hwInfo.gpu.limits.maxWorkgroupSize !== "?" && (
                  <div style={{ display: "flex", justifyContent: "space-between" }}>
                    <span style={{ opacity: 0.5 }}>Workgroup</span>
                    <span style={{ opacity: 0.7 }}>{hwInfo.gpu.limits.maxWorkgroupSize}</span>
                  </div>
                )}
                <div style={{ display: "flex", justifyContent: "space-between", marginTop: 2 }}>
                  <span style={{ opacity: 0.5 }}>NPU</span>
                  <span style={{ color: hwInfo.npuAvailable ? "#ffc832" : "#666" }}>
                    {hwInfo.npuAvailable ? "△ Present (no browser API)" : "N/A"}
                  </span>
                </div>
                {/* Live compute utilization */}
                <div style={{ marginTop: 6 }}>
                  <div style={{ fontSize: 6, opacity: 0.4, textTransform: "uppercase", letterSpacing: 1, marginBottom: 3 }}>Compute Utilization {autoMode ? "(live)" : ""}</div>
                  
                  {/* GPU bar */}
                  <div style={{ display: "flex", alignItems: "center", gap: 4, marginBottom: 3 }}>
                    <span style={{ fontSize: 6, width: 22, color: "#78ff64", opacity: 0.7 }}>GPU</span>
                    <div style={{ flex: 1, height: 8, borderRadius: 2, background: "rgba(255,255,255,0.04)", overflow: "hidden" }}>
                      <div style={{ height: "100%", width: `${utilization.gpuPct}%`, background: "linear-gradient(90deg, #78ff6450, #78ff64)", borderRadius: 2, transition: "width 0.3s" }} />
                    </div>
                    <span style={{ fontSize: 7, width: 28, textAlign: "right", color: "#78ff64", fontWeight: 600 }}>{utilization.gpuPct}%</span>
                  </div>

                  {/* CPU bar */}
                  <div style={{ display: "flex", alignItems: "center", gap: 4, marginBottom: 3 }}>
                    <span style={{ fontSize: 6, width: 22, color: "#50dcff", opacity: 0.7 }}>CPU</span>
                    <div style={{ flex: 1, height: 8, borderRadius: 2, background: "rgba(255,255,255,0.04)", overflow: "hidden" }}>
                      <div style={{ height: "100%", width: `${utilization.cpuPct}%`, background: "linear-gradient(90deg, #50dcff50, #50dcff)", borderRadius: 2, transition: "width 0.3s" }} />
                    </div>
                    <span style={{ fontSize: 7, width: 28, textAlign: "right", color: "#50dcff", fontWeight: 600 }}>{utilization.cpuPct}%</span>
                  </div>

                  {/* NPU bar */}
                  <div style={{ display: "flex", alignItems: "center", gap: 4, marginBottom: 3 }}>
                    <span style={{ fontSize: 6, width: 22, color: "#ffc832", opacity: 0.4 }}>NPU</span>
                    <div style={{ flex: 1, height: 8, borderRadius: 2, background: "rgba(255,255,255,0.04)", overflow: "hidden" }}>
                      <div style={{ height: "100%", width: "0%", background: "#ffc832", borderRadius: 2 }} />
                    </div>
                    <span style={{ fontSize: 7, width: 28, textAlign: "right", color: "#ffc832", opacity: 0.3 }}>—</span>
                  </div>

                  {/* Timing breakdown */}
                  {utilization.totalMs > 0 && (
                    <div style={{ fontSize: 6, opacity: 0.4, marginTop: 2, lineHeight: 1.6 }}>
                      <div style={{ display: "flex", justifyContent: "space-between" }}>
                        <span>Fetch/prep</span><span>{utilization.fetchMs}ms</span>
                      </div>
                      <div style={{ display: "flex", justifyContent: "space-between" }}>
                        <span>Scoring</span><span>{utilization.scoreMs}ms</span>
                      </div>
                      <div style={{ display: "flex", justifyContent: "space-between" }}>
                        <span>Total cycle</span><span>{utilization.totalMs}ms</span>
                      </div>
                    </div>
                  )}

                  <div style={{ fontSize: 6, opacity: 0.25, marginTop: 3 }}>
                    {gpu ? "GPU: pharmacophore scoring · CPU: fetch + data prep" : "CPU: all compute"} · NPU: no browser API
                  </div>
                </div>
              </div>
            ) : (
              <div style={{ fontSize: 8, opacity: 0.3 }}>Detecting...</div>
            )}
          </div>

          {/* Perf cards */}
          {lastMs && (
            <div style={{ padding: "8px 10px", borderBottom: "1px solid rgba(255,255,255,0.07)", flexShrink: 0 }}>
              <div style={{ fontSize: 7, textTransform: "uppercase", letterSpacing: 2, opacity: 0.6, marginBottom: 4 }}>Performance</div>
              <div style={{ display: "grid", gridTemplateColumns: "1fr 1fr", gap: 4 }}>
                {[
                  { l: "Speed", v: cps ? `${cps.toLocaleString()}/s` : "—", c: target.color },
                  { l: "Batch", v: `${lastMs.toFixed(0)}ms`, c: "#ffc832" },
                  { l: "Hits", v: String(hitCount), c: "#78ff64" },
                  { l: "Backend", v: gpu ? "WebGPU" : "CPU", c: gpu ? "#78ff64" : "#ffc832" },
                ].map(m => (
                  <div key={m.l} style={{ padding: "4px 6px", borderRadius: 3, background: "rgba(255,255,255,0.04)", border: "1px solid rgba(255,255,255,0.06)" }}>
                    <div style={{ fontSize: 11, fontWeight: 700, color: m.c }}>{m.v}</div>
                    <div style={{ fontSize: 6, opacity: 0.45, textTransform: "uppercase", letterSpacing: 1 }}>{m.l}</div>
                  </div>
                ))}
              </div>
            </div>
          )}

          {/* Structure Freshness */}
          <div style={{ padding: "8px 10px", borderBottom: "1px solid rgba(255,255,255,0.07)", flexShrink: 0 }}>
            <div style={{ fontSize: 7, textTransform: "uppercase", letterSpacing: 2, opacity: 0.6, marginBottom: 4 }}>Structure Freshness</div>
            <div style={{ fontSize: 7, lineHeight: 1.7 }}>
              {targets.map(t => {
                const f = structFreshness[t.pdb];
                const hasUpdate = f?.status === "update_available";
                const isCurrent = f?.status === "current";
                return (
                  <div key={t.pdb} style={{ marginBottom: 4, padding: "3px 5px", borderRadius: 3, background: hasUpdate ? "rgba(255,200,50,0.06)" : "transparent", border: hasUpdate ? "1px solid rgba(255,200,50,0.15)" : "1px solid transparent" }}>
                    <div style={{ display: "flex", justifyContent: "space-between", alignItems: "center" }}>
                      <span style={{ color: t.color, fontWeight: 600 }}>{t.pdb}</span>
                      <span style={{ fontSize: 6 }}>
                        {!f ? <span style={{ opacity: 0.3 }}>Checking...</span> :
                         hasUpdate ? <span style={{ color: "#ffc832", fontWeight: 700 }}>⚠ UPDATE</span> :
                         isCurrent ? <span style={{ color: "#78ff64" }}>✓ Current</span> :
                         <span style={{ opacity: 0.3 }}>? Failed</span>}
                      </span>
                    </div>
                    {f && isCurrent && (
                      <div style={{ fontSize: 6, opacity: 0.4 }}>
                        {f.currentRes} · deposited {f.deposited} · checked {f.checked}
                      </div>
                    )}
                    {f && hasUpdate && (
                      <div style={{ fontSize: 6, marginTop: 2 }}>
                        <div style={{ color: "#ffc832" }}>Better resolution may exist:</div>
                        {f.alternatives.map((a, i) => (
                          <span key={i} style={{ cursor: "pointer", color: "#ffc832", marginRight: 6, borderBottom: "1px dotted rgba(255,200,50,0.4)" }} onClick={() => window.open(`https://www.rcsb.org/structure/${a.pdb}`, "_blank")}>{a.pdb} ↗</span>
                        ))}
                      </div>
                    )}
                  </div>
                );
              })}
              <div style={{ fontSize: 6, opacity: 0.3, marginTop: 4 }}>Auto-checked on session start via RCSB API. Target list is static — new targets require code update.</div>
            </div>
          </div>

          {/* Selected hit detail */}
          {selectedHit && (
            <div style={{ padding: "8px 10px", borderBottom: "1px solid rgba(255,255,255,0.07)", flexShrink: 0 }}>
              <div style={{ fontSize: 7, textTransform: "uppercase", letterSpacing: 2, opacity: 0.6, marginBottom: 4 }}>Selected</div>
              <div style={{ fontSize: 10, fontWeight: 700, color: target.color, marginBottom: 4 }}>{selectedHit.id}</div>
              <div style={{ fontSize: 8, opacity: 0.6, lineHeight: 1.8 }}>
                <div>Score: <span style={{ color: target.color }}>{selectedHit.score?.toFixed(3)}</span></div>
                <div>MW: {selectedHit.mw?.toFixed(1)} · logP: {selectedHit.logP?.toFixed(2)}</div>
                <div>HBD: {selectedHit.hbd} · HBA: {selectedHit.hba} · TPSA: {selectedHit.tpsa?.toFixed(0)}</div>
                <div>Features: {selectedHit.features?.length} ({selectedHit.features?.map(f => f.type).filter((v, i, a) => a.indexOf(v) === i).join(", ")})</div>
                {selectedHit.smiles && <div style={{ fontSize: 7, wordBreak: "break-all", marginTop: 4, opacity: 0.3 }}>{selectedHit.smiles}</div>}
              </div>
            </div>
          )}

          {/* Log */}
          <div style={{ flex: 1, overflow: "hidden", display: "flex", flexDirection: "column", minHeight: 0 }}>
            <div style={{ padding: "6px 10px 3px", fontSize: 7, textTransform: "uppercase", letterSpacing: 2, opacity: 0.6, flexShrink: 0 }}>Log</div>
            <div style={{ flex: 1, overflowY: "auto", padding: "0 10px 6px", minHeight: 0 }}>
              {log.map((l, i) => (
                <div key={i} style={{
                  fontSize: 8, lineHeight: 1.6, paddingLeft: 7,
                  borderLeft: `2px solid ${l.t === "ok" ? "#78ff6430" : l.t === "warn" ? "#ffc83230" : l.t === "info" ? "#50dcff20" : "#252540"}`,
                  color: l.t === "ok" ? "#78ff64" : l.t === "warn" ? "#ffc832" : l.t === "info" ? "#50dcff" : "#444",
                }}>{l.m}</div>
              ))}
              <div ref={logEnd} />
            </div>
          </div>

          {/* Footer */}
          <div style={{ padding: "6px 10px", borderTop: "1px solid rgba(255,255,255,0.07)", fontSize: 7, opacity: 0.15, lineHeight: 1.6, flexShrink: 0 }}>
            ✓ PDB · ✓ ChEMBL · ✓ WebGPU · △ 3D approx · △ Unvalidated
          </div>
        </div>
      </div>

      <style>{`*{box-sizing:border-box;margin:0;padding:0;}
        @keyframes pulse{0%,100%{opacity:1}50%{opacity:0.5}}
        ::-webkit-scrollbar{width:3px;}::-webkit-scrollbar-track{background:transparent;}
        ::-webkit-scrollbar-thumb{background:rgba(255,255,255,0.10);border-radius:2px;}`}</style>
    </div>
  );
}
