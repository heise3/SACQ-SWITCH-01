"""Download documented public inputs into this isolated checkout; no analysis is run."""
from pathlib import Path
import json,hashlib,urllib.request
ROOT=Path(__file__).resolve().parents[1]
for r in json.loads((ROOT/'documentation/EXTERNAL_INPUTS.json').read_text()):
 p=ROOT/r['directory']/r['name'];p.parent.mkdir(parents=True,exist_ok=True)
 if not p.exists():
  tmp=p.with_suffix(p.suffix+'.partial');print('Downloading',r['name'],flush=True)
  req=urllib.request.Request(r['url'],headers={'User-Agent':'SACQ-SWITCH-01 reproducibility/1.3.2'})
  with urllib.request.urlopen(req,timeout=180) as src,tmp.open('wb') as dst:
   while b:=src.read(8*1024*1024):dst.write(b)
  tmp.replace(p)
 h=hashlib.sha256()
 with p.open('rb') as f:
  while b:=f.read(8*1024*1024):h.update(b)
 if h.hexdigest()!=r['sha256']:raise RuntimeError('Input hash differs from analyzed release: '+r['name'])
 print('Verified',r['name'],flush=True)
