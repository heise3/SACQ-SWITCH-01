from pathlib import Path
import hashlib,json
root=Path(__file__).resolve().parents[1]
m=json.loads((root/'SHA256SUMS.json').read_text())
for name,h in m.items():
 p=root/name
 assert p.exists() and hashlib.sha256(p.read_bytes()).hexdigest()==h,name
print('PASS:',len(m),'release files match SHA-256 manifest')
