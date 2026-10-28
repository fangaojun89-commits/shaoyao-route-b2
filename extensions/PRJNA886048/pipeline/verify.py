"""Verify exact payload bytes and independently check all retained statistics."""
from pathlib import Path
import argparse, hashlib, json, subprocess
p=argparse.ArgumentParser()
p.add_argument('--rscript',default='Rscript')
a=p.parse_args()
root=Path(__file__).resolve().parents[1]
m=json.loads((root/'MANIFEST.json').read_text(encoding='utf-8'))
for name,pin in m['files'].items():
    f=root/name
    if not f.is_file() or f.stat().st_size!=pin['bytes'] or hashlib.sha256(f.read_bytes()).hexdigest()!=pin['sha256']:
        raise SystemExit('Payload identity mismatch: '+name)
r=subprocess.run([a.rscript,'--vanilla',str(root/'pipeline/verify_tables.R'),str(root)],check=True,capture_output=True,text=True)
print(r.stdout,end='')
print(json.dumps({'status':'PASS_PAYLOAD_AND_CONTENT_READBACK','payloads':len(m['files']),'models_fitted':0,'raw_files_modified':0}))
