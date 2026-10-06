#!/usr/bin/env python3
import hashlib, json, subprocess, sys
from pathlib import Path
ROOT=Path(__file__).resolve().parent
manifest=json.loads((ROOT/'MANIFEST.json').read_text(encoding='utf-8'))
errors=[]
for rel,meta in manifest['files'].items():
    p=ROOT/rel
    if not p.is_file(): errors.append(f'missing: {rel}'); continue
    if p.stat().st_size!=meta['bytes']: errors.append(f'size: {rel}')
    if hashlib.sha256(p.read_bytes()).hexdigest()!=meta['sha256']: errors.append(f'hash: {rel}')
known=set(manifest['files'])|{'MANIFEST.json'}
actual={p.relative_to(ROOT).as_posix() for p in ROOT.rglob('*') if p.is_file() and '.git' not in p.parts and 'dist' not in p.parts}
for rel in sorted(actual-known): errors.append(f'unlisted: {rel}')
lock_path=ROOT/'awf/workflows/gromacs-organic-blend-zc/releases/1.0.0-rc.6/versions.lock'
if lock_path.is_file():
    lock=json.loads(lock_path.read_text())
    if lock.get('workflow_release')!='1.0.0-rc.6': errors.append('wrong workflow release')
else: errors.append('missing workflow lock')
for p in ROOT.rglob('*'):
    if p.is_dir() and p.name=='__pycache__': errors.append(f'cache: {p.relative_to(ROOT)}')
    if p.is_file() and p.suffix.lower() in {'.pyc','.log','.gjf'}: errors.append(f'excluded file type: {p.relative_to(ROOT)}')
if errors:
    print('FAIL'); print(chr(10).join(errors)); sys.exit(1)
tool=ROOT/'awf/components/tools/opv-library/releases/1.0.0-rc.2/source/opv_library.py'
library=ROOT/'awf/components/artifacts/opv-fragment-library/releases/1.0.0-rc.2/source/opv_fragment_library'
result=subprocess.run([sys.executable,str(tool),'--library',str(library),'verify'],text=True,capture_output=True)
if result.returncode:
    print('FAIL: library verification'); print(result.stdout); print(result.stderr); sys.exit(result.returncode)
report=json.loads(result.stdout)
if (report.get('fragment_records_checked'),report.get('molecules_checked'),report.get('checksums_checked'),report.get('status')) != (91,10,218,'PASS'):
    print('FAIL: unexpected library inventory'); print(result.stdout); sys.exit(1)
print(f"PASS: {len(manifest['files'])} files; opv-md {manifest['workflow_release']}; 91 fragments; 10 molecules; 218 checksums")
