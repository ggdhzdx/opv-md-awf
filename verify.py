#!/usr/bin/env python3
import ast, hashlib, json, re, subprocess, sys, tomllib
from pathlib import Path
ROOT=Path(__file__).resolve().parent; errors=[]
manifest=json.loads((ROOT/'MANIFEST.json').read_text(encoding='utf-8'))
for rel,meta in manifest['files'].items():
 p=ROOT/rel
 if not p.is_file(): errors.append(f'missing: {rel}'); continue
 data=p.read_bytes()
 if len(data)!=meta['bytes']: errors.append(f'size: {rel}')
 if hashlib.sha256(data).hexdigest()!=meta['sha256']: errors.append(f'hash: {rel}')
known=set(manifest['files'])|{'MANIFEST.json'}
actual={p.relative_to(ROOT).as_posix() for p in ROOT.rglob('*') if p.is_file() and '.git' not in p.parts and 'dist' not in p.parts}
for rel in sorted(actual-known): errors.append(f'unlisted: {rel}')
lock_path=ROOT/'awf/workflows/gromacs-organic-blend-zc/releases/1.0.0-rc.7/versions.lock'
expected={'artifact.opv-fragment-library':'1.0.0-rc.2','recipe.gromacs-organic-blend-production':'1.0.0-rc.2','recipe.opv-morphology-analysis':'1.0.0-rc.2','skill.analyze-opv-morphology':'1.0.0-rc.2','skill.run-gromacs-organic-blend':'1.0.0-rc.1','tool.opv-library':'1.0.0-rc.2'}
if lock_path.is_file():
 lock=json.loads(lock_path.read_text())
 if lock.get('workflow_release')!='1.0.0-rc.7': errors.append('wrong workflow release')
 if lock.get('components')!=expected: errors.append('wrong component lock')
else: errors.append('missing workflow lock')
analysis=ROOT/'awf/components/recipes/opv-morphology-analysis/releases/1.0.0-rc.2/source'
for rel in ('opv_morphology/bulkgraph_cli.py','opv_morphology/phase_separation.py','pyproject.toml'):
 if not (analysis/rel).is_file(): errors.append(f'missing analysis entry point: {rel}')
for p in ROOT.rglob('*'):
 if '.git' in p.parts or 'dist' in p.parts: continue
 if p.is_dir() and p.name=='__pycache__': errors.append(f'cache: {p.relative_to(ROOT)}')
 if not p.is_file(): continue
 if p.suffix.lower() in {'.pyc','.log','.gjf'}: errors.append(f'excluded file type: {p.relative_to(ROOT)}')
 if p.suffix=='.py':
  try: ast.parse(p.read_text(encoding='utf-8'),filename=str(p))
  except Exception as exc: errors.append(f'python syntax: {p.relative_to(ROOT)}: {exc}')
 if p.suffix=='.toml':
  try: tomllib.loads(p.read_text(encoding='utf-8'))
  except Exception as exc: errors.append(f'toml: {p.relative_to(ROOT)}: {exc}')
 if p.suffix.lower() in {'.json','.yaml'}:
  try: json.loads(p.read_text(encoding='utf-8'))
  except Exception as exc: errors.append(f'json: {p.relative_to(ROOT)}: {exc}')
 if p.name!='verify.py' and p.stat().st_size<5_000_000:
  try: text=p.read_text(encoding='utf-8')
  except UnicodeDecodeError: text=''
  if re.search(r'(/home/zc|/mnt/z|\\192\.168|gaowei|mmdata|20260920|D18_A-|qq\.com|Zc@zst|BEGIN (?:RSA |OPENSSH |EC )?PRIVATE KEY)',text,re.I): errors.append(f'private content: {p.relative_to(ROOT)}')
if errors: print('FAIL'); print(chr(10).join(errors)); sys.exit(1)
tool=ROOT/'awf/components/tools/opv-library/releases/1.0.0-rc.2/source/opv_library.py'
library=ROOT/'awf/components/artifacts/opv-fragment-library/releases/1.0.0-rc.2/source/opv_fragment_library'
result=subprocess.run([sys.executable,str(tool),'--library',str(library),'verify'],text=True,capture_output=True)
if result.returncode: print('FAIL: library verification'); print(result.stdout); print(result.stderr); sys.exit(result.returncode)
report=json.loads(result.stdout)
if (report.get('fragment_records_checked'),report.get('molecules_checked'),report.get('checksums_checked'),report.get('status'))!=(91,10,218,'PASS'): print('FAIL: unexpected library inventory'); print(result.stdout); sys.exit(1)
print(f"PASS: {len(manifest['files'])} files; opv-md {manifest['workflow_release']}; executable analysis; 91 fragments; 10 molecules; 218 checksums")
