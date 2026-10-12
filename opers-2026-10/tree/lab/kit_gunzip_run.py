"""KIT HELPER (new in v1.2.5, not an archived file): decompress an archived record that the kit ships gzipped (to keep every
kit file under 50 MB), verify it against the archived SHA-256 of the uncompressed file, write it at its archived path
inside the temporary copy of tree/, then run the archived script unchanged.
Usage (cwd = the archived script's directory): kit_gunzip_run.py <record.json.gz> <sha256 of record.json> <script.py> [args...]"""
import gzip, hashlib, os, runpy, shutil, sys

gz, want, script = sys.argv[1], sys.argv[2], sys.argv[3]
out = gz[:-3]
with gzip.open(gz, 'rb') as f, open(out, 'wb') as g:
    shutil.copyfileobj(f, g, 1 << 20)
h = hashlib.sha256()
with open(out, 'rb') as f:
    for b in iter(lambda: f.read(1 << 20), b''):
        h.update(b)
if h.hexdigest() != want:
    print(f'KIT RECORD HASH MISMATCH {out}: {h.hexdigest()} != {want}')
    sys.exit(3)
sys.argv = [script] + sys.argv[4:]
runpy.run_path(script, run_name='__main__')
