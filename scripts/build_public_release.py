"""Build a source-only ZIP from the reviewed release-files.txt allowlist."""
import argparse
import hashlib
import json
from pathlib import Path, PurePosixPath
import re
import zipfile


def sources(root):
    root = Path(root).resolve()
    entries = (root/'release-files.txt').read_text(encoding='utf-8').splitlines()
    selected = {}
    for entry in entries:
        if not entry or entry.startswith('#'):
            continue
        rel = PurePosixPath(entry)
        if rel.is_absolute() or '..' in rel.parts or '\\' in entry or ':' in entry:
            raise ValueError('Release paths must be relative POSIX paths')
        path = root.joinpath(*rel.parts)
        if path.is_symlink() or not path.resolve().is_relative_to(root):
            raise ValueError('Release sources must stay inside the project')
        if entry in selected or not path.is_file():
            raise ValueError('Duplicate or missing release source: '+entry)
        if path.suffix not in ('.md','.py','.txt','.json','.yaml','.yml') and path.name != '.gitignore':
            raise ValueError('Unsupported public file type: '+entry)
        data = path.read_bytes()
        content = data.decode('utf-8')
        if re.search(r'(?<!\w)[A-Za-z]:[\\/]',content):
            raise ValueError('Absolute Windows path in public source: '+entry)
        if re.search(r'-----BEGIN (?:RSA |EC |OPENSSH )?PRIVATE KEY-----',content):
            raise ValueError('Private key in public source: '+entry)
        selected[entry] = data
    if not selected:
        raise ValueError('Empty release allowlist')
    return selected


def build(root, output):
    payload = sources(root)
    output = Path(output)
    output.parent.mkdir(parents=True,exist_ok=True)
    manifest = {name:hashlib.sha256(data).hexdigest() for name,data in sorted(payload.items())}
    payload['RELEASE_MANIFEST.json'] = (json.dumps(manifest,indent=2)+'\n').encode('utf-8')
    with zipfile.ZipFile(output,'x',compression=zipfile.ZIP_DEFLATED) as archive:
        for name,data in sorted(payload.items()):
            info = zipfile.ZipInfo('radiomics-skill-harness/'+name,date_time=(2026,1,1,0,0,0))
            info.compress_type = zipfile.ZIP_DEFLATED
            info.external_attr = 0o100644 << 16
            archive.writestr(info,data)
    print(f'Packaged {len(payload)} files: {output}')
    return manifest


if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args = parser.parse_args()
    build(Path(__file__).resolve().parents[1],args.output)
