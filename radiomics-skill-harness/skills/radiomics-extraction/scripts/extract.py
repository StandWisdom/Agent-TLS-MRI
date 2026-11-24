"""Extract reviewed four-sequence MRI features; preserve failures in QC."""
import argparse
import hashlib
import json
import platform
from pathlib import Path
import shutil

import numpy as np
import pandas as pd
import SimpleITK as sitk
import yaml
import radiomics
from radiomics import featureextractor

SEQUENCES = ('ADC', 'T1C_ph1', 'T1C_ph3', 'T2')


def sha(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1048576), b''):
            h.update(block)
    return h.hexdigest()


def run(manifest, parameter_file, output):
    manifest, parameter_file, output = map(Path, (manifest, parameter_file, output))
    cases = json.loads(manifest.read_text(encoding='utf-8'))['cases']
    ids = [str(c['id']) for c in cases]
    if not cases or len(ids) != len(set(ids)) or any(not i.strip() for i in ids):
        raise ValueError('Require nonempty unique IDs')
    with parameter_file.open(encoding='utf-8-sig') as f:
        config = yaml.safe_load(f)
    config.setdefault('setting', {}).update(geometryTolerance=1e-5, correctMask=False)
    extractor = featureextractor.RadiomicsFeatureExtractor(config)
    output.mkdir(parents=True, exist_ok=False)
    shutil.copyfile(parameter_file, output/'requested_parameters.yaml')
    tables = {s: [] for s in SEQUENCES}
    qc = []
    for case in cases:
        for seq in SEQUENCES:
            record = dict(id=str(case['id']), sequence=seq, status='FAILED')
            try:
                if case.get('reviewed') is not True:
                    raise ValueError('Input not reviewed')
                pair = case['sequences'][seq]
                ip, mp = [(manifest.parent / pair[k]).resolve() for k in ('image', 'mask')]
                record.update(image_sha256=sha(ip), mask_sha256=sha(mp))
                image, mask = sitk.ReadImage(str(ip)), sitk.ReadImage(str(mp))
                if any(i.GetDimension() != 3 or i.GetNumberOfComponentsPerPixel() != 1 for i in (image, mask)):
                    raise ValueError('Require scalar 3D inputs')
                if image.GetSize() != mask.GetSize() or any(not np.allclose(getattr(image,k)(), getattr(mask,k)(), atol=1e-5, rtol=0) for k in ('GetSpacing','GetOrigin','GetDirection')):
                    raise ValueError('Geometry mismatch; return to MRI preparation')
                if not np.isfinite(sitk.GetArrayViewFromImage(image)).all():
                    raise ValueError('Nonfinite image; explicit preprocessing policy required')
                ma = sitk.GetArrayViewFromImage(mask)
                if not np.isfinite(ma).all() or (ma < 0).any() or not (ma == np.floor(ma)).all():
                    raise ValueError('Mask must contain nonnegative integer labels')
                values = extractor.execute(image, mask)
                features = {k: float(v) for k,v in values.items() if not k.startswith('diagnostics_')}
                if not features or not np.isfinite(list(features.values())).all():
                    raise ValueError('Empty/nonfinite extracted features')
                tables[seq].append(dict(ID=str(case['id']), **features))
                record.update(status='EXTRACTED', feature_count=len(features))
            except Exception as exc:
                record['reason'] = str(exc)
            qc.append(record)
    for seq, rows in tables.items():
        pd.DataFrame(rows if rows else {'ID': []}).to_excel(output/f'{seq}_features.xlsx', index=False)
    report = dict(python=platform.python_version(), pyradiomics=radiomics.__version__, numpy=np.__version__,
                  simpleitk=sitk.Version_VersionString(), yaml_sha256=sha(parameter_file),
                  manifest_sha256=sha(manifest), script_sha256=sha(__file__), effective_settings=extractor.settings,
                  image_types=extractor.enabledImagetypes, feature_classes=extractor.enabledFeatures, cases=qc)
    (output/'extraction_qc.json').write_text(json.dumps(report, indent=2, default=str), encoding='utf-8')
    return int(any(r['status'] == 'FAILED' for r in qc))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('manifest', type=Path)
    p.add_argument('--yaml', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    raise SystemExit(run(a.manifest, a.yaml, a.output))
