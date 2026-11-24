"""Legacy association screening followed by stability-aware feature reduction."""
import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from statistics_core import bh_fdr, vectorized_pearson

SEQUENCES = ('ADC', 'T1C_ph1', 'T1C_ph3', 'T2')
MARKERS = ('tls_area_total', 'tls_area_max', 'tls_area_mean', 'tls_area_std',
           'Agg_num', 'FL1_num', 'FL2_num', 'Total_num', 'Agg_flag', 'FL1_flag',
           'FL2_flag', 'TLS_012_flag')


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def table(path):
    frame = pd.read_excel(path, dtype={'ID': str})
    if 'ID' not in frame or frame.ID.isna().any() or frame.ID.str.strip().eq('').any() or frame.ID.duplicated().any():
        raise ValueError('Missing/empty/duplicate ID')
    return frame.set_index('ID')


def stability(x, y, full_r, args, seed):
    size = max(4, math.ceil(len(y)*args.fraction))
    counts = np.zeros(x.shape[1])
    directions = counts.copy()
    valid = 0
    rng = np.random.default_rng(seed)
    for _ in range(args.bootstrap):
        idx = rng.choice(len(y), size=size, replace=True)
        if np.std(y[idx]) == 0:
            continue
        r, p = vectorized_pearson(x[idx], y[idx])
        same = np.isfinite(r) & (np.sign(r) == np.sign(full_r))
        directions += same
        counts += same & (p < .05) & (np.abs(r) >= args.min_r)
        valid += 1
    return counts/max(valid, 1), directions/max(valid, 1), valid


def gates(r, p, q, freq, direction, args):
    legacy = (p < .05) & (np.abs(r) > .1)
    evidence = (q < args.q) | (freq >= args.frequency)
    if args.combine == 'and':
        evidence = (q < args.q) & (freq >= args.frequency)
    reduced = legacy & (np.abs(r) >= args.min_r) & (direction >= args.direction) & evidence
    return legacy, reduced


def run(args):
    if args.bootstrap < 1 or any(not 0 < x <= 1 for x in (args.fraction,args.q,args.min_r,args.frequency,args.direction)):
        raise ValueError('Invalid thresholds or bootstrap count')
    tls = table(args.tls)
    if not set(MARKERS).issubset(tls.columns):
        raise ValueError('TLS schema requires all 12 markers')
    args.output.mkdir(parents=True, exist_ok=False)
    report = dict(parameters={k:str(v) if isinstance(v,Path) else v for k,v in vars(args).items()},
                  tls_sha256=sha(args.tls), script_sha256=sha(__file__),
                  statistics_sha256=sha(Path(__file__).with_name('statistics_core.py')),
                  numpy=np.__version__, pandas=pd.__version__, sequences={})
    failed = False
    for si, seq in enumerate(SEQUENCES):
        audit = {}
        try:
            source = args.features/f'{seq}_features.xlsx'
            frame = table(source)
            names = [c for c in frame if c.startswith(('original_', 'wavelet-', 'log-'))]
            if not names:
                raise ValueError('No supported radiomics columns')
            ids = sorted(set(frame.index) & set(tls.index))
            audit.update(n_matched=len(ids), mri_only=sorted(set(frame.index)-set(tls.index)),
                         tls_only=sorted(set(tls.index)-set(frame.index)), source_sha256=sha(source))
            if len(ids) < 4:
                raise ValueError('Need at least four paired cases')
            x = frame.loc[ids,names].apply(pd.to_numeric, errors='raise').to_numpy(float)
            if not np.isfinite(x).all():
                raise ValueError('Nonfinite MRI features')
            rows = []
            for mi, marker in enumerate(MARKERS):
                y = pd.to_numeric(tls.loc[ids,marker], errors='raise').to_numpy(float)
                if not np.isfinite(y).all():
                    raise ValueError('Missing TLS values require an explicit analysis policy')
                r,p = vectorized_pearson(x,y)
                # Constant features remain in the audit, but not the finite BH family.
                q = bh_fdr(np.where(np.isfinite(r), p, np.nan))
                freq = direction = np.zeros(len(names))
                rounds = 0
                if args.stage == 'reduced':
                    freq,direction,rounds = stability(x,y,r,args,args.seed+si*10000+mi*100)
                legacy,reduced = gates(r,p,q,freq,direction,args)
                for j,name in enumerate(names):
                    rows.append(dict(Sequence=seq,TLSMarker=marker,Feature=name,N=len(ids),
                                     PearsonR=r[j],P=p[j],BHAdjustedP=q[j],
                                     LegacySelected=bool(legacy[j]),BootstrapFrequency=freq[j],
                                     DirectionConsistency=direction[j],ValidRounds=rounds,
                                     ReducedSelected=bool(reduced[j]) if args.stage=='reduced' else False))
            detail = pd.DataFrame(rows)
            column = 'LegacySelected' if args.stage=='legacy' else 'ReducedSelected'
            selected = sorted(set(detail.loc[detail[column],'Feature']))
            out = frame.loc[ids,selected].reset_index()
            with pd.ExcelWriter(args.output/f'{seq}_{args.stage}_features.xlsx', engine='openpyxl') as writer:
                out.to_excel(writer,sheet_name='Features',index=False)
                detail.to_excel(writer,sheet_name='AllTests',index=False)
                detail.groupby('TLSMarker')[column].sum().reset_index(name='AssociationCount').to_excel(writer,sheet_name='Counts',index=False)
            audit.update(status='COMPLETE',unique_features=len(selected))
        except Exception as exc:
            failed = True
            audit.update(status='FAILED',reason=str(exc))
        report['sequences'][seq] = audit
        print(seq,audit['status'])
    (args.output/'selection_qc.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    return int(failed)


def parser():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--features',type=Path,required=True,help='Original four-sequence feature tables, including all tests')
    p.add_argument('--tls',type=Path,required=True,help='Private local pathology workbook')
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--stage',choices=['legacy','reduced'],required=True)
    p.add_argument('--bootstrap',type=int,default=1000)
    p.add_argument('--fraction',type=float,default=.8)
    p.add_argument('--q',type=float,default=.1)
    p.add_argument('--min-r',type=float,default=.2)
    p.add_argument('--frequency',type=float,default=.65)
    p.add_argument('--direction',type=float,default=.8)
    p.add_argument('--combine',choices=['or','and'],default='or')
    p.add_argument('--seed',type=int,default=20260711)
    return p


if __name__=='__main__':
    raise SystemExit(run(parser().parse_args()))
