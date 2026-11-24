"""Local-only binary outcome modelling primitives shared with external validation."""
import hashlib
import importlib.metadata
import json
from pathlib import Path
import re

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.metrics import (roc_auc_score, average_precision_score, brier_score_loss,
                             confusion_matrix, cohen_kappa_score)
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler


def sha(path):
    h = hashlib.sha256()
    with open(path, 'rb') as handle:
        for block in iter(lambda: handle.read(1024*1024), b''):
            h.update(block)
    return h.hexdigest()


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False), encoding='utf-8')


def read_table(path):
    frame = pd.read_excel(path, dtype={'ID': str})
    if 'ID' not in frame or frame.empty:
        raise ValueError('Require nonempty patient table with ID')
    if frame.ID.isna().any() or frame.ID.str.strip().eq('').any() or frame.ID.duplicated().any():
        raise ValueError('Missing/empty/duplicate patient ID')
    if not frame.ID.eq(frame.ID.str.strip()).all():
        raise ValueError('IDs must not contain surrounding whitespace')
    return frame.set_index('ID').sort_index()


def read_signature(path):
    config = json.loads(Path(path).read_text(encoding='utf-8'))
    names = config['features']
    pattern = r'^(ADC|T1C_ph1|T1C_ph3|T2)::(original_|wavelet-|log-).+'
    if (config.get('origin') != 'independent_tls_discovery' or not names
            or len(names) != len(set(names))
            or any(not isinstance(n, str) or not re.fullmatch(pattern, n) for n in names)):
        raise ValueError('Require independent TLS signature and sequence-prefixed radiomics names')
    return names


def feature_matrix(frame, names, label=None):
    expected = set(names) | ({label} if label else set())
    if set(frame.columns) != expected:
        raise ValueError('Columns must exactly match the frozen feature manifest and optional label')
    x = frame.loc[:, names].apply(pd.to_numeric, errors='raise')
    if not np.isfinite(x.to_numpy(float)).all():
        raise ValueError('Missing/nonfinite features require a prespecified training-only policy')
    return x.to_numpy(float)


def labels(values, mode):
    y = pd.to_numeric(values, errors='raise').to_numpy(float)
    allowed = [0, 1, 2, 3] if mode == 'trg' else [0, 1]
    if not np.isfinite(y).all() or not np.isin(y, allowed).all():
        raise ValueError('Unexpected outcome coding; no rounding or coercion permitted')
    return (y != 0).astype(int) if mode == 'trg' else y.astype(int)


def disjoint(ids, other):
    if set(ids) & set(other):
        raise ValueError('Patient overlap across discovery/development/validation partitions')


def preprocessing(scale=False, pca=0):
    steps = [('scale', StandardScaler() if scale or pca else 'passthrough')]
    if pca:
        steps.append(('pca', PCA(n_components=pca, svd_solver='full')))
    return Pipeline(steps)


def probability(model, x, positive):
    classes = list(model.classes_)
    if set(classes) != {0, 1}:
        raise ValueError('Require fitted binary classes 0 and 1')
    p = np.asarray(model.predict_proba(x), float)[:, classes.index(positive)]
    if p.shape != (len(x),) or not np.isfinite(p).all() or ((p < 0) | (p > 1)).any():
        raise ValueError('Invalid predicted probabilities')
    return p


def metric_row(y, p, positive, threshold):
    truth = (np.asarray(y) == positive).astype(int)
    pred = (np.asarray(p) >= threshold).astype(int)
    tn, fp, fn, tp = confusion_matrix(truth, pred, labels=[0, 1]).ravel()
    ratio = lambda a, b: float(a/b) if b else None
    both = len(np.unique(truth)) == 2
    kappa = cohen_kappa_score(truth, pred) if both else float('nan')
    return dict(N=len(truth), AUC=float(roc_auc_score(truth, p)) if both else None,
                AUPRC=float(average_precision_score(truth, p)) if both else None,
                ACC=float((tp+tn)/len(truth)), Sensitivity=ratio(tp,tp+fn),
                Specificity=ratio(tn,tn+fp), PPV=ratio(tp,tp+fp), NPV=ratio(tn,tn+fn),
                F1=ratio(2*tp,2*tp+fp+fn), Brier=float(brier_score_loss(truth,p)),
                Kappa=float(kappa) if np.isfinite(kappa) else None,
                TN=int(tn), FP=int(fp), FN=int(fn), TP=int(tp))


def report_tables(y, p, positive, threshold):
    truth = (np.asarray(y) == positive).astype(int)
    bins = np.minimum((np.asarray(p)*10).astype(int), 9)
    calibration = pd.DataFrame({'Bin': bins, 'Probability': p, 'Outcome': truth}).groupby('Bin').agg(
        N=('Outcome','size'), MeanProbability=('Probability','mean'), ObservedRate=('Outcome','mean')).reset_index()
    dca = []
    for t in np.linspace(.05,.95,19):
        pred = p >= t
        nb = (np.sum(pred & (truth == 1))-np.sum(pred & (truth == 0))*t/(1-t))/len(truth)
        dca.append(dict(Threshold=float(t),Model=float(nb),TreatNone=0.,
                        TreatAll=float(truth.mean()-(1-truth.mean())*t/(1-t))))
    return pd.DataFrame([metric_row(y,p,positive,threshold)]), calibration, pd.DataFrame(dca)


def versions(backend):
    packages = ['numpy','pandas','scikit-learn','joblib','openpyxl']
    if backend == 'tabpfn':
        packages += ['tabpfn','torch']
    return {p: importlib.metadata.version(p) for p in packages}


def classifier(args, seed):
    if args.backend == 'logistic-smoke':
        from sklearn.linear_model import LogisticRegression
        return LogisticRegression(random_state=seed, max_iter=2000)
    from tabpfn import TabPFNClassifier
    return TabPFNClassifier(model_path=str(args.checkpoint.resolve()), device=args.device,
                            n_estimators=args.n_estimators, random_state=seed,
                            memory_saving_mode=True)
