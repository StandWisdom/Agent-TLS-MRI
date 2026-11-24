"""Fixed patient-level 70/30 split; freeze the 70%-trained model for external use."""
import argparse
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split

from model_common import (classifier, disjoint, feature_matrix, labels,
                          preprocessing, probability, read_signature, read_table,
                          report_tables, sha, versions, write_json)


def run(args):
    if not 0 < args.threshold < 1 or args.pca_components < 0 or args.n_estimators < 1:
        raise ValueError('Invalid frozen parameters')
    names = read_signature(args.signature)
    discovery = read_table(args.discovery_ids)
    frame = read_table(args.data)
    disjoint(frame.index, discovery.index)
    x = feature_matrix(frame, names, args.label)
    y = labels(frame[args.label], args.label_mode)
    if min(np.bincount(y, minlength=2)) < 2:
        raise ValueError('Stratification requires at least two patients per class')
    env = versions(args.backend)
    if args.backend == 'tabpfn':
        if not args.tabpfn_version or env['tabpfn'] != args.tabpfn_version:
            raise ValueError('Specify the exact installed TabPFN package version')
        if not args.checkpoint or not args.checkpoint.is_file():
            raise ValueError('Require an explicit local checkpoint; no automatic model download')
        from tabpfn.model_loading import save_fitted_tabpfn_model
    elif not all(i.startswith('synthetic') for i in list(frame.index)+list(discovery.index)):
        raise ValueError('logistic-smoke is restricted to synthetic IDs and software testing')
    args.output.mkdir(parents=True, exist_ok=False)
    write_json(args.output/'status.json', {'status':'INCOMPLETE'})
    train,val = train_test_split(np.arange(len(y)),test_size=.30, stratify=y,random_state=args.seed)
    if len(np.unique(y[train])) != 2 or len(np.unique(y[val])) != 2:
        raise ValueError('Both partitions must contain both outcome classes')
    if args.pca_components > min(x.shape[1], len(train)):
        raise ValueError('PCA count exceeds training dimensions')
    transform = preprocessing(args.scale,args.pca_components)
    model = classifier(args,args.seed)
    model.fit(transform.fit_transform(x[train],y[train]),y[train])
    probs = probability(model,transform.transform(x[val]),args.positive_label)
    predictions = pd.DataFrame(dict(ID=frame.index[val],TrueLabel=y[val],
                                    PositiveProbability=probs,
                                    PredictedLabel=np.where(probs>=args.threshold,args.positive_label,1-args.positive_label)))
    metrics,calibration,dca = report_tables(y[val],probs,args.positive_label,args.threshold)
    assignments = [dict(ID=frame.index[i],Role=role) for role,idx in [('train',train),('val',val)] for i in idx]
    with pd.ExcelWriter(args.output/'split_data_private.xlsx') as writer:
        frame.iloc[train].reset_index().to_excel(writer,sheet_name='Train70',index=False)
        frame.iloc[val].reset_index().to_excel(writer,sheet_name='Validation30',index=False)
    with pd.ExcelWriter(args.output/'internal_validation.xlsx') as writer:
        for name,df in [('HeldOutPredictions',predictions),('OverallMetrics',metrics),
                        ('Calibration',calibration),('DCA',dca),('Splits',pd.DataFrame(assignments))]:
            df.to_excel(writer,sheet_name=name,index=False)
    # Persist the SAME model used for internal validation; do not refit on the 30%.
    joblib.dump(transform,args.output/'preprocessor.joblib')
    model_name = 'model.tabpfn_fit' if args.backend=='tabpfn' else 'synthetic_model.joblib'
    if args.backend=='tabpfn':
        save_fitted_tabpfn_model(model,args.output/model_name)
    else:
        joblib.dump(model,args.output/model_name)
    contract = dict(status='COMPLETE',backend=args.backend,features=names,model_file=model_name,
                    positive_label=args.positive_label,label_mode=args.label_mode,label_column=args.label,
                    threshold=args.threshold,seed=args.seed,n_estimators=args.n_estimators,
                    scale=args.scale,pca_components=args.pca_components,versions=env,
                    development_ids=list(frame.index),discovery_ids=list(discovery.index),
                    training_ids=list(frame.index[train]),validation_ids=list(frame.index[val]),
                    input_sha256=sha(args.data),signature_sha256=sha(args.signature),
                    discovery_sha256=sha(args.discovery_ids),script_sha256=sha(__file__),
                    common_sha256=sha(Path(__file__).with_name('model_common.py')),
                    checkpoint_sha256=sha(args.checkpoint) if args.backend=='tabpfn' else None,
                    checkpoint_path=str(args.checkpoint.resolve()) if args.backend=='tabpfn' else None,
                    model_sha256=sha(args.output/model_name),
                    preprocessor_sha256=sha(args.output/'preprocessor.joblib'),
                    evaluation='stratified_70_train_30_validation_no_refit',
                    data_use='research' if args.backend=='tabpfn' else 'SYNTHETIC_SOFTWARE_TEST_ONLY')
    write_json(args.output/'private_contract.json',contract)
    write_json(args.output/'status.json',{'status':'COMPLETE'})
    return contract


def parser():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('data','signature','discovery-ids','output'):
        p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--label',default='TRG')
    p.add_argument('--label-mode',choices=['trg','binary'],required=True)
    p.add_argument('--positive-label',type=int,choices=[0,1],default=0)
    p.add_argument('--threshold',type=float,default=.5)
    p.add_argument('--seed',type=int,default=42)
    p.add_argument('--scale',action='store_true')
    p.add_argument('--pca-components',type=int,default=0)
    p.add_argument('--n-estimators',type=int,default=8)
    p.add_argument('--device',default='cuda')
    p.add_argument('--backend',choices=['tabpfn','logistic-smoke'],default='tabpfn')
    p.add_argument('--checkpoint',type=Path)
    p.add_argument('--tabpfn-version')
    return p


if __name__=='__main__':
    run(parser().parse_args())
