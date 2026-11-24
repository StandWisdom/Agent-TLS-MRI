"""Predict with a frozen 70%-trained model, then score separately after label release."""
import argparse
import json
from pathlib import Path
import sys

import joblib
import numpy as np
import pandas as pd

sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'internal-model-development/scripts'))
from model_common import (disjoint, feature_matrix, labels, probability, read_table,
                          report_tables, sha, versions, write_json)


def predict(args):
    if not args.trust_local_model:
        raise ValueError('Model deserialization requires --trust-local-model; never load untrusted files')
    contract_path = args.bundle/'private_contract.json'
    contract = json.loads(contract_path.read_text(encoding='utf-8'))
    if contract['status'] != 'COMPLETE' or contract['evaluation'] != 'stratified_70_train_30_validation_no_refit':
        raise ValueError('Require a complete, frozen 70/30 model bundle')
    if contract['versions'] != versions(contract['backend']):
        raise ValueError('Environment differs from development; restore the recorded environment')
    model_file = args.bundle/contract['model_file']
    if model_file.resolve().parent != args.bundle.resolve():
        raise ValueError('Model path must stay within the trusted bundle')
    if (sha(model_file) != contract['model_sha256'] or
            sha(args.bundle/'preprocessor.joblib') != contract['preprocessor_sha256']):
        raise ValueError('Frozen artifact checksum mismatch')
    frame = read_table(args.features)
    disjoint(frame.index,contract['development_ids']+contract['discovery_ids'])
    x = feature_matrix(frame,contract['features'])
    if contract['backend']=='logistic-smoke' and not all(i.startswith('synthetic') for i in frame.index):
        raise ValueError('Synthetic software test only')
    transform = joblib.load(args.bundle/'preprocessor.joblib')
    if contract['backend']=='tabpfn':
        checkpoint = Path(contract['checkpoint_path'])
        if not checkpoint.is_file() or sha(checkpoint) != contract['checkpoint_sha256']:
            raise ValueError('Original local checkpoint must remain available and unchanged')
        from tabpfn.model_loading import load_fitted_tabpfn_model
        model = load_fitted_tabpfn_model(model_file,device=args.device)
    else:
        model = joblib.load(model_file)
    p = probability(model,transform.transform(x),contract['positive_label'])
    args.output.mkdir(parents=True,exist_ok=False)
    result = pd.DataFrame(dict(ID=frame.index,PositiveProbability=p,
                               PredictedLabel=np.where(p>=contract['threshold'],contract['positive_label'],1-contract['positive_label'])))
    result.to_excel(args.output/'predictions_private.xlsx',index=False)
    write_json(args.output/'prediction_receipt.json',dict(
        status='FROZEN_PREDICTIONS',positive_label=contract['positive_label'],
        threshold=contract['threshold'],label_mode=contract['label_mode'],label_column=contract['label_column'],
        data_use=contract['data_use'],prediction_sha256=sha(args.output/'predictions_private.xlsx'),
        input_sha256=sha(args.features),bundle_contract_sha256=sha(contract_path),
        script_sha256=sha(__file__),N=len(result)))
    return result


def score(args):
    receipt = json.loads((args.predictions/'prediction_receipt.json').read_text(encoding='utf-8'))
    path = args.predictions/'predictions_private.xlsx'
    if receipt['status'] != 'FROZEN_PREDICTIONS' or sha(path) != receipt['prediction_sha256']:
        raise ValueError('Prediction file changed after freezing')
    pred = read_table(path)
    truth = read_table(args.labels)
    if set(truth.columns) != {receipt['label_column']} or set(truth.index) != set(pred.index):
        raise ValueError('Labels must contain exactly the frozen patient set and outcome column')
    y = labels(truth.loc[pred.index,receipt['label_column']],receipt['label_mode'])
    p = pred.PositiveProbability.to_numpy(float)
    if not np.isfinite(p).all() or ((p<0)|(p>1)).any():
        raise ValueError('Invalid frozen predictions')
    metrics,calibration,dca = report_tables(y,p,receipt['positive_label'],receipt['threshold'])
    args.output.mkdir(parents=True,exist_ok=False)
    with pd.ExcelWriter(args.output/'external_validation.xlsx') as writer:
        for name,df in [('Metrics',metrics),('Calibration',calibration),('DCA',dca)]:
            df.to_excel(writer,sheet_name=name,index=False)
    write_json(args.output/'evaluation_receipt.json',dict(
        status='COMPLETE',data_use=receipt['data_use'],labels_sha256=sha(args.labels),
        prediction_receipt_sha256=sha(args.predictions/'prediction_receipt.json'),
        threshold=receipt['threshold'],positive_label=receipt['positive_label'],script_sha256=sha(__file__)))
    return metrics


def parser():
    p = argparse.ArgumentParser(description=__doc__)
    commands = p.add_subparsers(dest='command',required=True)
    infer = commands.add_parser('predict')
    for field in ('bundle','features','output'):
        infer.add_argument('--'+field,type=Path,required=True)
    infer.add_argument('--device',default='cuda')
    infer.add_argument('--trust-local-model',action='store_true')
    evaluate = commands.add_parser('score')
    for field in ('predictions','labels','output'):
        evaluate.add_argument('--'+field,type=Path,required=True)
    return p


if __name__=='__main__':
    args = parser().parse_args()
    predict(args) if args.command=='predict' else score(args)
