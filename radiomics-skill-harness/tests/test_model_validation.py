"""Synthetic regression tests for split isolation, frozen inference and scoring."""
import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

import joblib
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT/'skills/internal-model-development/scripts'
sys.path.insert(0,str(SCRIPTS))
from model_common import labels, metric_row, read_table, read_signature, sha


def load(name, path):
    spec = importlib.util.spec_from_file_location(name,path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


develop = load('develop_model',SCRIPTS/'develop_model.py')
external = load('validate_external',ROOT/'skills/external-model-validation/scripts/validate_external.py')


class ModelValidationTests(unittest.TestCase):
    def prepare(self,root):
        rng = np.random.default_rng(31)
        self.names = ['ADC::original_firstorder_Mean','T2::original_firstorder_Mean']
        signature = {'origin':'independent_tls_discovery','features':self.names}
        (root/'signature.json').write_text(json.dumps(signature),encoding='utf-8')
        pd.DataFrame({'ID':['synthetic-discovery-1']}).to_excel(root/'discovery.xlsx',index=False)
        df = pd.DataFrame(rng.normal(size=(40,2)),columns=self.names)
        df.insert(0,'ID',[f'synthetic-dev-{i:03}' for i in range(40)])
        df['TRG'] = np.tile([0,1,0,3],10)
        df.to_excel(root/'dev.xlsx',index=False)
        ext = pd.DataFrame(rng.normal(10,2,size=(12,2)),columns=self.names)
        ext.insert(0,'ID',[f'synthetic-ext-{i:03}' for i in range(12)])
        ext.to_excel(root/'ext.xlsx',index=False)
        pd.DataFrame({'ID':ext.ID,'TRG':np.tile([0,1],6)}).to_excel(root/'labels.xlsx',index=False)
        args = develop.parser().parse_args([
            '--data',str(root/'dev.xlsx'),'--signature',str(root/'signature.json'),
            '--discovery-ids',str(root/'discovery.xlsx'),'--output',str(root/'bundle'),
            '--label-mode','trg','--backend','logistic-smoke','--scale','--pca-components','1'])
        return args

    def infer_args(self,root,output='pred'):
        return external.parser().parse_args(['predict','--bundle',str(root/'bundle'),
            '--features',str(root/'ext.xlsx'),'--output',str(root/output),'--trust-local-model','--device','cpu'])

    def test_70_30_train_only_and_no_refit(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            args = self.prepare(root)
            fit_sizes = []
            original = LogisticRegression.fit
            def spy(model,x,y,*a,**kw):
                fit_sizes.append(len(y))
                return original(model,x,y,*a,**kw)
            with patch.object(LogisticRegression,'fit',spy):
                contract = develop.run(args)
            self.assertEqual(fit_sizes,[28])
            self.assertEqual(len(contract['validation_ids']),12)
            self.assertFalse(set(contract['training_ids']) & set(contract['validation_ids']))
            frame = read_table(root/'dev.xlsx')
            transform = joblib.load(root/'bundle/preprocessor.joblib')
            expected = frame.loc[contract['training_ids'],self.names].mean().to_numpy()
            np.testing.assert_allclose(transform.named_steps['scale'].mean_,expected)
            self.assertFalse(np.allclose(expected,frame[self.names].mean()))
            results = pd.read_excel(root/'bundle/internal_validation.xlsx',sheet_name='HeldOutPredictions')
            self.assertEqual(set(results.ID),set(contract['validation_ids']))
            self.assertEqual(pd.read_excel(root/'bundle/split_data_private.xlsx',sheet_name='Train70').shape[0],28)

    def test_external_no_fit_and_separate_label_scoring(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            develop.run(self.prepare(root))
            before = sha(root/'bundle/preprocessor.joblib')
            with patch.object(LogisticRegression,'fit',side_effect=AssertionError('External fit prohibited')):
                first = external.predict(self.infer_args(root))
            self.assertEqual(before,sha(root/'bundle/preprocessor.joblib'))
            score = external.parser().parse_args(['score','--predictions',str(root/'pred'),
                '--labels',str(root/'labels.xlsx'),'--output',str(root/'score')])
            metrics = external.score(score)
            self.assertEqual(metrics.N.iloc[0],12)
            original_pred_hash = sha(root/'pred/predictions_private.xlsx')
            truth = pd.read_excel(root/'labels.xlsx')
            truth.TRG = 1-truth.TRG
            truth.to_excel(root/'labels.xlsx',index=False)
            second = external.predict(self.infer_args(root,'pred2'))
            np.testing.assert_array_equal(first.PositiveProbability,second.PositiveProbability)
            self.assertEqual(original_pred_hash,sha(root/'pred/predictions_private.xlsx'))

    def test_external_rejects_labels_overlap_missing_columns_and_tampering(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            develop.run(self.prepare(root))
            good = pd.read_excel(root/'ext.xlsx')
            for bad in (good.assign(TRG=0),good.drop(columns=self.names[0]),
                        good.assign(ID=['synthetic-dev-000']+list(good.ID[1:]))):
                bad.to_excel(root/'ext.xlsx',index=False)
                with self.assertRaises(ValueError):
                    external.predict(self.infer_args(root))
            good.to_excel(root/'ext.xlsx',index=False)
            external.predict(self.infer_args(root))
            predfile = root/'pred/predictions_private.xlsx'
            changed = pd.read_excel(predfile)
            changed.loc[0,'PositiveProbability'] = .123
            changed.to_excel(predfile,index=False)
            args = external.parser().parse_args(['score','--predictions',str(root/'pred'),
                '--labels',str(root/'labels.xlsx'),'--output',str(root/'score')])
            with self.assertRaises(ValueError):
                external.score(args)

    def test_discovery_overlap_and_target_in_signature_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            args = self.prepare(root)
            pd.DataFrame({'ID':['synthetic-dev-000']}).to_excel(root/'discovery.xlsx',index=False)
            with self.assertRaises(ValueError):
                develop.run(args)
            (root/'signature.json').write_text(json.dumps({'origin':'independent_tls_discovery','features':['TRG']}))
            with self.assertRaises(ValueError):
                read_signature(root/'signature.json')

    def test_labels_probability_direction_and_single_class(self):
        np.testing.assert_array_equal(labels(pd.Series([0,1,2,3]),'trg'),[0,1,1,1])
        with self.assertRaises(ValueError):
            labels(pd.Series([0,.5,1]),'binary')
        result = metric_row(np.array([0,1,0,1]),np.array([.9,.1,.8,.2]),0,.5)
        self.assertEqual(result['Kappa'],1.)
        self.assertEqual(result['AUC'],1.)
        self.assertIsNone(metric_row(np.array([0,0]),np.array([.8,.9]),0,.5)['AUC'])


if __name__=='__main__':
    unittest.main()
