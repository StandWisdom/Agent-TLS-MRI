"""Verify the upload archive contains reviewed sources only."""
import hashlib
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
import zipfile

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('release_builder',ROOT/'scripts/build_public_release.py')
builder = importlib.util.module_from_spec(spec)
spec.loader.exec_module(builder)


class ReleaseTests(unittest.TestCase):
    def test_archive_allowlist_hashes_and_no_overwrite(self):
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp)/'release.zip'
            manifest = builder.build(ROOT,output)
            with zipfile.ZipFile(output) as archive:
                self.assertIsNone(archive.testzip())
                prefix = 'radiomics-skill-harness/'
                self.assertEqual(set(archive.namelist()),
                                 {prefix+n for n in manifest} | {prefix+'RELEASE_MANIFEST.json'})
                saved = json.loads(archive.read(prefix+'RELEASE_MANIFEST.json'))
                self.assertEqual(saved,manifest)
                for name,digest in manifest.items():
                    self.assertEqual(hashlib.sha256(archive.read(prefix+name)).hexdigest(),digest)
                    self.assertNotIn('__pycache__',name)
                    self.assertNotIn('.model-test-env',name)
                    self.assertNotIn(Path(name).suffix,['.xlsx','.joblib','.nii','.ckpt'])
            with self.assertRaises(FileExistsError):
                builder.build(ROOT,output)

    def test_reject_path_escape_and_data_files(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            for entry in ('../secret.txt','/secret.txt','data.xlsx'):
                (root/'release-files.txt').write_text(entry,encoding='utf-8')
                (root/'data.xlsx').write_bytes(b'not-for-release')
                with self.assertRaises(ValueError):
                    builder.sources(root)


if __name__=='__main__':
    unittest.main()
