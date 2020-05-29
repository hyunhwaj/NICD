
import unittest
from unittest import mock
import argparse
import tempfile
import os
class Test_NICD(unittest.TestCase):
        
    
    def test_simple_run(self):
        try:
            from NICD.main import main
            with tempfile.TemporaryDirectory() as tmpdir:
                mock_args = {
                    "disease_a" : "C0524851",
                    "disease_b" : "C0520679",
                    "niter" : 0,
                    "random_seed": 123,
                    "outpath": tmpdir
                }
                main(**mock_args)
                self.assertTrue(os.path.isfile(f"{tmpdir}/permutation-test.csv"))
                self.assertTrue(os.path.isfile(f"{tmpdir}/edges.csv"))
                self.assertTrue(os.path.isfile(f"{tmpdir}/nodes.csv"))
        except:
            self.fail("calling NICD.main() failed!")