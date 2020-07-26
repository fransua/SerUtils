import unittest
from ..stats.correlate import fit_with_uncertainty, latex_formula
import numpy as np

class TestCorrelate(unittest.TestCase):
    def test_fit(self):
        x = [50, 52, 53, 54, 58, 60, 62, 64, 66, 67, 68, 70, 72, 74, 76, 55, 50, 45, 65]
        y = [25, 50, 55, 75, 80, 85, 50, 65, 85, 55, 45, 45, 50, 75, 95, 65, 50, 40, 45]
        x = np.array(x)
        y = np.array(y)
        func_string = "A*x+B"
        p_x, p_y, z, confs, preds, r2 = fit_with_uncertainty(x, y, func_string, use_odr=False)
        formula = latex_formula(func_string, z)
        self.assertEqual(formula, '0.77x+13')
        self.assertEqual(int(p_x.sum()), 6050)
        self.assertEqual(int(p_y.sum()), 5926)
        self.assertEqual(int(confs.sum()), 1231)
        self.assertEqual(int(preds.sum()), 4040)
        self.assertEqual(int(r2*10000), 1432)
        self.assertEqual(formula, '0.77x+13')

if __name__ == '__main__':
    unittest.main()
