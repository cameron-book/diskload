import unittest
from .context import diskload

ALPHA = 0.1
W = 1.0
THETA = 0.05794132044

class Agreement(unittest.TestCase):
    def setUp(self):
        pass

    @classmethod
    def setUpClass(cls):
        love = cls.love = diskload.love_numbers.read()
        cls.extrapolated = diskload.love_numbers.extrapolate(love, 400000)
        pass
    
    def test_uncompensated_approximation(self):
        love = Agreement.extrapolated
        uU, vU, gU = diskload.truncated( ALPHA, diskload.Compensation.UNCOMPENSATED, THETA, W, 400000, love )
        u, v, g = diskload.elliptic( ALPHA, diskload.Compensation.UNCOMPENSATED, THETA, W, 40000, love )
        self.assertTrue( abs( ( uU - u ) / u ) < 1e-6 )
        self.assertTrue( abs( ( vU - v ) / v ) < 1e-3 )  
        self.assertTrue( abs( ( gU - g ) / g ) < 1e-6 )
        pass
    
if __name__ == '__main__':
    unittest.main()
