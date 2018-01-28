import unittest
from .context import diskload

ALPHA = 0.1
W = 1.0
THETA = 0.05794132044
EPSILON = 1e-9

class Truncated(unittest.TestCase):
    def setUp(self):
        pass

    @classmethod
    def setUpClass(cls):
        love = cls.love = diskload.love_numbers.read()
        cls.extrapolated = diskload.love_numbers.extrapolate(love, 400000)
        pass
    
    def test_truncated_too_few_love_numbers(self):
        love = Truncated.love
        self.assertRaises(ValueError,
                          lambda: diskload.truncated( 0, diskload.Compensation.UNCOMPENSATED, 0, 0, len(love[0])*2, love ))
        pass

    def test_specific_input(self):
        love = Truncated.love
        u, v, g = diskload.truncated( ALPHA, diskload.Compensation.UNCOMPENSATED, THETA, W, 40000, love )
        self.assertTrue( abs( ( -2.0449592462e+00 - u ) / u ) < EPSILON ) 
        self.assertTrue( abs( ( -2.0687091115e-01 - v ) / v ) < EPSILON )
        self.assertTrue( abs( ( 4.2902044944e-01 - g ) / g ) < EPSILON )
        pass

    def test_four_hundred_thousand(self):
        love = Truncated.extrapolated        
        u, v, g = diskload.truncated( ALPHA, diskload.Compensation.UNCOMPENSATED, THETA, W, 400000, love )
        self.assertTrue( abs( ( -2.0446039399e+00 - u ) / u ) < EPSILON )
        self.assertTrue( abs( ( -2.0707066586e-01 - v ) / v ) < EPSILON )
        self.assertTrue( abs( ( 4.2896329466e-01  - g ) / g ) < EPSILON )
        pass

    def test_compensated_four_thousand(self):
        love = Truncated.love        
        u, v, g = diskload.truncated( ALPHA, diskload.Compensation.COMPENSATED, THETA, W, 4000, love )        
        self.assertTrue( abs( ( -2.0997389244e+00 - u ) / u ) < EPSILON )
        self.assertTrue( abs( ( -2.1688909314e-01 - v ) / v ) < EPSILON )
        self.assertTrue( abs( (  4.3742883583e-01 - g ) / g ) < EPSILON )
        pass
    
    def test_compensated_forty_thousand(self):
        love = Truncated.love        
        u, v, g = diskload.truncated( ALPHA, diskload.Compensation.COMPENSATED, THETA, W, 40000, love )                
        self.assertTrue( abs( ( -2.0448711591e+00 - u ) / u ) < EPSILON )
        self.assertTrue( abs( ( -2.0687106866e-01 - v ) / v ) < EPSILON )
        self.assertTrue( abs( (  4.2860575533e-01 - g ) / g ) < EPSILON )  
        pass
    
if __name__ == '__main__':
    unittest.main()
