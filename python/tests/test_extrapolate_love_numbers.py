import unittest
from .context import diskload

class ExtrapolateLoveNumbers(unittest.TestCase):
    def setUp(self):
        pass

    @classmethod
    def setUpClass(cls):
        love = cls.love = diskload.love_numbers.read()
        cls.oo = love.shape[1] - 1
        cls.nmax = cls.oo * 2
        cls.love = diskload.love_numbers.extrapolate(love,cls.nmax)
        pass

    def test_extrapolation_length(self):
        love = ExtrapolateLoveNumbers.love        
        self.assertEqual(len(love[0]), ExtrapolateLoveNumbers.nmax + 1, 'wrong size after resize')        
        pass

    def test_extrapolation_h(self):
        love = ExtrapolateLoveNumbers.love
        h = love[0]
        self.assertEqual(h[45000], h[40000] )
        pass
 
    def test_extrapolation_l(self):
        love = ExtrapolateLoveNumbers.love
        oo = ExtrapolateLoveNumbers.oo
        n = round((ExtrapolateLoveNumbers.oo + ExtrapolateLoveNumbers.nmax)/2)
        l = love[2]
        l_oo = l[oo]
        self.assertEqual(l[n], oo * l_oo / n )
        pass   

    def test_extrapolation_k(self):
        love = ExtrapolateLoveNumbers.love
        oo = ExtrapolateLoveNumbers.oo
        n = round((ExtrapolateLoveNumbers.oo + ExtrapolateLoveNumbers.nmax)/2)
        k = love[1]
        k_oo = k[oo]
        self.assertEqual(k[n], oo * k_oo / n )
        pass   

    def test_first_extrapolation_l(self):
        love = ExtrapolateLoveNumbers.love
        oo = ExtrapolateLoveNumbers.oo
        n = ExtrapolateLoveNumbers.oo + 1
        l = love[2]
        l_oo = l[oo]
        self.assertEqual(l[n], oo * l_oo / n )
        pass   

    def test_first_extrapolation_k(self):
        love = ExtrapolateLoveNumbers.love
        oo = ExtrapolateLoveNumbers.oo
        n = ExtrapolateLoveNumbers.oo + 1
        k = love[1]
        k_oo = k[oo]
        self.assertEqual(k[n], oo * k_oo / n )
        pass       

    def test_last_extrapolation_l(self):
        love = ExtrapolateLoveNumbers.love
        oo = ExtrapolateLoveNumbers.oo
        n = ExtrapolateLoveNumbers.nmax
        l = love[2]
        l_oo = l[oo]
        self.assertEqual(l[n], oo * l_oo / n )
        pass   

    def test_last_extrapolation_k(self):
        love = ExtrapolateLoveNumbers.love
        oo = ExtrapolateLoveNumbers.oo
        n = ExtrapolateLoveNumbers.nmax
        k = love[1]
        k_oo = k[oo]
        self.assertEqual(k[n], oo * k_oo / n )
        pass           
    
if __name__ == '__main__':
    unittest.main()
