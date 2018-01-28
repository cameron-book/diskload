import unittest
from .context import diskload

class ReadLoveNumbers(unittest.TestCase):
    def setUp(self):
        pass

    love = None
    
    @classmethod
    def setUpClass(cls):
        ReadLoveNumbers.love = diskload.love_numbers.read()        
        pass
    
    def test_right_length(self):
        love = ReadLoveNumbers.love
        self.assertEqual(len(love[0]), 40001)
        self.assertEqual(len(love[1]), 40001)
        self.assertEqual(len(love[2]), 40001)
        self.assertEqual(len(love), 3)
        pass

    def test_right_numbers(self):
        love = ReadLoveNumbers.love        
        self.assertEqual(love[0][0], -0.216)
        self.assertEqual(love[1][0], 0 )
        self.assertEqual(love[2][0], 0 )
        self.assertEqual(love[0][3901], -6.2144186090000 )
        self.assertEqual(love[1][3901], 0.0004851105826 )
        self.assertEqual(love[2][3901], -0.0007828382211 );
    
if __name__ == '__main__':
    unittest.main()
