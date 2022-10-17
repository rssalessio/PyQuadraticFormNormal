import unittest
from PyQuadraticFormNormal import davies_method


class TestDaviesMethod(unittest.TestCase):
    def test_returns(self):
        """ Test that the function returns the correct number of results """
        result = davies_method([2], [1], [0], [1], sigma=0)
        self.assertEqual(len(result), 1, "Should return only one result")
        result = davies_method([2, 1], [1], [0], [1], sigma=0)
        self.assertEqual(len(result), 2, "Should return two results")

    def test_value_chi_one(self):
        """ Test from a chi^2(1) """
        result = davies_method([2], [1], [0], [1], sigma=0)
        self.assertEqual(result[0], 0.84269977882736, "Incorrect result")

        result = davies_method([2], [1], [0], [1], sigma=0, accuracy=1e-6)
        self.assertEqual(result[0], 0.8427006907020171, "Incorrect result")

        result = davies_method([1], [1], [0], [1], sigma=0)
        self.assertEqual(result[0], 0.6826886153722266, "Incorrect result")

        result = davies_method([.5], [1], [0], [1], sigma=0)
        self.assertEqual(result[0], 0.5204991504663418, "Incorrect result")

        result = davies_method([.1], [1], [0], [1], sigma=0)
        self.assertEqual(result[0], 0.24816987371489302, "Incorrect result")

    def test_value_chi_two(self):
        """ Test from a chi^2(2) """
        result = davies_method([2], [1], [0], [2], sigma=0)
        self.assertEqual(result[0], 0.632119574889251, "Incorrect result")

        result = davies_method([2], [1, 1], [0, 0], [1, 1], sigma=0)
        self.assertEqual(result[0], 0.6321195748892511, "Incorrect result")

    def test_value_normal(self):
        """ Test from a normal distribution N(0, 1) """
        result = davies_method([1], [0], [0], [0], sigma=1)
        self.assertEqual(result[0], 0.8413445419723138, "Incorrect result")

        result = davies_method([-1], [0], [0], [0], sigma=1)
        self.assertEqual(result[0], 0.1586554580276862, "Incorrect result")

    def test_multiple_values(self):
        result = davies_method([1, 1, 1, 1], [1], [0], [2], sigma=0)
        
        for x in result:
            self.assertEqual(x, result[0], 'Values are not the same')

        result = davies_method([1, 1, 1, 1], [1, 1, 1], [0, 0, 0], [1, 1, 1], sigma=0)
        
        for x in result:
            self.assertEqual(x, result[0], 'Values are not the same')

    def test_huge_values(self):
        """ Test method with lot of parameters """
        N = 10000
        result = davies_method([9600], [1] * N, [0] * N, [1] * N, sigma=0)[0]
        self.assertEqual(result, 0.0021031643039009507, 'Value is not the same')

        result = davies_method([10600], [1] * N, [0] * N, [1] * N, sigma=0)[0]
        self.assertEqual(result, 0.9999844082691355, 'Value is not the same')

    def test_invalid_arguments(self):
        """ Test invalid arguments """
        self.assertRaises(ValueError, lambda:  davies_method([], [0], [0], [0]))
        self.assertRaises(ValueError, lambda:  davies_method([0], [0], [0], []))
        self.assertRaises(ValueError, lambda:  davies_method([0], [0,1], [0], [0]))
        self.assertRaises(ValueError, lambda:  davies_method([0], [0], [0], [0], accuracy=0))
        self.assertRaises(ValueError, lambda:  davies_method([0], [0], [0], [0], accuracy=-1))
        self.assertRaises(ValueError, lambda:  davies_method([1], [0,1], [0,1], [0,1,1]))
        self.assertRaises(IndexError, lambda:  davies_method([0], [0], [0], [0], limit=0))
        self.assertRaises(IndexError, lambda:  davies_method([0], [0], [0], [0], limit=-1))
        self.assertRaises(ValueError, lambda:  davies_method([0], [0], [0], [-1]))


    def test_generic(self):
        result = davies_method([1], [6, 3, 1], [0, 0, 0] , [1, 1, 1], limit=1000, accuracy = 0.0001, sigma=0)[0]
        self.assertAlmostEqual(result,  0.054213, places=6, msg='Value is not the same')

        result = davies_method([7], [6, 3, 1], [0, 0, 0] , [1, 1, 1], limit=1000, accuracy = 0.0001, sigma=0)[0]
        self.assertAlmostEqual(result,  0.493555, places=6, msg='Value is not the same')

        result = davies_method([20], [6, 3, 1], [0, 0, 0] , [1, 1, 1], limit=1000, accuracy = 0.0001, sigma=0)[0]
        self.assertAlmostEqual(result,  0.876027, places=6, msg='Value is not the same')

        result = davies_method([2], [6, 3, 1], [0, 0, 0] , [2, 2, 2], limit=1000, accuracy = 0.0001, sigma=0)[0]
        self.assertAlmostEqual(result,  0.006435, places=6, msg='Value is not the same')

        result = davies_method([20], [6, 3, 1], [0, 0, 0] , [2, 2, 2], limit=1000, accuracy = 0.0001, sigma=0)[0]
        self.assertAlmostEqual(result,  0.600208, places=6, msg='Value is not the same')

        result = davies_method([60], [6, 3, 1], [0, 0, 0] , [2, 2, 2], limit=1000, accuracy = 0.0001, sigma=0)[0]
        self.assertAlmostEqual(result,  0.983897, places=6, msg='Value is not the same')

        result = davies_method([10], [6, 3, 1], [0, 0, 0] , [6, 4, 2], limit=1000, accuracy = 0.0001, sigma=0)[0]
        self.assertAlmostEqual(result,  0.002697, places=6, msg='Value is not the same')

        result = davies_method([70], [7, 3, 7, 3], [6, 2, 6, 2] , [6, 2, 1, 1], limit=1000, accuracy = 0.0001, sigma=0)[0]
        self.assertAlmostEqual(result,  0.043679, places=6, msg='Value is not the same')
 
        result = davies_method([160], [7, 3, 7, 3], [6, 2, 6, 2] , [6, 2, 1, 1], limit=1000, accuracy = 0.0001, sigma=0)[0]
        self.assertAlmostEqual(result,  0.584765, places=6, msg='Value is not the same')
       
        result = davies_method([260], [7, 3, 7, 3], [6, 2, 6, 2] , [6, 2, 1, 1], limit=1000, accuracy = 0.0001, sigma=0)[0]
        self.assertAlmostEqual(result,  0.953774, places=6, msg='Value is not the same')

        result = davies_method([-40], [7, 3, -7, -3], [6, 2, 6, 2] , [6, 2, 1, 1], limit=1000, accuracy = 0.0001, sigma=0)[0]
        self.assertAlmostEqual(result,  0.078208, places=6, msg='Value is not the same')

        result = davies_method([40], [7, 3, -7, -3], [6, 2, 6, 2] , [6, 2, 1, 1], limit=1000, accuracy = 0.0001, sigma=0)[0]
        self.assertAlmostEqual(result,  0.522108, places=6, msg='Value is not the same')

        result = davies_method([140], [7, 3, -7, -3], [6, 2, 6, 2] , [6, 2, 1, 1], limit=1000, accuracy = 0.0001, sigma=0)[0]
        self.assertAlmostEqual(result,  0.960370, places=6, msg='Value is not the same')

        result = davies_method([120], [6, 3, 1, 6, 3, 1, 7, 3, 7, 3], [0, 0, 0, 0, 0, 0, 6, 2, 6, 2] , [6, 4, 2, 2, 4, 6, 6, 2, 1, 1], limit=1000, accuracy = 0.0001, sigma=0)[0]
        self.assertAlmostEqual(result,  0.015844, places=6, msg='Value is not the same')

        result = davies_method([240], [6, 3, 1, 6, 3, 1, 7, 3, 7, 3], [0, 0, 0, 0, 0, 0, 6, 2, 6, 2] , [6, 4, 2, 2, 4, 6, 6, 2, 1, 1], limit=1000, accuracy = 0.0001, sigma=0)[0]
        self.assertAlmostEqual(result,  0.573625, places=6, msg='Value is not the same')

        result = davies_method([400], [6, 3, 1, 6, 3, 1, 7, 3, 7, 3], [0, 0, 0, 0, 0, 0, 6, 2, 6, 2] , [6, 4, 2, 2, 4, 6, 6, 2, 1, 1], limit=1000, accuracy = 0.0001, sigma=0)[0]
        self.assertAlmostEqual(result,  0.988332, places=6, msg='Value is not the same')

        result = davies_method([5], [30,1], [0, 0] , [1, 10], limit=1000, accuracy = 0.0001, sigma=0)[0]
        self.assertAlmostEqual(result,  0.015392, places=6, msg='Value is not the same')

        result = davies_method([25], [30,1], [0, 0] , [1, 10], limit=1000, accuracy = 0.0001, sigma=0)[0]
        self.assertAlmostEqual(result,  0.510819, places=6, msg='Value is not the same')

        result = davies_method([100], [30,1], [0, 0] , [1, 10], limit=1000, accuracy = 0.0001, sigma=0)[0]
        self.assertAlmostEqual(result,  0.916340, places=6, msg='Value is not the same')

        result = davies_method([10], [30,1], [0, 0] , [1, 20], limit=1000, accuracy = 0.0001, sigma=0)[0]
        self.assertAlmostEqual(result,  0.004925, places=6, msg='Value is not the same')

        result = davies_method([40], [30,1], [0, 0] , [1, 20], limit=1000, accuracy = 0.0001, sigma=0)[0]
        self.assertAlmostEqual(result,  0.573251, places=6, msg='Value is not the same')

        result = davies_method([100], [30,1], [0, 0] , [1, 20], limit=1000, accuracy = 0.0001, sigma=0)[0]
        self.assertAlmostEqual(result,  0.896501, places=6, msg='Value is not the same')

        result = davies_method([20], [30,1], [0, 0] , [1, 30], limit=1000, accuracy = 0.0001, sigma=0)[0]
        self.assertAlmostEqual(result,  0.017101, places=6, msg='Value is not the same')

        result = davies_method([50], [30,1], [0, 0] , [1, 30], limit=1000, accuracy = 0.0001, sigma=0)[0]
        self.assertAlmostEqual(result,  0.566488, places=6, msg='Value is not the same')

        result = davies_method([100], [30,1], [0, 0] , [1, 30], limit=1000, accuracy = 0.0001, sigma=0)[0]
        self.assertAlmostEqual(result,  0.871323, places=6, msg='Value is not the same')

        result = davies_method([5], [30,1], [0, 0] , [1, 10], limit=1000, accuracy = 0.0001, sigma=0)[0]
        self.assertAlmostEqual(result,  0.015392401492801577, msg='Value is not the same')

        result = davies_method([50], [6, 3, 1], [0, 0, 0] , [6, 4, 2], limit=1000, accuracy = 0.0001, sigma=0)[0]
        self.assertAlmostEqual(result,  0.564753, places=6, msg='Value is not the same')

        result = davies_method([120], [6, 3, 1], [0, 0, 0] , [6, 4, 2], limit=1000, accuracy = 0.0001, sigma=0)[0]
        self.assertAlmostEqual(result,  0.991229, places=6, msg='Value is not the same')

        result = davies_method([10], [6, 3, 1], [0, 0, 0] , [2, 4, 6], limit=1000, accuracy = 0.0001, sigma=0)[0]
        self.assertAlmostEqual(result,  0.033357, places=6, msg='Value is not the same')
  
        result = davies_method([30], [6, 3, 1], [0, 0, 0] , [2, 4, 6], limit=1000, accuracy = 0.0001, sigma=0)[0]
        self.assertAlmostEqual(result,  0.580446, places=6, msg='Value is not the same')

        result = davies_method([80], [6, 3, 1], [0, 0, 0] , [2, 4, 6], limit=1000, accuracy = 0.0001, sigma=0)[0]
        self.assertAlmostEqual(result,  0.991283, places=6, msg='Value is not the same')

        result = davies_method([20], [7, 3], [6, 2] , [6, 2], limit=1000, accuracy = 0.0001, sigma=0)[0]
        self.assertAlmostEqual(result,  0.006125, places=6, msg='Value is not the same')

        result = davies_method([100], [7, 3], [6, 2] , [6, 2], limit=1000, accuracy = 0.0001, sigma=0)[0]
        self.assertAlmostEqual(result,  0.591339, places=6, msg='Value is not the same')

        result = davies_method([200], [7, 3], [6, 2] , [6, 2], limit=1000, accuracy = 0.0001, sigma=0)[0]
        self.assertAlmostEqual(result,  0.977914, places=6, msg='Value is not the same')

        result = davies_method([10], [7, 3], [6, 2] , [1, 1], limit=1000, accuracy = 0.0001, sigma=0)[0]
        self.assertAlmostEqual(result,  0.045126, places=6, msg='Value is not the same')

        result = davies_method([60], [7, 3], [6, 2] , [1, 1], limit=1000, accuracy = 0.0001, sigma=0)[0]
        self.assertAlmostEqual(result,  0.592431, places=6, msg='Value is not the same')

        result = davies_method([150], [7, 3], [6, 2] , [1, 1], limit=1000, accuracy = 0.0001, sigma=0)[0]
        self.assertAlmostEqual(result,  0.977648, places=6, msg='Value is not the same')

        result = davies_method([45], [6, 3, 1, 12, 6, 2], [0, 0, 0, 0, 0, 0] , [6, 4, 2, 2, 4, 6], limit=1000, accuracy = 0.0001, sigma=0)[0]
        self.assertAlmostEqual(result,  0.010950, places=6, msg='Value is not the same')

        result = davies_method([120], [6, 3, 1, 12, 6, 2], [0, 0, 0, 0, 0, 0] , [6, 4, 2, 2, 4, 6], limit=1000, accuracy = 0.0001, sigma=0)[0]
        self.assertAlmostEqual(result,  0.654735, places=6, msg='Value is not the same')

        result = davies_method([210], [6, 3, 1, 12, 6, 2], [0, 0, 0, 0, 0, 0] , [6, 4, 2, 2, 4, 6], limit=1000, accuracy = 0.0001, sigma=0)[0]
        self.assertAlmostEqual(result,  0.984606, places=6, msg='Value is not the same')

    




if __name__ == '__main__':
    unittest.main()