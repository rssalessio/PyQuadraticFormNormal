import unittest
from PyQuadraticFormNormal import davies_method


class TestDaviesMethod(unittest.TestCase):
    def test_basic(self):
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
        N = 10000
        result = davies_method([9600], [1] * N, [0] * N, [1] * N, sigma=0)[0]
        self.assertEqual(result, 0.0021031643039009507, 'Value is not the same')

        result = davies_method([10600], [1] * N, [0] * N, [1] * N, sigma=0)[0]
        self.assertEqual(result, 0.9999844082691355, 'Value is not the same')

if __name__ == '__main__':
    unittest.main()