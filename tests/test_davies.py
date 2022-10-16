import unittest
from QuadraticFormDistributions import davies_method


class TestDaviesMethod(unittest.TestCase):
    def test_basic(self):
        result, trace, fault = davies_method([2], [1], [0], [1], sigma=0)
        self.assertEqual(len(result), 1, "Should return only one result")
        self.assertEqual(fault, 0, "Fault needs to be zero")
        self.assertEqual(len(trace), 7, "Trace needs to have 7 elements")

        result, trace, fault = davies_method([2, 1], [1], [0], [1], sigma=0)
        self.assertEqual(len(result), 2, "Should return two results")
        self.assertEqual(fault, 0, "Fault needs to be zero")
        self.assertEqual(len(trace), 7, "Trace needs to have 7 elements")

    def test_value_chi_one(self):
        """ Test from a chi^2(1) """
        result, trace, fault = davies_method([2], [1], [0], [1], sigma=0)
        self.assertEqual(result[0], 0.8426997788273605, "Incorrect result")
        self.assertEqual(fault, 0, "Fault needs to be 0")

        result, trace, fault = davies_method([2], [1], [0], [1], sigma=0, accuracy=1e-6)
        self.assertEqual(result[0], 0.8427006907020175, "Incorrect result")
        self.assertEqual(fault, 0, "Fault needs to be 0")

        result, trace, fault = davies_method([1], [1], [0], [1], sigma=0)
        self.assertEqual(result[0], 0.6826886153722268, "Incorrect result")
        self.assertEqual(fault, 0, "Fault needs to be 0")

        result, trace, fault = davies_method([.5], [1], [0], [1], sigma=0)
        self.assertEqual(result[0], 0.5204991504663418, "Incorrect result")
        self.assertEqual(fault, 0, "Fault needs to be 0")

        result, trace, fault = davies_method([.1], [1], [0], [1], sigma=0)
        self.assertEqual(result[0], 0.24816987371489274, "Incorrect result")
        self.assertEqual(fault, 0, "Fault needs to be 0")

    def test_value_chi_two(self):
        """ Test from a chi^2(2) """
        result, trace, fault = davies_method([2], [1], [0], [2], sigma=0)
        self.assertEqual(result[0], 0.6321195748892512, "Incorrect result")
        self.assertEqual(fault, 0, "Fault needs to be 0")

        result, trace, fault = davies_method([2], [1, 1], [0, 0], [1, 1], sigma=0)
        self.assertEqual(result[0], 0.6321195748892512, "Incorrect result")
        self.assertEqual(fault, 0, "Fault needs to be 0")

    def test_value_normal(self):
        """ Test from a normal distribution N(0, 1) """
        result, trace, fault = davies_method([1], [0], [0], [0], sigma=1)
        self.assertEqual(result[0], 0.841344541972314, "Incorrect result")
        self.assertEqual(fault, 0, "Fault needs to be 0")

        result, trace, fault = davies_method([-1], [0], [0], [0], sigma=1)
        self.assertEqual(result[0], 0.15865545802768594, "Incorrect result")
        self.assertEqual(fault, 0, "Fault needs to be 0")


if __name__ == '__main__':
    unittest.main()