from generalized_chi_squared import davies_method
import numpy as np

print(1)
print(davies_method(np.zeros((2)), np.zeros((2)), np.zeros(2, dtype=np.int32)))
print('-----------')
print(2)
print(davies_method((0,0,0), (0,0,0), (0,0,0)))