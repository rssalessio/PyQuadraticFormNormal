import numpy as np
from compquadform import davies

print(1)
print(davies([0], np.zeros((2)), np.zeros((2)), np.zeros(2, dtype=np.int32)))
print('-----------')
print(2)
print(davies([0], (0,0,0), (0,0,0), (0,0,0)))