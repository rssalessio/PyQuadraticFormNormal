from multiprocessing.sharedctypes import Value
from PyQuadraticFormNormal import davies_method
from momentchi2 import hbe

print(hbe([1, 1], [1]))
print(davies_method([1,1,1], [1], [0], [2], sigma=0))
print(davies_method([1], [1], [0], [2], sigma=0)[0])
print(davies_method([1], [1], [0], [1], sigma=0)[0])
print(davies_method([1], [1, 1], [0, 0], [1, 1], sigma=0)[0])
print(davies_method([1], [1, 1,1], [0, 0,0], [1, 1,1], sigma=0)[0])
print(davies_method([1], [1, 1,1,1,1], [0, 0,0,0,0], [1,1, 1,1,1], sigma=0)[0])
print(davies_method([1], [1, 1,1,1,1,1], [0,0, 0,0,0,0], [1,1,1, 1,1,1], sigma=0)[0])
print(davies_method([1], [1, 1,1,1,1,1,1], [0,0,0, 0,0,0,0], [1,1,1,1, 1,1,1], sigma=0)[0])
