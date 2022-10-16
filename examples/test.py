from multiprocessing.sharedctypes import Value
from QuadraticFormDistributions import davies_method

#print(davies_method([0], [0], [0], [0]))
#print(davies_method([1], [1], [1], [1]))
print(davies_method([2], [1], [0], [1], sigma=0))
print(davies_method([2], [1], [0], [1], sigma=0, accuracy=1e-6))
print(davies_method([1], [1], [0], [1], sigma=0))
print(davies_method([1], [0], [0], [0], sigma=1))
# print(davies_method([.5], [1], [0], [1], sigma=0))
# print(davies_method([.1], [1], [0], [1], sigma=0))
# print(davies_method([0.1], [1], [0], [1]))
# print(davies_method([0.1], [10], [0], [1]))
try:
    print(davies_method([], [0], [0], [0]))
except ValueError as e:
    print(e)

try:
    print(davies_method([0], [0], [], [0]))
except ValueError as e:
    print(e)


try:
    print(davies_method([0], [0,1], [0], [0]))
except ValueError as e:
    print(e)

try:
    print(davies_method([0], [0], [0], []))
except ValueError as e:
    print(e)

try:
    print(davies_method([0], [0], [0], [0], accuracy=0))
except ValueError as e:
    print(e)
