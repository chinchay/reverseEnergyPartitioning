################################################################################
# python3
import os
print(os.getcwd())
os.chdir('/Users/chinchay/Documents/9_Git/reverseEnergyPartitioning/')
print(os.getcwd())

# import graphlab
################################################################################

x = [1.1, 2.1, 3.1, 4.1]
y = [162, 209, 221, 224]
from matplotlib import pyplot as plt
plt.plot(x, y, 'o')
