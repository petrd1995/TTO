import scipy as sc
import numpy as np

a = np.zeros(10,10)
a[0:2, 0:2] = np.ones((2, 2))
a[2:4, 0:2] = 2*np.ones((2, 2))
a[0:2, 2:4] = 3*np.ones((2, 2))
a[2:4, 2:4] = 4*np.ones((2, 2))
a
