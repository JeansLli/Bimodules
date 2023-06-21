import numpy as np
import gudhi as gd  
import pickle as pickle
from pylab import *
from sklearn.neighbors import KernelDensity
import pdb

fname = "../data/function_rips_GMM_points.txt"
band = 0.2


pointcloud = np.loadtxt(fname)

plt.scatter(pointcloud[:, 0], pointcloud[:, 1])
plt.show()

xval = np.arange(0, 10, 0.05)
yval = np.arange(0, 10, 0.05)
nx = len(xval)
ny = len(yval)



kde = KernelDensity(kernel = 'linear', bandwidth = band).fit(pointcloud)
positions = np.array([[u, v] for u in xval for v in yval])

filt_values = -kde.score_samples(X = positions)


cc_density_crater = gd.CubicalComplex(
    dimensions = [nx ,ny], 
    top_dimensional_cells = filt_values
)

print("type: ",type(cc_density_crater))
print("dimension: ",cc_density_crater.dimension())
print("num_simplices:",cc_density_crater.num_simplices())

BarCodes_Rips0 = cc_density_crater.persistence()



for i in range(len(BarCodes_Rips0)):
    print(BarCodes_Rips0[i])


