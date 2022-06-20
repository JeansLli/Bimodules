import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import argparse

def visualize_barcodes(m):
    x_range = m.shape[0]-1
    y_range = m.shape[1]-1

    data = np.ones((x_range,y_range)) * np.nan

    fig, ax = plt.subplots(1, 1, figsize=(5,5),tight_layout=True)
    for x in range(x_range+1):
        ax.axvline(x, lw=1, color='k', zorder=5, linestyle=":",) # Horizontal
    for y in range(y_range+1):
        ax.axhline(y, lw=1, color='k', zorder=5, linestyle=":",) # Ordinate

    for i in range(m.shape[0]):
        for j in range(m.shape[1]):
            for k in range(i,m.shape[0]):
                for l in range(j,m.shape[1]):
                    if(m[i,j,k,l]):
                        print("m(({0},{1}),({2},{3}))={4}".format(i,j,k,l,int(m[i,j,k,l])))
                    if(m[i,j,k,l]>0):
                        if(i==k and j==l):
                            ax.plot([i, k], [j, l],'bo', linestyle="-")
                        elif(i==k):
                            ax.plot([i+0.08, k+0.08], [j, l],'bv', linestyle="-")
                        elif(j==l):
                            ax.plot([i, k], [j+0.08, l+0.08],'bv',linestyle="-")
                        else:
                            ax.plot([i, k], [j, l],'bv', linestyle="-")
                    if(m[i,j,k,l]<0):
                        if(i==k and j==l):
                            ax.plot([i, k], [j, l],'ro', linestyle="-")
                        elif(i==k):
                            ax.plot([i-0.08, k-0.08], [j, l],'rv', linestyle="-")
                        elif(j==l):
                            ax.plot([i, k], [j-0.08, l-0.08],'rv', linestyle="-")
                        else:
                            ax.plot([i, k], [j, l], 'rv', linestyle="-")

    #ax.imshow(data, interpolation='none', extent=[-1,x_range+1, -1, y_range+1], zorder=0)
    #plt.axis('off')
    plt.show()
    #ax.axis('off')


def visualize_barcodes_lines_intensity(m):
    x_range = m.shape[0]-1
    y_range = m.shape[1]-1

    data = np.ones((x_range,y_range)) * np.nan

    fig, ax = plt.subplots(1, 1, figsize=(5,5),tight_layout=True)
    for x in range(x_range+2):
        ax.axvline(x, lw=1, color='k', zorder=5,linestyle=":",alpha=0.3) # Horizontal

    for y in range(y_range+2):
        ax.axhline(y, lw=1, color='k', zorder=5,linestyle=":",alpha=0.3) #) # Ordinate

    for i in range(m.shape[0]):
        for j in range(m.shape[1]):
            for k in range(i,m.shape[0]):
                for l in range(j,m.shape[1]):
                    if(m[i,j,k,l]):
                        print("m(({0},{1}),({2},{3}))={4}".format(i,j,k,l,int(m[i,j,k,l])))
                    if(m[i,j,k,l]>0 and i!=k and j!=l):
                        min_value = min(k-i,l-j)/min(x_range,y_range)
                        ax.plot([i, k], [j, l],'bv', linestyle="-",alpha=min_value)
                    if(m[i,j,k,l]<0 and i!=k and j!=l):
                        min_value = min(k-i,l-j)/min(x_range,y_range)
                        ax.plot([i, k], [j, l], 'rv', linestyle="-",alpha=min_value)

    #ax.imshow(data, interpolation='none', extent=[-1,x_range+1, -1, y_range+1], zorder=0)
    #plt.axis('off')
    plt.show()
    #ax.axis('off')

def visualize_barcodes_points_intensity(m):
    x_range = m.shape[0]-1
    y_range = m.shape[1]-1

    data = np.ones((x_range,y_range)) * np.nan

    fig, ax = plt.subplots(1, 1, figsize=(5,5),tight_layout=True)
    for x in range(-x_range-1,x_range+2):
        ax.axvline(x, lw=1, color='k', zorder=5,linestyle=":",alpha=0.3) # Horizontal
    ax.axvline(0, lw=1, color='k', zorder=5,linestyle="-") # Horizontal

    for y in range(-y_range-1,y_range+2):
        ax.axhline(y, lw=1, color='k', zorder=5,linestyle=":",alpha=0.3) #) # Ordinate
    ax.axhline(0, lw=1, color='k', zorder=5,linestyle="-") # Horizontal

    for i in range(m.shape[0]):
        for j in range(m.shape[1]):
            for k in range(i,m.shape[0]):
                for l in range(j,m.shape[1]):
                    if(m[i,j,k,l]):
                        print("m(({0},{1}),({2},{3}))={4}".format(i,j,k,l,int(m[i,j,k,l])))
                    if(m[i,j,k,l]>0 and i!=k and j!=l):
                        min_value = min(k-i,l-j)/min(x_range,y_range)
                        ax.plot([k-i], [l-j],'bo', linestyle="-",alpha=min_value)
                    if(m[i,j,k,l]<0 and i!=k and j!=l):
                        min_value = min(k-i,l-j)/min(x_range,y_range)
                        ax.plot([i-k], [j-l], 'ro', linestyle="-",alpha=min_value)

    #ax.imshow(data, interpolation='none', extent=[-1,x_range+1, -1, y_range+1], zorder=0)
    #plt.axis('off')
    plt.show()
    #ax.axis('off')


parser = argparse.ArgumentParser(description='test')
parser.add_argument('--input', type=str, help='input barcodes file')

args = parser.parse_args()

file_name = args.input
print("file_name=",file_name)


#file_name="./result/dim_0_barcodes_1025.txt"
f = open(file_name)
line = f.readline()
grid_size = line.split()
x_range = int(grid_size[0])
y_range = int(grid_size[1])
print("grid size is",x_range+1,"*",y_range+1)

m = np.array(np.zeros((x_range+1,y_range+1,x_range+1,y_range+1)))

line = f.readline()
while(line!=''):
    barcode=line.split()
    m[int(barcode[0]),int(barcode[1]),int(barcode[2]),int(barcode[3])] = int(barcode[4])
    line = f.readline()

f.close()

visualize_barcodes(m)
visualize_barcodes_lines_intensity(m)
visualize_barcodes_points_intensity(m)




