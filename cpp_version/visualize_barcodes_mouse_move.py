import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import argparse
import pdb
import math

parser = argparse.ArgumentParser(description='test')
parser.add_argument('--input', type=str, help='input barcodes file')

args = parser.parse_args()

file_name = args.input
print("file_name=",file_name)

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

x_range = m.shape[0]-1
y_range = m.shape[1]-1

m_max = -100000000
m_min = 100000000

for i in range(m.shape[0]):
    for j in range(m.shape[1]):
        for k in range(i,m.shape[0]):
            for l in range(j,m.shape[1]):
                if(i!=k and j!=l and m[i,j,k,l]<m_min):
                    m_min = m[i,j,k,l]
                if(i!=k and j!=l and m[i,j,k,l]>m_max):
                    m_max = m[i,j,k,l]

m_max_final = max(m_max, -m_min)

data = np.ones((x_range,y_range)) * np.nan

fig, ax = plt.subplots(1, 1, figsize=(5,5),tight_layout=True)
for x in range(x_range+2):
    ax.axvline(x, lw=1, color='k', zorder=5,linestyle=":",alpha=0.3) # Horizontal

for y in range(y_range+2):
    ax.axhline(y, lw=1, color='k', zorder=5,linestyle=":",alpha=0.3) #) # Ordinate

po_annotation = []

bbox1 = dict(boxstyle="round", fc='red', alpha=0.3)
bbox2 = dict(boxstyle="round", fc='blue', alpha=0.3)

for i in range(m.shape[0]):
    for j in range(m.shape[1]):
        for k in range(i,m.shape[0]):
            for l in range(j,m.shape[1]):
                if(m[i,j,k,l]):
                    print("m(({0},{1}),({2},{3}))={4}".format(i,j,k,l,int(m[i,j,k,l])))
                if(m[i,j,k,l]>0 and i!=k and j!=l):
                    min_value = min(k-i,l-j)/min(x_range,y_range)
                    wid_value = m[i,j,k,l] / m_max_final
                    ax.plot([i, k], [j, l],'bv', linestyle="-",alpha=min_value, linewidth=2*wid_value)
                    show_text = "m("+str(i)+","+str(j)+ ")(" + str(k)+","+str(l)+")="+str(int(m[i,j,k,l]))
                    annotation = plt.annotate(show_text,xy=((i+k)/2,(j+l)/2), xytext=(-15, 15), textcoords='offset points',bbox=bbox2, size=15)
                    annotation.set_visible(False)
                    po_annotation.append([[i,j,k,l,m[i,j,k,l]],annotation])

                if(m[i,j,k,l]<0 and i!=k and j!=l):
                    min_value = min(k-i,l-j)/min(x_range,y_range)
                    wid_value = -m[i,j,k,l] / m_max_final
                    ax.plot([i, k], [j, l], 'rv', linestyle="-",alpha=min_value, linewidth=2*wid_value)
                    show_text = "m("+str(i)+","+str(j)+ ")(" + str(k)+","+str(l)+")="+str(int(m[i,j,k,l]))
                    annotation = plt.annotate(show_text,xy=((i+k)/2,(j+l)/2), xytext=(-15, 15), textcoords='offset points',bbox=bbox1, size=15)
                    annotation.set_visible(False)
                    po_annotation.append([[i,j,k,l,m[i,j,k,l]],annotation])



def on_move(event):
    xdata = event.xdata
    ydata = event.ydata
    if(xdata==None or ydata==None):
        return
    for points, annotation in po_annotation:
        [i,j,k,l,m_value] = points
        vec1 = [xdata-i,ydata-j]
        vec2 = [k-i,l-j]
        vec1 = vec1/np.linalg.norm(vec1)
        vec2 = vec2/np.linalg.norm(vec2)

        ang = np.arccos(np.dot(vec1,vec2))
        ang = ang*180/math.pi
        
        if xdata>=i and xdata<=k and ydata>=j and ydata<=l and ang<1:
            annotation.set_visible(True)
        else:
            annotation.set_visible(False)

    plt.draw()


on_move_id = fig.canvas.mpl_connect('motion_notify_event', on_move)

plt.show()
