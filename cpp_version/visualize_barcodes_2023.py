# 
#                            Signed barcodes visualization 
#                          Code by Jingyi Li and Steve Oudot
#                         Copyright 2022, all rights reserved
#

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as pch
import mpl_interactions.ipyplot as iplt
import matplotlib
import argparse
import math
import pdb


def visualize_barcodes_segment(m,barcodes):
    x_range = m.shape[0]-1
    y_range = m.shape[1]-1

    m_max = -100000000
    m_min = 100000000
    for i in range(m.shape[0]):
        for j in range(m.shape[1]):
            for k in range(i,m.shape[0]):
                for l in range(j,m.shape[1]):
#                    if(i!=k and j!=l and m[i,j,k,l]<m_min):
                    if(m[i,j,k,l]<m_min):
                        m_min = m[i,j,k,l]
#                    if(i!=k and j!=l and m[i,j,k,l]>m_max):
                    if(m[i,j,k,l]>m_max):
                        m_max = m[i,j,k,l]

    print("m_max=",m_max)
    print("m_min=",m_min)

    data = np.ones((x_range,y_range)) * np.nan

    fig, ax = plt.subplots(1, 1, figsize=(5,5),tight_layout=True)
    for x in range(x_range+2):
        ax.axvline(x, lw=1, color='k', zorder=5,linestyle=":",alpha=0.3) # Horizontal

    for y in range(y_range+2):
        ax.axhline(y, lw=1, color='k', zorder=5,linestyle=":",alpha=0.3) #) # Ordinate


##########




    # draw grey scale for Hilbert funtion
    greys = np.zeros((x_range,y_range))
    for i in range(x_range):
        for j in range(y_range):
            grey_value = 0
            for barcode in barcodes:
                if i>=barcode[0] and i<barcode[2] and j>=barcode[1] and j<barcode[3]:
                    grey_value += barcode[4]
            greys[i][j]=grey_value


    grey_min = greys.min()
    grey_max = greys.max()
    print("grey_min=",grey_min)
    print("grey_max=",grey_max)
    for i in range(x_range):
        for j in range(y_range):
            ######## for x-y switch
            rect = pch.Rectangle(xy=(i,j),width=1, height=1, color=str(1-math.log(1+greys[i][j])/math.log(1+grey_max)/1.7))
#            rect = pch.Rectangle(xy=(i,j),width=1, height=1, color=str(0.5+0.5*math.exp(-greys[i][j])))
#            rect = pch.Rectangle(xy=(i,j),width=1, height=1, color=str(1-greys[i][j]/2/grey_max))
            ########
            ax.add_patch(rect)
            #plt.text(i+1/2,j+1/2,int(greys[j][i]))


        
    po_annotation = []
    bbox1 = dict(boxstyle="round", fc='red', alpha=0.3)
    bbox2 = dict(boxstyle="round", fc='blue', alpha=0.3)


    for i in range(m.shape[0]):
        for j in range(m.shape[1]):
            for k in range(i,m.shape[0]):
                for l in range(j,m.shape[1]):
#                    if(m[i,j,k,l]):
#                        print("m(({0},{1}),({2},{3}))={4}".format(i,j,k,l,int(m[i,j,k,l])))
                    alpha_value = (1+min(k-i,l-j))/(1+min(x_range,y_range))
#                    alpha_value = 1
                    wid_value = abs(m[i,j,k,l])

                    if(m[i,j,k,l]>0):
#                    if(m[i,j,k,l]>0 and i!=k and j!=l):

                        ########## for x-y switch
                        ax.plot([i, k], [j, l],'bv', linestyle="-",alpha=alpha_value,linewidth=1)
#                        ax.plot([i, k], [j, l],'bv', linestyle="-",alpha=alpha_value,linewidth=wid_value)
                        show_text = "m("+str(i)+","+str(j)+ ")(" + str(k)+","+str(l)+")="+str(int(m[i,j,k,l]))
                        
                        annotation = plt.annotate(show_text,xy=((i+k)/2,(j+l)/2), xytext=(-15, 15), textcoords='offset points',bbox=bbox2, size=15)
                        annotation.set_visible(False)
                        po_annotation.append([[i,j,k,l,m[i,j,k,l]],annotation])
                        ##########

                    elif(m[i,j,k,l]<0):
#                    elif(m[i,j,k,l]<0 and i!=k and j!=l):

                        ########## for x-y switch
                        ax.plot([i, k], [j, l], 'rv', linestyle="-",alpha=alpha_value,linewidth=1)
#                        ax.plot([i, k], [j, l], 'rv', linestyle="-",alpha=alpha_value,linewidth=wid_value)
                        show_text = "m("+str(i)+","+str(j)+ ")(" + str(k)+","+str(l)+")="+str(int(m[i,j,k,l]))
                        annotation = plt.annotate(show_text,xy=((i+k)/2,(j+l)/2), xytext=(-15, 15), textcoords='offset points',bbox=bbox1, size=15)
                        annotation.set_visible(False)
                        po_annotation.append([[i,j,k,l,m[i,j,k,l]],annotation])
    
#### for mouse interaction

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
            if vec2 != [0,0]:
                vec2 = vec2/np.linalg.norm(vec2)
            
            ang = np.arccos(np.dot(vec1,vec2))
            ang = ang*180/math.pi
            
            if xdata>=i and xdata<=k and ydata>=j and ydata<=l and ang<1:
                annotation.set_visible(True)
            else:
                annotation.set_visible(False)
                
        plt.draw()


    on_move_id = fig.canvas.mpl_connect('motion_notify_event', on_move)

### end mouse interaction


            
    def line_y(x,angle, offset):
        return math.tan(angle/180*math.pi)*x+offset

    def line_x(y,angle,offset):
        if angle == 0:
            return 0
        else:
            return (y-offset)/(math.tan(angle/180*math.pi))


    # define the function
    def segment(x, angle, offset):
        y = math.tan(angle/180*math.pi)*x+offset
#        return (x[y<=y_range], y[y<=y_range])
#        y[y>y_range] = y_range
        return y


    x = np.linspace(-0.5, x_range+0.5, 1000)


    
    def bars_b(angle, offset): 
        bars_bb=[]
        lines = {}
        for b in barcodes:
            if b[1] <= line_y(b[0],angle,offset):
                line_start = [b[0],line_y(b[0],angle,offset)]
                if b[3] <= line_y(b[2],angle,offset):
                    line_end = [line_x(b[3],angle,offset),b[3]]
                else:
                    line_end = [b[2],line_y(b[2],angle,offset)]

                if line_end[0]>=line_start[0] and line_end[1]>=line_start[1]:
                    if not (line_start[0], line_start[1], line_end[0], line_end[1]) in lines:
                        lines[(line_start[0], line_start[1], line_end[0], line_end[1])] = b[4]
                    else:
                        lines[(line_start[0], line_start[1], line_end[0], line_end[1])] += b[4]
                        if lines[(line_start[0], line_start[1], line_end[0], line_end[1])] == 0:
                            del lines[(line_start[0], line_start[1], line_end[0], line_end[1])]

            elif b[0] <= line_x(b[1],angle,offset):
                line_start = [line_x(b[1],angle,offset),b[1]]
                if b[3] <= line_y(b[2],angle,offset):
                    line_end = [line_x(b[3],angle,offset),b[3]]
                else:
                    line_end = [b[2],line_y(b[2],angle,offset)]

                if line_end[0]>=line_start[0] and line_end[1]>=line_start[1]:
                    if not (line_start[0], line_start[1], line_end[0], line_end[1]) in lines:
                        lines[(line_start[0], line_start[1], line_end[0], line_end[1])] = b[4]
                    else:
                        lines[(line_start[0], line_start[1], line_end[0], line_end[1])] += b[4]
                        if lines[(line_start[0], line_start[1], line_end[0], line_end[1])] == 0:
                            del lines[(line_start[0], line_start[1], line_end[0], line_end[1])]

        for l in lines:
            assert lines[l] > 0
            
        for l in lines:
            bars_bb.append(math.sqrt(l[0]*l[0] + (l[1]-offset)*(l[1]-offset)))
#            bars_bb.append(l[0])
        
        return bars_bb

    def bars_d(x,angle, offset): 
        if x==[]:
            return []
        bars_dd=[]
        lines = {}
        for b in barcodes:
            if b[1] <= line_y(b[0],angle,offset):
                line_start = [b[0],line_y(b[0],angle,offset)]
                if b[3] <= line_y(b[2],angle,offset):
                    line_end = [line_x(b[3],angle,offset),b[3]]
                else:
                    line_end = [b[2],line_y(b[2],angle,offset)]

                if line_end[0]>=line_start[0] and line_end[1]>=line_start[1]:
                    if not (line_start[0], line_start[1], line_end[0], line_end[1]) in lines:
                        lines[(line_start[0], line_start[1], line_end[0], line_end[1])] = b[4]
                    else:
                        lines[(line_start[0], line_start[1], line_end[0], line_end[1])] += b[4]
                        if lines[(line_start[0], line_start[1], line_end[0], line_end[1])] == 0:
                            del lines[(line_start[0], line_start[1], line_end[0], line_end[1])]
                        
            elif b[0] <= line_x(b[1],angle,offset):
                line_start = [line_x(b[1],angle,offset),b[1]]
                if b[3] <= line_y(b[2],angle,offset):
                    line_end = [line_x(b[3],angle,offset),b[3]]
                else:
                    line_end = [b[2],line_y(b[2],angle,offset)]

                if line_end[0]>=line_start[0] and line_end[1]>=line_start[1]:
                    if not (line_start[0], line_start[1], line_end[0], line_end[1]) in lines:
                        lines[(line_start[0], line_start[1], line_end[0], line_end[1])] = b[4]
                    else:
                        lines[(line_start[0], line_start[1], line_end[0], line_end[1])] += b[4]
                        if lines[(line_start[0], line_start[1], line_end[0], line_end[1])] == 0:
                            del lines[(line_start[0], line_start[1], line_end[0], line_end[1])]

        for l in lines:
            assert lines[l] > 0
            
        for l in lines:
            bars_dd.append(math.sqrt(l[2]*l[2] + (l[3]-offset)*(l[3]-offset)))
#            bars_dd.append(l[2])

        return bars_dd
    


    ### Set the range of angle and offset
    angle=np.linspace(0, 90, 91)
    offset=np.linspace(-y_range-0.5, y_range+0.5, 1+(int)(10*(2*y_range+1)))

#    print(type(x))
#    print(x.size)
#    print([item for item in x if segment(item,angle,offset) <= y_range])
    
    controls = iplt.plot(x, segment, angle=angle, offset=offset,ax=ax,color="green", xlim=[-0.5, x_range+0.5], ylim=[-0.5, y_range+0.5])
    fig2, ax2 = plt.subplots(1, 1, figsize=(5,5),tight_layout=True)
    ax2.plot([0, x_range+1],[0, x_range+1], linestyle="-")
    _ = iplt.scatter(bars_b, bars_d, controls=controls,ax=ax2, xlim='auto', ylim='auto')
    
    
    #ax.imshow(data, interpolation='none', extent=[-1,x_range+1, -1, y_range+1], zorder=0)
    #plt.axis('off')
    plt.show()
    #ax.axis('off')


parser = argparse.ArgumentParser(description='command-line parser')
parser.add_argument('input', type=str, help='input barcode file')

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
barcodes = []
while(line!=''):
    barcode=line.split()
    m[int(barcode[0]),int(barcode[1]),int(barcode[2]),int(barcode[3])] = int(barcode[4])
    barcodes.append([int(barcode[0]),int(barcode[1]),int(barcode[2]),int(barcode[3]),int(barcode[4])])
    line = f.readline()

f.close()

visualize_barcodes_segment(m,barcodes)



