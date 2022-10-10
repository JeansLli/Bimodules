import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as pch
import mpl_interactions.ipyplot as iplt
import matplotlib
import argparse
import math
import pdb

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


def visualize_barcodes_lines_width(m):
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

    print("m_max=",m_max)
    print("m_min=",m_min)

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
                        wid_value = m[i,j,k,l] / m_max
                        ax.plot([i, k], [j, l],'bv', linestyle="-",alpha=min_value,linewidth=m[i,j,k,l])
                        print("text position ",(i+k)/2, (j+l)/2)
                        plt.text((i+k)/2,(j+l)/2,str(m[i,j,k,l]))
                    if(m[i,j,k,l]<0 and i!=k and j!=l):
                        min_value = min(k-i,l-j)/min(x_range,y_range)
                        wid_value = m[i,j,k,l] / m_min
                        ax.plot([i, k], [j, l], 'rv', linestyle="-",alpha=min_value,linewidth=m[i,j,k,l])

    #ax.imshow(data, interpolation='none', extent=[-1,x_range+1, -1, y_range+1], zorder=0)
    #plt.axis('off')
    plt.show()
    #ax.axis('off')


def visualize_barcodes_grey(m,barcodes):
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

    print("m_max=",m_max)
    print("m_min=",m_min)

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
                        wid_value = m[i,j,k,l] / m_max

                        ########## for x-y switch
                        #ax.plot([i, k], [j, l],'bv', linestyle="-",alpha=min_value,linewidth=m[i,j,k,l])
                        ax.plot([j, l],[i, k], 'bv', linestyle="-",alpha=min_value,linewidth=m[i,j,k,l])
                        ##########

                    if(m[i,j,k,l]<0 and i!=k and j!=l):
                        min_value = min(k-i,l-j)/min(x_range,y_range)
                        wid_value = m[i,j,k,l] / m_min

                        ########## for x-y switch
                        #ax.plot([i, k], [j, l], 'rv', linestyle="-",alpha=min_value,linewidth=m[i,j,k,l])
                        ax.plot([j, l],[i, k], 'rv', linestyle="-",alpha=min_value,linewidth=m[i,j,k,l])
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


    grey_mim = greys.min()
    grey_max = greys.max()
    
    for i in range(x_range):
        for j in range(y_range):
            ######## for x-y switch
            #rect = pch.Rectangle(xy=(i,j),width=1, height=1, color=str(1-greys[i][j]/grey_max))
            rect = pch.Rectangle(xy=(i,j),width=1, height=1, color=str(1-greys[j][i]/grey_max))
            ########
            ax.add_patch(rect)
            plt.text(i+1/2,j+1/2,greys[j][i])


    # define the function
    def segment(x, angle, offset):
        y = math.tan(angle/180*math.pi)*x+offset
        y[y<=0] = 0
        y[y>=y_range] = y_range
        return y


    x = np.linspace(0, x_range, 10000)

    angle = np.linspace(0, 90)
    offset = np.linspace(-20, 20)

    #fig, ax = plt.subplots()

    controls = iplt.plot(x, segment, angle=angle, offset=offset)
                        

    #ax.imshow(data, interpolation='none', extent=[-1,x_range+1, -1, y_range+1], zorder=0)
    #plt.axis('off')
    plt.show()
    #ax.axis('off')


def visualize_barcodes_segment(m,barcodes):
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

    print("m_max=",m_max)
    print("m_min=",m_min)

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
                        wid_value = m[i,j,k,l] / m_max

                        ########## for x-y switch
                        #ax.plot([i, k], [j, l],'bv', linestyle="-",alpha=min_value,linewidth=m[i,j,k,l])
                        ax.plot([j, l],[i, k], 'bv', linestyle="-",alpha=min_value,linewidth=m[i,j,k,l])
                        ##########

                    if(m[i,j,k,l]<0 and i!=k and j!=l):
                        min_value = min(k-i,l-j)/min(x_range,y_range)
                        wid_value = m[i,j,k,l] / m_min

                        ########## for x-y switch
                        #ax.plot([i, k], [j, l], 'rv', linestyle="-",alpha=min_value,linewidth=m[i,j,k,l])
                        ax.plot([j, l],[i, k], 'rv', linestyle="-",alpha=min_value,linewidth=m[i,j,k,l])
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


    grey_mim = greys.min()
    grey_max = greys.max()
    print("grey_max=",grey_max)
    for i in range(x_range):
        for j in range(y_range):
            ######## for x-y switch
            #rect = pch.Rectangle(xy=(i,j),width=1, height=1, color=str(1-greys[i][j]/grey_max))
            rect = pch.Rectangle(xy=(i,j),width=1, height=1, color=str(1-greys[j][i]/grey_max))
            ########
            ax.add_patch(rect)
            plt.text(i+1/2,j+1/2,int(greys[j][i]))
    

    def line_y(x,angle, offset):
        return math.tan(angle/180*math.pi)*x+offset

    def line_x(y,angle,offset):
        return (y-offset)/(math.tan(angle/180*math.pi))






    # define the function
    def segment(x, angle, offset):
        y = math.tan(angle/180*math.pi)*x+offset
        y[y<=0] = 0
        y[y>=y_range] = y_range
        return y



    x = np.linspace(0, x_range, 10000)


    ######### project 2D to 1D. No, it's wrong now.
    #def bars_b(angle, offset): 
    #    bar_b=[]
    #    for b in barcodes:
    #        if b[2]>line_y(b[1],angle,offset) and b[0]<line_y(b[3],angle,offset):
    #            bar_b.append(b[1])
    #    print("bar_b=",bar_b)
    #    return bar_b
    #
    #def bars_d(x,angle, offset): 
    #    bar_d=[]
    #    for b in barcodes:
    #        print("angle=",angle)
    #        print("k=",math.tan(angle/180*math.pi))
    #        print("offset=",offset)
    #        print("barcode=",b)
    #        if b[2]>line_y(b[1],angle,offset) and b[0]<line_y(b[3],angle,offset):
    #            if(angle==0): 
    #                bar_d.append(b[3])
    #                print("1death=",b[3])
    #            else:
    #                bar_d.append(line_x(b[2],angle,offset))
    #                print("2death=",line_x(b[2],angle,offset))
    #    print("bar_d=",bar_d)
    #    return bar_d
    ####################
    def bars_b(angle, offset): 
        bars_bb=[]
        lines_pos = []
        lines_neg = []
        print("angle=",angle)
        print("k=",math.tan(angle/180*math.pi))
        print("offset=",offset)
        for b in barcodes:
            if b[0] <= line_y(b[1],angle,offset):
                line_start = [b[1],line_y(b[1],angle,offset)]
                if b[2] <= line_y(b[3],angle,offset):
                    line_end = [line_x(b[2],angle,offset),b[2]]
                else:
                    line_end = [b[3],line_y(b[3],angle,offset)]

                if line_end[1]>=line_start[1] and line_end[0]>=line_start[0]:
                    if b[4]>0:
                        lines_pos.append([line_start[0],line_start[1],line_end[0],line_end[1],b[4]])
                    else:
                        lines_neg.append([line_start[0],line_start[1],line_end[0],line_end[1],b[4]])

        for line_p in lines_pos:
            for line_n in lines_neg:

                if line_p[:4]== line_n[:4] and line_p[4]+line_n[4]<=0:
                    lines_pos.remove(line_p)
                    break

        for line_p in lines_pos:
            bars_bb.append(line_p[0])
        
        return bars_bb

    def bars_d(x,angle, offset): 
        if x==[]:
            return []
        bars_dd=[]
        lines_pos = []
        lines_neg = []
        for b in barcodes:
            if b[0] <= line_y(b[1],angle,offset):
                line_start = [b[1],line_y(b[1],angle,offset)]
                if b[2] <= line_y(b[3],angle,offset):
                    line_end = [line_x(b[2],angle,offset),b[2]]
                else:
                    line_end = [b[3],line_y(b[3],angle,offset)]
                if line_end[1]>=line_start[1] and line_end[0]>=line_start[0]:
                    if b[4]>0:
                        lines_pos.append([line_start[0],line_start[1],line_end[0],line_end[1],b[4]])
                    else:
                        lines_neg.append([line_start[0],line_start[1],line_end[0],line_end[1],b[4]])

        for line_p in lines_pos:
            for line_n in lines_neg:
                if line_p[:4]== line_n[:4] and (line_p[4]+line_n[4])<=0:
                    lines_pos.remove(line_p)
                    break

        for line_p in lines_pos:
            bars_dd.append(line_p[2])

        return bars_dd




    angle=np.linspace(0, 90)
    offset=np.linspace(-20, 20)

    controls = iplt.plot(x, segment, angle=angle, offset=offset,ax=ax,color="yellow")
    fig2, ax2 = plt.subplots(1, 1, figsize=(5,5),tight_layout=True)
    ax2.plot([0, 20],[0, 20], linestyle="-")
    _ = iplt.scatter(bars_b, bars_d, controls=controls,ax=ax2)
    

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
barcodes = []
while(line!=''):
    barcode=line.split()
    m[int(barcode[0]),int(barcode[1]),int(barcode[2]),int(barcode[3])] = int(barcode[4])
    barcodes.append([int(barcode[0]),int(barcode[1]),int(barcode[2]),int(barcode[3]),int(barcode[4])])
    line = f.readline()

f.close()

#visualize_barcodes(m)
#visualize_barcodes_lines_intensity(m)
#visualize_barcodes_points_intensity(m)
#visualize_barcodes_lines_width(m)
#visualize_barcodes_grey(m,barcodes)
visualize_barcodes_segment(m,barcodes)



