import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import re
import os
#Potential distribution's equation 
def Pot (h,c,e):
    return c*(h**e)

print("\n")

#Clean Plot style
sns.set()
sns.set_style("ticks")

cwd = os.getcwd()
indir = os.getcwd() + '/Output/Debug'
outdir = os.getcwd() + '/Output/Plots'
try: 
    os.mkdir(outdir)
except: 
    try: 
        os.chdir(outdir)
    except:
        OSError
        print("Output directory %s couldn't be created. \n" % outdir)
    else:
        print("Output directory %s already existed \n" % outdir)
        os.chdir(cwd)
    
try:        
            os.chdir(indir)
            nodes = pd.read_csv('Node_position.txt')
            indexes = list(nodes["index"])
            x_node = list(nodes["x"])
            y_node = list(nodes["y"])
            for iter in range(len(indexes)):
                print(iter)
                #interface = pd.read_csv('Interface_'+str(indexes[iter])+'.txt')
                boundary = pd.read_csv('boundary_stencil_'+str(indexes[iter])+'.txt')
                stencil = pd.read_csv('stencil_'+str(indexes[iter])+'.txt')
                fig, ax = plt.subplots()
                plt.title("Node "+str(indexes[iter]))
                ax.plot(boundary['x'], boundary['y'],color = 'red', label = 'boundary')
                ax.scatter(stencil['x'], stencil['y'], s = 50, color = 'black', marker = "x", label = "stencil")
                for i, txt in enumerate(stencil['index']):
                    ax.annotate(txt, (stencil['x'][i], stencil['y'][i]))
                #ax.scatter(interface['x'], interface['y'], s = 50, color = 'blue', marker = "x", label = "interface")
                #for i, txt in enumerate(interface['index']):
                    #ax.annotate(txt, (interface['x'][i], interface['y'][i]))
                ax.scatter(x_node[iter], y_node[iter], s = 50, color = 'cyan', marker = "x", label = "node")
                plt.legend()
                os.chdir(outdir)
                plt.savefig("Node_"+str(indexes[iter])+".png", dpi = 200)
                os.chdir(indir)
                plt.close()
    
except:
    print("Input directory %s couldn't be plotted. \n" %indir)
else:
    print("Input directory %s was succesfully plotted.\n" %indir)


