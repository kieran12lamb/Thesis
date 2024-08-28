import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches

class MutationVisualiser:

    def __init__(self, ref,data,metadata,Orfs):
        self.data = data
        self.metadata = metadata
        self.ref = ref
        self.Orfs = Orfs


    def lollipop(self):
        fig, axes = plt.subplots(figsize=(30, 5))
        palette = ["mediumorchid","orange", "brown", '#04BBEC','black','#E42824', "teal", "gold", "mediumblue",'grey','#A0CF63','#EEC4C4']
        start=self.Orfs.iloc[0].Start
        totalsize = 0
        axes.hlines(y=1, xmin=0, xmax=29903,color='#000000',linewidth=0.5)
        for i in range(len(self.Orfs.index)):
            size = (self.Orfs.iloc[i].End - self.Orfs.iloc[i].Start)
            rectangle = patches.Rectangle((start,0.5),width=size,height=1,facecolor=palette[i])
            axes.add_patch(rectangle)
            start =start+ size
            totalsize = totalsize+size
        axes.set_xlim(left = -(totalsize*0.05),right=totalsize+(totalsize*0.05))
        axes.set_ylim(0,2)
        for pos in self.data.POS:
            axes.axvline(x=pos, ymin=0.25, ymax=0.75,color='#2E5894',linewidth=0.5)
            
