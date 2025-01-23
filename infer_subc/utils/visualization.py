import napari
import itertools
import math

# matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import matplotlib.gridspec as gridspec

import numpy as np

# infer-subc
from infer_subc.utils.stats import create_contact, find_non_redundant_contacts

##################
## Color Constants
##################
ORANGE = '#FFA500'
BOPBLUE = '#20ADF8'
MAROON = '#800000'
WHITE = '#FFFFFF'
BLACK = '#000000'

##############################
## Determines Ideal Font Size
##############################
def find_size(width: float, text: str, base_size:float, fig, multiplier:float=1.00):
    base_size = int(base_size/multiplier)
    r = fig.canvas.get_renderer()
    for fs in list(range(base_size+1))[::-1]:
        t = plt.text(0.5, 0.5,text, fontsize=fs*multiplier)
        tw = t.get_window_extent(renderer=r).width
        if width >= tw:
            t.remove()
            return fs*multiplier, multiplier
        t.remove()
    return find_size(width, text, 10, fig, multiplier*0.1)

#####################
## Plots Overlaps
#####################
def plot_overlaps(orgs: str, splitter: str, organelle_segs: dict, scale: any, wspace: float=0.05, hspace: float=0.05, close: bool=True, holo: bool=True):
    viewer = napari.Viewer() 

    rw = (len(orgs.split(splitter))+1)
    rw = (math.comb(len(orgs.split(splitter)), 2)+1)
    fig = plt.figure(figsize=(3,rw))
    gs= gridspec.GridSpec(rw,3, wspace=wspace, hspace=hspace)

    for i in range(rw*3):
        ax = plt.subplot(gs[i])
        ax.set_aspect('equal')


    ##################################
    ## ORG SEGMENTATIONS
    ##################################
    for org in orgs.split(splitter):
        viewer.add_image(organelle_segs[org]>0, colormap=BOPBLUE, blending ='additive', scale=scale, visible=False, name=f"{org} Blue")
        viewer.add_image(organelle_segs[org]>0, colormap=ORANGE,  blending ='additive', scale=scale, visible=False, name=f"{org} Orange")

    ##################################
    ## 2nd ORDERS INVOLVED
    ##################################
    pos_overlaps = itertools.combinations(orgs.split(splitter), 2)
    overlaps = [splitter.join(inter) for inter in pos_overlaps]
    if holo:
        titles = ['Not In Higher Order', 'All', 'In Higher Order', 'Organelle A', 'Overlap', 'Organelle B']
        legend_ele = [Patch(facecolor=BOPBLUE, edgecolor=BLACK, label="Organelle A"),
                      Patch(facecolor=ORANGE, edgecolor=BLACK, label="Organelle B"),
                      Patch(facecolor=MAROON, edgecolor=BLACK, label="Highlighted Overlap"),
                      Patch(facecolor=WHITE, edgecolor=BLACK, label="Unhighlighted Overlap")]
        for row, inter in enumerate(overlaps):
            LOc = find_non_redundant_contacts(create_contact(inter, organelle_segs, splitter), 
                                            inter, organelle_segs, splitter)
            viewer.add_labels(LOc>0, colormap={1:MAROON}, blending='translucent', 
                            opacity=1.00, scale=scale,visible=False, name=f"{overlaps} Not In Higher Order")
            viewer.add_labels((np.invert(LOc>0)*create_contact(inter, organelle_segs, splitter))>0, 
                            colormap={1:MAROON}, blending='translucent', opacity=1.00, 
                            scale=scale,visible=False, name=f"{overlaps} In Higher Order")
            for org_num, org in enumerate(orgs.split(splitter)):
                if org == inter.split(splitter)[0]:
                    viewer.layers[(org_num*2)].visible = True
                elif org == inter.split(splitter)[1]:
                    viewer.layers[((org_num*2)+1)].visible = True
            highlight_Ac = viewer.export_figure()
            viewer.layers[((len(orgs.split(splitter))*2)+(row*2))].visible = True
            highlight_HOc = viewer.export_figure()
            viewer.layers[((len(orgs.split(splitter))*2)+(row*2))].visible = False
            viewer.layers[((len(orgs.split(splitter))*2)+(row*2)+1)].visible = True
            highlight_LOc = viewer.export_figure()

            # Leftmost Column
            fig.axes[((row*3)+3)].imshow(highlight_LOc, interpolation='nearest')

            # Middle Column
            fig.axes[((row*3)+4)].imshow(highlight_Ac, interpolation='nearest')

            # Rightmost Column
            fig.axes[((row*3)+5)].imshow(highlight_HOc, interpolation='nearest')

            for layer in viewer.layers:
                layer.visible= False
    else:
        titles = ['Organelle A', 'Overlap', 'Organelle B', 'Not In Higher Order', 'All', 'In Higher Order']
        legend_ele = [Patch(facecolor=BOPBLUE, edgecolor=BLACK, label="Organelle A"),
                      Patch(facecolor=ORANGE, edgecolor=BLACK, label="Organelle B"),
                      Patch(facecolor=WHITE, edgecolor=BLACK, label="Overlap")]
        for row, inter in enumerate(overlaps):
            for org_num, org in enumerate(orgs.split(splitter)):
                if org == inter.split(splitter)[0]:
                    viewer.layers[(org_num*2)].visible = True
                elif org == inter.split(splitter)[1]:
                    viewer.layers[((org_num*2)+1)].visible = True
            highlight_Ac = viewer.export_figure()

            for layer in viewer.layers:
                layer.visible= False
            
            for org_num, org in enumerate(orgs.split(splitter)):
                if org == inter.split(splitter)[0]:
                    viewer.layers[(org_num*2)].visible = True
                    highlight_LOc = viewer.export_figure()
                    viewer.layers[(org_num*2)].visible = False
                elif org == inter.split(splitter)[1]:
                    viewer.layers[((org_num*2)+1)].visible = True
                    highlight_HOc = viewer.export_figure()
                    viewer.layers[(org_num*2)].visible = False

            # Leftmost Column
            fig.axes[((row*3)+3)].imshow(highlight_LOc, interpolation='nearest')

            # Middle Column
            fig.axes[((row*3)+4)].imshow(highlight_Ac, interpolation='nearest')

            # Rightmost Column
            fig.axes[((row*3)+5)].imshow(highlight_HOc, interpolation='nearest')

            for layer in viewer.layers:
                layer.visible= False

    if close:
        viewer.close()

    ####################################
    # Determine Font Size
    ####################################
    titles += overlaps
    bbox = fig.axes[3].get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    width = bbox.width * fig.dpi
    fs = plt.rcParams['font.size']
    m = 1.00
    for title in titles:
        if m == 1.00:
            f, m = find_size(width, title, 10, fig, m)
        else:
            f, m = find_size(width, title, fs, fig, m)
        if f < fs:
            fs = f

    ##################################
    ## Row Naming
    ##################################
    for row, inter in enumerate(overlaps):
            fig.axes[((row*3)+3)].set_ylabel(f'{inter}', fontsize=fs)

    ##################################
    ## Column Naming
    ##################################
    fig.axes[3].set_title(titles[0], fontsize=fs)
    fig.axes[0].set_facecolor(('white',0.0))
    fig.axes[4].set_title(titles[1], fontsize=fs)
    fig.axes[1].set_facecolor(('white',0.0))
    fig.axes[5].set_title(titles[2], fontsize=fs)
    fig.axes[2].set_facecolor(('white',0.0))

    #################################
    ## Legend
    #################################
    fig.axes[1].legend(handles=legend_ele, frameon=False, loc='center', fontsize=fs)

    for ax in fig.axes:
        ax.spines['top'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_aspect('equal')
        
    return plt.show()