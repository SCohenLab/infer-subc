import napari
import itertools
import math

# matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import matplotlib.gridspec as gridspec

import numpy as np

# infer-subc
from infer_subc.utils.stats import create_overlap, find_non_redundant_overlaps

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

    #rw = (len(orgs.split(splitter))+1)
    rw = (math.comb(len(orgs.split(splitter)), 2)+1)
    fig = plt.figure(figsize=(4,rw))
    gs= gridspec.GridSpec(rw,4, wspace=wspace, hspace=hspace)

    for i in range(rw*4):
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
        titles = ['Merge', 'Overlap Only', 'Not In Higher Order', 'In Higher Order', 'Organelle A', 'Organelle B']
        legend_ele = [Patch(facecolor=BOPBLUE, edgecolor=BLACK, label="Organelle A"),
                      Patch(facecolor=ORANGE, edgecolor=BLACK, label="Organelle B"),
                      Patch(facecolor=MAROON, edgecolor=BLACK, label="Higher/Lower Order Overlap"),
                      Patch(facecolor=WHITE, edgecolor=BLACK, label="Overlap Region")]
        for row, inter in enumerate(overlaps):
            LOi = find_non_redundant_overlaps(create_overlap(inter, organelle_segs, splitter), 
                                            inter, organelle_segs, splitter)
            viewer.add_labels((create_overlap(inter, organelle_segs, splitter)>0),
                              colormap={1:WHITE}, blending='translucent', opacity=1.00, 
                              scale=scale,visible=False, name=f"{overlaps} Overlap Only")
            viewer.add_labels(LOi>0, colormap={1:MAROON}, blending='translucent', 
                              opacity=1.00, scale=scale,visible=False, name=f"{overlaps} Not In Higher Order")
            viewer.add_labels((np.invert(LOi>0)*create_overlap(inter, organelle_segs, splitter))>0, 
                              colormap={1:MAROON}, blending='translucent', opacity=1.00, 
                              scale=scale,visible=False, name=f"{overlaps} In Higher Order")

            for org_num, org in enumerate(orgs.split(splitter)):
                if org == inter.split(splitter)[0]:
                    viewer.layers[(org_num*2)].visible = True
                elif org == inter.split(splitter)[1]:
                    viewer.layers[((org_num*2)+1)].visible = True
            highlight_Ai = viewer.export_figure()
            for layer in viewer.layers:
                layer.visible= False
            viewer.layers[((len(orgs.split(splitter))*2)+(row*3))].visible = True
            highlight_O = viewer.export_figure()
            viewer.layers[((len(orgs.split(splitter))*2)+(row*3)+1)].visible = True
            highlight_HOi = viewer.export_figure()
            viewer.layers[((len(orgs.split(splitter))*2)+(row*3)+1)].visible = False
            viewer.layers[((len(orgs.split(splitter))*2)+(row*3)+2)].visible = True
            highlight_LOi = viewer.export_figure()

            # Column A (Merge)
            fig.axes[((row*4)+4)].imshow(highlight_Ai, interpolation='nearest')

            # Column B (Overlap Only)
            fig.axes[((row*4)+5)].imshow(highlight_O, interpolation='nearest')

            # Column C (Overlap + Not in Higher Order)
            fig.axes[((row*4)+6)].imshow(highlight_LOi, interpolation='nearest')

            # Column D (Overlap + In Higher Order)
            fig.axes[((row*4)+7)].imshow(highlight_HOi, interpolation='nearest')

            for layer in viewer.layers:
                layer.visible= False
    else:
        titles = ['Organelle A', 'Organelle B', 'Merge', 'Overlap Only', 'Not In Higher Order', 'All', 'In Higher Order']
        legend_ele = [Patch(facecolor=BOPBLUE, edgecolor=BLACK, label="Organelle A"),
                      Patch(facecolor=ORANGE, edgecolor=BLACK, label="Organelle B"),
                      Patch(facecolor=WHITE, edgecolor=BLACK, label="Overlap")]
        for row, inter in enumerate(overlaps):
            for org_num, org in enumerate(orgs.split(splitter)):
                if org == inter.split(splitter)[0]:
                    viewer.layers[(org_num*2)].visible = True
                elif org == inter.split(splitter)[1]:
                    viewer.layers[((org_num*2)+1)].visible = True
            highlight_Ai = viewer.export_figure()

            for layer in viewer.layers:
                layer.visible= False
            
            for org_num, org in enumerate(orgs.split(splitter)):
                if org == inter.split(splitter)[0]:
                    viewer.layers[(org_num*2)].visible = True
                    highlight_LOi = viewer.export_figure()
                    viewer.layers[(org_num*2)].visible = False
                elif org == inter.split(splitter)[1]:
                    viewer.layers[((org_num*2)+1)].visible = True
                    highlight_HOi = viewer.export_figure()
                    viewer.layers[((org_num*2)+1)].visible = False

            viewer.add_labels((create_overlap(inter, organelle_segs, splitter)>0),
                              colormap={1:WHITE}, blending='translucent', opacity=1.00, 
                              scale=scale,visible=True, name=f"{overlaps} Overlap Only")
            highlight_O = viewer.export_figure()

            # Column A (Org A)
            fig.axes[((row*4)+4)].imshow(highlight_LOi, interpolation='nearest')

            # Column B (Org B)
            fig.axes[((row*4)+5)].imshow(highlight_HOi, interpolation='nearest')

            # Column C (Merge)
            fig.axes[((row*4)+6)].imshow(highlight_Ai, interpolation='nearest')

            # Column D (Only Overlapping Region)
            fig.axes[((row*4)+7)].imshow(highlight_O, interpolation='nearest')

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
            fig.axes[((row*4)+4)].set_ylabel(f'{inter}', fontsize=fs)

    ##################################
    ## Column Naming
    ##################################
    fig.axes[4].set_title(titles[0], fontsize=fs)
    fig.axes[0].set_facecolor(('white',0.0))
    fig.axes[5].set_title(titles[1], fontsize=fs)
    fig.axes[1].set_facecolor(('white',0.0))
    fig.axes[6].set_title(titles[2], fontsize=fs)
    fig.axes[2].set_facecolor(('white',0.0))
    fig.axes[7].set_title(titles[3], fontsize=fs)
    fig.axes[3].set_facecolor(('white',0.0))

    #################################
    ## Legend
    #################################
    fig.axes[1].legend(handles=legend_ele, frameon=False, loc='center', fontsize=fs, bbox_to_anchor=(1, (0.5)))

    for ax in fig.axes:
        ax.spines['top'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_aspect('equal')
        
    return plt.show()