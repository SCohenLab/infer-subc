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
import operator

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
def find_size(width: float, text: str, base_size:float, fig, multiplier:float=1.00, weight = None):
    base_size = int(base_size/multiplier)
    r = fig.canvas.get_renderer()
    for fs in list(range(base_size+1))[::-1]:
        t = plt.text(0.5, 0.5,text, fontsize=fs*multiplier, fontweight=weight)
        tw = t.get_window_extent(renderer=r).width
        if width >= tw:
            t.remove()
            print(fs*multiplier)
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

def plot_n3_overlaps(orgs: str, splitter: str, organelle_segs: dict, scale: any, padding: float=0.05):
    # COLOR CONSTANTS
    ORG_A_COL = '#E69F00'
    ORG_B_COL = '#0072B2'
    ORG_C_COL = '#CC79A7'
    BC_MERGE  = '#56B4E9'
    AC_MERGE  = '#D55E00'
    AB_MERGE  = '#009E73'
    ABC_MERGE = '#FFFFFF'

    # Set all 2w colors to be the same color for visibility

    COLORS = [ORG_A_COL, ORG_B_COL, ORG_C_COL, AB_MERGE, AC_MERGE, BC_MERGE, ABC_MERGE]

    viewer = napari.Viewer() # Initializes the viewer
    plt_org = {}             # Dictionary of the images to plot
    axes = {}                # Dictionary of the axes to plot on

    # Generate all possible combinations of the organelles
    all_pos =[]
    for n in list(map(lambda x:x+2, (range(len(orgs.split(splitter))-1)))):
        all_pos += itertools.combinations(orgs.split(splitter), n)
    possib = [splitter.join(inter) for inter in all_pos]

    # Adds organelles and their combinations to the viewer, and exports the images to the plt_org dictionary
    for i, org in enumerate(orgs.split(splitter) + possib):
        print(org)
        if org in organelle_segs.keys():
            viewer.add_labels(organelle_segs[org]>0, name=org, colormap={1:COLORS[i]}, visible=True)
            plt_org[org] = viewer.export_figure()
        else:
            viewer.add_labels((create_overlap(org, organelle_segs, splitter)>0), name=org, colormap={1:COLORS[i]}, visible=True)
            for sub_org in (org.split(splitter) + [inter for inter in possib if all(sub in org for sub in inter.split(splitter))]):
                viewer.layers[sub_org].visible = True
            plt_org[org] = viewer.export_figure()
        for layer in viewer.layers:
            layer.visible = False
        if org == orgs:
            viewer.layers[org].visible = True
            plt_org[f"{org} Only"] = viewer.export_figure()
    
    # Set up the figure
    col = len(orgs.split(splitter))+1                               # Number of columns
    rw = (math.comb(len(orgs.split(splitter)), 2)+1)                # Number of rows

    fig = plt.figure(figsize=(col,rw))                              # Initialize the figure
    gs = fig.add_gridspec(rw, col, wspace=padding, hspace=padding)  # Initialize the grid spec

    large_plots = [(i+1)+(4*(x+1)) for i in range(3) for x in range(rw-1)] # Indices the large image is displayed at
    regular_plots = [i for i in range(rw*col) if i not in large_plots]     # Indices regular sized images are displayed


    for i, plot in zip(regular_plots, [f"{orgs} Only", f"{orgs.split(splitter)[0]}", f"{orgs.split(splitter)[1]}", 
                                       f"{orgs.split(splitter)[2]}"] + possib[:-1]):
        axes[plot] = fig.add_subplot(gs[i])
        axes[plot].imshow(plt_org[plot], interpolation='nearest')

    axes[orgs] = fig.add_subplot(gs[1:, 1:])
    axes[orgs].imshow(plt_org[orgs], interpolation='nearest')
 
    for ax in axes.values():
        ax.spines['top'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_aspect('equal')

    ######################################################################
    #                    Col A        Col B        Col C        Col D    #
    #                 "A+B+C Only"   "Org A"      "Org B"      "Org C"   #
    # r1                (A+B+C)      (org a)      (org b)      (org c)   #
    #                                                                    #
    # r2   "A+B merge"  (a + b)   M                                      #
    #                             E                                      #
    # r3   "A+C merge"  (a + c)   R               (a+b+c)                #
    #                             G                                      #
    # r4   "B+C merge"  (b + c)   E                                      #
    ######################################################################

    #(a+b+c) plot is ~3x height and width of others
    # math: (a+b+c) width = (width a) + (width b) + (width c) + <padding> + <padding>
    # math: (a+b+c) height = (height a) + (height b) + (height c) + <padding> + <padding>
    # tldr: (a+b+c) plot is 3x + 2(padding) for both height and width
    # NOTE: unsure if it is actually 4x padding, must test to see

    # Color A: #E69F00 |   Color = Orange   | Org A
    # Color B: #56B4E9 |  Color = Sky Blue  | B + C
    # Color C: #009E73 | Color = Teal Green | A + B
    # Color D: #0072B2 |    Color = Blue    | Org B
    # Color E: #D55E00 | Color = Vermillion | A + C
    # Color F: #CC79A7 |    Color = Pink    | Org C
    # Color G: #000000 |    Color = Black   | Background
    # Color H: #FFFFFF |    Color = White   | A + B + C
    #viewer.close()
    return plt.show()

def crop_to_square(img):
    yy,xx,cc = img.shape
    if yy > xx:
        bounding = (xx, xx, cc)
    elif xx > yy:
        bounding = (yy, yy, cc)
    else:
        return img
    start = tuple(map(lambda a, da: a//2-da//2, img.shape, bounding))
    end = tuple(map(operator.add, start, bounding))
    slices = tuple(map(slice, start, end))
    return img[slices]


def plot_n_overlaps(orgs: str, splitter: str, organelle_segs: dict, scale: any, padding: float=0.05, close_viewer: bool=False, view='2D'):
    """
    Plots nth order overlaps of organelles in a grid format. Only works for up to 7 organelles.

    Parameters
    ----------

    Returns
    -------
    plt
        Figure displaying the nth dimensional image alongside the lower order interactions

    """
    # ORG_A_COL = '#E69F00'
    # ORG_B_COL = '#0072B2'
    # ORG_C_COL = '#CC79A7'
    # ORG_D_COL = '#009E73'
    # ORG_E_COL = '#F0E442'
    # ORG_F_COL = '#D55E00'
    # ORG_G_COL = '#56B4E9'
    # HO_MERGE = '#FFFFFF'
    # LO_MERGE = '#999999'

    ORG_A_COL = "#FFAE00"
    ORG_B_COL = "#00A0FC"
    ORG_C_COL = "#FF93CE"
    ORG_D_COL = "#01CE97"
    ORG_E_COL = '#F0E442'
    ORG_F_COL = "#FF6F01"
    ORG_G_COL = "#60C5FF"
    HO_MERGE = '#FFFFFF'
    LO_MERGE = '#999999'
    ORG_COLORS = [ORG_A_COL, ORG_B_COL, ORG_C_COL, ORG_D_COL, ORG_E_COL, ORG_F_COL, ORG_G_COL]
    viewer = napari.Viewer()# Initializes the viewer
    viewer.window.resize(width=800, height=800)
    plt_org = {}            # Dictionary of the images to plot
    axes = {}               # Dictionary of the axes to plot on
    counts = []             # List of number of plots per section
    y_titles = {}           # Dictionary of titles used vertically
    x_titles = {}           # Dictionary of titles used horizontally

    if view == '3D':
        dim = 3
        op = 0.5
    elif view == '2D':
        dim = 2
        op = 1.0
    elif type(view) == int:
        dim = 2
        op = 1.0
    else:
        raise ValueError("view must be either '2D' or '3D'")

    # Generate all possible combinations of the organelles
    all_pos =[]
    for n in list(map(lambda x:x+2, (range(len(orgs.split(splitter))-1)))):
        all_pos += itertools.combinations(orgs.split(splitter), n)
        counts.append(math.comb(len(orgs.split(splitter)), n))
    possib = [splitter.join(inter) for inter in all_pos]

    # Determine number of each order of LO overlaps
    counts = counts[:-1] + [len(orgs.split(splitter))+1]

    grid_width = math.lcm(*counts)
    grid_height = sum(grid_width/h for h in counts[:-1]) + ((grid_width/(len(orgs.split(splitter))+1))*(len(orgs.split(splitter)))) # height of large plot
    
    # Determine the minimum gridspec area
    min_spec = grid_width//max(counts)

    # Adds organelles and their combinations to the viewer, and exports the images to the plt_org dictionary
    for i, org in enumerate(orgs.split(splitter) + possib):
        if org in organelle_segs.keys():
            viewer.add_labels(organelle_segs[org]>0, name=f"{org}_LO", scale=scale, opacity=op, colormap={0:None,1:ORG_COLORS[i]}, visible=True)
            viewer.dims.ndisplay = dim
            if type(view) == int:
                viewer.dims.set_point(0, view)
            plt_org[org] = crop_to_square(viewer.screenshot(canvas_only=True))
        else:
            viewer.add_labels((create_overlap(org, organelle_segs, splitter)>0), name=f"{org}_LO", scale=scale, opacity=1.0, colormap={0:None,1:LO_MERGE}, visible=False)
            viewer.add_labels((create_overlap(org, organelle_segs, splitter)>0), name=f"{org}_HO", scale=scale, opacity=1.0, colormap={0:None,1:HO_MERGE}, visible=True)
            viewer.dims.ndisplay = dim
            if type(view) == int:
                viewer.dims.set_point(0, view)
            plt_org[f"{org}_ol"] = crop_to_square(viewer.screenshot(canvas_only=True))
            for sub_org in (org.split(splitter) + [inter for inter in possib if (all(sub in org for sub in inter.split(splitter)) and inter != org)]):
                viewer.layers[f"{sub_org}_LO"].visible = True
            plt_org[org] = crop_to_square(viewer.screenshot(canvas_only=True))
        for layer in viewer.layers:
            layer.visible = False
    
    if close_viewer:
        viewer.close()

    fig = plt.figure(figsize=((2*grid_width),grid_height)) 

    # Initialize the grid spec
    if grid_width > grid_height:
        # fig = plt.figure(figsize=((2*grid_width/grid_height)*12,12)) 
        gs = fig.add_gridspec(int(grid_height), int(2*grid_width), wspace=(padding*2*(grid_width/grid_height)), hspace=padding)
    elif grid_height > grid_width:
        # fig = plt.figure(figsize=((2)*12,(grid_height/grid_width)*12)) 
        gs = fig.add_gridspec(int(grid_height), int(2*grid_width), wspace=padding, hspace=(padding*2*(grid_height/grid_width)))  
    else:
        # fig = plt.figure(figsize=((2)*12,12)) 
        gs = fig.add_gridspec(int(grid_height), int(2*grid_width), wspace=padding, hspace=padding)

    # Merge Plots
    for n in list(map(lambda x:x+2, (range(len(orgs.split(splitter))-1)))): # n = overlap order number
        if n == (len(orgs.split(splitter))):
            # Single organelle 
            axes[orgs] = fig.add_subplot(gs[int(grid_height - (len(orgs.split(splitter))*(grid_width/(len(orgs.split(splitter))+1)))):,
                                            int(grid_width/(len(orgs.split(splitter))+1)):int(grid_width)])
            axes[f"{orgs}_ol"] = fig.add_subplot(gs[int(grid_height - (len(orgs.split(splitter))*(grid_width/(len(orgs.split(splitter))+1)))):,
                                                    int(grid_width):int((2*grid_width)-(grid_width/(len(orgs.split(splitter))+1)))])
            y = int(sum((grid_width/counts[j-1]) for j in range(1, len(orgs.split(splitter)))))
            x_titles[orgs] = orgs
            x_titles["orgs"] = "Organelles"
            for i, org in enumerate(orgs.split(splitter)):
                axes[org] = fig.add_subplot(gs[int((y-(grid_width/(len(orgs.split(splitter))+1))+(i*(grid_width/(len(orgs.split(splitter))+1))))):int(y+((i)*(grid_width/(len(orgs.split(splitter))+1)))),
                                               int(0):int(grid_width/(len(orgs.split(splitter))+1))])
                y_titles[org] = org

        else:
            y = int(sum((grid_width/counts[j-1]) for j in range(1, n)))
            for i, org in enumerate([splitter.join(inter) for inter in itertools.combinations(orgs.split(splitter), n)]):   # i = overlap image number within the order
                axes[org] = fig.add_subplot(gs[int(y - (grid_width/counts[n-2])):y,
                                            int(((2*i)*(grid_width/counts[n-2]))):int(((2*i*(grid_width/counts[n-2]))+(grid_width/counts[n-2])))])
                axes[f"{org}_ol"] = fig.add_subplot(gs[int(y - (grid_width/counts[n-2])):y,
                                                       int(((2*i)*(grid_width/counts[n-2]))+(grid_width/counts[n-2])):int(((2*i)*((grid_width/counts[n-2]))+(2*(grid_width/counts[n-2]))))])
                x_titles[org] = org
                if i == 0:
                    y_titles[org] = f"Interaction Order: {n}"
    
    # combine dictionaries of titles
    all_titles = {"Key":"Legend"}
    all_titles.update(x_titles)
    all_titles.update(y_titles)

    #######################
    # Determine Font Size #
    #######################
    # bbox = axes[list(axes.keys())[0]].get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    # width = bbox.width * fig.dpi
    width = min([axes[ax].get_window_extent().transformed(fig.dpi_scale_trans.inverted()).width * fig.dpi for ax in axes.keys()])
    fs = grid_width
    m = 1.00
    for title in all_titles.keys():
        if m == 1.00:
            f, m = find_size(width, all_titles[title], grid_width*2, fig, m)
        else:
            f, m = find_size(width, all_titles[title], fs, fig, m)
        if f < fs:
            fs = f
    
    ######################
    # Assigning Y Titles #
    ###################### 
    for title in y_titles.keys():
        axes[title].set_ylabel(y_titles[title], size=fs)

    ########################
    # Cleaning Image Plots #
    ########################
    for org in axes.keys():
        axes[org].spines['top'].set_visible(False)
        axes[org].spines['bottom'].set_visible(False)
        axes[org].spines['right'].set_visible(False)
        axes[org].spines['left'].set_visible(False)
        axes[org].set_xticks([])
        axes[org].set_yticks([])
        axes[org].set_aspect('equal')
        axes[org].imshow(plt_org[org], interpolation='nearest')
    
    #####################
    # Setting Up Legend #
    #####################
    axes["key"] = fig.add_subplot(gs[int(grid_height - (len(orgs.split(splitter))*(grid_width/(len(orgs.split(splitter))+1)))):,
                                     int((2*grid_width)-(grid_width/(len(orgs.split(splitter))+1))):])
    axes["key"].set_xticks([])
    axes["key"].set_yticks([])
    axes["key"].spines['top'].set_visible(False)
    axes["key"].spines['bottom'].set_visible(False)
    axes["key"].spines['right'].set_visible(False)
    axes["key"].spines['left'].set_visible(False)
    legend_ele = [Patch(facecolor=ORG_COLORS[i], edgecolor="#000000", label=org) for i, org in enumerate(orgs.split(splitter))]
    legend_ele += [Patch(facecolor=HO_MERGE, edgecolor="#000000", label="Nth Order Overlap"), 
                   Patch(facecolor=LO_MERGE, edgecolor="#000000", label="Lower Order Overlap")]
    axes["key"].legend(handles=legend_ele, fontsize=((3*fs)//4), loc='right', 
                       frameon=False, mode='expand', title="Legend", title_fontsize=fs)
    print("Done")
    return fig