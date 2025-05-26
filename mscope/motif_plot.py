import sys
import matplotlib
from matplotlib import pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection, LineCollection, QuadMesh

import scipy.cluster.hierarchy as shc
from scipy.spatial.distance import squareform

import numpy as np
import copy

def heatmap(cfg, seq_order, grouped_motif_seq, sequence_lengths, dim_reduction, motif_counts, cmap, ax, cbar_ax, cbar_kws= None,cbar_sb_ax=None):
    assert len(seq_order) == len(grouped_motif_seq), "Number of sequences in seq_order and grouped_motif_seq should be the same"
    assert len(seq_order) == len(sequence_lengths), "Number of sequences in seq_order and sequence_lengths should be the same"
   
    edgecolor = cfg.edgecolor
    linewidth = cfg.linewidth
    #linewidth = 0.0001

    #to add singlebase edges, set singlebase_edges to True in the Config section in the 'motifscope' file
    if cfg.singlebase_edges:
        singlebase_linewidth = linewidth
        singlebase_edgecolor = edgecolor
    else:
        singlebase_linewidth = 0
        singlebase_edgecolor = 'face'
    


    #handling single bp colors
    single_motifs = ['A', 'C', 'G', 'T']
    if 'N' in motif_counts:
        single_motifs.append('N')                 # keeps the same logic you had

    used_single_motifs = [x for x in motif_counts if x in single_motifs]
    # ── 2. build ONE integer mapping and ONE colormap ─────────────────────────────
    motif_to_idx = {m:i for i, m in enumerate(single_motifs)}   # A→0, C→1, …

    # greys for A-T, yellow for N – same colours you chose
    base_colours = [
        (0.35, 0.35, 0.35),   # A
        (0.55, 0.55, 0.55),   # C
        (0.75, 0.75, 0.75),   # G
        (0.90, 0.90, 0.90)    # T
    ]
    if 'N' in single_motifs:
        base_colours.append('yellow')
    

    singlebase_used_cmap = matplotlib.colors.ListedColormap(base_colours)
    #convert_color = {nuc: singlebase_used_cmap(idx) for nuc, idx in motif_to_idx.items()}

    #show single value after the decimal point
    sblabels = [
    f"{base} ({motif_counts[base] / float(len(grouped_motif_seq)):.1f})"
    for base in single_motifs          # keeps A-C-G-T order
    if base in motif_counts            # skip bases that don’t occur
    ]
    #sblabels = [f"{key} ({motif_counts[key] / float(len(grouped_motif_seq)):.1f})" for key in single_motifs[::-1]]
    if len(sblabels) > 0:
        single_bp_color(cfg, sblabels, cbar_sb_ax, singlebase_used_cmap)

    convert_color = {nuc: singlebase_used_cmap(idx) for nuc, idx in motif_to_idx.items()}
    convert_color.update({nuc.lower(): rgba for nuc, rgba in convert_color.items()})

    nmotifs = dim_reduction['motif'].nunique()
    umotifs = dim_reduction['motif'].tolist()

    # Create the motif colormap
    cmap1 = plt.get_cmap('YlGnBu_r')
    motif_colors = cmap(np.linspace(0, 1, nmotifs))
    motif_cmap = matplotlib.colors.ListedColormap(motif_colors)

    # Normalize for the motif part of the combined colormap
    norm = matplotlib.colors.BoundaryNorm(np.arange(len(umotifs) + 1) - 0.5, len(umotifs))

    colorcache = {}
    def cache_color(value):
        if value not in colorcache:
            colorcache[value] = motif_cmap(norm(value))
        return colorcache[value]
 
    singlebase_rectangles = [] 
    rectangles = []
    motif_sep_lines = []
    for seq in seq_order:
        motifs = grouped_motif_seq[seq]
        ypos = seq_order.index(seq)
        for mp in motifs:
            start, end, motif, count, color_value = mp.score_items()
            
            if mp.singlebase: #single bases
                s = len(motif)
                if s == 1: #fast path
                    color = convert_color[motif]
                    singlebase_rectangles.append(Rectangle((start, ypos), 1, 1, fill=True, lw = singlebase_linewidth, edgecolor=edgecolor, facecolor=color, clip_on=False))
                else: #use quadmesh for faster rendering
                    xcolors = [convert_color[b] for b in motif]
                    xpos = np.linspace(start, end, len(motif) + 1)
                    coords = np.zeros((2, len(xpos), 2))
                    coords[0,:,0] = xpos
                    coords[1,:,0] = xpos
                    coords[0,:,1] = ypos
                    coords[1,:,1] = ypos + 1

                    #quadmesh_coords = np.array([[[x, ypos], [x, ypos + 1], [x + 1, ypos + 1], [x + 1, ypos]] for x in xpos])

                    q = QuadMesh(coords, facecolors=xcolors, rasterized=False, lw=singlebase_linewidth, edgecolor=singlebase_edgecolor)
                    ax.add_collection(q)
            else: 
                color = cache_color(color_value)
                if not mp.leftborder or not mp.rightborder:
                    rectangles.append(Rectangle((start, ypos), end - start, 1, fill=True,  lw=0, facecolor=color, clip_on=False))
                    if mp.leftborder:
                        motif_sep_lines.append(((start, ypos), (start, ypos + 1)))
                    if mp.rightborder:
                        motif_sep_lines.append(((end, ypos), (end, ypos + 1)))
                else:
                    rectangles.append(Rectangle((start, ypos), end - start, 1, fill=True, edgecolor=edgecolor, lw=linewidth, facecolor=color, clip_on=False))
                if count > 0:
                    xstart = start + mp.first_motif_offset
                    xcount = int(np.ceil(count - mp.first_motif_offset / len(mp.motif)))
                    xend = xstart + len(mp.motif) * xcount
                    sep_xpos = np.linspace(xstart, xend, int(xcount) + 1)[1:-1]
                    motif_sep_lines.extend([((x, ypos), (x, ypos + 1)) for x in sep_xpos])

    p = PatchCollection(singlebase_rectangles, match_original = True)
    ax.add_collection(p)

    p = PatchCollection(rectangles, match_original = True)            
    ax.add_collection(p)
    l = LineCollection(motif_sep_lines, color=edgecolor, lw=linewidth)
    ax.add_collection(l)
    ax.set_ylim(0, len(grouped_motif_seq))
    
    ypos = np.arange(len(seq_order))
    rectangles = [Rectangle((0, y), sequence_lengths[seq], 1, fill=False, edgecolor=edgecolor, lw=linewidth, facecolor=edgecolor, clip_on=False) for seq,y in zip(seq_order,ypos)]
    p = PatchCollection(rectangles, match_original = True)
    ax.add_collection(p)

    cb = plt.colorbar(matplotlib.cm.ScalarMappable(cmap=motif_cmap, norm=norm), cax=cbar_ax, ticks=np.arange(len(umotifs)), spacing='uniform', orientation='vertical')
    #ulabels = [f"{motif} ({motif_counts[motif] / float(len(grouped_motif_seq)):.1f})" for motif in umotifs]
    ulabels = [f"{motif} ({motif_counts.get(motif, 0) / float(len(grouped_motif_seq)):.1f})" for motif in umotifs]

    cbar_ax.set_yticklabels(ulabels)
    cbar_ax.tick_params(labelsize=cfg.cbar_fontsize)
    cbar_ax.title.set_fontsize(cfg.cbar_fontsize)

    cb.outline.set_linewidth(0.1)
    cb.outline.set_edgecolor('black')
    return cb



def pop_heatmap(cfg, data, poplabels, ax, cbar, cmap, clabel):
    data = np.array(data)
    
    bounds = np.arange(0, len(poplabels) + 1) - 0.5
    norm = matplotlib.colors.BoundaryNorm(bounds, cmap.N)

    ypos = np.arange(data.shape[0] + 1)
    coords = np.zeros((len(data) + 1, 2, 2))
    coords[:,0,1] = ypos
    coords[:,1,1] = ypos
    coords[:,0,0] = 0
    coords[:,1,0] = 1
    #q = QuadMesh(coords, facecolors=cmap(norm(data)), rasterized=False, edgecolor='silver', lw=0.1)
    q = QuadMesh(coords, facecolors=cmap(norm(data)), rasterized=False, edgecolor='black', lw=0.1)
    ax.add_collection(q)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, len(data))
    ax.set_ylabel("")
    ax.set_yticks([])
    ax.set_xticks([0.5])
    ax.set_xticklabels([clabel],rotation=90,fontsize=cfg.heatmap_labels_xfontsize)


    cb = plt.colorbar(matplotlib.cm.ScalarMappable(cmap=cmap, norm=norm), cax=cbar, ticks=np.arange(len(poplabels)), spacing='uniform', orientation='vertical')
    cbar.set_yticklabels(poplabels)
    cb.outline.set_linewidth(0.1)
    cb.outline.set_edgecolor('black')
    cbar.tick_params(labelsize=cfg.cbar_fontsize)
    cbar.title.set_fontsize(cfg.cbar_fontsize)

def single_bp_color(cfg, sblabels, cbar, cmap):
    bounds = np.arange(0, len(sblabels) + 1) - 0.5
    norm = matplotlib.colors.BoundaryNorm(bounds, cmap.N)

    cb = plt.colorbar(matplotlib.cm.ScalarMappable(cmap=cmap, norm=norm), cax=cbar, ticks=np.arange(len(sblabels)), spacing='uniform', orientation='vertical')
    cbar.set_yticklabels(sblabels)
    cb.outline.set_linewidth(0.1)
    cb.outline.set_edgecolor('black')
    cbar.tick_params(labelsize=cfg.cbar_fontsize)
    cbar.title.set_fontsize(cfg.cbar_fontsize)

def count_motifs(grouped_motif_seq):
    motif_counts = {}
    for seq, motifs in grouped_motif_seq.items():
        for mp in motifs:
            if mp.singlebase:
                for letter in mp.motif:
                    if letter not in motif_counts:
                        motif_counts[letter] = 0
                    motif_counts[letter] += 1
            else:                
                motif = mp.motif
                if motif not in motif_counts:
                    motif_counts[motif] = 0
                motif_counts[motif] += mp.count
    return motif_counts

class MotifPlot:
    def __init__(self, cfg, title, file):
        self.plot_dendrogram = False
        self.plot_classes = False
        self.cfg = cfg
        self.figtitle = title
        self.figfile = file
        self.max_classlbl_length = self.max_seqlbl_length = self.max_motif_length = 0

    def load_data(self, grouped_motif_seq, sequence_lengths, dim_reduction):
        self.grouped_motif_seq = grouped_motif_seq
        self.dim_reduction = dim_reduction
        
        self.sequence_lengths = sequence_lengths
        self.max_seq_length = max(self.sequence_lengths.values())
        self.nsamples = len(grouped_motif_seq)
        self.motif_counts = count_motifs(self.grouped_motif_seq)
        
        self.max_seqlbl_length = max([len(e) for e in self._sequence_id_to_label(list(self.grouped_motif_seq.keys()))])
        self.max_motif_length = max([len(e) for e in self.dim_reduction['motif']])


    def enable_dendrogram(self, seq_distance_df):
        assert len(self.grouped_motif_seq) == len(seq_distance_df), "Number of sequences in grouped_motif_seq and seq_distance_df should be the same"
        self.plot_dendrogram = True
        self.seq_distance_df = seq_distance_df

    def enable_classes(self, classes_df):
        self.plot_classes = True
        self.classes_df = classes_df
        #get name of second column, used as label
        self.classes_label = classes_df.columns[1]
        self.max_classlbl_length = max([len(e) for e in classes_df[self.classes_label]])
       
        self.nclasses = len(classes_df[self.classes_label].unique())

    def create_figure(self):
        fig, axes = self._prepare_fig()

        if self.plot_dendrogram:
            seq_order_list = self._plot_dendrogram(fig, axes)
        else:
            #seq_order_list = np.arange(len(self.grouped_motif_seq))
            seq_order_list = list(self.grouped_motif_seq.keys())

        self._plot_heatmap(fig, axes, seq_order_list)

        if self.plot_classes:
            self._plot_classes(fig, axes, seq_order_list)

        axes['heatmap'].title.set_text(self.figtitle)
        axes['heatmap'].title.set_fontsize(self.cfg.title_fontsize)
        sys.stderr.write("Saving figure...\n"); sys.stderr.flush()
        plt.savefig(self.figfile, bbox_inches = "tight")

    def _prepare_fig(self):       

        COLORBAR_TITLE_VSPACE = 0.6
        CBAR_UNIT = 1.0
        nsbmotifs = len([x for x in self.motif_counts if x in ["A", "T", "C", "G"]])
        nmotifs = len(self.motif_counts) - nsbmotifs

        #nmotifs = len(self.motif_counts) - 4 #correct for single base motifs
        heatmap_height = 0.5 * self.nsamples

        #adapt CBAR_UNIT
        CBAR_UNIT =  min(1.0, max(0.5, 0.5 * (float(self.nsamples) / nmotifs)))

        #calculate vertical space needed for colorbars if we put them in one column
        vspace_needed = CBAR_UNIT * len(self.motif_counts) + COLORBAR_TITLE_VSPACE
        if self.plot_classes:
            vspace_needed += CBAR_UNIT * self.nclasses + COLORBAR_TITLE_VSPACE


        #determine if we need two columns for colorbars
        col_destination = {'motifs':1, 'single_bp':1, 'classes':1} #default: all colorbars in col1
        if heatmap_height < vspace_needed:
            ncols = 2
            #calculate minimum vspace needed for each column
            if self.plot_classes:
                vspace_needed1 = max(nmotifs * CBAR_UNIT, CBAR_UNIT * (nsbmotifs + self.nclasses) + COLORBAR_TITLE_VSPACE) #motifs in col1 , single bp and classes in col2
                vspace_needed2 = max(len(self.motif_counts) * CBAR_UNIT + COLORBAR_TITLE_VSPACE, CBAR_UNIT * self.nclasses) #motifs and single bp in col1, classes in col2
                if vspace_needed1 < vspace_needed2:
                    col_destination['single_bp'] = 2
                col_destination['classes'] = 2
                vspace_needed = min(vspace_needed1, vspace_needed2)
            else:
                vspace_needed = CBAR_UNIT * max(nmotifs, nsbmotifs)
                col_destination['single_bp'] = 2
        else:
            ncols = 1
            
        #determine space needed for labels of colorbars
        lblspace = {'motifs': (self.max_motif_length + 6) * (self.cfg.cbar_fontsize / 72.0), 
                    'classes': (self.max_classlbl_length + 6) * (self.cfg.cbar_fontsize / 72.0),
                    'single_bp': (1 + 6) * (self.cfg.cbar_fontsize / 72.0)
                    }
        col_lblspace = [0] * ncols
        for key, value in col_destination.items():
            col_lblspace[value - 1] = max(col_lblspace[value - 1], lblspace[key])


        #determine figure height (max of heatmap height and colorbar height)
        height = max(heatmap_height, vspace_needed) * 1.2
       
       
        #determine figure width
        COLORBAR_SPACE = [(1 + c) for c in col_lblspace]
        CLASSES_SPACE=0.35
        DENDROGRAM_SPACE=4
        SEQLBL_SPACE = (self.max_seqlbl_length + 3) * (self.cfg.heatmap_labels_yfontsize / 72.0)
        width = sum(COLORBAR_SPACE) + max(0.015 * self.max_seq_length, 20)
        width += SEQLBL_SPACE

        width += CLASSES_SPACE if self.plot_classes else 0
        width += DENDROGRAM_SPACE if self.plot_dendrogram else 0
        width = min(width, 120)

       
        #split up figure in 4 parts: dendrogram, classes, heatmap, colorbar
        width_ratios = []
        heatmap_width = width - sum(COLORBAR_SPACE)
        if self.plot_dendrogram:
            width_ratios.append(DENDROGRAM_SPACE/width)
            heatmap_width -= DENDROGRAM_SPACE
        if self.plot_classes:
            width_ratios.append(CLASSES_SPACE/width)
            heatmap_width -= CLASSES_SPACE
        width_ratios.append(heatmap_width / width)
        width_ratios.append(SEQLBL_SPACE / width)
        for c in COLORBAR_SPACE:
            width_ratios.append(c / width)
        

        fig = plt.figure(figsize=(width, height), dpi = self.cfg.dpi)

        if heatmap_height < height:
            #reduce height of heatmap compared to colorbars if needed (many motifs, few sequences), create two rows in that case.
            spec = fig.add_gridspec(ncols=len(width_ratios), nrows=2, width_ratios=width_ratios, height_ratios=[heatmap_height / height, 1 - heatmap_height / height], wspace=0.02)
        else:
            spec = fig.add_gridspec(ncols=len(width_ratios), nrows=1, width_ratios=width_ratios, wspace=0.02)
        
        axes = {}
        cur_pos = 0
        if self.plot_dendrogram:
            axes['dendrogram'] = fig.add_subplot(spec[0, cur_pos])
            cur_pos += 1
        if self.plot_classes:
            axes['classes'] = fig.add_subplot(spec[0, cur_pos])
            cur_pos += 1
        axes['heatmap'] = fig.add_subplot(spec[0, cur_pos])


        colorbar_specs = [spec[:,-(i+1)] for i in range(ncols)][::-1]

        colorbar_xpos = [s.get_position(fig).x0 for s in colorbar_specs]

            
        colorbar_pos = spec[:,-ncols:].get_position(fig)

        #print(heatmap_pos)
        
        max_height = (colorbar_pos.y1 - colorbar_pos.y0) * height
        max_width = (colorbar_specs[0].get_position(fig).x1 - colorbar_specs[-1].get_position(fig).x0) * width
        colorbar_axes_width = 1.0 / width

        sep_dist_colorbar_y = COLORBAR_TITLE_VSPACE / max_height
        column1_max = column2_max = colorbar_pos.y1
      
        #add color bar for motifs
        cbar_height = min(CBAR_UNIT * nmotifs, max_height)
        #print(cbar_height, height, nmotifs, max_height)

        

        axes['cbar_heatmap'] = fig.add_axes([colorbar_xpos[0], column1_max  - cbar_height / height, colorbar_axes_width, cbar_height / height], title="motifs")
        column1_max = axes['cbar_heatmap'].get_position().y0 - sep_dist_colorbar_y

        #add color bar for singe bp. 
        #cbar_sb_height = min(1.0 * 4, max_height)
        cbar_sb_height = min(len([x for x in self.motif_counts if x in ["A", "T", "C", "G"]]) * CBAR_UNIT, max_height)

        # see if we can fit the single bp motif colorbar below the motif colorbar
        if col_destination['single_bp'] == 1:
            axes['cbar_single_bp'] = fig.add_axes([colorbar_xpos[0], column1_max - cbar_sb_height / height, colorbar_axes_width, cbar_sb_height / height], title="single bp\nmotifs")
            column1_max = axes['cbar_single_bp'].get_position().y0 - sep_dist_colorbar_y
        else:
            #put it next to it.
            axes['cbar_single_bp'] = fig.add_axes([colorbar_xpos[1], column2_max - cbar_sb_height / height, colorbar_axes_width, cbar_sb_height / height], title="single bp\nmotifs")
            column2_max = axes['cbar_single_bp'].get_position().y0 - sep_dist_colorbar_y
        

        if self.plot_classes:
            cbar_cls_height = min(CBAR_UNIT * self.nclasses, max_height)
            if col_destination['classes'] == 1:
                axes['cbar_classes'] = fig.add_axes([colorbar_xpos[0], column1_max - cbar_cls_height / height, colorbar_axes_width, cbar_cls_height / height], title=self.classes_label)
            else:
                axes['cbar_classes'] = fig.add_axes([colorbar_xpos[1], column2_max - cbar_cls_height / height, colorbar_axes_width, cbar_cls_height / height], title=self.classes_label)

        return fig, axes


    def _plot_dendrogram(self, fig, axes):
        clusters = shc.linkage(squareform(self.seq_distance_df), method='average', metric="euclidean", optimal_ordering=True)
        d = shc.dendrogram(Z = clusters, ax = axes['dendrogram'], labels = self.seq_distance_df.index, orientation = "left")
        axes['dendrogram'].tick_params(right=False, left = False, top=False, bottom=False, labelright=False, labelleft=False,labeltop=False)
        axes['dendrogram'].set_xticks([])

        seq_order_list = list(d["ivl"])
        return seq_order_list


    def _sequence_id_to_sample(self, seq_order_list):
        return [seq.split("#")[0] for seq in seq_order_list]

    def _sequence_id_to_label(self, seq_order_list):
        res = []
        for seq in seq_order_list:
            parts = seq.split("#")
            if len(parts) >= 2:
                res.append(f"{parts[0]}:{parts[1]}")
            elif len(parts) == 1:
                res.append(parts[0])
        return res                

    def _plot_heatmap(self, fig, axes, seq_order_list):
        ax = axes['heatmap']

        sequences = self._sequence_id_to_label(seq_order_list)
       
        if self.cfg.args.embed_motif_method == "random":
            g = len(self.dim_reduction.motif.unique())
            colors = np.random.rand(g, 3)
            # Create a colormap from the random colors
            cmap = matplotlib.colors.ListedColormap(colors)
        else:
            cmap = copy.copy(plt.get_cmap('YlGnBu_r'))
            #cmap = copy.copy(plt.get_cmap('Spectral'))
            #cmap = copy.copy(plt.get_cmap('YlGnBu'))
            cmap.set_over('none')
            #scale colorbar to the range of the dimension reduction


        heatmap(self.cfg, seq_order_list, self.grouped_motif_seq, self.sequence_lengths, self.dim_reduction, 
                        cmap=cmap, ax=ax, cbar_ax = axes['cbar_heatmap'], 
                        cbar_sb_ax = axes['cbar_single_bp'],
                        motif_counts = self.motif_counts)

        ax.set(ylabel="")

        pos_step_sizes = [1, 2, 5, 10, 25, 50, 100, 250, 500, 1000, 2500, 5000, 10000]
        pos_cur = 0
        while (self.max_seq_length / pos_step_sizes[pos_cur]) > 10.0 and pos_cur < len(pos_step_sizes) - 1:
            pos_cur += 1
        step_size = pos_step_sizes[pos_cur]
        xticks_positions = np.arange(0, ((self.max_seq_length//step_size) + 1) * step_size + 1, step_size)
        xticks_labels = [str(x) for x in xticks_positions]
        ax.tick_params(right=True, left = False, top=False, labelright=True, labelleft=False, labeltop=False, rotation=0)
        ax.set_xticks(xticks_positions + 0.5)
        ax.set_xticklabels(xticks_labels, rotation=90, fontsize=self.cfg.heatmap_labels_xfontsize)
        ax.tick_params(axis ='x', which ='major')
        
        ax.set_yticks(np.arange(self.nsamples) + 0.5)
        ax.set_yticklabels(sequences, fontsize=self.cfg.heatmap_labels_yfontsize)
        ax.set_xlabel("Sequence position", fontsize=self.cfg.heatmap_labels_xfontsize)
       
        
        #row_index_map = {label: idx for idx, label in enumerate(seq_order_list)}


    def _plot_classes(self, fig, axes, seq_order_list):
        #Provide population information
        ax = axes['classes']
        cbaxes = axes['cbar_classes']

        samples = self._sequence_id_to_sample(seq_order_list)
        sequences = self._sequence_id_to_label(seq_order_list)

        group_dict = {sample: group for sample, group in zip(self.classes_df[self.classes_df.columns[0]], self.classes_df[self.classes_label]) if sample in samples}
        grouplabels = []
        for seq in seq_order_list:
            sample_id = seq.split('#')[0]
            if sample_id in group_dict:
                grouplabels.append(group_dict[sample_id])
            else:
                grouplabels.append(np.nan)
        unique_groups = sorted(list(set(grouplabels)))
        groupvalues = [unique_groups.index(group) if group in unique_groups else np.nan for group in grouplabels]

        if len(unique_groups) <= 10:
            pop_colormap = plt.get_cmap('tab10')
        else:
            if len(unique_groups) > 20:
                sys.stderr.write(f'WARNING: {len(unique_groups)} groups detected. Population colormap provides only 20 colors. Some groups will have the same color.\n'); sys.stderr.flush()
            pop_colormap = plt.get_cmap('tab20c')

        pop_heatmap(self.cfg, groupvalues, unique_groups, ax, cbaxes, pop_colormap, self.classes_label)



