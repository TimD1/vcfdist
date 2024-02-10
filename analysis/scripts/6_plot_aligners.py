import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams.update(mpl.rcParamsDefault)
mpl.rcParams['text.usetex'] = True

fig, axes = plt.subplots(1,2, figsize=(10,5))
sr_tools = ['BWA', 'Sentieon BWA', 'DRAGEN', 'GRAF', 'Novoalign']
sr_counts = [27, 5, 1, 1, 1]
sr_fracs = [x/sum(sr_counts) for x in sr_counts]
sr_colors = ["C7", "C7", "C7", "C7", "C4"]

lr_tools = [r'minimap2 \texttt{map-ont}', r'minimap2 \texttt{map-pb}',
        'pbmm2', 'winnowmap']
lr_counts = [9, 9, 14, 2]
lr_fracs = [x/sum(lr_counts) for x in lr_counts]
lr_colors = ["C2", "C3", "C6", "C2"]


axes[0].pie(sr_counts, labels=sr_tools, colors=sr_colors,
        autopct = lambda p: f'{(p/100)*sum(sr_counts):.0f}')
axes[1].pie(lr_counts, labels=lr_tools, colors=lr_colors,
        autopct=lambda p: f'{(p/100)*sum(lr_counts):.0f}')
axes[0].set_title(r"\textbf{\Large Short Read Aligners}")
axes[1].set_title(r"\textbf{\Large Long Read Aligners}")
plt.tight_layout()

fig.savefig('img/6_aligners.png')
