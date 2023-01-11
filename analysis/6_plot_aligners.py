import numpy as np
# import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
# matplotlib.rcParams.update({'font.size': 14})

fig, axes = plt.subplots(1,2, figsize=(10,5))
sr_tools = ['BWA', 'Sentieon BWA', 'DRAGEN', 'GRAF', 'Novoalign']
sr_counts = [27, 5, 1, 1, 1]
sr_fracs = [x/sum(sr_counts) for x in sr_counts]

lr_tools = ['minimap2', 'pbmm2', 'winnowmap']
lr_counts = [18, 14, 2]
lr_fracs = [x/sum(lr_counts) for x in lr_counts]


axes[0].pie(sr_counts, labels=sr_tools, 
        autopct = lambda p: f'{(p/100)*sum(sr_counts):.0f}')
axes[1].pie(lr_counts, labels=lr_tools,
        autopct=lambda p: f'{(p/100)*sum(lr_counts):.0f}')
axes[0].set_title("Short Read Aligners")
axes[1].set_title("Long Read Aligners")
plt.tight_layout()

fig.savefig('img/6_aligners.png')
