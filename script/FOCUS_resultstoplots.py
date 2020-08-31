#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
plt.switch_backend('agg')
import glob

#data reading
file= glob.glob("*chr_all*")
dat = pd.read_csv(file[0], header=0,sep="\t", engine="python")
df = pd.DataFrame({"ID": dat["mol_name"],
                   "CHR": dat["chrom"],
                "CHR:pos": 100000000000*(dat["chrom"]) + dat["tx_start"],
                  "Z": dat["twas_z"],
                   "pip": dat["pip"]})
df = df.dropna(how='any',axis=0)
df2 = pd.DataFrame({"ID": dat["mol_name"],
                   "CHR": dat["chrom"],
                "CHR:pos": 100000000000*(dat["chrom"]) + dat["tx_start"],
                  "Z": dat["twas_z"],
                   "pip": dat["pip"]})
df2 = df2.dropna(how='any',axis=0)

#plot
df.CHR = df.CHR.astype('category')
df = df.sort_values('CHR:pos')
df['ind'] = range(len(df))
df_grouped = df.groupby(('CHR'))

df2.CHR= df2.CHR.astype('category')
df2 = df2.sort_values('CHR:pos')
df2['ind'] = range(len(df2))
df2_grouped = df2.groupby(('CHR'))
fig = plt.figure(figsize=(15,9),dpi=200)
ax = fig.add_subplot(111)
x_labels = []
x_labels_pos = []
for num, (name, group) in enumerate(df_grouped):
    group.plot(kind='scatter', x='ind', y='Z',color="darkgrey", ax=ax, s=10)
    x_labels.append(name)
    x_labels_pos.append((group['ind'].iloc[-1] - (group['ind'].iloc[-1] - group['ind'].iloc[0])/2))
    ax.set_xticks(x_labels_pos)
    ax.set_xticklabels(x_labels)
    ax.set_xlim([0, len(df)])
    plt.ylim(min(df2.Z)-2, max(df2.Z)+2)
    ax.set_xlabel('Chromosome')
    ax.set_ylabel('TWAS Z')
    xmin, xmax = 0, len(df)


fname = file[0]+".png"
fname2 = fname.replace(".tsv","")
plt.savefig(fname2, dpi = 64, facecolor = "white", tight_layout = True)

