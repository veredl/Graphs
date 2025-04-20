
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from scipy import stats
from scipy.stats import ttest_ind

sns.set(rc={"figure.dpi":300, 'savefig.dpi':300})
sns.set_context('notebook')
sns.set_style("ticks")

z_score_max = 2


"""--------------------area--------------------"""

area_xx_minus = pd.read_csv('area_xx_minus.txt')
area_xx_plus = pd.read_csv('area_xx_plus.txt')
area_yy_minus = pd.read_csv('area_yy_minus.txt')
area_yy_plus = pd.read_csv('area_yy_plus.txt')

area_xx_minus = area_xx_minus[(np.abs(stats.zscore(area_xx_minus)) < z_score_max).all(axis=1)]
area_xx_plus = area_xx_plus[(np.abs(stats.zscore(area_xx_plus)) < z_score_max).all(axis=1)]
area_yy_minus = area_yy_minus[(np.abs(stats.zscore(area_yy_minus)) < z_score_max).all(axis=1)]
area_yy_plus = area_yy_plus[(np.abs(stats.zscore(area_yy_plus)) < z_score_max).all(axis=1)]

area_all_Relative_territory = len(area_xx_minus) + len(area_xx_plus) + len(area_yy_minus) + len(area_yy_plus)

df_area = pd.DataFrame(np.zeros((area_all_Relative_territory,3)), columns = ["Chromosome number", 'Cell type', 'Relative territory'])

df_area.loc[:(len(area_xx_minus) + len(area_xx_plus) -1),"Chromosome number"] = 'Chr. #xx'
df_area.loc[len(area_xx_minus) + len(area_xx_plus):,"Chromosome number"] = 'Chr. #yy'

df_area.loc[:len(area_xx_plus) -1 ,'Cell type'] = 'WT'
df_area.loc[len(area_xx_plus):len(area_xx_minus) + len(area_xx_plus) -1,'Cell type'] = 'KO'
df_area.loc[len(area_xx_minus) + len(area_xx_plus):len(area_xx_minus) + len(area_xx_plus) + len(area_yy_plus) -1 ,'Cell type'] = 'WT'
df_area.loc[len(area_xx_minus) + len(area_xx_plus) + len(area_yy_plus):,'Cell type'] = 'KO'

df_area['Relative territory'] = area_xx_plus.append(area_xx_minus, ignore_index = True).append(area_yy_plus, ignore_index = True).append(area_yy_minus, ignore_index = True)



plt.figure(figsize= (7,7))
sns.set_theme(style="white", palette=None)
sns.set_theme(style="ticks",font="arial",font_scale=2.7)

Violin_area = sns.violinplot(x="Chromosome number", y="Relative territory", hue="Cell type", data=df_area, palette="muted",)
Violin_area.set(xlabel=None)
Violin_area.set_ylim(ymax=0.49)
ticks = [0,0.1,0.2,0.3,0.4]
Violin_area.set_yticks(ticks)
Violin_area.legend(frameon=False)


#statistics
# area_order_grouped = df_area.groupby(by=["Chromosome number", 'Cell type'])

# area_order_grouped.mean()['Relative territory']
# area_order_grouped.std()['Relative territory']

# area_order_grouped.median()['Relative territory']

# area_order_grouped.quantile(0.25)['Relative territory']
# area_order_grouped.quantile(0.75)['Relative territory']

# area_order_grouped.min()['Relative territory']
# area_order_grouped.max()['Relative territory']

# ttest_ind(area_xx_minus, area_xx_plus)
# ttest_ind(area_yy_minus, area_yy_plus)


"""--------------------circ--------------------"""

circ_xx_minus = pd.read_csv('circ_xx_minus.txt')
circ_xx_plus = pd.read_csv('circ_xx_plus.txt')
circ_yy_minus = pd.read_csv('circ_yy_minus.txt')
circ_yy_plus = pd.read_csv('circ_yy_plus.txt')

circ_xx_minus = circ_xx_minus[(np.abs(stats.zscore(circ_xx_minus)) < z_score_max).all(axis=1)]
circ_xx_plus = circ_xx_plus[(np.abs(stats.zscore(circ_xx_plus)) < z_score_max).all(axis=1)]
circ_yy_minus = circ_yy_minus[(np.abs(stats.zscore(circ_yy_minus)) < z_score_max).all(axis=1)]
circ_yy_plus = circ_yy_plus[(np.abs(stats.zscore(circ_yy_plus)) < z_score_max).all(axis=1)]

circ_all_Circularity = len(circ_xx_minus) + len(circ_xx_plus) + len(circ_yy_minus) + len(circ_yy_plus)

df_circ = pd.DataFrame(np.zeros((circ_all_Circularity,3)), columns = ["Chromosome number", 'Cell type', 'Circularity'])

df_circ.loc[:(len(circ_xx_minus) + len(circ_xx_plus) -1),"Chromosome number"] = 'Chr. #xx'
df_circ.loc[len(circ_xx_minus) + len(circ_xx_plus):,"Chromosome number"] = 'Chr. #yy'

df_circ.loc[:len(circ_xx_plus) -1 ,'Cell type'] = 'WT'
df_circ.loc[len(circ_xx_plus):len(circ_xx_minus) + len(circ_xx_plus) -1,'Cell type'] = 'KO'
df_circ.loc[len(circ_xx_minus) + len(circ_xx_plus):len(circ_xx_minus) + len(circ_xx_plus) + len(circ_yy_plus) -1 ,'Cell type'] = 'WT'
df_circ.loc[len(circ_xx_minus) + len(circ_xx_plus) + len(circ_yy_plus):,'Cell type'] = 'KO'

df_circ['Circularity'] = circ_xx_plus.append(circ_xx_minus, ignore_index = True).append(circ_yy_plus, ignore_index = True).append(circ_yy_minus, ignore_index = True)



plt.figure(figsize= (7,7))
sns.set_theme(style="white", palette=None)
sns.set_theme(style="ticks",font="arial",font_scale=2.7)

Violin_circ = sns.violinplot(x="Chromosome number", y="Circularity", hue="Cell type", data=df_circ, palette="muted")
Violin_circ.legend(frameon=False)
Violin_circ.set(xlabel=None)
# ticks = [1,2,3,4,5]
# Violin_circ.set_yticks(ticks)


#statistics
# circ_order_grouped = df_circ.groupby(by=["Chromosome number", 'Cell type'])

# circ_order_grouped.mean()['Circularity']
# circ_order_grouped.std()['Circularity']

# circ_order_grouped.median()['Circularity']

# circ_order_grouped.quantile(0.25)['Circularity']
# circ_order_grouped.quantile(0.75)['Circularity']

# circ_order_grouped.min()['Circularity']
# circ_order_grouped.max()['Circularity']

# ttest_ind(circ_xx_minus, circ_xx_plus)
# ttest_ind(circ_yy_minus, circ_yy_plus)


"""--------------------coloc--------------------"""

coloc_xx_yy_minus = pd.read_csv('coloc_xx_yy_minus.txt')
coloc_xx_yy_plus = pd.read_csv('coloc_xx_yy_plus.txt')

coloc_xx_yy_minus = coloc_xx_yy_minus[(np.abs(stats.zscore(coloc_xx_yy_minus)) < z_score_max).all(axis=1)]
coloc_xx_yy_plus = coloc_xx_yy_plus[(np.abs(stats.zscore(coloc_xx_yy_plus)) < z_score_max).all(axis=1)]

coloc_all_Relative_colocalization = len(coloc_xx_yy_minus) + len(coloc_xx_yy_plus) 
df_coloc = pd.DataFrame(np.zeros((coloc_all_Relative_colocalization,2)), columns = ['Cell type', 'Relative colocalization'])

df_coloc.loc[:(len(coloc_xx_yy_plus) -1),'Cell type'] = 'WT'
df_coloc.loc[len(coloc_xx_yy_plus):,'Cell type'] = 'KO'

df_coloc['Relative colocalization'] = coloc_xx_yy_plus.append(coloc_xx_yy_minus, ignore_index = True)



plt.figure(figsize= (7,7))
sns.set_theme(style="white", palette=None)
sns.set_theme(style="ticks",font="arial",font_scale=2.7)

Violin_coloc = sns.violinplot(x="Cell type", y="Relative colocalization", data=df_coloc, palette="muted")
# ticks = [0,0.05,0.1,0.15,0.2]
# Violin_coloc.set_yticks(ticks)

                    
#statistics
# coloc_order_grouped = df_coloc.groupby(by=['Cell type'])

# coloc_order_grouped.mean()['Relative colocalization']
# coloc_order_grouped.std()['Relative colocalization']

# coloc_order_grouped.median()['Relative colocalization']

# coloc_order_grouped.quantile(0.25)['Relative colocalization']
# coloc_order_grouped.quantile(0.75)['Relative colocalization']

# coloc_order_grouped.min()['Relative colocalization']
# coloc_order_grouped.max()['Relative colocalization']

# ttest_ind(coloc_xx_yy_minus, coloc_xx_yy_plus)






#subplot

sns.set_theme(style="ticks",font="arial",font_scale=3.3)

fig, axes = plt.subplots(1, 3, figsize=(30, 10))
fig.tight_layout()    

aa = sns.violinplot(ax=axes[0], x="Chromosome number", y="Relative territory", hue="Cell type", data=df_area, palette="muted") 
aa.legend(frameon=False, loc='upper left')
aa.set_ylim(ymax=0.49)
aa.set(xlabel=None)

bb = sns.violinplot(ax=axes[1], x="Chromosome number", y="Circularity", hue="Cell type", data=df_circ, palette="muted")
bb.get_legend().remove()
bb.set(xlabel=None)

cc = sns.violinplot(ax=axes[2], x="Cell type", y="Relative colocalization", data=df_coloc, palette="muted")

