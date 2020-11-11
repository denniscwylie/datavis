import plotnine as gg
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pandas.plotting import scatter_matrix
import seaborn as sns

plt.ion()
sns.set_style('white')

def ggplot(data, mapping, *args, **kwargs):
    out = gg.ggplot(data, mapping, *args, **kwargs)
    out += gg.theme_bw()
    return out

def rt(f, sep='\t', index_col=0, header=0, *args, **kwargs):
    return pd.read_csv(f, sep=sep, index_col=index_col, header=header,
                       *args, **kwargs)

logTpm = rt('gse75386_logtpm_filtered.tsv.gz')
trxAnnot = rt('Mus_musculus_GRCm38_82_TranscriptMap.tsv.gz')
trxAnnot = trxAnnot.reindex(logTpm.index)
annot = rt('gse75386_processed_annot.tsv')
annot.head()

cellType = annot['type']
simpleType = pd.Series({
    'CA1 cholecystokinin cell' : "Cck",
    'CA1 parvalbumin cell' : "Pvalb",
    'CA1 pyramidal cell' : "Pyramidal"
})[cellType]
simpleType.index = cellType.index
simpleType.value_counts()

gse75386 = pd.DataFrame({
    'class' : simpleType[logTpm.columns],
    'Pvalb' : logTpm.loc['ENSMUST00000005860'],
    'Cck' : logTpm.loc['ENSMUST00000035120'],
    'Gad1' : logTpm.loc['ENSMUST00000140478']
}, index = logTpm.columns)
gse75386.head()


## -----------------------------------------------------------------
## GSE75386 stripchart example
## -----------------------------------------------------------------
ggstrip = ggplot(
    data = gse75386,
    mapping = gg.aes(
        x = 'Gad1',
        y = 'class'
    )
)
ggstrip += gg.geom_point()
print(ggstrip)
# ggstrip.save('gse75386_gad1_stripchart_bw.pdf', format='pdf',
#              height=1, width=6)

## can also use seaborn for strip plotting...
plt.close()
# plt.figure(figsize=(6, 1))
sns.stripplot(data=gse75386, y='class', x='Gad1', color='black')
# plt.savefig('gse75386_gad1_stripchart_bw.pdf',
#             format='pdf', bbox_inches='tight')


## -----------------------------------------------------------------
## GSE75386 overplotted bars
## -----------------------------------------------------------------
plt.close()
ggbar = ggplot(gse75386, gg.aes(x='class', y='Gad1'))
ggbar += gg.geom_bar(alpha=0.1, position='identity', stat='identity')
ggbar += gg.coord_flip()
print(ggbar)
# ggbar.save('gse75386_gad1_barchart_id.pdf', format='pdf',
#            height=1, width=6)


## -----------------------------------------------------------------
## GSE75386 mean bars + SE lines
## -----------------------------------------------------------------
plt.close()
## use pandas functionality to compute stat transformations
gse75386means = gse75386[['class', 'Gad1']]\
                .groupby('class').agg(np.mean).iloc[:, 0]
gse75386ses = gse75386[['class', 'Gad1']]\
              .groupby('class').agg(lambda x: x.std() / np.sqrt(len(x)))\
              .iloc[:, 0]
gse75386stats = pd.DataFrame({'Gad1 (Mean)' : gse75386means,
                              'SE' : gse75386ses,
                              'ymin' : gse75386means - gse75386ses,
                              'ymax' : gse75386means + gse75386ses,
                              'class' : gse75386means.index.values})
ggbarse = ggplot(gse75386stats, gg.aes(x='class', y='Gad1 (Mean)')) +\
          gg.geom_bar(alpha=0.6, stat='identity') +\
          gg.geom_errorbar(mapping=gg.aes(ymin='ymin', ymax='ymax'), width=0.0001) +\
          gg.coord_flip()
print(ggbarse)
# ggbarse.save('gse75386_gad1_barchart_stat.pdf', format='pdf',
#              height=1, width=6)

## mean bars +/- standard error using seaborn
plt.close()
# plt.figure(figsize=(6, 1))
sns.barplot(data=gse75386, y='class', x='Gad1', color='slategray', ci=68)
# plt.savefig('gse75386_gad1_barchart_stat.pdf',
#             format='pdf', bbox_inches='tight')


## -----------------------------------------------------------------
## GSE75386 boxplot + stripchart
## -----------------------------------------------------------------
plt.close()
ggbox = ggplot(gse75386, gg.aes(x='class', y='Gad1')) +\
        gg.geom_boxplot(stat='boxplot', outlier_size=0.0001) +\
        gg.geom_point(alpha=0.5) +\
        gg.coord_flip()
print(ggbox)
# ggbox.save('gse75386_gad1_boxplot.pdf', format='pdf', height=1, width=6)

plt.close()
# plt.figure(figsize=(6, 1))
sns.boxplot(data=gse75386, y='class', x='Gad1', color='white')
sns.stripplot(data=gse75386, y='class', x='Gad1', color='black')
# plt.savefig('gse75386_gad1_boxplot.pdf',
#             format='pdf', bbox_inches='tight')


## -----------------------------------------------------------------
## GSE75386 scatterplot
## -----------------------------------------------------------------
plt.close()
ggscat = ggplot(
    gse75386,
    gg.aes(x='Gad1', y='Cck', color='class')
)
ggscat += gg.geom_point(alpha=0.75)
ggscat += gg.scale_color_manual(
        values=['darkslategray', 'goldenrod', 'lightseagreen'])
print(ggscat)
# ggscat.save('gse75386_cck_vs_gad1.pdf', format='pdf',
#             height=5, width=7)

def binarize(x, column, brk):
    out = pd.Series(['low ' + column]*x.shape[0], index=x.index)
    out.loc[x[column] > brk] = 'high ' + column
    return out

gse75386['Pvalb (cut)'] = binarize(gse75386, 'Pvalb', 5)
gse75386['Gad1 (cut)'] = binarize(gse75386, 'Gad1', 6)
gse75386.head()

plt.close()
ggscat = ggplot(
    gse75386,
    gg.aes(x='Gad1', y='Cck', color='class', size='Pvalb (cut)')
)
ggscat += gg.geom_point(alpha=0.75)
ggscat += gg.scale_color_manual(
        values=['darkslategray', 'goldenrod', 'lightseagreen'])
print(ggscat)
# ggscat.save('gse75386_cck_vs_gad1_sized_by_pvalb.pdf', format='pdf',
#             height=5, width=7)


## -----------------------------------------------------------------
## GSE75386 scatterplot + text layer
## -----------------------------------------------------------------
gse75386['odd'] = annot.loc[logTpm.columns, 'title']
## Pyramidal cells with low Gad1 and low Pvalb are not odd
gse75386.loc[(gse75386['class'] == 'Pyramidal') &
             (gse75386['Gad1 (cut)'] == 'low Gad1') &
             (gse75386['Pvalb (cut)'] == 'low Pvalb'),
             'odd'] = ''
## Pvalb cells with high Gad1 and high Pvalb are not odd
gse75386.loc[(gse75386['class'] == 'Pvalb') &
             (gse75386['Gad1 (cut)'] == 'high Gad1') &
             (gse75386['Pvalb (cut)'] == 'high Pvalb'),
             'odd'] = ''
## Cck cells with high Gad1 and low Pvalb are not odd
gse75386.loc[(gse75386['class'] == 'Cck') &
             (gse75386['Gad1 (cut)'] == 'high Gad1') &
             (gse75386['Pvalb (cut)'] == 'low Pvalb'),
             'odd'] = ''

plt.close()
ggscat = ggplot(
    gse75386,
    gg.aes(x='Gad1', y='Cck', color='class',
           size='Pvalb (cut)', label='odd')
)
ggscat += gg.scale_color_manual(
        values=['darkslategray', 'goldenrod', 'lightseagreen'])
ggscat += gg.geom_point(alpha=0.75)
ggscat += gg.geom_text(size=15, nudge_y= 0.65, show_legend=False)
ggscat += gg.scale_size_manual(values=[4., 1.5])
print(ggscat)
# ggscat.save('gse75386_cck_vs_gad1_sized_by_pvalb_odds_labeled.pdf',
#             format='pdf', height=5, width=7)

## alternately can generate similar scatterplot using seaborn
plt.close()
plt.figure(figsize=(5, 7))
p = sns.lmplot(data=gse75386.sort_index(), x='Gad1', y='Cck', hue='class',
               palette={'Cck' : 'darkslategray',
                        'Pvalb' : 'goldenrod',
                        'Pyramidal' : 'lightseagreen'},
               scatter_kws={'alpha': 0.75},
               legend=False, fit_reg=False)
# plt.savefig('gse75386_cck_vs_gad1.pdf',
#             format='pdf', bbox_inches='tight')

## adding text to seaborn plot is a bit more painful than for ggplot...
for i in range(gse75386.shape[0]):
    p.fig.text(0.12 + 0.8*gse75386['Gad1'].iloc[i] / gse75386['Gad1'].max(),
               0.15 + 0.78*gse75386['Cck'].iloc[i] / gse75386['Cck'].max(),
               gse75386['odd'].iloc[i])

# plt.savefig('gse75386_cck_vs_gad1_odds_labeled.pdf',
#             format='pdf', bbox_inches='tight')


## -----------------------------------------------------------------
## minard plotting
## -----------------------------------------------------------------
troops = rt('minard-troops.tsv', index_col=None)
cities = rt('minard-cities.tsv', index_col=None)
troops.head()

plt.close()
ggtroops = ggplot(troops, gg.aes('long', 'lat'))
ggtroops += gg.geom_path(gg.aes(
    size = 'survivors',
    color = 'direction',
    group = 'group'
))
print(ggtroops)
# ggtroops.save('ggplot_minard_troops.pdf', format='pdf',
#               height=4, width=12)

plt.close()
ggboth = ggtroops + gg.geom_text(
    gg.aes(label = 'city'),
    size = 12,
    data = cities
)
print(ggboth)
# ggboth.save('ggplot_minard_both.pdf', format='pdf',
#             height=4, width=12)

plt.close()
ggboth += gg.scale_size(
    range = [1, 10],
    breaks = [1e5, 2e5, 3e5]
)
ggboth = ggboth + gg.scale_color_manual(values = ["#d2b48c", "black"])
ggboth = ggboth + gg.xlab('') + gg.ylab('')
print(ggboth)
# ggboth.save('ggplot_minard_both_formatted.pdf', format='pdf',
#             height=4, width=12)


## -----------------------------------------------------------------
## Small multiples and facetting
## -----------------------------------------------------------------
anscombe = rt('anscombe_orig.tsv')
anscombe = pd.DataFrame({
    'x' : pd.Series(list(anscombe['x0'])*3 + list(anscombe['x4'])).values,
    'y' : pd.concat([anscombe['y1'], anscombe['y2'],
                     anscombe['y3'], anscombe['y4']]).values,
    'set' : 'set' + pd.Series(['1'] * anscombe.shape[0] +
                              ['2'] * anscombe.shape[0] +
                              ['3'] * anscombe.shape[0] +
                              ['4'] * anscombe.shape[0]).values
})
anscombe.head()

plt.close()
ggo = gg.ggplot(anscombe, gg.aes(x='x', y='y')) +\
      gg.facet_wrap('~ set') +\
      gg.geom_point() +\
      gg.theme_bw()
print(ggo)
# ggo.save('anscombe_points.pdf', format='pdf', height=5, width=5)

plt.close()
ggo += gg.stat_smooth(method='lm')
print(ggo)
## ggo.save('anscombe_lm.pdf', format='pdf', height=5, width=5)

## seaborn's lmplot function often useful in same situations
## one would want stat_smooth in R with ggplot2
plt.close()
sns.lmplot(data=anscombe, x='x', y='y', col='set')

plt.close()
sns.lmplot(data=anscombe, x='x', y='y', col='set', robust=True, ci=None)

plt.close()
sns.lmplot(data=anscombe, x='x', y='y', col='set', lowess=True)


## -----------------------------------------------------------------
## GSE75386 scatterplot matrix (a.k.a. pairs plot)
##-----------------------------------------------------------------
plt.close()
## for pairs plot / scatterplot matrix can either use seaborn:
sns.pairplot(gse75386[['Gad1', 'Pvalb', 'Cck', 'class']],
             hue='class',
             palette={'Cck' : 'darkslategray',
                      'Pvalb' : 'goldenrod',
                      'Pyramidal' : 'lightseagreen'})

plt.close()
## or pandas own scatter_matrix function:
scatter_matrix(gse75386[['Gad1', 'Pvalb', 'Cck', 'class']])
## neither one includes the categorical variable class, though


## -----------------------------------------------------------------
## clustered heatmap
##-----------------------------------------------------------------
plt.close()
## seaborn's clustermap function is similar to R's pheatmap
theGenes = [
    'Npy',
    'Cacna1d',
    'Hcn1',
    'Erbb4',
    'Gad1',
    'Pvalb',
    'Slc17a8',
    'Kcna1',
    'Bcl11b',
    'Chrm1',
    'Calb1',
    'Gabra1',
    'Cck',
    'S100a10',
    'Vip'
]
theGeneData = logTpm.loc[trxAnnot.loc[logTpm.index, 'gene_name'].isin(theGenes)]
## remove duplicate transcripts for same gene...
theGeneData = theGeneData[~theGeneData.isin([
    'ENSMUST00000094934',
    'ENSMUST00000141336'
])]
## use gene_name instead of ensembl transcript id to identify genes
theGeneData.index = trxAnnot.loc[theGeneData.index, 'gene_name']
heatmapData = theGeneData.subtract(theGeneData.mean(axis=1), axis=0)
heatmapColors = simpleType.loc[heatmapData.columns].copy()
heatmapColors.loc[heatmapColors == 'Cck'] = 'darkslategray'
heatmapColors.loc[heatmapColors == 'Pvalb'] = 'goldenrod'
heatmapColors.loc[heatmapColors == 'Pyramidal'] = 'lightseagreen'
sns.clustermap(heatmapData, method='average', col_colors=heatmapColors)
