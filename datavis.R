library(dplyr)     ## Hadley Wickham's grammar of data manipulation package
library(ggplot2)   ## Hadley Wickham's grammar of graphics plotting package
library(GGally)    ## ggplot version of pairs plot
#library(ggrepel)   ## non-overlapping geom_text labels
library(MASS)
library(pheatmap)  ## 'pretty' clustered heatmap
#library(robust)    ## robust regression tools
library(scales)    ## useful for creating/adjusting ggplot scales

ggplot = function(..., bw=TRUE, grid=FALSE) {
    ## change defaults for ggplot: white background, no grid
    out = ggplot2::ggplot(...)
    if (bw) {out = out + theme_bw()}
    if (!grid) {
        out = out + theme(
            panel.grid.minor = element_blank(),
            panel.grid.major = element_blank()
        )
    }
    invisible(out)
}

rt = function(..., sep='\t', row.names=1, header=TRUE, check.names=FALSE) {
    ## read.table with different defaults and shorter name
    return(read.table(..., sep=sep, row.names=row.names, header=header,
           check.names=check.names))
}

logTpm = rt('gse75386_logtpm_filtered.tsv.gz')
trxAnnot = rt('Mus_musculus_GRCm38_82_TranscriptMap.tsv.gz')
trxAnnot = trxAnnot[rownames(logTpm), ]
annot = droplevels(rt('gse75386_processed_annot.tsv')[colnames(logTpm), ])
head(annot)

type = structure(annot$type, names=rownames(annot))
simpleType = factor(c(
    `CA1 cholecystokinin cell` = "Cck",
    `CA1 parvalbumin cell` = "Pvalb",
    `CA1 pyramidal cell` = "Pyramidal"
))[as.character(type)]
names(simpleType) = names(type)
table(simpleType)

gse75386 = data.frame(
    row.names = colnames(logTpm),
    class = simpleType[colnames(logTpm)],
    Pvalb = as.numeric(logTpm['ENSMUST00000005860', ]),
    Cck = as.numeric(logTpm['ENSMUST00000035120', ]),
    Gad1 = as.numeric(logTpm['ENSMUST00000140478', ]),
    dummy = ''
)
head(gse75386)


## -----------------------------------------------------------------
## GSE75386 stripchart example
## -----------------------------------------------------------------
ggstrip = ggplot(
    data = gse75386,
    mapping = aes(
        x = Gad1,
        y = class
    )
)
ggstrip = ggstrip + geom_point()
## pdf('gse75386_gad1_stripchart_bw.pdf', h=1, w=6)
print(ggstrip)
## garbage = dev.off()


## -----------------------------------------------------------------
## GSE75386 overplotted bars
## -----------------------------------------------------------------
ggbar = ggplot(gse75386, aes(x=class, y=Gad1))
ggbar = ggbar + geom_bar(alpha=0.1,
                         position='identity', stat='identity')
ggbar = ggbar + coord_flip()
## pdf('gse75386_gad1_barchart_id.pdf', h=1, w=6)
print(ggbar)
## garbage = dev.off()


## -----------------------------------------------------------------
## GSE75386 mean bars + SE lines
## -----------------------------------------------------------------
## use dplyr functionality to compute stat transformations
gse75386stats = gse75386 %>%
                group_by(class) %>%
                summarize(
                    `Gad1 (Mean)` = mean(Gad1),
                    SE = sd(Gad1) / sqrt(length(Gad1))
                )
ggbarse = ggplot(gse75386stats, aes(x=class, y=`Gad1 (Mean)`))
ggbarse = ggbarse + geom_bar(alpha=0.6, stat='identity')
ggbarse = ggbarse + geom_errorbar(aes(ymin=`Gad1 (Mean)` - SE,
                                      ymax=`Gad1 (Mean)` + SE),
                                  width=0)
ggbarse = ggbarse + coord_flip()
## pdf('gse75386_gad1_barchart_stat.pdf', h=1, w=6)
print(ggbarse)
## garbage = dev.off()


## -----------------------------------------------------------------
## GSE75386 boxplot + stripchart
## -----------------------------------------------------------------
ggbox = ggplot(gse75386, aes(x=class, y=Gad1))
ggbox = ggbox + geom_boxplot(stat='boxplot',
                             outlier.size=0)
ggbox = ggbox + geom_point(alpha=0.5)
ggbox = ggbox + coord_flip()
## pdf('gse75386_gad1_boxplot.pdf', h=1, w=6)
print(ggbox)
## garbage = dev.off()


## -----------------------------------------------------------------
## GSE75386 scatterplot
## -----------------------------------------------------------------
ggscat = ggplot(
    gse75386,
    aes(x=Gad1, y=Cck, color=class)
)
ggscat = ggscat + geom_point(alpha=0.75)
ggscat = ggscat + scale_color_manual(
values=c('darkslategray', 'goldenrod', 'lightseagreen'))
## pdf('gse75386_cck_vs_gad1.pdf', h=4, w=5.5)
print(ggscat)
## garbage = dev.off()

binarize = function(x, column, brk) {
    out = cut(x[ , column], breaks=c(-Inf, brk, Inf))
    levels(out) = paste(c('low', 'high'), column)
    names(out) = rownames(x)
    return(out)
}
gse75386$'Pvalb (cut)' = binarize(gse75386, 'Pvalb', 5)
gse75386$'Gad1 (cut)' = binarize(gse75386, 'Gad1', 6)
head(gse75386)

ggscat = ggplot(
    gse75386,
    aes(x=Gad1, y=Cck, color=class, size=`Pvalb (cut)`)
)
ggscat = ggscat + geom_point(alpha=0.75)
ggscat = ggscat + scale_color_manual(
        values=c('darkslategray', 'goldenrod', 'lightseagreen'))
ggscat = ggscat + scale_size_manual(values=c(2, 4))
## pdf('gse75386_cck_vs_gad1_sized_by_pvalb.pdf', h=4, w=5.6)
print(ggscat)
## garbage = dev.off()


## -----------------------------------------------------------------
## GSE75386 scatterplot + text layer
## -----------------------------------------------------------------
gse75386$odd = annot[colnames(logTpm), 'title']
## Pyramidal cells with low Gad1 and low Pvalb are not odd
gse75386[gse75386$class == 'Pyramidal' &
         gse75386$'Gad1 (cut)' == 'low Gad1' &
         gse75386$'Pvalb (cut)' == 'low Pvalb',
         'odd'] = NA
## Pvalb cells with high Gad1 and high Pvalb are not odd
gse75386[gse75386$class == 'Pvalb' &
         gse75386$'Gad1 (cut)' == 'high Gad1' &
         gse75386$'Pvalb (cut)' == 'high Pvalb',
         'odd'] = NA
## Cck cells with high Gad1 and low Pvalb are not odd
gse75386[gse75386$class == 'Cck' &
         gse75386$'Gad1 (cut)' == 'high Gad1' &
         gse75386$'Pvalb (cut)' == 'low Pvalb',
         'odd'] = NA

ggscat = ggplot(
    gse75386,
    aes(x=Gad1, y=Cck, color=class, size=`Pvalb (cut)`)
)
ggscat = ggscat + geom_point(alpha=0.75)
ggscat = ggscat + scale_color_manual(
        values=c('darkslategray', 'goldenrod', 'lightseagreen'))
ggscat = ggscat + scale_size_manual(values=c(2, 4))
## pdf('gse75386_cck_vs_gad1_sized_by_pvalb.pdf', h=4, w=5.6)
print(ggscat)
## garbage = dev.off()

ggscat = ggscat + geom_text(
    aes(label=odd),
    vjust = -0.85,
    size = 3,
    show.legend = FALSE
)
## pdf('gse75386_cck_vs_gad1_sized_by_pvalb_odds_labeled.pdf', h=4, w=5.6)
print(ggscat)
## garbage = dev.off()


## -----------------------------------------------------------------
## minard plotting
## -----------------------------------------------------------------
troops = rt('minard-troops.tsv', row.names=NULL)
cities = rt('minard-cities.tsv', row.names=NULL)
head(troops)

ggtroops = ggplot(troops, aes(long, lat))
ggtroops = ggtroops + geom_path(aes(
    size = survivors,
    color = direction,
    group = group
))
## pdf('ggplot_minard_troops.pdf', h=4, w=12)
print(ggtroops)
## garbage = dev.off()

ggboth = ggtroops + geom_text(
    aes(label = city),
    size = 4,
    data = cities
)
## pdf('ggplot_minard_both.pdf', h=4, w=12)
print(ggboth)
## garbage = dev.off()

ggboth = ggboth + scale_size(
    range = c(1, 10),
    breaks = c(1, 2, 3) * 10^5,
    labels = comma(c(1, 2, 3) * 10^5)
)
ggboth = ggboth + scale_color_manual(values = c("#d2b48c","black"))
ggboth = ggboth + xlab(NULL) + ylab(NULL)
## pdf('ggplot_minard_both_formatted.pdf', h=4, w=12)
print(ggboth)
## garbage = dev.off()


## -----------------------------------------------------------------
## Small multiples and facetting
## -----------------------------------------------------------------
anscombe = rt('anscombe_orig.tsv')
anscombe = data.frame(
    x = c(rep(anscombe$x0, 3), anscombe$x4),
    y = c(anscombe$y1, anscombe$y2, anscombe$y3, anscombe$y4),
    set = paste('set', c(rep('1', nrow(anscombe)), rep('2', nrow(anscombe)),
    rep('3', nrow(anscombe)), rep('4', nrow(anscombe))))
)
head(anscombe)

ggo = ggplot(anscombe, aes(x=x, y=y))
ggo = ggo + facet_wrap(~ set)
ggo = ggo + geom_point()
## ggo = ggo + ylim(3.1, 15.25)
## pdf('anscombe_points.pdf', h=5, w=5)
print(ggo)
## garbage = dev.off()

ggo = ggo + stat_smooth(method=lm)
## pdf('anscombe_lm.pdf', h=5, w=5)
print(ggo)
## garbage = dev.off()

lmDescription = function(x, y=NULL, method=lm, digits=2) {
    if (length(y) == 0) {y=x$y; x=x$x}
    lmo = method(y ~ x)
    return(paste0(
        'y = ',
        round(coef(lmo)[2], digits),
        '*x + ',
        round(coef(lmo)[1], digits),
        ' + e \n',
        'Var[e] = ',
        round(summary(lmo)$sigma^2, digits)
    ))
}
lmDescriptions = function(data, method=lm, digits=2) {
    lmDescs = lapply(
        split(data[ , c('x', 'y')], data$set),
        FUN = lmDescription,
        method = method
    )
    return(data.frame(
        x = 13.5,
        y = 5.25,
        set = names(lmDescs),
        text = as.character(unlist(lmDescs)),
        stringsAsFactors = FALSE
    ))
}

ggo = ggplot(anscombe, aes(x=x, y=y))
ggo = ggo + geom_point()
ggo = ggo + geom_text(aes(label=text), data=lmDescriptions(anscombe))
ggo = ggo + stat_smooth(method=lm)
ggo = ggo + facet_wrap(~ set)
## pdf('anscombe_lm_text.pdf', h=5, w=5)
print(ggo)
## garbage = dev.off()

ggo = ggplot(anscombe, aes(x=x, y=y))
ggo = ggo + geom_point()
ggo = ggo + geom_text(aes(label=text),
data=lmDescriptions(anscombe, rlm))
ggo = ggo + stat_smooth(method=rlm)
ggo = ggo + facet_wrap(~ set)
## pdf('anscombe_rlm.pdf', h=5, w=5)
print(ggo)
## garbage = dev.off()

ggo = ggplot(anscombe, aes(x=x, y=y))
ggo = ggo + geom_point()
ggo = ggo + stat_smooth(method=loess)
ggo = ggo + facet_wrap(~ set)
## pdf('anscombe_loess.pdf', h=5, w=5)
print(ggo)
## garbage = dev.off()

wrappedHistogram = function(data, mapping, ...) {
    ggobj = eval(parse(text=paste0(
            'ggplot(data, aes(x=`', as.character(mapping[['x']]), '`))')))
    if ("y" %in% names(mapping)) {
        ggobj = eval(parse(text=paste0(
            'ggobj + facet_grid(`',
            as.character(mapping[['y']]),
            '` ~ .)'
        )))
    }
    ggobj = ggobj + geom_histogram(aes(fill=class, y=..density..), binwidth=2)
    ggobj = ggobj + scale_fill_manual(
            values=c('darkslategray', 'goldenrod', 'lightseagreen'))
    return(ggobj)
}
wrappedBar = function(data, mapping, ...) {
    ggobj = ggplot(data, mapping)
    ggobj = ggobj + geom_bar(aes(fill=class))
    ggobj = ggobj + scale_fill_manual(
            values=c('darkslategray', 'goldenrod', 'lightseagreen'))
    return(ggobj)
}
wrappedBox = function(data, mapping, ...) {
    ggobj = ggplot(data, mapping)
    ggobj = ggobj + geom_boxplot(aes(color=class), outlier.size=1)
    ggobj = ggobj + scale_color_manual(
            values=c('darkslategray', 'goldenrod', 'lightseagreen'))
    return(ggobj)
}
wrappedLoess = function(data, mapping,
        method='rlm', method.args=list(deg=1), span=1.5, ...) {
    ggobj = ggplot(data, mapping)
    ggobj = ggobj + geom_point(size=1, alpha=0.75, aes(color=class))
    ggobj = ggobj + scale_color_manual(
            values=c('darkslategray', 'goldenrod', 'lightseagreen'))
    ggobj = ggobj + geom_smooth(aes(color=class),
            method=method, method.args=method.args, span=span, se=FALSE, ...)
    return(ggobj)
}


## -----------------------------------------------------------------
## GSE75386 scatterplot matrix (a.k.a. pairs plot)
##-----------------------------------------------------------------
## pdf('gse75386_pairs.pdf', h=5, w=5)
ggpairs(
    gse75386[ , c('Gad1', 'Pvalb', 'Cck', 'class')],
    diag = list(
        continuous = wrappedHistogram,
        discrete = wrappedBar
    ),
    lower = list(
        continuous = wrappedLoess,
        combo = wrappedHistogram
    ),
    upper = list(
        combo = wrappedBox
    )
) + theme_bw() + theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
)
## garbage = dev.off()


## -----------------------------------------------------------------
## clustered heatmap
##-----------------------------------------------------------------
theGenes = c(
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
)
theGeneData = logTpm[trxAnnot[rownames(logTpm), 'gene_name'] %in% theGenes, ]
## remove duplicate transcripts for same gene...
theGeneData = theGeneData[!rownames(theGeneData) %in% c(
    'ENSMUST00000094934',
    'ENSMUST00000141336'
), ]
## use gene_name instead of ensembl transcript id to identify genes
rownames(theGeneData) = trxAnnot[rownames(theGeneData), 'gene_name']
heatmapData = sweep(theGeneData, 1, rowMeans(theGeneData), `-`)
## pdf('gse75386_int_gene_heatmap.pdf', h=3.25, w=12, onefile=FALSE)
pheatmap(
    heatmapData,
    annotation_col = data.frame(
        row.names = colnames(heatmapData),
        type = simpleType[colnames(heatmapData)]
    ),
    annotation_colors = list(type=c(
        Cck = 'darkslategray',
        Pvalb = 'goldenrod',
        Pyramidal = 'lightseagreen'
    )),
    cluster_method = 'mcquitty',
    show_colnames = FALSE
)
## garbage = dev.off()
