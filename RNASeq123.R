# RNASeq123 9/8/20

## ----setup, message=FALSE, echo = FALSE-------------------------------------------------
library(BiocStyle)
library(knitr)
options(digits=3)
options(width=90)

## ----setup2, message=FALSE, eval=TRUE---------------------------------------------------
library(limma)
library(Glimma)
library(edgeR)
library(Mus.musculus)

## ----downloadData, eval=TRUE------------------------------------------------------------
url <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE63310&format=file"
utils::download.file(url, destfile="GSE63310_RAW.tar", mode="wb") 
utils::untar("GSE63310_RAW.tar", exdir = ".")
files <- c("GSM1545535_10_6_5_11.txt", "GSM1545536_9_6_5_11.txt", "GSM1545538_purep53.txt",
           "GSM1545539_JMS8-2.txt", "GSM1545540_JMS8-3.txt", "GSM1545541_JMS8-4.txt",
           "GSM1545542_JMS8-5.txt", "GSM1545544_JMS9-P7c.txt", "GSM1545545_JMS9-P8c.txt")
for(i in paste(files, ".gz", sep=""))
  R.utils::gunzip(i, overwrite=TRUE)

## ----import1----------------------------------------------------------------------------
files <- c("GSM1545535_10_6_5_11.txt", "GSM1545536_9_6_5_11.txt", 
           "GSM1545538_purep53.txt", "GSM1545539_JMS8-2.txt", 
           "GSM1545540_JMS8-3.txt", "GSM1545541_JMS8-4.txt", 
           "GSM1545542_JMS8-5.txt", "GSM1545544_JMS9-P7c.txt", 
           "GSM1545545_JMS9-P8c.txt")
read.delim(files[1], nrow=5)

## ----import2----------------------------------------------------------------------------
x <- readDGE(files, columns=c(1,3))
class(x)
dim(x)

## ----annotatesamples--------------------------------------------------------------------
samplenames <- substring(colnames(x), 12, nchar(colnames(x)))
samplenames
colnames(x) <- samplenames
group <- as.factor(c("LP", "ML", "Basal", "Basal", "ML", "LP", 
                     "Basal", "ML", "LP"))
x$samples$group <- group
lane <- as.factor(rep(c("L004","L006","L008"), c(3,4,2)))
x$samples$lane <- lane
x$samples

## ----annotategenes, message=FALSE-------------------------------------------------------
geneid <- rownames(x)
genes <- select(Mus.musculus, keys=geneid, columns=c("SYMBOL", "TXCHROM"), 
                keytype="ENTREZID")
head(genes)

## ----removedups-------------------------------------------------------------------------
genes <- genes[!duplicated(genes$ENTREZID),]

## ----assigngeneanno---------------------------------------------------------------------
x$genes <- genes
x

## ----cpm--------------------------------------------------------------------------------
cpm <- cpm(x)
lcpm <- cpm(x, log=TRUE, prior.count=2)

## ----lcpm-------------------------------------------------------------------------------
L <- mean(x$samples$lib.size) * 1e-6
M <- median(x$samples$lib.size) * 1e-6
c(L, M)
summary(lcpm)

## ----zeroes-----------------------------------------------------------------------------
table(rowSums(x$counts==0)==9)

## ----filter-----------------------------------------------------------------------------
keep.exprs <- filterByExpr(x, group=group)
x <- x[keep.exprs,, keep.lib.sizes=FALSE]
dim(x)

## ----filterplot1, fig.height=4, fig.width=8, fig.cap="æ¯ä¸ªæ ·æœ¬è¿‡æ»¤å‰çš„åŽŸå§‹æ•°æ®ï¼ˆAï¼‰å’Œè¿‡æ»¤åŽï¼ˆBï¼‰çš„æ•°æ®çš„log-CPMå€¼å¯†åº¦ã€‚ç«–ç›´è™šçº¿æ ‡å‡ºäº†è¿‡æ»¤æ­¥éª¤ä¸­æ‰€ç”¨é˜ˆå€¼ï¼ˆç›¸å½“äºŽCPMå€¼ä¸ºçº¦0.2ï¼‰ã€‚"----
lcpm.cutoff <- log2(10/M + 2/L)
library(RColorBrewer)
nsamples <- ncol(x)
col <- brewer.pal(nsamples, "Paired")
par(mfrow=c(1,2))
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
title(main="A. Raw data", xlab="Log-cpm")
abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", samplenames, text.col=col, bty="n")
lcpm <- cpm(x, log=TRUE)
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
title(main="B. Filtered data", xlab="Log-cpm")
abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", samplenames, text.col=col, bty="n")

## ----normalize--------------------------------------------------------------------------
x <- calcNormFactors(x, method = "TMM")
x$samples$norm.factors

## ----normalizemodifieddata--------------------------------------------------------------
x2 <- x
x2$samples$norm.factors <- 1
x2$counts[,1] <- ceiling(x2$counts[,1]*0.05)
x2$counts[,2] <- x2$counts[,2]*5

## ----plotmodifieddata, fig.height=4, fig.width=8, fig.cap="æ ·ä¾‹æ•°æ®ï¼šlog-CPMå€¼çš„ç®±çº¿å›¾å±•ç¤ºäº†æœªç»å½’ä¸€åŒ–çš„æ•°æ®ï¼ˆAï¼‰åŠå½’ä¸€åŒ–åŽçš„æ•°æ®ï¼ˆBï¼‰ä¸­æ¯ä¸ªæ ·æœ¬çš„è¡¨è¾¾åˆ†å¸ƒã€‚æ•°æ®é›†ç»è¿‡è°ƒæ•´ï¼Œæ ·æœ¬1å’Œ2ä¸­çš„è¡¨è¾¾è®¡æ•°åˆ†åˆ«è¢«ç¼©æ”¾åˆ°å…¶åŽŸå§‹å€¼çš„5%å’Œ500%ã€‚"----
par(mfrow=c(1,2))
lcpm <- cpm(x2, log=TRUE)
boxplot(lcpm, las=2, col=col, main="")
title(main="A. Example: Unnormalised data",ylab="Log-cpm")
x2 <- calcNormFactors(x2)  
x2$samples$norm.factors
lcpm <- cpm(x2, log=TRUE)
boxplot(lcpm, las=2, col=col, main="")
title(main="B. Example: Normalised data",ylab="Log-cpm")

## ----MDS1, fig.height=4, fig.width=8, fig.cap="log-CPMå€¼åœ¨ç»´åº¦1å’Œ2çš„MDSå›¾ï¼Œä»¥æ ·å“åˆ†ç»„ä¸Šè‰²å¹¶æ ‡è®°ï¼ˆAï¼‰å’Œç»´åº¦3å’Œ4çš„MDSå›¾ï¼Œä»¥æµ‹åºé“ä¸Šè‰²å¹¶æ ‡è®°ï¼ˆBï¼‰ã€‚å›¾ä¸­çš„è·ç¦»å¯¹åº”äºŽæœ€ä¸»è¦çš„å€æ•°å˜åŒ–ï¼ˆfold changeï¼‰ï¼Œé»˜è®¤æƒ…å†µä¸‹ä¹Ÿå°±æ˜¯å‰500ä¸ªåœ¨æ¯å¯¹æ ·å“ä¹‹é—´å·®å¼‚æœ€å¤§çš„åŸºå› çš„å¹³å‡ï¼ˆå‡æ–¹æ ¹ï¼‰log2å€æ•°å˜åŒ–ã€‚"----
lcpm <- cpm(x, log=TRUE)
par(mfrow=c(1,2))
col.group <- group
levels(col.group) <-  brewer.pal(nlevels(col.group), "Set1")
col.group <- as.character(col.group)
col.lane <- lane
levels(col.lane) <-  brewer.pal(nlevels(col.lane), "Set2")
col.lane <- as.character(col.lane)
plotMDS(lcpm, labels=group, col=col.group)
title(main="A. Sample groups")
plotMDS(lcpm, labels=lane, col=col.lane, dim=c(3,4))
title(main="B. Sequencing lanes")

## ----GlimmaMDSplot----------------------------------------------------------------------
glMDSPlot(lcpm, labels=paste(group, lane, sep="_"), 
          groups=x$samples[,c(2,5)], launch=FALSE)

## ----design-----------------------------------------------------------------------------
design <- model.matrix(~0+group+lane)
colnames(design) <- gsub("group", "", colnames(design))
design

## ----contrasts--------------------------------------------------------------------------
contr.matrix <- makeContrasts(
  BasalvsLP = Basal-LP, 
  BasalvsML = Basal - ML, 
  LPvsML = LP - ML, 
  levels = colnames(design))
contr.matrix

## ----voom, fig.height=4, fig.width=8, fig.cap="å›¾ä¸­ç»˜åˆ¶äº†æ¯ä¸ªåŸºå› çš„å‡å€¼ï¼ˆxè½´ï¼‰å’Œæ–¹å·®ï¼ˆyè½´ï¼‰ï¼Œæ˜¾ç¤ºäº†åœ¨è¯¥æ•°æ®ä¸Šä½¿ç”¨`voom`å‰å®ƒä»¬ä¹‹é—´çš„ç›¸å…³æ€§ï¼ˆå·¦ï¼‰ï¼Œä»¥åŠå½“è¿ç”¨`voom`çš„ç²¾ç¡®æƒé‡åŽè¿™ç§è¶‹åŠ¿æ˜¯å¦‚ä½•æ¶ˆé™¤çš„ï¼ˆå³ï¼‰ã€‚å·¦ä¾§çš„å›¾æ˜¯ä½¿ç”¨`voom`å‡½æ•°ç»˜åˆ¶çš„ï¼Œå®ƒä¸ºè¿›è¡Œlog-CPMè½¬æ¢åŽçš„æ•°æ®æ‹Ÿåˆçº¿æ€§æ¨¡åž‹ä»Žè€Œæå–æ®‹å·®æ–¹å·®ã€‚ç„¶åŽï¼Œå¯¹æ–¹å·®å–å¹³æ–¹æ ¹ï¼ˆæˆ–å¯¹æ ‡å‡†å·®å–å¹³æ–¹æ ¹ï¼‰ï¼Œå¹¶ç›¸å¯¹æ¯ä¸ªåŸºå› çš„å¹³å‡è¡¨è¾¾ä½œå›¾ã€‚å‡å€¼é€šè¿‡å¹³å‡è®¡æ•°åŠ ä¸Š2å†è¿›è¡Œlog2è½¬æ¢è®¡ç®—å¾—åˆ°ã€‚å³ä¾§çš„å›¾ä½¿ç”¨`plotSA`ç»˜åˆ¶äº†log2æ®‹å·®æ ‡å‡†å·®ä¸Žlog-CPMå‡å€¼çš„å…³ç³»ã€‚å¹³å‡log2æ®‹å·®æ ‡å‡†å·®ç”±æ°´å¹³è“çº¿æ ‡å‡ºã€‚åœ¨è¿™ä¸¤å¹…å›¾ä¸­ï¼Œæ¯ä¸ªé»‘ç‚¹è¡¨ç¤ºä¸€ä¸ªåŸºå› ï¼Œçº¢çº¿ä¸ºå¯¹è¿™äº›ç‚¹çš„æ‹Ÿåˆã€‚"----
par(mfrow=c(1,2))
v <- voom(x, design, plot=TRUE)
v
vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
plotSA(efit, main="Final model: Mean-variance trend")

## ----decidetests------------------------------------------------------------------------
summary(decideTests(efit))

## ----treat------------------------------------------------------------------------------
tfit <- treat(vfit, lfc=1)
dt <- decideTests(tfit)
summary(dt)

## ----venn, fig.height=6, fig.width=6, fig.cap="éŸ¦æ©å›¾å±•ç¤ºäº†ä»…basalå’ŒLPï¼ˆå·¦ï¼‰ã€ä»…basalå’ŒMLï¼ˆå³ï¼‰çš„å¯¹æ¯”çš„DEåŸºå› æ•°é‡ï¼Œè¿˜æœ‰ä¸¤ç§å¯¹æ¯”ä¸­å…±åŒçš„DEåŸºå› æ•°é‡ï¼ˆä¸­ï¼‰ã€‚åœ¨ä»»ä½•å¯¹æ¯”ä¸­å‡ä¸å·®å¼‚è¡¨è¾¾çš„åŸºå› æ•°é‡æ ‡äºŽå³ä¸‹ã€‚"----
de.common <- which(dt[,1]!=0 & dt[,2]!=0)
length(de.common)
head(tfit$genes$SYMBOL[de.common], n=20)
vennDiagram(dt[,1:2], circle.col=c("turquoise", "salmon"))
write.fit(tfit, dt, file="results.txt")

## ----toptables--------------------------------------------------------------------------
basal.vs.lp <- topTreat(tfit, coef=1, n=Inf)
basal.vs.ml <- topTreat(tfit, coef=2, n=Inf)
head(basal.vs.lp)
head(basal.vs.ml)

## ----MDplot, fig.keep='none'------------------------------------------------------------
plotMD(tfit, column=1, status=dt[,1], main=colnames(tfit)[1], 
       xlim=c(-8,13))

## ----GlimmaMDplot-----------------------------------------------------------------------
glMDPlot(tfit, coef=1, status=dt, main=colnames(tfit)[1],
         side.main="ENTREZID", counts=lcpm, groups=group, launch=FALSE)

## ----heatmap, fig.height=8, fig.width=5, fig.cap="åœ¨basalå’ŒLPçš„å¯¹æ¯”ä¸­å‰100ä¸ªDEåŸºå› log-CPMå€¼çš„çƒ­å›¾ã€‚ç»è¿‡ç¼©æ”¾è°ƒæ•´åŽï¼Œæ¯ä¸ªåŸºå› ï¼ˆæ¯è¡Œï¼‰çš„è¡¨è¾¾å‡å€¼ä¸º0ï¼Œå¹¶ä¸”æ ‡å‡†å·®ä¸º1ã€‚ç»™å®šåŸºå› ç›¸å¯¹é«˜è¡¨è¾¾çš„æ ·æœ¬è¢«æ ‡è®°ä¸ºçº¢è‰²ï¼Œç›¸å¯¹ä½Žè¡¨è¾¾çš„æ ·æœ¬è¢«æ ‡è®°ä¸ºè“è‰²ã€‚æµ…è‰²å’Œç™½è‰²ä»£è¡¨ä¸­ç­‰è¡¨è¾¾æ°´å¹³çš„åŸºå› ã€‚æ ·æœ¬å’ŒåŸºå› å·²é€šè¿‡åˆ†å±‚èšç±»çš„æ–¹æ³•é‡æ–°æŽ’åºã€‚å›¾ä¸­æ˜¾ç¤ºæœ‰æ ·æœ¬èšç±»çš„æ ‘çŠ¶å›¾ã€‚", message=FALSE----
library(gplots)
basal.vs.lp.topgenes <- basal.vs.lp$ENTREZID[1:100]
i <- which(v$genes$ENTREZID %in% basal.vs.lp.topgenes)
mycol <- colorpanel(1000,"blue","white","red")
heatmap.2(lcpm[i,], scale="row",
          labRow=v$genes$SYMBOL[i], labCol=group, 
          col=mycol, trace="none", density.info="none", 
          margin=c(8,6), lhei=c(2,10), dendrogram="column")

## ----camera-----------------------------------------------------------------------------
load(system.file("extdata", "mouse_c2_v5p1.rda", package = "RNAseq123"))
idx <- ids2indices(Mm.c2,id=rownames(v))
cam.BasalvsLP <- camera(v,idx,design,contrast=contr.matrix[,1])
head(cam.BasalvsLP,5)
cam.BasalvsML <- camera(v,idx,design,contrast=contr.matrix[,2])
head(cam.BasalvsML,5)
cam.LPvsML <- camera(v,idx,design,contrast=contr.matrix[,3])
head(cam.LPvsML,5)

## ----barcodeplot, fig.height=6, fig.width=6, fig.cap="`LIM_MAMMARY_LUMINAL_MATURE_UP` ï¼ˆçº¢è‰²æ¡å½¢ï¼Œå›¾è¡¨ä¸Šæ–¹ï¼‰å’Œ`LIM_MAMMARY_LUMINAL_MATURE_DN`ï¼ˆè“è‰²æ¡å½¢ï¼Œå›¾è¡¨ä¸‹æ–¹ï¼‰åŸºå› é›†åœ¨LPå’ŒMLçš„å¯¹æ¯”ä¸­çš„æ¡ç å›¾ï¼Œæ¯ä¸ªåŸºå› é›†éƒ½æœ‰ä¸€æ¡å¯Œé›†çº¿å±•ç¤ºäº†ç«–ç›´æ¡å½¢åœ¨å›¾è¡¨æ¯éƒ¨åˆ†çš„ç›¸å¯¹å¯Œé›†ç¨‹åº¦ã€‚Limç­‰äººçš„å®žéªŒ[@Lim:BreastCancerRes:2010]éžå¸¸ç±»ä¼¼äºŽæˆ‘ä»¬çš„ï¼Œç”¨äº†ç›¸åŒçš„åˆ†é€‰æ–¹å¼æ¥èŽ·å–ä¸åŒçš„ç»†èƒžç¾¤ï¼Œåªæ˜¯ä»–ä»¬ä½¿ç”¨çš„æ˜¯å¾®é˜µåˆ—è€Œä¸æ˜¯RNA-seqæ¥æµ‹å®šåŸºå› è¡¨è¾¾ã€‚éœ€è¦æ³¨æ„çš„æ˜¯ï¼Œä¸Šè°ƒåŸºå› é›†å‘ç”Ÿä¸‹è°ƒè€Œä¸‹è°ƒåŸºå› é›†å‘ç”Ÿä¸Šè°ƒçš„é€†ç›¸å…³æ€§æ¥è‡ªäºŽå¯¹æ¯”çš„è®¾å®šæ–¹å¼ï¼ˆLPç›¸æ¯”äºŽMLï¼‰ï¼Œå¦‚æžœå°†å…¶å¯¹è°ƒï¼Œæ–¹å‘æ€§å°†ä¼šå»åˆã€‚"----
barcodeplot(efit$t[,3], index=idx$LIM_MAMMARY_LUMINAL_MATURE_UP, 
            index2=idx$LIM_MAMMARY_LUMINAL_MATURE_DN, main="LPvsML")

## ----softwareinfo-----------------------------------------------------------------------
sessionInfo()