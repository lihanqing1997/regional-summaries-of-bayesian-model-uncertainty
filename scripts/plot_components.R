# Reproduce Figure 3 from the accompanying CSV; no simulations are run.
args <- commandArgs(trailingOnly=FALSE)
script <- sub('^--file=', '', args[grep('^--file=', args)])
work <- dirname(normalizePath(script, winslash='/', mustWork=TRUE))
d <- read.csv(file.path(work,'../results/bma_components/figure3_data.csv'))
stopifnot(nrow(d)==200, length(unique(d$job))==40)
levels <- c('h_adaptive','j_equal','j_adaptive','budget_50000','budget_100000')
cols <- c('#D7E1E6','#D7E1E6','#388799','#89B9C3','#C0DAD8')
pdf(file.path(Sys.getenv('REBUILD_OUTPUT'),'figures','bma_confirmatory.pdf'),width=7,height=3.55,
    family='Helvetica',title='Interval refinement, allocation, and computation',
    useDingbats=FALSE)
par(mfrow=c(1,2),mar=c(3.5,4,2.6,.5),mgp=c(2.6,.7,0),tcl=-.25,
    cex=.85,bty='l',las=1)
set.seed(61001)
for(i in 1:2) {
  case <- c('stress_quartets_p20_b005','stress_quartets_p20_b020')[i]
  g <- d[d$case_id==case,]
  g$plot_condition <- factor(g$plot_condition,levels=levels)
  stopifnot(all(table(g$plot_condition)==20))
  boxplot(ratio~plot_condition,g,outline=FALSE,col=cols,border='#48616D',
          ylab=if(i==1)'Complexity / optimum' else '', xlab='', xaxt='n',
          ylim=c(.7,max(d$ratio)*1.06),cex.axis=.9,
          main=paste0('(',letters[i],') Between-group correlation ',c('.05','.20')[i]),
          cex.main=1.0,boxwex=.55)
  abline(h=axTicks(2),col='#E1E5E8',lwd=.5)
  boxplot(ratio~plot_condition,g,outline=FALSE,col=cols,border='#48616D',
          add=TRUE,axes=FALSE,boxwex=.55)
  stripchart(ratio~plot_condition,g,method='jitter',jitter=.14,vertical=TRUE,
             pch=16,cex=.55,col=adjustcolor('#203E4B',alpha.f=.55),add=TRUE)
  abline(h=1,lty=2,col='#A34F31',lwd=1)
  axis(1,at=1:5,labels=FALSE,tick=FALSE)
  mtext(expression(H[N],tilde(H)[N],tilde(H)[N],'50,000','100,000'),side=1,at=1:5,line=.6,cex=.76)
  mtext(c('adaptive','equal','adaptive','120 s','180 s'),side=1,at=1:5,line=1.6,cex=.76)
}
invisible(dev.off())
cat('Rendered two panels from 200 recorded outcomes on 40 paired datasets.\n')
