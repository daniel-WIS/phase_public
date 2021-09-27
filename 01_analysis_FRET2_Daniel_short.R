source("/media/elmicro/danielo/scripts_phase/MicroscopeToolBoxFRET_Daniel.R")
library(nnls)

get.all.cells = function(data){

    RES = c()
    WELLS = names(data[[1]])
    
    for (each.well in WELLS){

        if(length(RES)==0){
            RES = data[[1]][[each.well]]
            RES$well = each.well
        } else {
            RES.tmp = data[[1]][[each.well]]

            if(nrow(RES.tmp)>0){
                RES.tmp$well = each.well
                RES = rbind(RES, RES.tmp)
            }
        }
    }
    return(RES)
}
get.cell.from.wells = function(data, WELLS){

    RES = c()
    
    for (each.well in WELLS){

        if(length(RES)==0){
            RES = data[[1]][[each.well]]
        } else {
            RES = rbind(RES, data[[1]][[each.well]])
        }
    }
    return(RES)
}

##source("https://raw.githubusercontent.com/elevylab/RToolBox2/main/MicroscopeToolBox.R")

#args = commandArgs(trailingOnly = TRUE)

### Script has to be executed within a directory
current.dir = "/home/d/210824_FRET2_pairs_GFP/"

### Returns the same path withou the last directory (the plate) 
### screen.dir = gsub(pattern="^(.+/)[^/]+$", x=current.dir, replacement="\\1", perl=TRUE)
### Return only the last directory (the plate)
### subscreen.dir = gsub(pattern="^.*/([^/]+)$", x=current.dir, replacement="\\1", perl=TRUE)

output.dir = paste0(current.dir, "results/")

design.FRET = microscope.get.design(
    F=paste0(current.dir), 
    D=c("comp"),
    #pdef needs to have the same layout as: /media/elmicro/danielo/210611_FRET/210611_FRET_GFP_2_pairs/pdef_FRET_revision_pairs.csv
    PDEF=paste0(current.dir,"pdef_FRET2_pairs_no_D51_Y54_Y55A.csv"),
    FORMAT=384,
    OUT=paste0(current.dir,"images_output"),
    CHANELS=c("BFP", "GFP", "GFP2"),
    MEASURE.COL = "_int_b5",
    KEY=c(),
    DIR.res = output.dir 
)

data.FRET = microscope.load.data(design.FRET) 
uscope.count.cells(data.FRET)

design.FRET    = uscope.process.estimate.background(data.FRET, design.FRET)
data.FRET     = uscope.process.reorder(data.FRET, design=design.FRET)
data.FRET     = uscope.process.remove.first.pic(data.FRET, N=1)
#data.FRET     = uscope.process.rename.undefined(data.FRET)
data.FRET2     = uscope.process.remove.background.per.cell(data.FRET, design.FRET)
data.FRET2     = uscope.process.add.meanArea(data.FRET2)
data.FRET2     = uscope.process.remove.small(data.FRET2, MIN.size=800,MAX.size=2500)
data.FRET2     = uscope.process.add.ncells(data = data.FRET2)
#data.FRET2     = uscope.process.centroid(data.FRET2, X.center=504, Y.center=504)
data.FRET2     = uscope.process.centroid(data.FRET2, X.center=1304, Y.center=1324)
data.FRET2     = uscope.process.remove.centroid.outliers(data.FRET2, min=0,max=250)
uscope.count.cells(data.FRET2)

##### First I need to calibrate the crosstalk of mTurquoise 
par(mfcol=c(1,1))
LOG="xy"
XMAX = 300
YMAX = 300
NCELLS = 1000
cells = get.all.cells(data.FRET2)
cells = cells[cells$GFP2_int_bg < 115,]

cells.m = merge(cells, design.FRET$PDEF[,c("well","BFP","GFP","affinity")], by="well")
cells.m$affinity = -log10(as.numeric(cells.m$affinity))
my.col = rep("grey",nrow(cells.m))
my.col[ which(cells.m$affinity <= 5) ] = "red"
#my.col[ which(cells.m$affinity > 5 & cells.m$affinity < 7) ] = "orange"
#my.col[ which(cells.m$affinity >= 7 & cells.m$affinity < 8) ] = "green"
#my.col[ which(cells.m$affinity > 8 & cells.m$affinity < 10)] = "purple"
my.col[ which(cells.m$affinity >= 10)] = "blue"
#nTQ = mTQ[1:NCELLS,]
X = cells.m$BFP_int_b5+1  ## RAW PROT
Y = cells.m$GFP_int_b5+1 ## Basal "FRET" signal expected when there is no FRET
Z = cells.m$GFP2_int_b5+1

plot(X, Y, main="BFP-BFP polluted", log=LOG, ylab="GFP (FRET)", xlab="BFP", xlim=c(1,3000),ylim=c(1,10000), col = my.col)
plot(X, Z, main="BFP-GFP2", log=LOG, ylab="BFP-crosstalk", xlab="BFP", xlim=c(1,3000),ylim=c(1,60000), col = my.col)
plot(Y, Z, main="GFP2", log=LOG, ylab="GFP2-crosstalk", xlab="FRET", xlim=c(1,1000),ylim=c(1,60000), col = my.col)

plot(X[Z<50], Y[Z<50], main="BFP-BFP not-polluted", log=LOG, ylab="GFP (FRET)", xlab="BFP", xlim=c(1,3000),ylim=c(1,10000), col = my.col)

fit.fret.no.acceptor = lm(Y[Z<50]~X[Z<50])

SEL = X<15
myX = Z[SEL]
myY = Y[SEL]
fit.fret.no.bfp = lm(myY ~ myX + I(myX^0.5))
fit.fret.no.bfp.lin = lm(myY ~ myX)
myx = c(10,50,100,500,1000,5000,10000,60000)
myy = (myx^0.5*fit.fret.no.bfp$coefficients[3])+ (myx*fit.fret.no.bfp$coefficients[2]) + fit.fret.no.bfp$coefficients[1]
myy.lin = (myx*fit.fret.no.bfp.lin$coefficients[2]) + fit.fret.no.bfp.lin$coefficients[1]
plot(myX, myY, main="BFP-BFP not-polluted", log=LOG, ylab="GFP2 @ 488", xlab="GFP2 @ 405", col=my.col[SEL])
lines(myx, myy,col="orange", lwd=2)
lines(myx, myy.lin,col="red", lwd=2)


grid()
fit1 = lm( Y[Z<100] ~ X[Z<100])
abline(fit1)

#### DEFINES PHENO, which give WELLS given some affinity
E9.var1 = "E9-mNG"
E9.var2 = "E9-10aaL-mNG"
PHENO = list()
pheno.sel = design.FRET$PDEF$GFP == E9.var1 # | design.FRET$PDEF$GFP == E9.var2

all.pheno = names(table(design.FRET$PDEF$BFP[pheno.sel]))
all.aff = rep(0, length(all.pheno))
for(i in 1:length(all.pheno)){
    PHENO[[i]] = design.FRET$PDEF$well[ which(design.FRET$PDEF$BFP==all.pheno[i] & pheno.sel)]
    MUT = design.FRET$PDEF$BFP[ which(design.FRET$PDEF$BFP==all.pheno[i] & pheno.sel)]
    AFF = design.FRET$PDEF$affinity[which(design.FRET$PDEF$BFP==all.pheno[i] & pheno.sel)]
    all.aff[i] = paste(AFF[1],"---",MUT[1])
}
names(PHENO) = all.aff

AFFS = matrix(ncol=3, nrow=length(PHENO))
rownames(AFFS) = names(PHENO)


pdf(paste0(current.dir,"/FRET_ALL_plot_GFP2_bigScreen.pdf"), width=12, height=6)
for(i in 1:10){
    SELECT = Z > 250 & Z < 2000 & X < 100 & ((cells.m$BFP_int_b0/cells.m$BFP_int_b5) < 3) & cells.m$well %in% PHENO[[i]]
    OBS.FRET.1 = cells.m$GFP_int_b5[SELECT]-(10 + fit.fret.no.bfp[[1]][2]*cells.m$GFP2_int_b5[SELECT] + fit.fret.no.bfp[[1]][3]*(cells.m$GFP2_int_b5[SELECT]^0.5))
    EXP.FRET.1 = fit.fret.no.acceptor[[1]][1]+ (fit.fret.no.acceptor[[1]][2]*cells.m$BFP_int_b5[SELECT])
    my.col.sel = "orange"#my.col[SELECT]
    plot(100+EXP.FRET.1, 100+OBS.FRET.1, log=LOG, xlim=c(5,5000), ylim=c(5,5000),xlab="Expected FRET Baseline", ylab="Observed FRET Signal", cex=0.5, col=my.col.sel, main=names(PHENO)[i])
    abline(a=0, b=1, col="dark red")


    SELECT = Z > 2500 & X < 100 & ((cells.m$BFP_int_b0/cells.m$BFP_int_b5) < 3) & cells.m$well %in% PHENO[[i]]
    OBS.FRET.1 = cells.m$GFP_int_b5[SELECT]-(fit.fret.no.bfp.lin[[1]][2]*cells.m$GFP2_int_b5[SELECT])
    EXP.FRET.1 = fit.fret.no.acceptor[[1]][1]+ (fit.fret.no.acceptor[[1]][2]*cells.m$BFP_int_b5[SELECT])
    my.col.sel = "green"#my.col[SELECT]
    points(100+EXP.FRET.1, 100+OBS.FRET.1, cex=0.5, col=my.col.sel)

    grid()
}
dev.off()



####
#### Now we check if the FRET signal is dependent on GFP2 expression
#### 
#x11()
par(mfcol=c(2,3))

pdf(paste0(current.dir,"/FRET_ALL_GFP2_excessBFP.pdf"), width=12, height=6)
for (i in 1:length(PHENO)){
    FRET.1 = get.cell.from.wells(data.FRET2, PHENO[[i]])
    #FRET.1 = FRET.1[FRET.1$GFP2_int_b5 < 30000,]
                                        #FRET.1 = FRET.1[FRET.1$BFP_int_b5 < 60,]

    FRET.1 = FRET.1[FRET.1$GFP2_int_bg < 115,]
    FRET.1 = FRET.1[ (FRET.1$BFP_int_b0/FRET.1$BFP_int_b5) < 2 & FRET.1$BFP_int_b5 > 60 & FRET.1$BFP_int_b5 < 100,]
    print(names(PHENO)[i])

    OBS.FRET.1 = FRET.1$GFP_int_b3-(fit.fret.no.bfp.lin[[1]][2]*FRET.1$GFP2_int_b3)
    EXP.FRET.1 = fit.fret.no.acceptor[[1]][1]+ (fit.fret.no.acceptor[[1]][2]*FRET.1$BFP_int_b3)
    ratio = OBS.FRET.1/EXP.FRET.1
    reg.lo = loess(ratio ~ GFP, data = data.frame(ratio = ratio, GFP=FRET.1$GFP2_int_b5))
    GFP.exp = seq(1,5000,by=10)
    ratio.exp = predict(reg.lo, newdata = data.frame(GFP=GFP.exp))
                                        # AFFS[i,]=sapply(SP,median)    
    plot(FRET.1$GFP2_int_b5,ratio , log="x", main= paste0(names(PHENO)[i], " ncells=",length(ratio)), varwidth=TRUE, ylim=c(0.2,2),xlim=c(10,20000), xlab="E9 Abundance", ylab="FRET Signal")
    lines(GFP.exp,ratio.exp, lwd=2, col="orange")
    
    my.Y = ratio
    my.X = log2(FRET.1$GFP2_int_b5)
    #fitmodel = nls(my.Y ~ d + (a / (1 + exp(-b * (my.X-c)))), start=list(a=1 , b=1, c = 8, d=1), lower=list(a=0.1, b=0.1,c=4, d=0.7), upper=list(a=2,b=3,c=15, d=1.1), algorithm="port")
    #fitmodel = nls(my.Y ~ 1 + (a / (1 + exp(-b * (my.X-c)))), start=list(a=1 , b=1, c = 8), lower=list(a=0.1, b=0.1,c=4), upper=list(a=2,b=3,c=15), algorithm="port")
    fitmodel = nls(my.Y ~ 1 + (0.7 / (1 + exp(-b * (my.X-c)))), start=list(b=1, c = 8), lower=list(b=0.1,c=4), upper=list(b=3,c=15), algorithm="port")

    my.pars = fitmodel$m$getPars()

    AFFS[i,] = c(my.pars, -log10(as.numeric(design.FRET$PDEF$affinity[which(design.FRET$PDEF$BFP==all.pheno[i] & pheno.sel)][1])))
    
    ratio.exp.2 = predict(fitmodel, newdata = data.frame(my.X = log2(seq(1,5000,by=10))))
    #plot(my.X, my.Y, col="green", lwd=2)
    lines( seq(1,5000,by=10), ratio.exp.2, col="green", lwd=2)
                                        #boxplot(split(ratio,FRET.1$AB), log="y", main= names(PHENO)[i], varwidth=TRUE, ylim=c(0.2,10), xlab="E9 Abundance", ylab="FRET Signal")
    abline(h=1, col="dark red")
    grid()
}
#colnames(AFFS) = names(SP)
plot(AFFS[,2],AFFS[,3], cex=2)
text(AFFS[,2],AFFS[,3], lab=names(PHENO))
dev.off()
 
plot(EXP.FRET.1, OBS.FRET.1, log="xy", main=paste0("mTq2-10aaL-Im2wt + E9-mNG (", nrow(FRET.1)," cells)"))#, cex = log10(FRET.1$GFP2_int_b5))
abline(a=0, b=1, col="dark red")









### So, given some BFP measured we expect (-78+2.4*GFP2) GFP2 signal when there's no FRET
##### Second I need to calibrate the crosstalk of mNG
mNG = get.cell.from.wells(data.FRET2, WELLS=c("A04","A05","A06"))
X = mNG$GFP2_int_b5  ## RAW PROT
Y = mNG$GFP_int_b5 ## Basal "FRET" signal expected when there is no FRET
#plot(Y, X, main="GFP2 Alone", log=LOG, ylab="BASAL FRET SIGNAL EXPECTED", xlab="GFP2")#, xlim=c(0,XMAX),ylim=c(0,YMAX), col = "grey")
#grid()
#fit1 = lm( X ~ Y)
#abline(fit1)
### If we keep only cells where signal is <30K, then no correction is needed.
### So, given some BFP measured we expect (-78+2.4*GFP2) GFP2 signal when there's no FRET

######
###### Now the FRET
######
LOG=""
FRET.1 = get.cell.from.wells(data.FRET2, c("A07","A08","A09","A10","A11","A12"))
FRET.1 = FRET.1[FRET.1$GFP2_int_b5 < 20000,]
print(dim(FRET.1))
FRET.1 = FRET.1[1:NCELLS,]
OBS.FRET.1 = FRET.1$GFP_int_b5
EXP.FRET.1 = -78+ 2.4*FRET.1$BFP_int_b5

#ratio = OBS.FRET.1/EXP.FRET.1
plot(EXP.FRET.1, OBS.FRET.1, log=LOG, xlim=c(1,XMAX), ylim=c(1,YMAX), main=paste0("FRET with mTq2-Im2wt + E9-mNG (", nrow(FRET.1)," cells)"))
abline(a=0, b=1, col="dark red")

FRET.1 = get.cell.from.wells(data.FRET2, c("A13","A14","A15","A16","A17","A18"))
FRET.1 = FRET.1[FRET.1$GFP2_int_b5 < 20000,]
print(dim(FRET.1))
FRET.1 = FRET.1[1:NCELLS,]
OBS.FRET.1 = FRET.1$GFP_int_b5
EXP.FRET.1 = -78+ 2.4*FRET.1$BFP_int_b5
#ratio = OBS.FRET.1/EXP.FRET.1
plot(EXP.FRET.1, OBS.FRET.1, log=LOG, xlim=c(1,XMAX), ylim=c(1,YMAX), main=paste0("mTq2-10aaL-Im2wt + E9-mNG (", nrow(FRET.1)," cells)"))#, cex = log10(FRET.1$GFP2_int_b5))
abline(a=0, b=1, col="dark red")

FRET.1 = get.cell.from.wells(data.FRET2, c("A19","A20","A21","A22","A23","A24"))
FRET.1 = FRET.1[FRET.1$GFP2_int_b5 < 20000,]
print(dim(FRET.1))
FRET.1 = FRET.1[1:NCELLS,]
OBS.FRET.1 = FRET.1$GFP_int_b5
EXP.FRET.1 = -78+ 2.4*FRET.1$BFP_int_b5
#ratio = OBS.FRET.1/EXP.FRET.1
plot(EXP.FRET.1, OBS.FRET.1, log=LOG, xlim=c(1,XMAX), ylim=c(1,YMAX), main=paste0("sf-mTq2-OX-Im2 + E9-mNG (", nrow(FRET.1)," cells)") )
abline(a=0, b=1, col="dark red")

FRET.1 = get.cell.from.wells(data.FRET2, c("B01","B02","B03","B04","B05","B06"))
FRET.1 = FRET.1[FRET.1$GFP2_int_b5 < 20000,]
print(dim(FRET.1))
FRET.1 = FRET.1[1:NCELLS,]
OBS.FRET.1 = FRET.1$GFP_int_b5
EXP.FRET.1 = -78+ 2.4*FRET.1$BFP_int_b5
#ratio = OBS.FRET.1/EXP.FRET.1
plot(EXP.FRET.1, OBS.FRET.1, log=LOG, xlim=c(1,XMAX), ylim=c(1,YMAX), main=paste0("mTq2-Im2wt E9-10aaL-mNG (", nrow(FRET.1)," cells)") )
abline(a=0, b=1, col="dark red")

FRET.1 = get.cell.from.wells(data.FRET2, c("B07","B08","B09","B10","B11","B12"))
FRET.1 = FRET.1[FRET.1$GFP2_int_b5 < 20000,]
print(dim(FRET.1))
FRET.1 = FRET.1[1:NCELLS,]
OBS.FRET.1 = FRET.1$GFP_int_b5
EXP.FRET.1 = -78+ 2.4*FRET.1$BFP_int_b5
#ratio = OBS.FRET.1/EXP.FRET.1
plot(EXP.FRET.1, OBS.FRET.1, log=LOG, xlim=c(1,XMAX), ylim=c(1,YMAX), main=paste0("mTq2-10aaL-Im2wt E9-10aaL-mNG (", nrow(FRET.1)," cells)") )
abline(a=0, b=1, col="dark red")

FRET.1 = get.cell.from.wells(data.FRET2, c("B13","B14","B15","B16","B17","B18"))
FRET.1 = FRET.1[FRET.1$GFP2_int_b5 < 20000,]
print(dim(FRET.1))
FRET.1 = FRET.1[1:NCELLS,]
OBS.FRET.1 = FRET.1$GFP_int_b5
EXP.FRET.1 = -78+ 2.4*FRET.1$BFP_int_b5
#ratio = OBS.FRET.1/EXP.FRET.1
plot(EXP.FRET.1, OBS.FRET.1, log=LOG, xlim=c(1,XMAX), ylim=c(1,YMAX), main=paste0("sf-mTq2-OX-Im2wt E9-10aaL-mNG (", nrow(FRET.1)," cells)") )
abline(a=0, b=1, col="dark red")
