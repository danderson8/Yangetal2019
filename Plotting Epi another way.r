
All_data <- read.csv("AllData.csv")	#reads in the data file - name inside the quotation marks needs to match the filename. First column must be time.

attach(All_data)							# Assigns the first line - i.e. the head of all the columns - to be the variable name for all the vlues in those columns

#dhc[which(dhc <= 0.1)] = 0.1

pxm[which(pxm <= 0.1)] = 0.1

pte[which(ptm <= 0.1)] = 0.1

pxe[which(pxe <= 0.1)] = 0.1

pte[which(pte <= 0.1)] = 0.1


# DHC_conv <- log(dhc/(0.58*1300)*1000000/60, base=10)

PXM_conv <- log(pxm, base=10) #/(0.58*18300)*1000000/60, base=10)

PTM_conv <- log(ptm, base=10) #/(0.58*18300)*1000000/60, base=10)

PXE_conv <- log(pxe, base=10) #/(0.58*18300)*1000000/60, base=10)

PTE_conv <- log(pte, base=10) #/(0.58*18300)*1000000/60, base=10)

n <- length(PTM_conv)



# Now calculating the effect of each mutation in all possible backgrounds

y <- array(NA, dim=c(n/3))

z <- array(NA, dim=c(n/3))

yy <- array(NA, dim=c(n/3))

zz <- array(NA, dim=c(n/3))

for(j in 1:(n/3)) {
    #    x[j] <- (DHC_conv[(j*3)]+DHC_conv[(j*3-1)]+DHC_conv[(j*3-2)])/3
    y[j] <- (PXM_conv[(j*3)]+PXM_conv[(j*3-1)]+PXM_conv[(j*3-2)])/3
    z[j] <- (PTM_conv[(j*3)]+PTM_conv[(j*3-1)]+PTM_conv[(j*3-2)])/3
    yy[j] <- (PXE_conv[(j*3)]+PXE_conv[(j*3-1)]+PXE_conv[(j*3-2)])/3
    zz[j] <- (PTE_conv[(j*3)]+PTE_conv[(j*3-1)]+PTE_conv[(j*3-2)])/3
}


traj_pxm <- array(NA, dim=c(3,3,3,3,3))

traj_ptm <- array(NA, dim=c(3,3,3,3,3))

traj_pxe <- array(NA, dim=c(3,3,3,3,3))

traj_pte <- array(NA, dim=c(3,3,3,3,3))

for (i in 1:32) {
    # traj_dhc[pos72[i*3]+2, pos193[i*3]+2, pos258[i*3]+2, pos271[i*3]+2, pos273[i*3]+2] = x[i]
        traj_pxm[pos72[i*3]+2, pos193[i*3]+2, pos258[i*3]+2, pos271[i*3]+2, pos273[i*3]+2] = y[i]
        traj_ptm[pos72[i*3]+2, pos193[i*3]+2, pos258[i*3]+2, pos271[i*3]+2, pos273[i*3]+2] = z[i]
        traj_pxe[pos72[i*3]+2, pos193[i*3]+2, pos258[i*3]+2, pos271[i*3]+2, pos273[i*3]+2] = yy[i]
        traj_pte[pos72[i*3]+2, pos193[i*3]+2, pos258[i*3]+2, pos271[i*3]+2, pos273[i*3]+2] = zz[i]
}

#p72_dhc <- array(NA, dim=c(16))

p72_pxm <- array(NA, dim=c(16))

p72_ptm <- array(NA, dim=c(16))

p72_pxe <- array(NA, dim=c(16))

p72_pte <- array(NA, dim=c(16))

# Computing effects


p72_pxm[1] <- traj_pxm[3,1,1,1,1] - traj_pxm[1,1,1,1,1]
p72_pxm[2] <- traj_pxm[3,3,1,1,1] - traj_pxm[1,3,1,1,1]
p72_pxm[3] <- traj_pxm[3,1,3,1,1] - traj_pxm[1,1,3,1,1]
p72_pxm[4] <- traj_pxm[3,1,1,3,1] - traj_pxm[1,1,1,3,1]
p72_pxm[5] <- traj_pxm[3,1,1,1,3] - traj_pxm[1,1,1,1,3]
p72_pxm[6] <- traj_pxm[3,3,3,1,1] - traj_pxm[1,3,3,1,1]
p72_pxm[7] <- traj_pxm[3,3,1,3,1] - traj_pxm[1,3,1,3,1]
p72_pxm[8] <- traj_pxm[3,3,1,1,3] - traj_pxm[1,3,1,1,3]
p72_pxm[9] <- traj_pxm[3,1,3,3,1] - traj_pxm[1,1,3,3,1]
p72_pxm[10] <- traj_pxm[3,1,3,1,3] - traj_pxm[1,1,3,1,3]
p72_pxm[11] <- traj_pxm[3,1,1,3,3] - traj_pxm[1,1,1,3,3]
p72_pxm[12] <- traj_pxm[3,3,3,3,1] - traj_pxm[1,3,3,3,1]
p72_pxm[13] <- traj_pxm[3,3,3,1,3] - traj_pxm[1,3,3,1,3]
p72_pxm[14] <- traj_pxm[3,3,1,3,3] - traj_pxm[1,3,1,3,3]
p72_pxm[15] <- traj_pxm[3,1,3,3,3] - traj_pxm[1,1,3,3,3]
p72_pxm[16] <- traj_pxm[3,3,3,3,3] - traj_pxm[1,3,3,3,3]

p72_ptm[1] <- traj_ptm[3,1,1,1,1] - traj_ptm[1,1,1,1,1]
p72_ptm[2] <- traj_ptm[3,3,1,1,1] - traj_ptm[1,3,1,1,1]
p72_ptm[3] <- traj_ptm[3,1,3,1,1] - traj_ptm[1,1,3,1,1]
p72_ptm[4] <- traj_ptm[3,1,1,3,1] - traj_ptm[1,1,1,3,1]
p72_ptm[5] <- traj_ptm[3,1,1,1,3] - traj_ptm[1,1,1,1,3]
p72_ptm[6] <- traj_ptm[3,3,3,1,1] - traj_ptm[1,3,3,1,1]
p72_ptm[7] <- traj_ptm[3,3,1,3,1] - traj_ptm[1,3,1,3,1]
p72_ptm[8] <- traj_ptm[3,3,1,1,3] - traj_ptm[1,3,1,1,3]
p72_ptm[9] <- traj_ptm[3,1,3,3,1] - traj_ptm[1,1,3,3,1]
p72_ptm[10] <- traj_ptm[3,1,3,1,3] - traj_ptm[1,1,3,1,3]
p72_ptm[11] <- traj_ptm[3,1,1,3,3] - traj_ptm[1,1,1,3,3]
p72_ptm[12] <- traj_ptm[3,3,3,3,1] - traj_ptm[1,3,3,3,1]
p72_ptm[13] <- traj_ptm[3,3,3,1,3] - traj_ptm[1,3,3,1,3]
p72_ptm[14] <- traj_ptm[3,3,1,3,3] - traj_ptm[1,3,1,3,3]
p72_ptm[15] <- traj_ptm[3,1,3,3,3] - traj_ptm[1,1,3,3,3]
p72_ptm[16] <- traj_ptm[3,3,3,3,3] - traj_ptm[1,3,3,3,3]

p72_pxe[1] <- traj_pxe[3,1,1,1,1] - traj_pxe[1,1,1,1,1]
p72_pxe[2] <- traj_pxe[3,3,1,1,1] - traj_pxe[1,3,1,1,1]
p72_pxe[3] <- traj_pxe[3,1,3,1,1] - traj_pxe[1,1,3,1,1]
p72_pxe[4] <- traj_pxe[3,1,1,3,1] - traj_pxe[1,1,1,3,1]
p72_pxe[5] <- traj_pxe[3,1,1,1,3] - traj_pxe[1,1,1,1,3]
p72_pxe[6] <- traj_pxe[3,3,3,1,1] - traj_pxe[1,3,3,1,1]
p72_pxe[7] <- traj_pxe[3,3,1,3,1] - traj_pxe[1,3,1,3,1]
p72_pxe[8] <- traj_pxe[3,3,1,1,3] - traj_pxe[1,3,1,1,3]
p72_pxe[9] <- traj_pxe[3,1,3,3,1] - traj_pxe[1,1,3,3,1]
p72_pxe[10] <- traj_pxe[3,1,3,1,3] - traj_pxe[1,1,3,1,3]
p72_pxe[11] <- traj_pxe[3,1,1,3,3] - traj_pxe[1,1,1,3,3]
p72_pxe[12] <- traj_pxe[3,3,3,3,1] - traj_pxe[1,3,3,3,1]
p72_pxe[13] <- traj_pxe[3,3,3,1,3] - traj_pxe[1,3,3,1,3]
p72_pxe[14] <- traj_pxe[3,3,1,3,3] - traj_pxe[1,3,1,3,3]
p72_pxe[15] <- traj_pxe[3,1,3,3,3] - traj_pxe[1,1,3,3,3]
p72_pxe[16] <- traj_pxe[3,3,3,3,3] - traj_pxe[1,3,3,3,3]

p72_pte[1] <- traj_pte[3,1,1,1,1] - traj_pte[1,1,1,1,1]
p72_pte[2] <- traj_pte[3,3,1,1,1] - traj_pte[1,3,1,1,1]
p72_pte[3] <- traj_pte[3,1,3,1,1] - traj_pte[1,1,3,1,1]
p72_pte[4] <- traj_pte[3,1,1,3,1] - traj_pte[1,1,1,3,1]
p72_pte[5] <- traj_pte[3,1,1,1,3] - traj_pte[1,1,1,1,3]
p72_pte[6] <- traj_pte[3,3,3,1,1] - traj_pte[1,3,3,1,1]
p72_pte[7] <- traj_pte[3,3,1,3,1] - traj_pte[1,3,1,3,1]
p72_pte[8] <- traj_pte[3,3,1,1,3] - traj_pte[1,3,1,1,3]
p72_pte[9] <- traj_pte[3,1,3,3,1] - traj_pte[1,1,3,3,1]
p72_pte[10] <- traj_pte[3,1,3,1,3] - traj_pte[1,1,3,1,3]
p72_pte[11] <- traj_pte[3,1,1,3,3] - traj_pte[1,1,1,3,3]
p72_pte[12] <- traj_pte[3,3,3,3,1] - traj_pte[1,3,3,3,1]
p72_pte[13] <- traj_pte[3,3,3,1,3] - traj_pte[1,3,3,1,3]
p72_pte[14] <- traj_pte[3,3,1,3,3] - traj_pte[1,3,1,3,3]
p72_pte[15] <- traj_pte[3,1,3,3,3] - traj_pte[1,1,3,3,3]
p72_pte[16] <- traj_pte[3,3,3,3,3] - traj_pte[1,3,3,3,3]

write("Position 72 Effects", file="Pos72.csv")
write("genotype, dhc, px, pt, pom, pte", file="Pos72.csv", append=TRUE)
write(paste("xaaaa", p72_pxm[1],p72_ptm[1],p72_pxe[1],p72_pte[1], sep=","), file="Pos72.csv", append=TRUE)
write(paste("xdaaa", p72_pxm[2],p72_ptm[2],p72_pxe[2],p72_pte[2], sep=","), file="Pos72.csv", append=TRUE)
write(paste("xadaa", p72_pxm[3],p72_ptm[3],p72_pxe[3],p72_pte[3], sep=","), file="Pos72.csv", append=TRUE)
write(paste("xaada", p72_pxm[4],p72_ptm[4],p72_pxe[4],p72_pte[4], sep=","), file="Pos72.csv", append=TRUE)
write(paste("xaaad", p72_pxm[5],p72_ptm[5],p72_pxe[5],p72_pte[5], sep=","), file="Pos72.csv", append=TRUE)
write(paste("xddaa", p72_pxm[6],p72_ptm[6],p72_pxe[6],p72_pte[6], sep=","), file="Pos72.csv", append=TRUE)
write(paste("xdada", p72_pxm[7],p72_ptm[7],p72_pxe[7],p72_pte[7], sep=","), file="Pos72.csv", append=TRUE)
write(paste("xdaad", p72_pxm[8],p72_ptm[8],p72_pxe[8],p72_pte[8], sep=","), file="Pos72.csv", append=TRUE)
write(paste("xadda", p72_pxm[9],p72_ptm[9],p72_pxe[9],p72_pte[9], sep=","), file="Pos72.csv", append=TRUE)
write(paste("xadad", p72_pxm[10],p72_ptm[10],p72_pxe[10],p72_pte[10], sep=","), file="Pos72.csv", append=TRUE)
write(paste("xaadd", p72_pxm[11],p72_ptm[11],p72_pxe[11],p72_pte[11], sep=","), file="Pos72.csv", append=TRUE)
write(paste("xddda", p72_pxm[12],p72_ptm[12],p72_pxe[12],p72_pte[12], sep=","), file="Pos72.csv", append=TRUE)
write(paste("xddad", p72_pxm[13],p72_ptm[13],p72_pxe[13],p72_pte[13], sep=","), file="Pos72.csv", append=TRUE)
write(paste("xdadd", p72_pxm[14],p72_ptm[14],p72_pxe[14],p72_pte[14], sep=","), file="Pos72.csv", append=TRUE)
write(paste("xaddd", p72_pxm[15],p72_ptm[15],p72_pxe[15],p72_pte[15], sep=","), file="Pos72.csv", append=TRUE)
write(paste("xdddd", p72_pxm[16],p72_ptm[16],p72_pxe[16],p72_pte[16], sep=","), file="Pos72.csv", append=TRUE)


pdf(file=paste("72_pxm.pdf"))

barplot(p72_pxm)

dev.off()

pdf(file=paste("72_ptm.pdf"))

barplot(p72_ptm)

dev.off()

pdf(file=paste("72_pxe.pdf"))

barplot(p72_pxe)

dev.off()

pdf(file=paste("72_pte.pdf"))

barplot(p72_pte)

dev.off()



pdf(file=paste("72_pxmvspxe", "pdf", sep="."))

pxPXE_72 <- lm(p72_pxm~p72_pxe)
pxPXE_r2 = summary(pxPXE_72)$r.squared
pxPXE_m = coefficients(pxPXE_72)[[2]]

plot(p72_pxe, p72_pxm)

abline(pxPXE_72)

mtext(paste("m =", round(pxPXE_m, 3)), side=3, line=-1, at=0, cex=0.8, adj=0)
mtext(paste("r2 =", round(pxPXE_r2, 3)), side=3, line=-2, at=0, cex=0.8, adj=0)
pxPXE_p = summary(pxPXE_72)$coefficients[,4]
mtext(paste("p =", round(pxPXE_p, 3)), side=3, line=-3, at=0, cex=0.8, adj=0)

dev.off()


pdf(file=paste("72_pxmvsptm", "pdf", sep="."))

pxPTM_72 <- lm(p72_pxm~p72_ptm)
pxPTM_r2 = summary(pxPTM_72)$r.squared
pxPTM_m = coefficients(pxPTM_72)[[2]]

plot(p72_ptm, p72_pxm)

abline(pxPTM_72)

mtext(paste("m =", round(pxPTM_m, 3)), side=3, line=-1, at=0, cex=0.8, adj=0)
mtext(paste("r2 =", round(pxPTM_r2, 3)), side=3, line=-2, at=0, cex=0.8, adj=0)
pxPTM_p = summary(pxPTM_72)$coefficients[,4]
mtext(paste("p =", round(pxPTM_p, 3)), side=3, line=-3, at=0, cex=0.8, adj=0)

dev.off()


pdf(file=paste("72_pxmvspte", "pdf", sep="."))

pxpte_72 <- lm(p72_pxm~p72_pte)
pxpte_r2 = summary(pxpte_72)$r.squared
pxpte_m = coefficients(pxpte_72)[[2]]

plot(p72_pte, p72_pxm)

abline(pxpte_72)

mtext(paste("m =", round(pxpte_m, 3)), side=3, line=-1, at=0, cex=0.8, adj=0)
mtext(paste("r2 =", round(pxpte_r2, 3)), side=3, line=-2, at=0, cex=0.8, adj=0)
pxpte_p = summary(pxpte_72)$coefficients[,4]
mtext(paste("p =", round(pxpte_p, 3)), side=3, line=-3, at=0, cex=0.8, adj=0)

dev.off()


pdf(file=paste("72_ptmvspxe", "pdf", sep="."))

ptPXE_72 <- lm(p72_ptm~p72_pxe)
ptPXE_r2 = summary(ptPXE_72)$r.squared
ptPXE_m = coefficients(ptPXE_72)[[2]]

plot(p72_pxe, p72_ptm)

abline(ptPXE_72)

mtext(paste("m =", round(ptPXE_m, 3)), side=3, line=-1, at=0, cex=0.8, adj=0)
mtext(paste("r2 =", round(ptPXE_r2, 3)), side=3, line=-2, at=0, cex=0.8, adj=0)
ptPXE_p = summary(ptPXE_72)$coefficients[,4]
mtext(paste("p =", round(ptPXE_p, 3)), side=3, line=-3, at=0, cex=0.8, adj=0)

dev.off()


pdf(file=paste("72_ptmvspte", "pdf", sep="."))

ptpte_72 <- lm(p72_ptm~p72_pte)
ptpte_r2 = summary(ptpte_72)$r.squared
ptpte_m = coefficients(ptpte_72)[[2]]

plot(p72_pte, p72_ptm)

abline(ptpte_72)

mtext(paste("m =", round(ptpte_m, 3)), side=3, line=-1, at=0, cex=0.8, adj=0)
mtext(paste("r2 =", round(ptpte_r2, 3)), side=3, line=-2, at=0, cex=0.8, adj=0)
ptpte_p = summary(ptpte_72)$coefficients[,4]
mtext(paste("p =", round(ptpte_p, 3)), side=3, line=-3, at=0, cex=0.8, adj=0)

dev.off()

pdf(file=paste("72_pxevspte", "pdf", sep="."))

pompte_72 <- lm(p72_pxe~p72_pte)
pompte_r2 = summary(pompte_72)$r.squared
pompte_m = coefficients(pompte_72)[[2]]

plot(p72_pte, p72_pxe)

abline(pompte_72)

mtext(paste("m =", round(pompte_m, 3)), side=3, line=-1, at=0, cex=0.8, adj=0)
mtext(paste("r2 =", round(pompte_r2, 3)), side=3, line=-2, at=0, cex=0.8, adj=0)
pompte_p = summary(pompte_72)$coefficients[,4]
mtext(paste("p =", round(pompte_p, 3)), side=3, line=-3, at=0, cex=0.8, adj=0)

dev.off()


p193_pxm <- array(NA, dim=c(16))

p193_ptm <- array(NA, dim=c(16))

p193_pxe <- array(NA, dim=c(16))

p193_pte <- array(NA, dim=c(16))



# Computing effects


p193_pxm[1] <- traj_pxm[1,3,1,1,1] - traj_pxm[1,1,1,1,1]
p193_pxm[2] <- traj_pxm[3,3,1,1,1] - traj_pxm[3,1,1,1,1]
p193_pxm[3] <- traj_pxm[1,3,3,1,1] - traj_pxm[1,1,3,1,1]
p193_pxm[4] <- traj_pxm[1,3,1,3,1] - traj_pxm[1,1,1,3,1]
p193_pxm[5] <- traj_pxm[1,3,1,1,3] - traj_pxm[1,1,1,1,3]
p193_pxm[6] <- traj_pxm[3,3,3,1,1] - traj_pxm[3,1,3,1,1]
p193_pxm[7] <- traj_pxm[3,3,1,3,1] - traj_pxm[3,1,1,3,1]
p193_pxm[8] <- traj_pxm[3,3,1,1,3] - traj_pxm[3,1,1,1,3]
p193_pxm[9] <- traj_pxm[1,3,3,3,1] - traj_pxm[1,1,3,3,1]
p193_pxm[10] <- traj_pxm[1,3,3,1,3] - traj_pxm[1,1,3,1,3]
p193_pxm[11] <- traj_pxm[1,3,1,3,3] - traj_pxm[1,1,1,3,3]
p193_pxm[12] <- traj_pxm[3,3,3,3,1] - traj_pxm[3,1,3,3,1]
p193_pxm[13] <- traj_pxm[3,3,3,1,3] - traj_pxm[3,1,3,1,3]
p193_pxm[14] <- traj_pxm[3,3,1,3,3] - traj_pxm[3,1,1,3,3]
p193_pxm[15] <- traj_pxm[1,3,3,3,3] - traj_pxm[1,1,3,3,3]
p193_pxm[16] <- traj_pxm[3,3,3,3,3] - traj_pxm[3,1,3,3,3]

p193_ptm[1] <- traj_ptm[1,3,1,1,1] - traj_ptm[1,1,1,1,1]
p193_ptm[2] <- traj_ptm[3,3,1,1,1] - traj_ptm[3,1,1,1,1]
p193_ptm[3] <- traj_ptm[1,3,3,1,1] - traj_ptm[1,1,3,1,1]
p193_ptm[4] <- traj_ptm[1,3,1,3,1] - traj_ptm[1,1,1,3,1]
p193_ptm[5] <- traj_ptm[1,3,1,1,3] - traj_ptm[1,1,1,1,3]
p193_ptm[6] <- traj_ptm[3,3,3,1,1] - traj_ptm[3,1,3,1,1]
p193_ptm[7] <- traj_ptm[3,3,1,3,1] - traj_ptm[3,1,1,3,1]
p193_ptm[8] <- traj_ptm[3,3,1,1,3] - traj_ptm[3,1,1,1,3]
p193_ptm[9] <- traj_ptm[1,3,3,3,1] - traj_ptm[1,1,3,3,1]
p193_ptm[10] <- traj_ptm[1,3,3,1,3] - traj_ptm[1,1,3,1,3]
p193_ptm[11] <- traj_ptm[1,3,1,3,3] - traj_ptm[1,1,1,3,3]
p193_ptm[12] <- traj_ptm[3,3,3,3,1] - traj_ptm[3,1,3,3,1]
p193_ptm[13] <- traj_ptm[3,3,3,1,3] - traj_ptm[3,1,3,1,3]
p193_ptm[14] <- traj_ptm[3,3,1,3,3] - traj_ptm[3,1,1,3,3]
p193_ptm[15] <- traj_ptm[1,3,3,3,3] - traj_ptm[1,1,3,3,3]
p193_ptm[16] <- traj_ptm[3,3,3,3,3] - traj_ptm[3,1,3,3,3]

p193_pxe[1] <- traj_pxe[1,3,1,1,1] - traj_pxe[1,1,1,1,1]
p193_pxe[2] <- traj_pxe[3,3,1,1,1] - traj_pxe[3,1,1,1,1]
p193_pxe[3] <- traj_pxe[1,3,3,1,1] - traj_pxe[1,1,3,1,1]
p193_pxe[4] <- traj_pxe[1,3,1,3,1] - traj_pxe[1,1,1,3,1]
p193_pxe[5] <- traj_pxe[1,3,1,1,3] - traj_pxe[1,1,1,1,3]
p193_pxe[6] <- traj_pxe[3,3,3,1,1] - traj_pxe[3,1,3,1,1]
p193_pxe[7] <- traj_pxe[3,3,1,3,1] - traj_pxe[3,1,1,3,1]
p193_pxe[8] <- traj_pxe[3,3,1,1,3] - traj_pxe[3,1,1,1,3]
p193_pxe[9] <- traj_pxe[1,3,3,3,1] - traj_pxe[1,1,3,3,1]
p193_pxe[10] <- traj_pxe[1,3,3,1,3] - traj_pxe[1,1,3,1,3]
p193_pxe[11] <- traj_pxe[1,3,1,3,3] - traj_pxe[1,1,1,3,3]
p193_pxe[12] <- traj_pxe[3,3,3,3,1] - traj_pxe[3,1,3,3,1]
p193_pxe[13] <- traj_pxe[3,3,3,1,3] - traj_pxe[3,1,3,1,3]
p193_pxe[14] <- traj_pxe[3,3,1,3,3] - traj_pxe[3,1,1,3,3]
p193_pxe[15] <- traj_pxe[1,3,3,3,3] - traj_pxe[1,1,3,3,3]
p193_pxe[16] <- traj_pxe[3,3,3,3,3] - traj_pxe[3,1,3,3,3]

p193_pte[1] <- traj_pte[1,3,1,1,1] - traj_pte[1,1,1,1,1]
p193_pte[2] <- traj_pte[3,3,1,1,1] - traj_pte[3,1,1,1,1]
p193_pte[3] <- traj_pte[1,3,3,1,1] - traj_pte[1,1,3,1,1]
p193_pte[4] <- traj_pte[1,3,1,3,1] - traj_pte[1,1,1,3,1]
p193_pte[5] <- traj_pte[1,3,1,1,3] - traj_pte[1,1,1,1,3]
p193_pte[6] <- traj_pte[3,3,3,1,1] - traj_pte[3,1,3,1,1]
p193_pte[7] <- traj_pte[3,3,1,3,1] - traj_pte[3,1,1,3,1]
p193_pte[8] <- traj_pte[3,3,1,1,3] - traj_pte[3,1,1,1,3]
p193_pte[9] <- traj_pte[1,3,3,3,1] - traj_pte[1,1,3,3,1]
p193_pte[10] <- traj_pte[1,3,3,1,3] - traj_pte[1,1,3,1,3]
p193_pte[11] <- traj_pte[1,3,1,3,3] - traj_pte[1,1,1,3,3]
p193_pte[12] <- traj_pte[3,3,3,3,1] - traj_pte[3,1,3,3,1]
p193_pte[13] <- traj_pte[3,3,3,1,3] - traj_pte[3,1,3,1,3]
p193_pte[14] <- traj_pte[3,3,1,3,3] - traj_pte[3,1,1,3,3]
p193_pte[15] <- traj_pte[1,3,3,3,3] - traj_pte[1,1,3,3,3]
p193_pte[16] <- traj_pte[3,3,3,3,3] - traj_pte[3,1,3,3,3]

write("Position 193 Effects", file="Pos193.csv")
write("genotype, dhc, px, pt, pom, pte", file="Pos193.csv", append=TRUE)
write(paste("axaaa", p193_pxm[1],p193_ptm[1],p193_pxe[1],p193_pte[1], sep=","), file="Pos193.csv", append=TRUE)
write(paste("dxaaa", p193_pxm[2],p193_ptm[2],p193_pxe[2],p193_pte[2], sep=","), file="Pos193.csv", append=TRUE)
write(paste("axdaa", p193_pxm[3],p193_ptm[3],p193_pxe[3],p193_pte[3], sep=","), file="Pos193.csv", append=TRUE)
write(paste("axada", p193_pxm[4],p193_ptm[4],p193_pxe[4],p193_pte[4], sep=","), file="Pos193.csv", append=TRUE)
write(paste("axaad", p193_pxm[5],p193_ptm[5],p193_pxe[5],p193_pte[5], sep=","), file="Pos193.csv", append=TRUE)
write(paste("dxdaa", p193_pxm[6],p193_ptm[6],p193_pxe[6],p193_pte[6], sep=","), file="Pos193.csv", append=TRUE)
write(paste("dxada", p193_pxm[7],p193_ptm[7],p193_pxe[7],p193_pte[7], sep=","), file="Pos193.csv", append=TRUE)
write(paste("dxaad", p193_pxm[8],p193_ptm[8],p193_pxe[8],p193_pte[8], sep=","), file="Pos193.csv", append=TRUE)
write(paste("axdda", p193_pxm[9],p193_ptm[9],p193_pxe[9],p193_pte[9], sep=","), file="Pos193.csv", append=TRUE)
write(paste("axdad", p193_pxm[10],p193_ptm[10],p193_pxe[10],p193_pte[10], sep=","), file="Pos193.csv", append=TRUE)
write(paste("axadd", p193_pxm[11],p193_ptm[11],p193_pxe[11],p193_pte[11], sep=","), file="Pos193.csv", append=TRUE)
write(paste("dxdda", p193_pxm[12],p193_ptm[12],p193_pxe[12],p193_pte[12], sep=","), file="Pos193.csv", append=TRUE)
write(paste("dxdad", p193_pxm[13],p193_ptm[13],p193_pxe[13],p193_pte[13], sep=","), file="Pos193.csv", append=TRUE)
write(paste("dxadd", p193_pxm[14],p193_ptm[14],p193_pxe[14],p193_pte[14], sep=","), file="Pos193.csv", append=TRUE)
write(paste("axddd", p193_pxm[15],p193_ptm[15],p193_pxe[15],p193_pte[15], sep=","), file="Pos193.csv", append=TRUE)
write(paste("dxddd", p193_pxm[16],p193_ptm[16],p193_pxe[16],p193_pte[16], sep=","), file="Pos193.csv", append=TRUE)



pdf(file=paste("193_pxm.pdf"))

barplot(p193_pxm)

dev.off()

pdf(file=paste("193_ptm.pdf"))

barplot(p193_ptm)

dev.off()

pdf(file=paste("193_pxe.pdf"))

barplot(p193_pxe)

dev.off()

pdf(file=paste("193_pte.pdf"))

barplot(p193_pte)

dev.off()


pdf(file=paste("193_pxmvspxe", "pdf", sep="."))

pxPXE_193 <- lm(p193_pxm~p193_pxe)
pxPXE_r2 = summary(pxPXE_193)$r.squared
pxPXE_m = coefficients(pxPXE_193)[[2]]

plot(p193_pxe, p193_pxm)

abline(pxPXE_193)

mtext(paste("m =", round(pxPXE_m, 3)), side=3, line=-1, at=0, cex=0.8, adj=0)
mtext(paste("r2 =", round(pxPXE_r2, 3)), side=3, line=-2, at=0, cex=0.8, adj=0)
pxPXE_p = summary(pxPXE_193)$coefficients[,4]
mtext(paste("p =", round(pxPXE_p, 3)), side=3, line=-3, at=0, cex=0.8, adj=0)

dev.off()


pdf(file=paste("193_pxmvsptm", "pdf", sep="."))

pxPTM_193 <- lm(p193_pxm~p193_ptm)
pxPTM_r2 = summary(pxPTM_193)$r.squared
pxPTM_m = coefficients(pxPTM_193)[[2]]

plot(p193_ptm, p193_pxm)

abline(pxPTM_193)

mtext(paste("m =", round(pxPTM_m, 3)), side=3, line=-1, at=0, cex=0.8, adj=0)
mtext(paste("r2 =", round(pxPTM_r2, 3)), side=3, line=-2, at=0, cex=0.8, adj=0)
pxPTM_p = summary(pxPTM_193)$coefficients[,4]
mtext(paste("p =", round(pxPTM_p, 3)), side=3, line=-3, at=0, cex=0.8, adj=0)

dev.off()


pdf(file=paste("193_pxmvspte", "pdf", sep="."))

pxpte_193 <- lm(p193_px~p193_pte)
pxpte_r2 = summary(pxpte_193)$r.squared
pxpte_m = coefficients(pxpte_193)[[2]]

plot(p193_pte, p193_pxm)

abline(pxpte_193)

mtext(paste("m =", round(pxpte_m, 3)), side=3, line=-1, at=0, cex=0.8, adj=0)
mtext(paste("r2 =", round(pxpte_r2, 3)), side=3, line=-2, at=0, cex=0.8, adj=0)
pxpte_p = summary(pxpte_193)$coefficients[,4]
mtext(paste("p =", round(pxpte_p, 3)), side=3, line=-3, at=0, cex=0.8, adj=0)

dev.off()


pdf(file=paste("193_ptmvspxe", "pdf", sep="."))

ptPXE_193 <- lm(p193_ptm~p193_pxe)
ptPXE_r2 = summary(ptPXE_193)$r.squared
ptPXE_m = coefficients(ptPXE_193)[[2]]

plot(p193_pxe, p193_ptm)

abline(ptPXE_193)

mtext(paste("m =", round(ptPXE_m, 3)), side=3, line=-1, at=0, cex=0.8, adj=0)
mtext(paste("r2 =", round(ptPXE_r2, 3)), side=3, line=-2, at=0, cex=0.8, adj=0)
ptPXE_p = summary(ptPXE_193)$coefficients[,4]
mtext(paste("p =", round(ptPXE_p, 3)), side=3, line=-3, at=0, cex=0.8, adj=0)

dev.off()


pdf(file=paste("193_ptmvspte", "pdf", sep="."))

ptpte_193 <- lm(p193_ptm~p193_pte)
ptpte_r2 = summary(ptpte_193)$r.squared
ptpte_m = coefficients(ptpte_193)[[2]]

plot(p193_pte, p193_ptm)

abline(ptpte_193)

mtext(paste("m =", round(ptpte_m, 3)), side=3, line=-1, at=0, cex=0.8, adj=0)
mtext(paste("r2 =", round(ptpte_r2, 3)), side=3, line=-2, at=0, cex=0.8, adj=0)
ptpte_p = summary(ptpte_193)$coefficients[,4]
mtext(paste("p =", round(ptpte_p, 3)), side=3, line=-3, at=0, cex=0.8, adj=0)

dev.off()

pdf(file=paste("193_pxevspte", "pdf", sep="."))

pompte_193 <- lm(p193_pxe~p193_pte)
pompte_r2 = summary(pompte_193)$r.squared
pompte_m = coefficients(pompte_193)[[2]]

plot(p193_pte, p193_pxe)

abline(pompte_193)

mtext(paste("m =", round(pompte_m, 3)), side=3, line=-1, at=0, cex=0.8, adj=0)
mtext(paste("r2 =", round(pompte_r2, 3)), side=3, line=-2, at=0, cex=0.8, adj=0)
pompte_p = summary(pompte_193)$coefficients[,4]
mtext(paste("p =", round(pompte_p, 3)), side=3, line=-3, at=0, cex=0.8, adj=0)

dev.off()




p258_pxm <- array(NA, dim=c(16))

p258_ptm <- array(NA, dim=c(16))

p258_pxe <- array(NA, dim=c(16))

p258_pte <- array(NA, dim=c(16))

# Computing effects


p258_pxm[1] <- traj_pxm[1,1,3,1,1] - traj_pxm[1,1,1,1,1]
p258_pxm[2] <- traj_pxm[3,1,3,1,1] - traj_pxm[3,1,1,1,1]
p258_pxm[3] <- traj_pxm[1,3,3,1,1] - traj_pxm[1,3,1,1,1]
p258_pxm[4] <- traj_pxm[1,1,3,3,1] - traj_pxm[1,1,1,3,1]
p258_pxm[5] <- traj_pxm[1,1,3,1,3] - traj_pxm[1,1,1,1,3]
p258_pxm[6] <- traj_pxm[3,3,3,1,1] - traj_pxm[3,3,1,1,1]
p258_pxm[7] <- traj_pxm[3,1,3,3,1] - traj_pxm[3,1,1,3,1]
p258_pxm[8] <- traj_pxm[3,1,3,1,3] - traj_pxm[3,1,1,1,3]
p258_pxm[9] <- traj_pxm[1,3,3,3,1] - traj_pxm[1,3,1,3,1]
p258_pxm[10] <- traj_pxm[1,3,3,1,3] - traj_pxm[1,3,1,1,3]
p258_pxm[11] <- traj_pxm[1,1,3,3,3] - traj_pxm[1,1,1,3,3]
p258_pxm[12] <- traj_pxm[3,3,3,3,1] - traj_pxm[3,3,1,3,1]
p258_pxm[13] <- traj_pxm[3,3,3,1,3] - traj_pxm[3,3,1,1,3]
p258_pxm[14] <- traj_pxm[3,1,3,3,3] - traj_pxm[3,1,1,3,3]
p258_pxm[15] <- traj_pxm[1,3,3,3,3] - traj_pxm[1,3,1,3,3]
p258_pxm[16] <- traj_pxm[3,3,3,3,3] - traj_pxm[3,3,1,3,3]

p258_ptm[1] <- traj_ptm[1,1,3,1,1] - traj_ptm[1,1,1,1,1]
p258_ptm[2] <- traj_ptm[3,1,3,1,1] - traj_ptm[3,1,1,1,1]
p258_ptm[3] <- traj_ptm[1,3,3,1,1] - traj_ptm[1,3,1,1,1]
p258_ptm[4] <- traj_ptm[1,1,3,3,1] - traj_ptm[1,1,1,3,1]
p258_ptm[5] <- traj_ptm[1,1,3,1,3] - traj_ptm[1,1,1,1,3]
p258_ptm[6] <- traj_ptm[3,3,3,1,1] - traj_ptm[3,3,1,1,1]
p258_ptm[7] <- traj_ptm[3,1,3,3,1] - traj_ptm[3,1,1,3,1]
p258_ptm[8] <- traj_ptm[3,1,3,1,3] - traj_ptm[3,1,1,1,3]
p258_ptm[9] <- traj_ptm[1,3,3,3,1] - traj_ptm[1,3,1,3,1]
p258_ptm[10] <- traj_ptm[1,3,3,1,3] - traj_ptm[1,3,1,1,3]
p258_ptm[11] <- traj_ptm[1,1,3,3,3] - traj_ptm[1,1,1,3,3]
p258_ptm[12] <- traj_ptm[3,3,3,3,1] - traj_ptm[3,3,1,3,1]
p258_ptm[13] <- traj_ptm[3,3,3,1,3] - traj_ptm[3,3,1,1,3]
p258_ptm[14] <- traj_ptm[3,1,3,3,3] - traj_ptm[3,1,1,3,3]
p258_ptm[15] <- traj_ptm[1,3,3,3,3] - traj_ptm[1,3,1,3,3]
p258_ptm[16] <- traj_ptm[3,3,3,3,3] - traj_ptm[3,3,1,3,3]

p258_pxe[1] <- traj_pxe[1,1,3,1,1] - traj_pxe[1,1,1,1,1]
p258_pxe[2] <- traj_pxe[3,1,3,1,1] - traj_pxe[3,1,1,1,1]
p258_pxe[3] <- traj_pxe[1,3,3,1,1] - traj_pxe[1,3,1,1,1]
p258_pxe[4] <- traj_pxe[1,1,3,3,1] - traj_pxe[1,1,1,3,1]
p258_pxe[5] <- traj_pxe[1,1,3,1,3] - traj_pxe[1,1,1,1,3]
p258_pxe[6] <- traj_pxe[3,3,3,1,1] - traj_pxe[3,3,1,1,1]
p258_pxe[7] <- traj_pxe[3,1,3,3,1] - traj_pxe[3,1,1,3,1]
p258_pxe[8] <- traj_pxe[3,1,3,1,3] - traj_pxe[3,1,1,1,3]
p258_pxe[9] <- traj_pxe[1,3,3,3,1] - traj_pxe[1,3,1,3,1]
p258_pxe[10] <- traj_pxe[1,3,3,1,3] - traj_pxe[1,3,1,1,3]
p258_pxe[11] <- traj_pxe[1,1,3,3,3] - traj_pxe[1,1,1,3,3]
p258_pxe[12] <- traj_pxe[3,3,3,3,1] - traj_pxe[3,3,1,3,1]
p258_pxe[13] <- traj_pxe[3,3,3,1,3] - traj_pxe[3,3,1,1,3]
p258_pxe[14] <- traj_pxe[3,1,3,3,3] - traj_pxe[3,1,1,3,3]
p258_pxe[15] <- traj_pxe[1,3,3,3,3] - traj_pxe[1,3,1,3,3]
p258_pxe[16] <- traj_pxe[3,3,3,3,3] - traj_pxe[3,3,1,3,3]

p258_pte[1] <- traj_pte[1,1,3,1,1] - traj_pte[1,1,1,1,1]
p258_pte[2] <- traj_pte[3,1,3,1,1] - traj_pte[3,1,1,1,1]
p258_pte[3] <- traj_pte[1,3,3,1,1] - traj_pte[1,3,1,1,1]
p258_pte[4] <- traj_pte[1,1,3,3,1] - traj_pte[1,1,1,3,1]
p258_pte[5] <- traj_pte[1,1,3,1,3] - traj_pte[1,1,1,1,3]
p258_pte[6] <- traj_pte[3,3,3,1,1] - traj_pte[3,3,1,1,1]
p258_pte[7] <- traj_pte[3,1,3,3,1] - traj_pte[3,1,1,3,1]
p258_pte[8] <- traj_pte[3,1,3,1,3] - traj_pte[3,1,1,1,3]
p258_pte[9] <- traj_pte[1,3,3,3,1] - traj_pte[1,3,1,3,1]
p258_pte[10] <- traj_pte[1,3,3,1,3] - traj_pte[1,3,1,1,3]
p258_pte[11] <- traj_pte[1,1,3,3,3] - traj_pte[1,1,1,3,3]
p258_pte[12] <- traj_pte[3,3,3,3,1] - traj_pte[3,3,1,3,1]
p258_pte[13] <- traj_pte[3,3,3,1,3] - traj_pte[3,3,1,1,3]
p258_pte[14] <- traj_pte[3,1,3,3,3] - traj_pte[3,1,1,3,3]
p258_pte[15] <- traj_pte[1,3,3,3,3] - traj_pte[1,3,1,3,3]
p258_pte[16] <- traj_pte[3,3,3,3,3] - traj_pte[3,3,1,3,3]

write("Position 258 Effects", file="Pos258.csv")
write("genotype, pxm, ptm, pxe, pte", file="Pos258.csv", append=TRUE)
write(paste("aaxaa", p258_pxm[1],p258_ptm[1],p258_pxe[1],p258_pte[1], sep=","), file="Pos258.csv", append=TRUE)
write(paste("daxaa", p258_pxm[2],p258_ptm[2],p258_pxe[2],p258_pte[2], sep=","), file="Pos258.csv", append=TRUE)
write(paste("adxaa", p258_pxm[3],p258_ptm[3],p258_pxe[3],p258_pte[3], sep=","), file="Pos258.csv", append=TRUE)
write(paste("aaxda", p258_pxm[4],p258_ptm[4],p258_pxe[4],p258_pte[4], sep=","), file="Pos258.csv", append=TRUE)
write(paste("aaxad", p258_pxm[5],p258_ptm[5],p258_pxe[5],p258_pte[5], sep=","), file="Pos258.csv", append=TRUE)
write(paste("ddxaa", p258_pxm[6],p258_ptm[6],p258_pxe[6],p258_pte[6], sep=","), file="Pos258.csv", append=TRUE)
write(paste("daxda", p258_pxm[7],p258_ptm[7],p258_pxe[7],p258_pte[7], sep=","), file="Pos258.csv", append=TRUE)
write(paste("daxad", p258_pxm[8],p258_ptm[8],p258_pxe[8],p258_pte[8], sep=","), file="Pos258.csv", append=TRUE)
write(paste("adxda", p258_pxm[9],p258_ptm[9],p258_pxe[9],p258_pte[9], sep=","), file="Pos258.csv", append=TRUE)
write(paste("adxad", p258_pxm[10],p258_ptm[10],p258_pxe[10],p258_pte[10], sep=","), file="Pos258.csv", append=TRUE)
write(paste("aaxdd", p258_pxm[11],p258_ptm[11],p258_pxe[11],p258_pte[11], sep=","), file="Pos258.csv", append=TRUE)
write(paste("ddxda", p258_pxm[12],p258_ptm[12],p258_pxe[12],p258_pte[12], sep=","), file="Pos258.csv", append=TRUE)
write(paste("ddxad", p258_pxm[13],p258_ptm[13],p258_pxe[13],p258_pte[13], sep=","), file="Pos258.csv", append=TRUE)
write(paste("daxdd", p258_pxm[14],p258_ptm[14],p258_pxe[14],p258_pte[14], sep=","), file="Pos258.csv", append=TRUE)
write(paste("adxdd", p258_pxm[15],p258_ptm[15],p258_pxe[15],p258_pte[15], sep=","), file="Pos258.csv", append=TRUE)
write(paste("ddxdd", p258_pxm[16],p258_ptm[16],p258_pxe[16],p258_pte[16], sep=","), file="Pos258.csv", append=TRUE)



pdf(file=paste("258_pxm.pdf"))

barplot(p258_pxm)

dev.off()

pdf(file=paste("258_ptm.pdf"))

barplot(p258_ptm)

dev.off()

pdf(file=paste("258_pxe.pdf"))

barplot(p258_pxe)

dev.off()

pdf(file=paste("258_pte.pdf"))

barplot(p258_pte)

dev.off()


pdf(file=paste("258_pxmvspxe", "pdf", sep="."))

pxPXE_258 <- lm(p258_pxm~p258_pxe)
pxPXE_r2 = summary(pxPXE_258)$r.squared
pxPXE_m = coefficients(pxPXE_258)[[2]]

plot(p258_pxe, p258_pxm)

abline(pxPXE_258)

mtext(paste("m =", round(pxPXE_m, 3)), side=3, line=-1, at=0, cex=0.8, adj=0)
mtext(paste("r2 =", round(pxPXE_r2, 3)), side=3, line=-2, at=0, cex=0.8, adj=0)
pxPXE_p = summary(pxPXE_258)$coefficients[,4]
mtext(paste("p =", round(pxPXE_p, 3)), side=3, line=-3, at=0, cex=0.8, adj=0)

dev.off()


pdf(file=paste("258_pxmvsptm", "pdf", sep="."))

pxPTM_258 <- lm(p258_pxm~p258_ptm)
pxPTM_r2 = summary(pxPTM_258)$r.squared
pxPTM_m = coefficients(pxPTM_258)[[2]]

plot(p258_ptm, p258_pxm)

abline(pxPTM_258)

mtext(paste("m =", round(pxPTM_m, 3)), side=3, line=-1, at=0, cex=0.8, adj=0)
mtext(paste("r2 =", round(pxPTM_r2, 3)), side=3, line=-2, at=0, cex=0.8, adj=0)
pxPTM_p = summary(pxPTM_258)$coefficients[,4]
mtext(paste("p =", round(pxPTM_p, 3)), side=3, line=-3, at=0, cex=0.8, adj=0)

dev.off()


pdf(file=paste("258_pxmvspte", "pdf", sep="."))

pxpte_258 <- lm(p258_pxm~p258_pte)
pxpte_r2 = summary(pxpte_258)$r.squared
pxpte_m = coefficients(pxpte_258)[[2]]

plot(p258_pte, p258_pxm)

abline(pxpte_258)

mtext(paste("m =", round(pxpte_m, 3)), side=3, line=-1, at=0, cex=0.8, adj=0)
mtext(paste("r2 =", round(pxpte_r2, 3)), side=3, line=-2, at=0, cex=0.8, adj=0)
pxpte_p = summary(pxpte_258)$coefficients[,4]
mtext(paste("p =", round(pxpte_p, 3)), side=3, line=-3, at=0, cex=0.8, adj=0)

dev.off()


pdf(file=paste("258_ptmvspxe", "pdf", sep="."))

ptPXE_258 <- lm(p258_ptm~p258_pxe)
ptPXE_r2 = summary(ptPXE_258)$r.squared
ptPXE_m = coefficients(ptPXE_258)[[2]]

plot(p258_pxe, p258_ptm)

abline(ptPXE_258)

mtext(paste("m =", round(ptPXE_m, 3)), side=3, line=-1, at=0, cex=0.8, adj=0)
mtext(paste("r2 =", round(ptPXE_r2, 3)), side=3, line=-2, at=0, cex=0.8, adj=0)
ptPXE_p = summary(ptPXE_258)$coefficients[,4]
mtext(paste("p =", round(ptPXE_p, 3)), side=3, line=-3, at=0, cex=0.8, adj=0)

dev.off()


pdf(file=paste("258_ptmvspte", "pdf", sep="."))

ptpte_258 <- lm(p258_ptm~p258_pte)
ptpte_r2 = summary(ptpte_258)$r.squared
ptpte_m = coefficients(ptpte_258)[[2]]

plot(p258_pte, p258_ptm)

abline(ptpte_258)

mtext(paste("m =", round(ptpte_m, 3)), side=3, line=-1, at=0, cex=0.8, adj=0)
mtext(paste("r2 =", round(ptpte_r2, 3)), side=3, line=-2, at=0, cex=0.8, adj=0)
ptpte_p = summary(ptpte_258)$coefficients[,4]
mtext(paste("p =", round(ptpte_p, 3)), side=3, line=-3, at=0, cex=0.8, adj=0)

dev.off()

pdf(file=paste("258_pxevspte", "pdf", sep="."))

pompte_258 <- lm(p258_pxe~p258_pte)
pompte_r2 = summary(pompte_258)$r.squared
pompte_m = coefficients(pompte_258)[[2]]

plot(p258_pte, p258_pxe)

abline(pompte_258)

mtext(paste("m =", round(pompte_m, 3)), side=3, line=-1, at=0, cex=0.8, adj=0)
mtext(paste("r2 =", round(pompte_r2, 3)), side=3, line=-2, at=0, cex=0.8, adj=0)
pompte_p = summary(pompte_258)$coefficients[,4]
mtext(paste("p =", round(pompte_p, 3)), side=3, line=-3, at=0, cex=0.8, adj=0)

dev.off()



p271_pxm <- array(NA, dim=c(16))

p271_ptm <- array(NA, dim=c(16))

p271_pxe <- array(NA, dim=c(16))

p271_pte <- array(NA, dim=c(16))

# Computing effects

p271_pxm[1] <- traj_pxm[1,1,1,3,1] - traj_pxm[1,1,1,1,1]
p271_pxm[2] <- traj_pxm[3,1,1,3,1] - traj_pxm[3,1,1,1,1]
p271_pxm[3] <- traj_pxm[1,3,1,3,1] - traj_pxm[1,3,1,1,1]
p271_pxm[4] <- traj_pxm[1,1,3,3,1] - traj_pxm[1,1,3,1,1]
p271_pxm[5] <- traj_pxm[1,1,1,3,3] - traj_pxm[1,1,1,1,3]
p271_pxm[6] <- traj_pxm[3,3,1,3,1] - traj_pxm[3,3,1,1,1]
p271_pxm[7] <- traj_pxm[3,1,3,3,1] - traj_pxm[3,1,3,1,1]
p271_pxm[8] <- traj_pxm[3,1,1,3,3] - traj_pxm[3,1,1,1,3]
p271_pxm[9] <- traj_pxm[1,3,3,3,1] - traj_pxm[1,3,3,1,1]
p271_pxm[10] <- traj_pxm[1,3,1,3,3] - traj_pxm[1,3,1,1,3]
p271_pxm[11] <- traj_pxm[1,1,3,3,3] - traj_pxm[1,1,3,1,3]
p271_pxm[12] <- traj_pxm[3,3,3,3,1] - traj_pxm[3,3,3,1,1]
p271_pxm[13] <- traj_pxm[3,3,1,3,3] - traj_pxm[3,3,1,1,3]
p271_pxm[14] <- traj_pxm[3,1,3,3,3] - traj_pxm[3,1,3,1,3]
p271_pxm[15] <- traj_pxm[1,3,3,3,3] - traj_pxm[1,3,3,1,3]
p271_pxm[16] <- traj_pxm[3,3,3,3,3] - traj_pxm[3,3,3,1,3]

p271_ptm[1] <- traj_ptm[1,1,1,3,1] - traj_ptm[1,1,1,1,1]
p271_ptm[2] <- traj_ptm[3,1,1,3,1] - traj_ptm[3,1,1,1,1]
p271_ptm[3] <- traj_ptm[1,3,1,3,1] - traj_ptm[1,3,1,1,1]
p271_ptm[4] <- traj_ptm[1,1,3,3,1] - traj_ptm[1,1,3,1,1]
p271_ptm[5] <- traj_ptm[1,1,1,3,3] - traj_ptm[1,1,1,1,3]
p271_ptm[6] <- traj_ptm[3,3,1,3,1] - traj_ptm[3,3,1,1,1]
p271_ptm[7] <- traj_ptm[3,1,3,3,1] - traj_ptm[3,1,3,1,1]
p271_ptm[8] <- traj_ptm[3,1,1,3,3] - traj_ptm[3,1,1,1,3]
p271_ptm[9] <- traj_ptm[1,3,3,3,1] - traj_ptm[1,3,3,1,1]
p271_ptm[10] <- traj_ptm[1,3,1,3,3] - traj_ptm[1,3,1,1,3]
p271_ptm[11] <- traj_ptm[1,1,3,3,3] - traj_ptm[1,1,3,1,3]
p271_ptm[12] <- traj_ptm[3,3,3,3,1] - traj_ptm[3,3,3,1,1]
p271_ptm[13] <- traj_ptm[3,3,1,3,3] - traj_ptm[3,3,1,1,3]
p271_ptm[14] <- traj_ptm[3,1,3,3,3] - traj_ptm[3,1,3,1,3]
p271_ptm[15] <- traj_ptm[1,3,3,3,3] - traj_ptm[1,3,3,1,3]
p271_ptm[16] <- traj_ptm[3,3,3,3,3] - traj_ptm[3,3,3,1,3]

p271_pxe[1] <- traj_pxe[1,1,1,3,1] - traj_pxe[1,1,1,1,1]
p271_pxe[2] <- traj_pxe[3,1,1,3,1] - traj_pxe[3,1,1,1,1]
p271_pxe[3] <- traj_pxe[1,3,1,3,1] - traj_pxe[1,3,1,1,1]
p271_pxe[4] <- traj_pxe[1,1,3,3,1] - traj_pxe[1,1,3,1,1]
p271_pxe[5] <- traj_pxe[1,1,1,3,3] - traj_pxe[1,1,1,1,3]
p271_pxe[6] <- traj_pxe[3,3,1,3,1] - traj_pxe[3,3,1,1,1]
p271_pxe[7] <- traj_pxe[3,1,3,3,1] - traj_pxe[3,1,3,1,1]
p271_pxe[8] <- traj_pxe[3,1,1,3,3] - traj_pxe[3,1,1,1,3]
p271_pxe[9] <- traj_pxe[1,3,3,3,1] - traj_pxe[1,3,3,1,1]
p271_pxe[10] <- traj_pxe[1,3,1,3,3] - traj_pxe[1,3,1,1,3]
p271_pxe[11] <- traj_pxe[1,1,3,3,3] - traj_pxe[1,1,3,1,3]
p271_pxe[12] <- traj_pxe[3,3,3,3,1] - traj_pxe[3,3,3,1,1]
p271_pxe[13] <- traj_pxe[3,3,1,3,3] - traj_pxe[3,3,1,1,3]
p271_pxe[14] <- traj_pxe[3,1,3,3,3] - traj_pxe[3,1,3,1,3]
p271_pxe[15] <- traj_pxe[1,3,3,3,3] - traj_pxe[1,3,3,1,3]
p271_pxe[16] <- traj_pxe[3,3,3,3,3] - traj_pxe[3,3,3,1,3]

p271_pte[1] <- traj_pte[1,1,1,3,1] - traj_pte[1,1,1,1,1]
p271_pte[2] <- traj_pte[3,1,1,3,1] - traj_pte[3,1,1,1,1]
p271_pte[3] <- traj_pte[1,3,1,3,1] - traj_pte[1,3,1,1,1]
p271_pte[4] <- traj_pte[1,1,3,3,1] - traj_pte[1,1,3,1,1]
p271_pte[5] <- traj_pte[1,1,1,3,3] - traj_pte[1,1,1,1,3]
p271_pte[6] <- traj_pte[3,3,1,3,1] - traj_pte[3,3,1,1,1]
p271_pte[7] <- traj_pte[3,1,3,3,1] - traj_pte[3,1,3,1,1]
p271_pte[8] <- traj_pte[3,1,1,3,3] - traj_pte[3,1,1,1,3]
p271_pte[9] <- traj_pte[1,3,3,3,1] - traj_pte[1,3,3,1,1]
p271_pte[10] <- traj_pte[1,3,1,3,3] - traj_pte[1,3,1,1,3]
p271_pte[11] <- traj_pte[1,1,3,3,3] - traj_pte[1,1,3,1,3]
p271_pte[12] <- traj_pte[3,3,3,3,1] - traj_pte[3,3,3,1,1]
p271_pte[13] <- traj_pte[3,3,1,3,3] - traj_pte[3,3,1,1,3]
p271_pte[14] <- traj_pte[3,1,3,3,3] - traj_pte[3,1,3,1,3]
p271_pte[15] <- traj_pte[1,3,3,3,3] - traj_pte[1,3,3,1,3]
p271_pte[16] <- traj_pte[3,3,3,3,3] - traj_pte[3,3,3,1,3]

write("Position 271 Effects", file="Pos271.csv")
write("genotype, pxm, ptm, pxe, pte", file="Pos271.csv", append=TRUE)
write(paste("aaaxa", p271_pxm[1],p271_ptm[1],p271_pxe[1],p271_pte[1], sep=","), file="Pos271.csv", append=TRUE)
write(paste("daaxa", p271_pxm[2],p271_ptm[2],p271_pxe[2],p271_pte[2], sep=","), file="Pos271.csv", append=TRUE)
write(paste("adaxa", p271_pxm[3],p271_ptm[3],p271_pxe[3],p271_pte[3], sep=","), file="Pos271.csv", append=TRUE)
write(paste("aadxa", p271_pxm[4],p271_ptm[4],p271_pxe[4],p271_pte[4], sep=","), file="Pos271.csv", append=TRUE)
write(paste("aaaxd", p271_pxm[5],p271_ptm[5],p271_pxe[5],p271_pte[5], sep=","), file="Pos271.csv", append=TRUE)
write(paste("ddaxa", p271_pxm[6],p271_ptm[6],p271_pxe[6],p271_pte[6], sep=","), file="Pos271.csv", append=TRUE)
write(paste("dadxa", p271_pxm[7],p271_ptm[7],p271_pxe[7],p271_pte[7], sep=","), file="Pos271.csv", append=TRUE)
write(paste("daaxd", p271_pxm[8],p271_ptm[8],p271_pxe[8],p271_pte[8], sep=","), file="Pos271.csv", append=TRUE)
write(paste("addxa", p271_pxm[9],p271_ptm[9],p271_pxe[9],p271_pte[9], sep=","), file="Pos271.csv", append=TRUE)
write(paste("adaxd", p271_pxm[10],p271_ptm[10],p271_pxe[10],p271_pte[10], sep=","), file="Pos271.csv", append=TRUE)
write(paste("aadxd", p271_pxm[11],p271_ptm[11],p271_pxe[11],p271_pte[11], sep=","), file="Pos271.csv", append=TRUE)
write(paste("dddxa", p271_pxm[12],p271_ptm[12],p271_pxe[12],p271_pte[12], sep=","), file="Pos271.csv", append=TRUE)
write(paste("ddaxd", p271_pxm[13],p271_ptm[13],p271_pxe[13],p271_pte[13], sep=","), file="Pos271.csv", append=TRUE)
write(paste("dadxd", p271_pxm[14],p271_ptm[14],p271_pxe[14],p271_pte[14], sep=","), file="Pos271.csv", append=TRUE)
write(paste("addxd", p271_pxm[15],p271_ptm[15],p271_pxe[15],p271_pte[15], sep=","), file="Pos271.csv", append=TRUE)
write(paste("dddxd", p271_pxm[16],p271_ptm[16],p271_pxe[16],p271_pte[16], sep=","), file="Pos271.csv", append=TRUE)



pdf(file=paste("271_pxm.pdf"))

barplot(p271_pxm)

dev.off()

pdf(file=paste("271_ptm.pdf"))

barplot(p271_ptm)

dev.off()

pdf(file=paste("271_pxe.pdf"))

barplot(p271_pxe)

dev.off()

pdf(file=paste("271_pte.pdf"))

barplot(p271_pte)

dev.off()


pdf(file=paste("271_pxmvspxe", "pdf", sep="."))

pxPXE_271 <- lm(p271_pxm~p271_pxe)
pxPXE_r2 = summary(pxPXE_271)$r.squared
pxPXE_m = coefficients(pxPXE_271)[[2]]

plot(p271_pxe, p271_pxm)

abline(pxPXE_271)

mtext(paste("m =", round(pxPXE_m, 3)), side=3, line=-1, at=0, cex=0.8, adj=0)
mtext(paste("r2 =", round(pxPXE_r2, 3)), side=3, line=-2, at=0, cex=0.8, adj=0)
pxPXE_p = summary(pxPXE_271)$coefficients[,4]
mtext(paste("p =", round(pxPXE_p, 3)), side=3, line=-3, at=0, cex=0.8, adj=0)

dev.off()


pdf(file=paste("271_pxmvsptm", "pdf", sep="."))

pxPTM_271 <- lm(p271_pxm~p271_ptm)
pxPTM_r2 = summary(pxPTM_271)$r.squared
pxPTM_m = coefficients(pxPTM_271)[[2]]

plot(p271_ptm, p271_pxm)

abline(pxPTM_271)

mtext(paste("m =", round(pxPTM_m, 3)), side=3, line=-1, at=0, cex=0.8, adj=0)
mtext(paste("r2 =", round(pxPTM_r2, 3)), side=3, line=-2, at=0, cex=0.8, adj=0)
pxPTM_p = summary(pxPTM_271)$coefficients[,4]
mtext(paste("p =", round(pxPTM_p, 3)), side=3, line=-3, at=0, cex=0.8, adj=0)

dev.off()


pdf(file=paste("271_pxmvspte", "pdf", sep="."))

pxpte_271 <- lm(p271_pxm~p271_pte)
pxpte_r2 = summary(pxpte_271)$r.squared
pxpte_m = coefficients(pxpte_271)[[2]]

plot(p271_pte, p271_pxm)

abline(pxpte_271)

mtext(paste("m =", round(pxpte_m, 3)), side=3, line=-1, at=0, cex=0.8, adj=0)
mtext(paste("r2 =", round(pxpte_r2, 3)), side=3, line=-2, at=0, cex=0.8, adj=0)
pxpte_p = summary(pxpte_271)$coefficients[,4]
mtext(paste("p =", round(pxpte_p, 3)), side=3, line=-3, at=0, cex=0.8, adj=0)

dev.off()


pdf(file=paste("271_ptmvspxe", "pdf", sep="."))

ptPXE_271 <- lm(p271_ptm~p271_pxe)
ptPXE_r2 = summary(ptPXE_271)$r.squared
ptPXE_m = coefficients(ptPXE_271)[[2]]

plot(p271_pxe, p271_ptm)

abline(ptPXE_271)

mtext(paste("m =", round(ptPXE_m, 3)), side=3, line=-1, at=0, cex=0.8, adj=0)
mtext(paste("r2 =", round(ptPXE_r2, 3)), side=3, line=-2, at=0, cex=0.8, adj=0)
ptPXE_p = summary(ptPXE_271)$coefficients[,4]
mtext(paste("p =", round(ptPXE_p, 3)), side=3, line=-3, at=0, cex=0.8, adj=0)

dev.off()


pdf(file=paste("271_ptmvspte", "pdf", sep="."))

ptpte_271 <- lm(p271_ptm~p271_pte)
ptpte_r2 = summary(ptpte_271)$r.squared
ptpte_m = coefficients(ptpte_271)[[2]]

plot(p271_pte, p271_ptm)

abline(ptpte_271)

mtext(paste("m =", round(ptpte_m, 3)), side=3, line=-1, at=0, cex=0.8, adj=0)
mtext(paste("r2 =", round(ptpte_r2, 3)), side=3, line=-2, at=0, cex=0.8, adj=0)
ptpte_p = summary(ptpte_271)$coefficients[,4]
mtext(paste("p =", round(ptpte_p, 3)), side=3, line=-3, at=0, cex=0.8, adj=0)

dev.off()

pdf(file=paste("271_pxevspte", "pdf", sep="."))

pompte_271 <- lm(p271_pxe~p271_pte)
pompte_r2 = summary(pompte_271)$r.squared
pompte_m = coefficients(pompte_271)[[2]]

plot(p271_pte, p271_pxe)

abline(pompte_271)

mtext(paste("m =", round(pompte_m, 3)), side=3, line=-1, at=0, cex=0.8, adj=0)
mtext(paste("r2 =", round(pompte_r2, 3)), side=3, line=-2, at=0, cex=0.8, adj=0)
pompte_p = summary(pompte_271)$coefficients[,4]
mtext(paste("p =", round(pompte_p, 3)), side=3, line=-3, at=0, cex=0.8, adj=0)

dev.off()




#p273_dhc <- array(NA, dim=c(16))

p273_pxm <- array(NA, dim=c(16))

p273_ptm <- array(NA, dim=c(16))

p273_pxe <- array(NA, dim=c(16))

p273_pte <- array(NA, dim=c(16))

# Computing effects

p273_pxm[1] <- traj_pxm[1,1,1,1,3] - traj_pxm[1,1,1,1,1]
p273_pxm[2] <- traj_pxm[3,1,1,1,3] - traj_pxm[3,1,1,1,1]
p273_pxm[3] <- traj_pxm[1,3,1,1,3] - traj_pxm[1,3,1,1,1]
p273_pxm[4] <- traj_pxm[1,1,3,1,3] - traj_pxm[1,1,3,1,1]
p273_pxm[5] <- traj_pxm[1,1,1,3,3] - traj_pxm[1,1,1,3,1]
p273_pxm[6] <- traj_pxm[3,3,1,1,3] - traj_pxm[3,3,1,1,1]
p273_pxm[7] <- traj_pxm[3,1,3,1,3] - traj_pxm[3,1,3,1,1]
p273_pxm[8] <- traj_pxm[3,1,1,3,3] - traj_pxm[3,1,1,3,1]
p273_pxm[9] <- traj_pxm[1,3,3,1,3] - traj_pxm[1,3,3,1,1]
p273_pxm[10] <- traj_pxm[1,3,1,3,3] - traj_pxm[1,3,1,3,1]
p273_pxm[11] <- traj_pxm[1,1,3,3,3] - traj_pxm[1,1,3,3,1]
p273_pxm[12] <- traj_pxm[3,3,3,1,3] - traj_pxm[3,3,3,1,1]
p273_pxm[13] <- traj_pxm[3,3,1,3,3] - traj_pxm[3,3,1,3,1]
p273_pxm[14] <- traj_pxm[3,1,3,3,3] - traj_pxm[3,1,3,3,1]
p273_pxm[15] <- traj_pxm[1,3,3,3,3] - traj_pxm[1,3,3,3,1]
p273_pxm[16] <- traj_pxm[3,3,3,3,3] - traj_pxm[3,3,3,3,1]

p273_ptm[1] <- traj_ptm[1,1,1,1,3] - traj_ptm[1,1,1,1,1]
p273_ptm[2] <- traj_ptm[3,1,1,1,3] - traj_ptm[3,1,1,1,1]
p273_ptm[3] <- traj_ptm[1,3,1,1,3] - traj_ptm[1,3,1,1,1]
p273_ptm[4] <- traj_ptm[1,1,3,1,3] - traj_ptm[1,1,3,1,1]
p273_ptm[5] <- traj_ptm[1,1,1,3,3] - traj_ptm[1,1,1,3,1]
p273_ptm[6] <- traj_ptm[3,3,1,1,3] - traj_ptm[3,3,1,1,1]
p273_ptm[7] <- traj_ptm[3,1,3,1,3] - traj_ptm[3,1,3,1,1]
p273_ptm[8] <- traj_ptm[3,1,1,3,3] - traj_ptm[3,1,1,3,1]
p273_ptm[9] <- traj_ptm[1,3,3,1,3] - traj_ptm[1,3,3,1,1]
p273_ptm[10] <- traj_ptm[1,3,1,3,3] - traj_ptm[1,3,1,3,1]
p273_ptm[11] <- traj_ptm[1,1,3,3,3] - traj_ptm[1,1,3,3,1]
p273_ptm[12] <- traj_ptm[3,3,3,1,3] - traj_ptm[3,3,3,1,1]
p273_ptm[13] <- traj_ptm[3,3,1,3,3] - traj_ptm[3,3,1,3,1]
p273_ptm[14] <- traj_ptm[3,1,3,3,3] - traj_ptm[3,1,3,3,1]
p273_ptm[15] <- traj_ptm[1,3,3,3,3] - traj_ptm[1,3,3,3,1]
p273_ptm[16] <- traj_ptm[3,3,3,3,3] - traj_ptm[3,3,3,3,1]

p273_pxe[1] <- traj_pxe[1,1,1,1,3] - traj_pxe[1,1,1,1,1]
p273_pxe[2] <- traj_pxe[3,1,1,1,3] - traj_pxe[3,1,1,1,1]
p273_pxe[3] <- traj_pxe[1,3,1,1,3] - traj_pxe[1,3,1,1,1]
p273_pxe[4] <- traj_pxe[1,1,3,1,3] - traj_pxe[1,1,3,1,1]
p273_pxe[5] <- traj_pxe[1,1,1,3,3] - traj_pxe[1,1,1,3,1]
p273_pxe[6] <- traj_pxe[3,3,1,1,3] - traj_pxe[3,3,1,1,1]
p273_pxe[7] <- traj_pxe[3,1,3,1,3] - traj_pxe[3,1,3,1,1]
p273_pxe[8] <- traj_pxe[3,1,1,3,3] - traj_pxe[3,1,1,3,1]
p273_pxe[9] <- traj_pxe[1,3,3,1,3] - traj_pxe[1,3,3,1,1]
p273_pxe[10] <- traj_pxe[1,3,1,3,3] - traj_pxe[1,3,1,3,1]
p273_pxe[11] <- traj_pxe[1,1,3,3,3] - traj_pxe[1,1,3,3,1]
p273_pxe[12] <- traj_pxe[3,3,3,1,3] - traj_pxe[3,3,3,1,1]
p273_pxe[13] <- traj_pxe[3,3,1,3,3] - traj_pxe[3,3,1,3,1]
p273_pxe[14] <- traj_pxe[3,1,3,3,3] - traj_pxe[3,1,3,3,1]
p273_pxe[15] <- traj_pxe[1,3,3,3,3] - traj_pxe[1,3,3,3,1]
p273_pxe[16] <- traj_pxe[3,3,3,3,3] - traj_pxe[3,3,3,3,1]

p273_pte[1] <- traj_pte[1,1,1,1,3] - traj_pte[1,1,1,1,1]
p273_pte[2] <- traj_pte[3,1,1,1,3] - traj_pte[3,1,1,1,1]
p273_pte[3] <- traj_pte[1,3,1,1,3] - traj_pte[1,3,1,1,1]
p273_pte[4] <- traj_pte[1,1,3,1,3] - traj_pte[1,1,3,1,1]
p273_pte[5] <- traj_pte[1,1,1,3,3] - traj_pte[1,1,1,3,1]
p273_pte[6] <- traj_pte[3,3,1,1,3] - traj_pte[3,3,1,1,1]
p273_pte[7] <- traj_pte[3,1,3,1,3] - traj_pte[3,1,3,1,1]
p273_pte[8] <- traj_pte[3,1,1,3,3] - traj_pte[3,1,1,3,1]
p273_pte[9] <- traj_pte[1,3,3,1,3] - traj_pte[1,3,3,1,1]
p273_pte[10] <- traj_pte[1,3,1,3,3] - traj_pte[1,3,1,3,1]
p273_pte[11] <- traj_pte[1,1,3,3,3] - traj_pte[1,1,3,3,1]
p273_pte[12] <- traj_pte[3,3,3,1,3] - traj_pte[3,3,3,1,1]
p273_pte[13] <- traj_pte[3,3,1,3,3] - traj_pte[3,3,1,3,1]
p273_pte[14] <- traj_pte[3,1,3,3,3] - traj_pte[3,1,3,3,1]
p273_pte[15] <- traj_pte[1,3,3,3,3] - traj_pte[1,3,3,3,1]
p273_pte[16] <- traj_pte[3,3,3,3,3] - traj_pte[3,3,3,3,1]

write("Position 273 Effects", file="Pos273.csv")
write("genotype, dhc, px, pt, pom, pte", file="Pos273.csv", append=TRUE)
write(paste("aaaax", p273_dhc[1],p273_pxm[1],p273_ptm[1],p273_pxe[1],p273_pte[1], sep=","), file="Pos273.csv", append=TRUE)
write(paste("daaax", p273_dhc[2],p273_pxm[2],p273_ptm[2],p273_pxe[2],p273_pte[2], sep=","), file="Pos273.csv", append=TRUE)
write(paste("adaax", p273_dhc[3],p273_pxm[3],p273_ptm[3],p273_pxe[3],p273_pte[3], sep=","), file="Pos273.csv", append=TRUE)
write(paste("aadax", p273_dhc[4],p273_pxm[4],p273_ptm[4],p273_pxe[4],p273_pte[4], sep=","), file="Pos273.csv", append=TRUE)
write(paste("aaadx", p273_dhc[5],p273_pxm[5],p273_ptm[5],p273_pxe[5],p273_pte[5], sep=","), file="Pos273.csv", append=TRUE)
write(paste("ddaax", p273_dhc[6],p273_pxm[6],p273_ptm[6],p273_pxe[6],p273_pte[6], sep=","), file="Pos273.csv", append=TRUE)
write(paste("dadax", p273_dhc[7],p273_pxm[7],p273_ptm[7],p273_pxe[7],p273_pte[7], sep=","), file="Pos273.csv", append=TRUE)
write(paste("daadx", p273_dhc[8],p273_pxm[8],p273_ptm[8],p273_pxe[8],p273_pte[8], sep=","), file="Pos273.csv", append=TRUE)
write(paste("addax", p273_dhc[9],p273_pxm[9],p273_ptm[9],p273_pxe[9],p273_pte[9], sep=","), file="Pos273.csv", append=TRUE)
write(paste("adadx", p273_dhc[10],p273_pxm[10],p273_ptm[10],p273_pxe[10],p273_pte[10], sep=","), file="Pos273.csv", append=TRUE)
write(paste("aaddx", p273_dhc[11],p273_pxm[11],p273_ptm[11],p273_pxe[11],p273_pte[11], sep=","), file="Pos273.csv", append=TRUE)
write(paste("dddax", p273_dhc[12],p273_pxm[12],p273_ptm[12],p273_pxe[12],p273_pte[12], sep=","), file="Pos273.csv", append=TRUE)
write(paste("ddadx", p273_dhc[13],p273_pxm[13],p273_ptm[13],p273_pxe[13],p273_pte[13], sep=","), file="Pos273.csv", append=TRUE)
write(paste("daddx", p273_dhc[14],p273_pxm[14],p273_ptm[14],p273_pxe[14],p273_pte[14], sep=","), file="Pos273.csv", append=TRUE)
write(paste("adddx", p273_dhc[15],p273_pxm[15],p273_ptm[15],p273_pxe[15],p273_pte[15], sep=","), file="Pos273.csv", append=TRUE)
write(paste("ddddx", p273_dhc[16],p273_pxm[16],p273_ptm[16],p273_pxe[16],p273_pte[16], sep=","), file="Pos273.csv", append=TRUE)


pdf(file=paste("273_pxm.pdf"))

barplot(p273_pxm)

dev.off()

pdf(file=paste("273_ptm.pdf"))

barplot(p273_ptm)

dev.off()

pdf(file=paste("273_pxe.pdf"))

barplot(p273_pxe)

dev.off()

pdf(file=paste("273_pte.pdf"))

barplot(p273_pte)

dev.off()


pdf(file=paste("273_pxmvspxe", "pdf", sep="."))

pxPXE_273 <- lm(p273_pxm~p273_pxe)
pxPXE_r2 = summary(pxPXE_273)$r.squared
pxPXE_m = coefficients(pxPXE_273)[[2]]

plot(p273_pxe, p273_pxm)

abline(pxPXE_273)

mtext(paste("m =", round(pxPXE_m, 3)), side=3, line=-1, at=0, cex=0.8, adj=0)
mtext(paste("r2 =", round(pxPXE_r2, 3)), side=3, line=-2, at=0, cex=0.8, adj=0)
pxPXE_p = summary(pxPXE_273)$coefficients[,4]
mtext(paste("p =", round(pxPXE_p, 3)), side=3, line=-3, at=0, cex=0.8, adj=0)

dev.off()


pdf(file=paste("273_pxmvsptm", "pdf", sep="."))

pxPTM_273 <- lm(p273_pxm~p273_ptm)
pxPTM_r2 = summary(pxPTM_273)$r.squared
pxPTM_m = coefficients(pxPTM_273)[[2]]

plot(p273_pt, p273_px)

abline(pxPTM_273)

mtext(paste("m =", round(pxPTM_m, 3)), side=3, line=-1, at=0, cex=0.8, adj=0)
mtext(paste("r2 =", round(pxPTM_r2, 3)), side=3, line=-2, at=0, cex=0.8, adj=0)
pxPTM_p = summary(pxPTM_273)$coefficients[,4]
mtext(paste("p =", round(pxPTM_p, 3)), side=3, line=-3, at=0, cex=0.8, adj=0)

dev.off()


pdf(file=paste("273_pxmvspte", "pdf", sep="."))

pxpte_273 <- lm(p273_pxm~p273_pte)
pxpte_r2 = summary(pxpte_273)$r.squared
pxpte_m = coefficients(pxpte_273)[[2]]

plot(p273_pte, p273_pxm)

abline(pxpte_273)

mtext(paste("m =", round(pxpte_m, 3)), side=3, line=-1, at=0, cex=0.8, adj=0)
mtext(paste("r2 =", round(pxpte_r2, 3)), side=3, line=-2, at=0, cex=0.8, adj=0)
pxpte_p = summary(pxpte_273)$coefficients[,4]
mtext(paste("p =", round(pxpte_p, 3)), side=3, line=-3, at=0, cex=0.8, adj=0)

dev.off()


pdf(file=paste("273_ptmvspxe", "pdf", sep="."))

ptPXE_273 <- lm(p273_ptm~p273_pxe)
ptPXE_r2 = summary(ptPXE_273)$r.squared
ptPXE_m = coefficients(ptPXE_273)[[2]]

plot(p273_pxe, p273_ptm)

abline(ptPXE_273)

mtext(paste("m =", round(ptPXE_m, 3)), side=3, line=-1, at=0, cex=0.8, adj=0)
mtext(paste("r2 =", round(ptPXE_r2, 3)), side=3, line=-2, at=0, cex=0.8, adj=0)
ptPXE_p = summary(ptPXE_273)$coefficients[,4]
mtext(paste("p =", round(ptPXE_p, 3)), side=3, line=-3, at=0, cex=0.8, adj=0)

dev.off()


pdf(file=paste("273_ptmvspte", "pdf", sep="."))

ptpte_273 <- lm(p273_ptm~p273_pte)
ptpte_r2 = summary(ptpte_273)$r.squared
ptpte_m = coefficients(ptpte_273)[[2]]

plot(p273_pte, p273_ptm)

abline(ptpte_273)

mtext(paste("m =", round(ptpte_m, 3)), side=3, line=-1, at=0, cex=0.8, adj=0)
mtext(paste("r2 =", round(ptpte_r2, 3)), side=3, line=-2, at=0, cex=0.8, adj=0)
ptpte_p = summary(ptpte_273)$coefficients[,4]
mtext(paste("p =", round(ptpte_p, 3)), side=3, line=-3, at=0, cex=0.8, adj=0)

dev.off()

pdf(file=paste("273_pxevspte", "pdf", sep="."))

pompte_273 <- lm(p273_pxe~p273_pte)
pompte_r2 = summary(pompte_273)$r.squared
pompte_m = coefficients(pompte_273)[[2]]

plot(p273_pte, p273_pxe)

abline(pompte_273)

mtext(paste("m =", round(pompte_m, 3)), side=3, line=-1, at=0, cex=0.8, adj=0)
mtext(paste("r2 =", round(pompte_r2, 3)), side=3, line=-2, at=0, cex=0.8, adj=0)
pompte_p = summary(pompte_273)$coefficients[,4]
mtext(paste("p =", round(pompte_p, 3)), side=3, line=-3, at=0, cex=0.8, adj=0)

dev.off()


#dhc_all <- array(NA, c(80))

PXMall <- array(NA, c(80))

PTM_all <- array(NA, c(80))

PXE_all <- array(NA, c(80))

pte_all <- array(NA, c(80))


for (k in 1:16) {
    
    #  dhc_all[k] <- p72_dhc[k]
    
    PXM_all[k] <- p72_pxm[k]
    
    PTM_all[k] <- p72_ptm[k]
    
    PXE_all[k] <- p72_pxe[k]
    
    pte_all[k] <- p72_pte[k]
    
}

for (l in 1:16) {
    
    # dhc_all[l+16] <- p193_dhc[l]
    
    PXM_all[l+16] <- p193_pxm[l]
    
    PTM_all[l+16] <- p193_ptm[l]
    
    PXE_all[l+16] <- p193_pxe[l]
    
    pte_all[l+16] <- p193_pte[l]
    
}

for (m in 1:16) {
    
    #  dhc_all[m+32] <- p258_dhc[m]
    
    PXM_all[m+32] <- p258_pxm[m]
    
    PTM_all[m+32] <- p258_ptm[m]
    
    PXE_all[m+32] <- p258_pxe[m]
     
    pte_all[m+32] <- p258_pte[m]
    
}


for (n in 1:16) {
    
    #    dhc_all[n+48] <- p271_dhc[n]
    
    PXM_all[m+48] <- p271_pxm[n]
    
    PTM_all[m+48] <- p271_ptm[n]
    
    PXE_all[m+48] <- p271_pxe[n]
    
    pte_all[m+48] <- p271_pte[n]
    
}


for (o in 1:16) {
    
    #    dhc_all[o+64] <- p273_dhc[o]
    
    PXM_all[m+64] <- p273_pxm[o]
    
    PTM_all[m+64] <- p273_ptm[o]
    
    PXE_all[m+64] <- p273_pxe[o]
    
    pte_all[m+64] <- p273_pte[o]
    
}



pdf(file=paste("all_pxmvspxe", "pdf", sep="."))

pxPXE_all <- lm(PXMall~PXE_all)
pxPXE_r2 = summary(pxPXE_all)$r.squared
pxPXE_m = coefficients(pxPXE_all)[[2]]

plot(PXE_all, PXMall)

abline(pxPXE_all)

mtext(paste("m =", round(pxPXE_m, 3)), side=3, line=-1, at=0, cex=0.8, adj=0)
mtext(paste("r2 =", round(pxPXE_r2, 3)), side=3, line=-2, at=0, cex=0.8, adj=0)
pxPXE_p = summary(pxPXE_all)$coefficients[,4]
mtext(paste("p =", round(pxPXE_p, 3)), side=3, line=-3, at=0, cex=0.8, adj=0)

dev.off()


pdf(file=paste("all_pxmvsptm", "pdf", sep="."))

pxPTM_all <- lm(PXMall~PTM_all)
pxPTM_r2 = summary(pxPTM_all)$r.squared
pxPTM_m = coefficients(pxPTM_all)[[2]]

plot(PTM_all, PXMall)

abline(pxPTM_all)

mtext(paste("m =", round(pxPTM_m, 3)), side=3, line=-1, at=0, cex=0.8, adj=0)
mtext(paste("r2 =", round(pxPTM_r2, 3)), side=3, line=-2, at=0, cex=0.8, adj=0)
pxPTM_p = summary(pxPTM_all)$coefficients[,4]
mtext(paste("p =", round(pxPTM_p, 3)), side=3, line=-3, at=0, cex=0.8, adj=0)

dev.off()


pdf(file=paste("all_pxmvspte", "pdf", sep="."))

pxpte_all <- lm(PXMall~pte_all)
pxpte_r2 = summary(pxpte_all)$r.squared
pxpte_m = coefficients(pxpte_all)[[2]]

plot(pte_all, PXMall)

abline(pxpte_all)

mtext(paste("m =", round(pxpte_m, 3)), side=3, line=-1, at=0, cex=0.8, adj=0)
mtext(paste("r2 =", round(pxpte_r2, 3)), side=3, line=-2, at=0, cex=0.8, adj=0)
pxpte_p = summary(pxpte_all)$coefficients[,4]
mtext(paste("p =", round(pxpte_p, 3)), side=3, line=-3, at=0, cex=0.8, adj=0)

dev.off()


pdf(file=paste("all_ptmvspxe", "pdf", sep="."))

ptPXE_all <- lm(PTM_all~PXE_all)
ptPXE_r2 = summary(ptPXE_all)$r.squared
ptPXE_m = coefficients(ptPXE_all)[[2]]

plot(PXE_all, PTM_all)

abline(ptPXE_all)

mtext(paste("m =", round(ptPXE_m, 3)), side=3, line=-1, at=0, cex=0.8, adj=0)
mtext(paste("r2 =", round(ptPXE_r2, 3)), side=3, line=-2, at=0, cex=0.8, adj=0)
ptPXE_p = summary(ptPXE_all)$coefficients[,4]
mtext(paste("p =", round(ptPXE_p, 3)), side=3, line=-3, at=0, cex=0.8, adj=0)

dev.off()


pdf(file=paste("all_ptmvspte", "pdf", sep="."))

ptpte_all <- lm(PTM_all~pte_all)
ptpte_r2 = summary(ptpte_all)$r.squared
ptpte_m = coefficients(ptpte_all)[[2]]

plot(pte_all, PTM_all)

abline(ptpte_all)

mtext(paste("m =", round(ptpte_m, 3)), side=3, line=-1, at=0, cex=0.8, adj=0)
mtext(paste("r2 =", round(ptpte_r2, 3)), side=3, line=-2, at=0, cex=0.8, adj=0)
ptpte_p = summary(ptpte_all)$coefficients[,4]
mtext(paste("p =", round(ptpte_p, 3)), side=3, line=-3, at=0, cex=0.8, adj=0)

dev.off()

pdf(file=paste("all_pxevspte", "pdf", sep="."))

pompte_all <- lm(PXE_all~pte_all)
pompte_r2 = summary(pompte_all)$r.squared
pompte_m = coefficients(pompte_all)[[2]]

plot(pte_all, PXE_all)

abline(pompte_all)

mtext(paste("m =", round(pompte_m, 3)), side=3, line=-1, at=0, cex=0.8, adj=0)
mtext(paste("r2 =", round(pompte_r2, 3)), side=3, line=-2, at=0, cex=0.8, adj=0)
pompte_p = summary(pompte_all)$coefficients[,4]
mtext(paste("p =", round(pompte_p, 3)), side=3, line=-3, at=0, cex=0.8, adj=0)

dev.off()





