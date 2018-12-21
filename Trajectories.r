plot_data <- read.csv("Coordinates.csv")

attach(plot_data)

#write("Avgs", file="Average_values.csv")
#write(xcoord, file="Average_values.csv", append = TRUE)

# Now calculating the trajectory that takes the most positive steps, regardless of genotype or environmental change

traj_array_R1 <- array(NA, dim=c(3,3,3,3,3,8))
traj_array_R2 <- array(NA, dim=c(3,3,3,3,3,8))
traj_array_R3 <- array(NA, dim=c(3,3,3,3,3,8))

for (i in 1:256) {
        traj_array_R1[pos72[i]+2, pos193[i]+2, pos258[i]+2, pos271[i]+2, pos273[i]+2, metal[i]] = R1[i]
        traj_array_R2[pos72[i]+2, pos193[i]+2, pos258[i]+2, pos271[i]+2, pos273[i]+2, metal[i]] = R2[i]
        traj_array_R3[pos72[i]+2, pos193[i]+2, pos258[i]+2, pos271[i]+2, pos273[i]+2, metal[i]] = R3[i]
}



# write("pos72, pos193, pos258, pos271, pos273, metal", file="TrajFile_Ni.txt")

# write("1, 1, 1, 1, 1, 8", file="TrajFile_Ni.txt", append=TRUE)

# Computing Trajectory

for (h in 1:8) {
    endpoints <- array(0, dim=c(3,3,3,3,3))
    write("pos72, pos193, pos258, pos271, pos273", file=paste("TrajFile", h, "txt", sep="."))


    for (i in c(1,3)) {
        for (j in c(1,3)) {
            for (k in c(1,3)) {
                for (l in c(1,3)) {
                    for (m in c(1,3)) {
                
                        point <- c(i, j, k, l, m, h)
                        
                        write(c(point[1], point[2], point[3], point[4], point[5]), file=paste("TrajFile", h, "txt", sep="."), sep=",", ncol=6, append = TRUE)
                    
                        if (point[1]==1) {
                            test_72 <- c(point[1]+2, point[2], point[3], point[4], point[5], point[6])
                        }
                        if (point[1]==3) {
                            test_72 <- c(point[1]-2, point[2], point[3], point[4], point[5], point[6])
                        }
                            
                        diff_72_1 <- traj_array_R1[test_72[1], test_72[2], test_72[3], test_72[4], test_72[5], point[6]] - traj_array_R1[point[1], point[2], point[3], point[4], point[5], point[6]]
                        diff_72_2 <- traj_array_R1[test_72[1], test_72[2], test_72[3], test_72[4], test_72[5], point[6]] - traj_array_R2[point[1], point[2], point[3], point[4], point[5], point[6]]
                        diff_72_3 <- traj_array_R1[test_72[1], test_72[2], test_72[3], test_72[4], test_72[5], point[6]] - traj_array_R3[point[1], point[2], point[3], point[4], point[5], point[6]]
                        diff_72_4 <- traj_array_R2[test_72[1], test_72[2], test_72[3], test_72[4], test_72[5], point[6]] - traj_array_R1[point[1], point[2], point[3], point[4], point[5], point[6]]
                        diff_72_5 <- traj_array_R2[test_72[1], test_72[2], test_72[3], test_72[4], test_72[5], point[6]] - traj_array_R2[point[1], point[2], point[3], point[4], point[5], point[6]]
                        diff_72_6 <- traj_array_R2[test_72[1], test_72[2], test_72[3], test_72[4], test_72[5], point[6]] - traj_array_R3[point[1], point[2], point[3], point[4], point[5], point[6]]
                        diff_72_7 <- traj_array_R3[test_72[1], test_72[2], test_72[3], test_72[4], test_72[5], point[6]] - traj_array_R1[point[1], point[2], point[3], point[4], point[5], point[6]]
                        diff_72_8 <- traj_array_R3[test_72[1], test_72[2], test_72[3], test_72[4], test_72[5], point[6]] - traj_array_R2[point[1], point[2], point[3], point[4], point[5], point[6]]
                        diff_72_9 <- traj_array_R3[test_72[1], test_72[2], test_72[3], test_72[4], test_72[5], point[6]] - traj_array_R3[point[1], point[2], point[3], point[4], point[5], point[6]]
                            
                        passable <- sum(c(diff_72_1, diff_72_2, diff_72_3, diff_72_4, diff_72_5, diff_72_6, diff_72_7, diff_72_8, diff_72_9)>0)
                            
                        assess_72 = passable/9
                    
                        write(paste(c("Position 72:", assess_72),sep=","), file=paste("TrajFile", h, "txt", sep="."), sep=",", ncol=2, append=TRUE)
                    
                    
                        if (point[2]==1) {
                            test_193 <- c(point[1], point[2]+2, point[3], point[4], point[5], point[6])
                        }
                        if (point[2]==3) {
                            test_193 <- c(point[1], point[2]-2, point[3], point[4], point[5], point[6])
                        }
                    
                        diff_193_1 <- traj_array_R1[test_193[1], test_193[2], test_193[3], test_193[4], test_193[5], point[6]] - traj_array_R1[point[1], point[2], point[3], point[4], point[5], point[6]]
                        diff_193_2 <- traj_array_R1[test_193[1], test_193[2], test_193[3], test_193[4], test_193[5], point[6]] - traj_array_R2[point[1], point[2], point[3], point[4], point[5], point[6]]
                        diff_193_3 <- traj_array_R1[test_193[1], test_193[2], test_193[3], test_193[4], test_193[5], point[6]] - traj_array_R3[point[1], point[2], point[3], point[4], point[5], point[6]]
                        diff_193_4 <- traj_array_R2[test_193[1], test_193[2], test_193[3], test_193[4], test_193[5], point[6]] - traj_array_R1[point[1], point[2], point[3], point[4], point[5], point[6]]
                        diff_193_5 <- traj_array_R2[test_193[1], test_193[2], test_193[3], test_193[4], test_193[5], point[6]] - traj_array_R2[point[1], point[2], point[3], point[4], point[5], point[6]]
                        diff_193_6 <- traj_array_R2[test_193[1], test_193[2], test_193[3], test_193[4], test_193[5], point[6]] - traj_array_R3[point[1], point[2], point[3], point[4], point[5], point[6]]
                        diff_193_7 <- traj_array_R3[test_193[1], test_193[2], test_193[3], test_193[4], test_193[5], point[6]] - traj_array_R1[point[1], point[2], point[3], point[4], point[5], point[6]]
                        diff_193_8 <- traj_array_R3[test_193[1], test_193[2], test_193[3], test_193[4], test_193[5], point[6]] - traj_array_R2[point[1], point[2], point[3], point[4], point[5], point[6]]
                        diff_193_9 <- traj_array_R3[test_193[1], test_193[2], test_193[3], test_193[4], test_193[5], point[6]] - traj_array_R3[point[1], point[2], point[3], point[4], point[5], point[6]]
                    
                        passable <- sum(c(diff_193_1, diff_193_2, diff_193_3, diff_193_4, diff_193_5, diff_193_6, diff_193_7, diff_193_8, diff_193_9)>0)
                    
                        assess_193 = passable/9
                    
                        write(paste(c("Position 193:", assess_193),sep=","), file=paste("TrajFile", h, "txt", sep="."), sep=",", ncol=2, append=TRUE)
                    
                    
                    
                        if (point[3]==1) {
                            test_258 <- c(point[1], point[2], point[3]+2, point[4], point[5], point[6])
                        }
                        if (point[3]==3) {
                            test_258 <- c(point[1], point[2], point[3]-2, point[4], point[5], point[6])
                        }
                    
                        diff_258_1 <- traj_array_R1[test_258[1], test_258[2], test_258[3], test_258[4], test_258[5], point[6]] - traj_array_R1[point[1], point[2], point[3], point[4], point[5], point[6]]
                        diff_258_2 <- traj_array_R1[test_258[1], test_258[2], test_258[3], test_258[4], test_258[5], point[6]] - traj_array_R2[point[1], point[2], point[3], point[4], point[5], point[6]]
                        diff_258_3 <- traj_array_R1[test_258[1], test_258[2], test_258[3], test_258[4], test_258[5], point[6]] - traj_array_R3[point[1], point[2], point[3], point[4], point[5], point[6]]
                        diff_258_4 <- traj_array_R2[test_258[1], test_258[2], test_258[3], test_258[4], test_258[5], point[6]] - traj_array_R1[point[1], point[2], point[3], point[4], point[5], point[6]]
                        diff_258_5 <- traj_array_R2[test_258[1], test_258[2], test_258[3], test_258[4], test_258[5], point[6]] - traj_array_R2[point[1], point[2], point[3], point[4], point[5], point[6]]
                        diff_258_6 <- traj_array_R2[test_258[1], test_258[2], test_258[3], test_258[4], test_258[5], point[6]] - traj_array_R3[point[1], point[2], point[3], point[4], point[5], point[6]]
                        diff_258_7 <- traj_array_R3[test_258[1], test_258[2], test_258[3], test_258[4], test_258[5], point[6]] - traj_array_R1[point[1], point[2], point[3], point[4], point[5], point[6]]
                        diff_258_8 <- traj_array_R3[test_258[1], test_258[2], test_258[3], test_258[4], test_258[5], point[6]] - traj_array_R2[point[1], point[2], point[3], point[4], point[5], point[6]]
                        diff_258_9 <- traj_array_R3[test_258[1], test_258[2], test_258[3], test_258[4], test_258[5], point[6]] - traj_array_R3[point[1], point[2], point[3], point[4], point[5], point[6]]
                    
                        passable <- sum(c(diff_258_1, diff_258_2, diff_258_3, diff_258_4, diff_258_5, diff_258_6, diff_258_7, diff_258_8, diff_258_9)>0)
                    
                        assess_258 = passable/9
                    
                        write(paste(c("Position 258:", assess_258),sep=","), file=paste("TrajFile", h, "txt", sep="."), sep=",", ncol=2, append=TRUE)
                    
                    
                    
                        if (point[4]==1) {
                            test_271 <- c(point[1], point[2], point[3], point[4]+2, point[5], point[6])
                        }
                        if (point[4]==3) {
                            test_271 <- c(point[1], point[2], point[3], point[4]-2, point[5], point[6])
                        }
                    
                        diff_271_1 <- traj_array_R1[test_271[1], test_271[2], test_271[3], test_271[4], test_271[5], point[6]] - traj_array_R1[point[1], point[2], point[3], point[4], point[5], point[6]]
                        diff_271_2 <- traj_array_R1[test_271[1], test_271[2], test_271[3], test_271[4], test_271[5], point[6]] - traj_array_R2[point[1], point[2], point[3], point[4], point[5], point[6]]
                        diff_271_3 <- traj_array_R1[test_271[1], test_271[2], test_271[3], test_271[4], test_271[5], point[6]] - traj_array_R3[point[1], point[2], point[3], point[4], point[5], point[6]]
                        diff_271_4 <- traj_array_R2[test_271[1], test_271[2], test_271[3], test_271[4], test_271[5], point[6]] - traj_array_R1[point[1], point[2], point[3], point[4], point[5], point[6]]
                        diff_271_5 <- traj_array_R2[test_271[1], test_271[2], test_271[3], test_271[4], test_271[5], point[6]] - traj_array_R2[point[1], point[2], point[3], point[4], point[5], point[6]]
                        diff_271_6 <- traj_array_R2[test_271[1], test_271[2], test_271[3], test_271[4], test_271[5], point[6]] - traj_array_R3[point[1], point[2], point[3], point[4], point[5], point[6]]
                        diff_271_7 <- traj_array_R3[test_271[1], test_271[2], test_271[3], test_271[4], test_271[5], point[6]] - traj_array_R1[point[1], point[2], point[3], point[4], point[5], point[6]]
                        diff_271_8 <- traj_array_R3[test_271[1], test_271[2], test_271[3], test_271[4], test_271[5], point[6]] - traj_array_R2[point[1], point[2], point[3], point[4], point[5], point[6]]
                        diff_271_9 <- traj_array_R3[test_271[1], test_271[2], test_271[3], test_271[4], test_271[5], point[6]] - traj_array_R3[point[1], point[2], point[3], point[4], point[5], point[6]]
                    
                        passable <- sum(c(diff_271_1, diff_271_2, diff_271_3, diff_271_4, diff_271_5, diff_271_6, diff_271_7, diff_271_8, diff_271_9)>0)
                    
                        assess_271 = passable/9
                    
                        write(paste(c("Position 271:", assess_271),sep=","), file=paste("TrajFile", h, "txt", sep="."), sep=",", ncol=2, append=TRUE)
                    
                    
                    
                        if (point[5]==1) {
                            test_273 <- c(point[1], point[2], point[3], point[4], point[5]+2, point[6])
                        }
                        if (point[5]==3) {
                            test_273 <- c(point[1], point[2], point[3], point[4], point[5]-2, point[6])
                        }
                    
                        diff_273_1 <- traj_array_R1[test_273[1], test_273[2], test_273[3], test_273[4], test_273[5], point[6]] - traj_array_R1[point[1], point[2], point[3], point[4], point[5], point[6]]
                        diff_273_2 <- traj_array_R1[test_273[1], test_273[2], test_273[3], test_273[4], test_273[5], point[6]] - traj_array_R2[point[1], point[2], point[3], point[4], point[5], point[6]]
                        diff_273_3 <- traj_array_R1[test_273[1], test_273[2], test_273[3], test_273[4], test_273[5], point[6]] - traj_array_R3[point[1], point[2], point[3], point[4], point[5], point[6]]
                        diff_273_4 <- traj_array_R2[test_273[1], test_273[2], test_273[3], test_273[4], test_273[5], point[6]] - traj_array_R1[point[1], point[2], point[3], point[4], point[5], point[6]]
                        diff_273_5 <- traj_array_R2[test_273[1], test_273[2], test_273[3], test_273[4], test_273[5], point[6]] - traj_array_R2[point[1], point[2], point[3], point[4], point[5], point[6]]
                        diff_273_6 <- traj_array_R2[test_273[1], test_273[2], test_273[3], test_273[4], test_273[5], point[6]] - traj_array_R3[point[1], point[2], point[3], point[4], point[5], point[6]]
                        diff_273_7 <- traj_array_R3[test_273[1], test_273[2], test_273[3], test_273[4], test_273[5], point[6]] - traj_array_R1[point[1], point[2], point[3], point[4], point[5], point[6]]
                        diff_273_8 <- traj_array_R3[test_273[1], test_273[2], test_273[3], test_273[4], test_273[5], point[6]] - traj_array_R2[point[1], point[2], point[3], point[4], point[5], point[6]]
                        diff_273_9 <- traj_array_R3[test_273[1], test_273[2], test_273[3], test_273[4], test_273[5], point[6]] - traj_array_R3[point[1], point[2], point[3], point[4], point[5], point[6]]
                    
                        passable <- sum(c(diff_273_1, diff_273_2, diff_273_3, diff_273_4, diff_273_5, diff_273_6, diff_273_7, diff_273_8, diff_273_9)>0)
                    
                        assess_273 = passable/9
                    
                        write(paste(c("Position 273:", assess_273),sep=","), file=paste("TrajFile", h, "txt", sep="."), sep=",", ncol=2, append=TRUE)
                    
                    
                    
                    
                    }
                }
            }
        }
    }
}