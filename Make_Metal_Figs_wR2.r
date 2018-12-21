data <- read.csv("PTM_Terms.csv")	#reads in the data file - name inside the quotation marks needs to match the filename. First column must be time.

attach(data)							# Assigns the first line - i.e. the head of all the columns - to be the variable name for all the vlues in those columns

var_names <- colnames(data)
var_names <- var_names[2:length(var_names)] # This saves all the variable names, excluding the first column, which should be "terms"

data2 <- read.csv("PTM_R2_Terms.csv")

attach(data2)

n = length(Terms)

colours <- c("blue", "blue", "blue", "blue", "blue", "green", "green", "green", "green", "green", "green", "green", "green", "green", "green", "red", "red", "red", "red", "red", "red", "red", "red", "red", "red", "orange", "orange", "orange", "orange", "orange", "purple")

for(i in var_names) {					# Here we start the "loop" that will run through each column, and therefore analyze all of the datasets in the file separately
	x1 <- as.numeric(data[[i]])
    
    width_var <- paste(i,"R2", sep="")
    
    width <- as.numeric(data2[[width_var]])
    
    # print(width)
  
    x1 <- as.numeric(data[[i]])
    
    width_var <- paste(i,"R2",sep="")
    
    width <- as.numeric(data2[[width_var]])
  
    x1_main <- x1[1:5]
    
    width_main <- width[1:5]

    x1_second <- x1[6:15]

    width_second <- c(width[7]/10,width[7]/10,width[7]/10,width[7]/10,width[7]/10,width[7]/10,width[7]/10,width[7]/10,width[7]/10,width[7]/10)

    x1_third <- x1[16:25]

    width_third <- c(width[8]/10,width[8]/10,width[8]/10,width[8]/10,width[8]/10,width[8]/10,width[8]/10,width[8]/10,width[8]/10,width[8]/10)

    x1_fourth <- x1[26:30]

    width_fourth <- c(width[9]/5,width[9]/5,width[9]/5,width[9]/5,width[9]/5)

    x1_fifth <- x1[31]

    width_fifth <- width[10]

    width_fig <- c(width_main, width_second, width_third, width_fourth, width_fifth)
    
    print(width_fig)
    
    pdf(file=paste(i,"Effects","pdf", sep="."))	# This begins the output that will go into a pdf figure, which will be saved when we are done.
    
    op <- par(mar = c(10,4,4,2) + 0.1)
    
    barplot(x1, width_fig, col=colours, names.arg=Terms, space=0, las=2, ylab="Effect on activity (log(nM/s))", main=paste("Genetic Effects on", i, sep=" "))
    par(op)
    
    dev.off()
    
    pdf(file=paste(i,"Main Effects","pdf", sep="."))
    
    op <- par(mar = c(10,4,4,2) + 0.1)
    
    plottitle = paste(i,"main effects", sep=" ")
    
    barplot(x1_main, width_main, col=colours[1:5], names.arg=Terms[1:5], las=2, space=0, ylab="Effect on activity (log(nM/s)", main=plottitle)
    par(op)
    
    dev.off()
    
    pdf(file=paste(i,"Second Order","pdf", sep="."))
    
    op <- par(mar = c(10,4,4,2) + 0.1)
    
    plottitle = paste("Second Order Epistatic Effects")
    
    barplot(x1_second, width_second, col=colours[6:15], names.arg=Terms[6:15], las=2, space=0, ylab="Effect on activity (log(nM/s)", main=plottitle)
    par(op)
    
    dev.off()
    
    pdf(file=paste(i,"Third Order Effects","pdf", sep="."))
    
    op <- par(mar = c(10,4,4,2) + 0.1)
    
    plottitle = paste("Third Order Effects")
    
    barplot(x1_third, width_third, col=colours[16:25], names.arg=Terms[16:25], las=2, space=0, ylab="Effect on activity (log(nM/s)", main=plottitle)
    par(op)
    
    dev.off()
    
    pdf(file=paste(i,"Fourth Order Effects","pdf", sep="."))
    
    op <- par(mar = c(10,4,4,2) + 0.1)
    
    plottitle = paste("Fourth Order Effects")
    
    barplot(x1_fourth, width_fourth, col=colours[26:30], names.arg=Terms[26:30], las=2, space=0, ylab="Effect on activity (log(nM/s)", main=plottitle)
    par(op)
    
    dev.off()
    
    pdf(file=paste(i,"Fifth Order Effects","pdf", sep="."))
    
    op <- par(mar = c(10,4,4,2) + 0.1)
    
    plottitle = paste("Fifth Order Effects", sep=" ")
    
    barplot(x1_fifth, width_fifth, col=colours[31], names.arg=Terms[31], space=0, ylab="Effect on activity (log(nM/s)", main=plottitle)
    par(op)
    
    dev.off()
    
 
}

