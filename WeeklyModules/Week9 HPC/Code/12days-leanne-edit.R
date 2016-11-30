# written by james rosindell james@rosindell.org Imperial college london released open source under an MIT license

#clear namespace
rm(list=ls())
graphics.off()

# vector of presents
presents_classic <- c("Windows computer crashing","Repositories","Shell scripts","Databases","Processors","Gits a pushing","Samraats singing","Laptops running","Printers printing","Errors flashing","Coding CMEEs","Months of Coding")

# define function
twelvedays <- function(presents=presents_classic) {
    # to store the total
    total1 <- 0
    num_days <- length(presents)
    # loop over all days
    for (days in 1:num_days){
        print("                                    ")
        part1 <- paste ("On the",days)
        part2 <- ""
        # an if statement to work out if its 'th' 'nd' or 'rd' to add on the end of the number
        if ((days%%100 > 10)&&(days%%100 <20)){
            part2 <- "th"
        } else if (days%%10 == 1) {
            part2 <- "st"
        } else if (days%%10 == 2) {
            part2 <- "nd"
        } else if (days%%10 == 3) {
            part2 <- "rd"
        } else {
            part2 <- "th"
        }
        part3 <- " day of CODE-mas my Samraat git-pushed to me"
        # print first line of song
        print(paste(part1,part2,part3,sep = ""))
        # loop over all presents starting high working down
        for (pres_i in days:1) {
            total1 <- total1 + pres_i # add to total
            if (pres_i ==1 ) {
                # we need to say a not a number
                if (days != 1) {
                    # we're giving a list so need to say and
                    print (paste ("and a",presents[pres_i]))
                } else {
                    print (paste ("a",presents[pres_i]))
                }
            } else {
                if (pres_i == 2) {
                    # send to last entry followed by and so no comma needed
                    print (paste (pres_i," ",presents[pres_i],sep=""))
                }
                else {
                    # comma needed for proper list
                    print (paste (pres_i," ",presents[pres_i],",",sep=""))
                }
            }
        }
    }
    print("                                    ")
    print("total number of bugs is:")
    print (paste (total1,"by summing"))
    # check the total from formula (x+2)(x+1)x/6
    print (paste ((num_days+2)*(num_days+1)*(num_days)/6,"by formula"))
}

twelvedays()