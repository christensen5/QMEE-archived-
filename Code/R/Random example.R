# written by james rosindell james@rosindell.org Imperial college london released open source under an MIT license

# Assign random speciation rates to class
CMEE_2014 <- c(
"Fatima Alafifi",
"Louise Archer",
"Dominic Barker",
"Holly Black",
"Sarah Bunney",
"Vignesh Chundru",
"Paul Cronin",
"Fiona Deacon",
"Leonie Gough",
"Robert Hohan",
"Sean Jordan",
"Michael Massam",
"Dermott McMorrough",
"Samuel North",
"Matishalin Patel",
"Sharmila Rana",
"Thomas Smallwood",
"Thomas Smith",
"Matthew Speight",
"Samuel Thompson",
"Leah Trigg",
"Benjamin Watkinson-Powell")

choose_student <- function(class) {
    print(sample(class,1))
}

choose_student_2 <- function(class,seedin = 1) {
    set.seed(seedin)
    print(sample(class,1))
}

choose_student_3 <- function(class,seedin=-1) {
    if (seedin <= 0){
        set.seed(floor(proc.time()[3]*1000))
    }
    else {
        set.seed(seedin)
    }
    print(sample(class,1))
}

assign_student_number_1 <- function(class,min,max,seedin=-1) {
    
    if (seedin <= 0){
        set.seed(floor(proc.time()[3]*1000))
    }
    else {
        set.seed(seedin)
    }
    number <- runif(length(class))*(max-min)+min
    
    return(cbind(number,class))
}

assign_student_number_2 <- function(class,min,max,sigfig=4,seedin=-1) {
    if (seedin <= 0){
        set.seed(floor(proc.time()[3]*1000))
    }
    else {
        set.seed(seedin)
    }
    number <- signif(runif(length(class))*(max-min)+min,sigfig)
   
    return(cbind(number,class))
}

assign_student_number_3 <- function(class,seedin=-1,min=0.002,max=0.007,sigfig=4,unique=TRUE) {
    if (seedin <= 0){
        set.seed(floor(proc.time()[3]*1000))
    }
    else {
        set.seed(seedin)
    }
    speciation_values <- signif(runif(length(class))*(max-min)+min,sigfig)
    if (unique){
        while(length(unique(speciation_values)) < length(class)){
            speciation_values <- signif(runif(length(class))*(max-min)+min,sigfig)
        }
    }
    return(cbind(speciation_values,class))
}
