# written by james rosindell james@rosindell.org Imperial college london released open source under an MIT license
#clear namespace
rm(list=ls())
graphics.off()

for (i in 1:12)
{
    print(paste("On day",i," of christmas my true love sent to me"))
    for (j in i:1) print (paste (j," ",c("Drummers Drumming","Pipers Piping","Lords A-leaping","Ladies Dancing","Maids A-milking","Swans A-swimming","Geese A-laying","Golden Rings","Calling Birds","French Hens","Turtle Doves","Partridge in a pear tree")[13-j]))
}
