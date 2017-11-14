import re

# Read the file
f = open('../Data/blackbirds.txt', 'r')
text = f.read()
f.close()

# remove \t\n and put a space in:
text = text.replace('\t',' ')
text = text.replace('\n',' ')

# note that there are "strange characters" (these are accents and
# non-ascii symbols) because we don't care for them, first transform
# to ASCII:
text = text.decode('ascii', 'ignore')

# Now write a regular expression my_reg that captures # the Kingdom, 
# Phylum and Species name for each species and prints it out neatly:
my_reg = r'(Kingdom\s\w*).*?(Phylum\s\w*).*?(Species\s[\w\s]*)'
result = re.findall(my_reg, text)

# Print result out nicely
for i in range(len(result)):
    print (str(result[i][0]))
    print (str(result[i][1]))
    print (str(result[i][2]))
    print("\n")
