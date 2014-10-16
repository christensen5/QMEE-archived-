import re

# Read the file
f = open('../Data/blackbirds.txt', 'r')
text = f.read()
f.close()

# remove \t\n and put a space in
text = text.replace('\t',' ')
text = text.replace('\n',' ')

# note that there are "strange characters" (these are accents and
# non-ascii symbols) because we don't care for them, let's transform
# to ASCII
text = text.decode('ascii', 'ignore')

## Now write a regular expression that captures 
## the Kingdom, Phylum and Species name for each species

my_reg = ??????

# such that re.findall(my_reg, text) should  return
#[(u'Animalia', u'Chordata', u'Euphagus carolinus'),
# (u'Animalia', u'Chordata', u'Euphagus cyanocephalus'),
# (u'Animalia', u'Chordata', u'Turdus boulboul'),
# (u'Animalia', u'Chordata', u'Agelaius assimilis')]
