import sys

#band = 7
#sex = 6 
#pair info = 18
#mate id = 20

out = open(sys.argv[2],'w')

header = []
data = []

for line in open(sys.argv[1], 'r'):
	dat = line.rstrip()
	if dat.split()[0] == "date":
		header.append(dat)
	else:
		data.append(dat)
	
for h in header:
	out.write(str(h)+'\t'+str(h)+'\n')

mates = []

with open(sys.argv[1],'r') as bands:
	next(bands)
	for line in bands:
		if str(line.split('\t')[18]) != "NA":
			band = line.split('\t')[7]
			mate = line.split('\t')[20]
			for dat in data:
				if dat.split('\t')[7] == mate:
					mate_dat = str(dat)
			mates.append(line.rstrip()+'\t'+mate_dat.rstrip())

# Make sure the exact pair (unless from two years) is NOT in the list already.
# Checks that first mate in pair is male, which accounts for doubles.

mate1_list = []
mate2_list = []

for m in mates:
	mate1 = m.split('\t')[7]
	mate1_sex = m.split('\t')[6]
	mate2 = m.split('\t')[79]
	mate2_sex = m.split('\t')[78]
	if str(mate1_sex) == "m" or str(mate1_sex) =="m?":
		out.write(str(m)+'\n')

# OLDER DOUBLE CHECK (KEEPING JUST IN CASE)	
# 		mate1_list.append(mate1)
# 	else:
# 		mate2_list.append(mate2)
# 	
# 	if mate1 not in mate2_list:
# 		if mate2 not in mate1_list:
# 			out.write(str(m)+'\n')
