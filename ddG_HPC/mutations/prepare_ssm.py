import os
import glob

#sequence = "AKGDPHVLLTTSAGNIELELDKQKAPVSVQNFVDYVNSGFYNNTTFHRVIPGFMIQGGGFTEQMQQKKPNPPIKNEADNGLRNTRGTIAMARTADKDSATSQFFINVADNAFLDHGQRDFGYAVFGKVVKGMDVADKISQVPTHDVGPYQNVPSKPVVILSATVLP"
mut=['L1', 'Y2', 'I3', 'Q4', 'W5', 'L6', 'K7', 'D8']
aa_list = ['M', 'F','G', 'A', 'I', 'P', 'S']

mutations = []

for aa in range(0, len(mut)):
	pos = mut[aa]
	AA = mut[aa][0]
#	pos = AA + str(aa+1)
	if not os.path.exists(pos):
		os.makedirs(pos)
	for m in aa_list:
		if m != AA:
			mutation = pos+m
			mutations.append(mutation)
			if not os.path.exists(pos+"/"+mutation):
				os.makedirs(pos+"/"+mutation)
			with open(pos+"/"+mutation+"/mutfile",'w') as mut_file:
				mut_file.write("total 1\n")
				mut_file.write("1\n")
				mut_file.write("%s %d %s\n" %(AA, int(mut[aa][1:]), m))
with open("task.list", "w") as task_list:
	for i in mutations:
		task_list.write(i+"\n")
