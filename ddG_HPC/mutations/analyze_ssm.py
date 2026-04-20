import os
import glob
import numpy as np

#sequence = "MVTFHTNHGDIVIKTFDDKAPETVKNFLDYCREGFYNNTIFHRVINGFMIQGGGFEPGMKQKATKEPIKNEANNGLKNTRGTLAMARTQAPHSATAQFFINVVDNDFLNFSGESLQGWGYCVFAEVVDGMDEVDKIKGVATGRSGMHQDVPKEDVIIESVTVSE"
#mut=["L118","L120","L2","L4","L60","L62","V119","V121","V122","V124","V3","V5","V6","V61","V63","V64","V66","V8"]
mut=["C101","C114","C484","C572"]
#aa_list = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']
aa_list = ["G","P","E","D","R","K","H","Q","N","T","S","Y","W","F","M","C","I","L","V","A"]

ala_scan = {}
ssm = np.zeros([20, len(mut)], dtype=float)

for aa in range(0, len(mut)):
	pos = mut[aa]
	AA = mut[aa][0]
	if AA != "A":
		with open(pos+"/"+pos+"A/mutfile.ddg", 'r') as mutfile:
			n_WT = 0
			n_MUT = 0
			score_WT = 0
			score_MUT = 0
			for line in mutfile:
				if "WT_" in line:
					score = float(line.split()[3])
					n_WT += 1
					score_WT += score
				elif "MUT_" in line:
					score = float(line.split()[3])
					n_MUT += 1
					score_MUT += score
			score_WT = score_WT/n_WT
			score_MUT = score_MUT/n_MUT
			ddG = score_MUT - score_WT
			ala_scan[mut[aa][1:]] = ddG
	else:
		ddG = 0
		ala_scan[mut[aa][1:]] = ddG

	for i in range(0, len(aa_list)):
		amino_acid = aa_list[i]
		if AA != amino_acid:
			mutation = mut[aa]+amino_acid
			n_WT = 0
			n_MUT = 0
			score_WT = 0
			score_MUT = 0
			with open(pos+"/"+mutation+"/mutfile.ddg", 'r') as mutfile:
				for line in mutfile:
					if "WT_" in line:
						score = float(line.split()[3])
						n_WT += 1
						score_WT += score
					elif "MUT_" in line:
						score = float(line.split()[3])
						n_MUT += 1
						score_MUT += score
			score_WT = score_WT/n_WT
			score_MUT = score_MUT/n_MUT
			ddG = score_MUT - score_WT
			ssm[i,aa] = float(ddG)
print(ssm)
print(ala_scan)
with open("TprC_N_AF2.pdb", "r") as in_pdb:
	with open("CB001_Ala_scan.pdb", "w") as out_pdb:
		for line in in_pdb:
			if line.startswith("ATOM "):
				vals = line.split()
				if vals[5] in ala_scan.values():
					ddG = ala_scan[vals[5]]
				else:
					ddG = 0.0
				line_ddg = line[0:61] + '{0: <8}'.format(str("{:.2f}".format(ddG))) + line[66:]
				print(line_ddg)
				out_pdb.write(line_ddg)


np.savetxt("SSM_ddg.csv", ssm, delimiter=",")
