from glob import glob
import os

accession_to_organism = {}
with open("accession_to_info.tsv") as INFO:
    next(INFO)
    for line in INFO:
        line = line.strip().split('\t')
        sci_name = line[3]
        acc = line[0].split('.')[0]
        sci_name = sci_name.replace('-', '_')
        sci_name = sci_name.replace('/', '')
        accession_to_organism[acc] = sci_name

# print(accession_to_organism)

for old_sig in glob('sigs/*sig'):
    acc = old_sig.split('/')[-1].split('.fna')[0].split('.')[0]
    if acc not in accession_to_organism:
        print(f"[ERROR] acc({acc}) not found!")
        continue
    organism = accession_to_organism[acc].replace(' ', '_')
    new_name = old_sig.replace(acc, f"{acc}_{organism}")
    os.rename(old_sig, new_name)