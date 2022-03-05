import sys
import os
from glob import glob
import screed
from tqdm import tqdm

genomes_dir = sys.argv[1]
smash_taxonomy_file = sys.argv[2]

all_genomes = glob(f"{genomes_dir}/*gz")
accession_to_ident = {}
print("parsing genomes")
for genome in tqdm(all_genomes):
    for record in screed.open(genome):
        record_name = record.name
        ident, *remainder = record_name.split(' ', 1)
        accession = os.path.basename(genome).replace(".fna.gz", '').split('.')[0]
        accession_to_ident[accession] = ident
        break

print("writing fixed tax info")
new_file = os.path.join(os.path.dirname(smash_taxonomy_file), "fixed_" + os.path.basename(smash_taxonomy_file))
with open(smash_taxonomy_file) as tax, open(new_file, 'w') as NEW:
    NEW.write(next(tax))
    for line in tax:
        accession = line.split(',')[0]
        line = line.replace(accession, accession_to_ident[accession])
        NEW.write(line)
