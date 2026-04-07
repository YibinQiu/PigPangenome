import os

input_path = os.path.expanduser('~/pangenome/genofile/clean/merge_SVtransMissing_retain.vcf')

# Recessive model (0/1 -> 0/0)
with open(input_path, 'r') as infile, open('merge_SVtransMissing_retain_recessive.vcf', 'w') as outfile:
    for line in infile:
        if line.startswith('#'):
            outfile.write(line)
        else:
            outfile.write(line.replace('0/1', '0/0'))

# Dominant model (0/1 -> 1/1)
with open(input_path, 'r') as infile, open('merge_SVtransMissing_retain_dominant.vcf', 'w') as outfile:
    for line in infile:
        if line.startswith('#'):
            outfile.write(line)
        else:
            outfile.write(line.replace('0/1', '1/1'))
