import glob
from Bio import SeqIO

"""
Note that this code will fail on phages with non-standard genetic codes!
"""
############################################################################
###Constants
############################################################################
fasta_dir = '../../Data/model_data/phage_data_nmicro2017/phage_fasta_files/'
min_prot_length = 40
############################################################################
############################################################################

for fasta_file in glob.glob(fasta_dir+'*.fasta')[:]:
    print(fasta_file)
    nt_record = SeqIO.read(fasta_file, 'fasta')
    db_id = fasta_file.split('/')[-1].split('.fasta')[0]
    prots = []
    for i in [0, 1, 2]:###Each of the three reading frames
        ###Translate the regular strand
        tempy = nt_record[i:].translate()
        ###Split according to stop codons
        seq = str(tempy.seq).split('*')
        ###Append the sequences if they are longer than the defined minimum
        for j in seq:
            if len(j) >= min_prot_length:
                prots.append(j)
        ###Repeat for the reverse complement
        tempy = nt_record.reverse_complement()[i:].translate()
        seq = str(tempy.seq).split('*')
        for j in seq:
            if len(j) >= min_prot_length:
                prots.append(j)
    ###Write out the results in the same directory
    with open(fasta_file.replace('.fasta', '_6frame.faa'), 'w') as outfile:
        for i, seq in enumerate(prots):
            outfile.write('>{}_{}\n{}\n'.format(db_id, i, seq))
