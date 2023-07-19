from sp_approx import sp_approx
from helpers.utils import parse_fasta_multiple, create_score_matrix
'''
This script is for the second part of the project, where we are to align 
the 8 sequences of length 200 in brca1-testseqs.fasta as close to optimum 
as possible (the sequences corresponds to first 200 nucleotides in the mRNA 
from the BRCA1 gene in cow, wolf, chicken, human, macaque, mouse, 
chimpanzee, and rat).
'''

def write_to_file(alignments, names, filename):
    with open(filename, 'w') as f_out:
        for i in range(0, len(alignments)):
            f_out.write('>' + names[i] + '\n')
            f_out.write(alignments[i] + '\n')

if __name__ == "__main__":
    print("Beep boop!")
    
    seqs, names = parse_fasta_multiple('testdata/brca1-testseqs.fasta')
    score_matrix = create_score_matrix('testdata/score_matrix.txt')

    print("Aligning the sequences corresponds to first 200 nucleotides in the mRNA from\nthe BRCA1 gene in cow, wolf, chicken, human, macaque, mouse, chimpanzee, and rat...")
    cost, alignment  = sp_approx(seqs, score_matrix, 5)

    new_alignment = list(map(list, zip(*alignment)))

    holder_list=[]
    for list in new_alignment:
        holder_list.append(''.join(list))
    
    opfile = 'brca1-alignment.fasta'
    write_to_file(holder_list, names, opfile)   
    print("Done! Written to the file,", opfile, ". :-)\nBye!\n")
