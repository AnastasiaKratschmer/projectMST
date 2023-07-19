from sp_approx import sp_approx, compute_cost
from sp_exact_3 import sp_exact_3, backtrack_3seq
from helpers.utils import get_cost_3, parse_fasta_multiple, create_score_matrix

'''
Script used for alignments  needed in the presentation
'''


def write_to_file(alignments, names, filename):
    with open(filename, 'w') as f_out:
        for i in range(0, len(alignments)):
            f_out.write('>' + names[i] + '\n')
            f_out.write(alignments[i] + '\n')

def align_and_write_to_file(seqs, names, opfile):
    score_matrix = create_score_matrix('testdata/score_matrix.txt')
    cost, alignment  = sp_approx(seqs, score_matrix, 5)
    print("Cost:", cost)

    new_alignment = list(map(list, zip(*alignment))) # NOTE for some reason python does not recognise its OWN list() function hera

    holder_list=[]
    for list in new_alignment:
        holder_list.append(''.join(list))
    
    write_to_file(holder_list, names, opfile)

def exact3_align_and_write_to_file(seqs, names, opfile):
    score_matrix = create_score_matrix('testdata/score_matrix.txt')
    D = sp_exact_3(seqs[0], seqs[1], seqs[2], score_matrix, 5)
    cost  = get_cost_3(D)
    print("Cost:", cost)

    alignment = backtrack_3seq(seqs[0], seqs[1], seqs[2], D, score_matrix, 5)

    holder_list=[]
    for list in alignment:
        holder_list.append(''.join(list))
    write_to_file(holder_list, names, opfile)

if __name__ == "__main__":
    print("Beep boop!")

    print("Seqs 1-3 (exact):")
    seqs, names = parse_fasta_multiple('testdata/brca1-testseqs.fasta')
    exact3_align_and_write_to_file(seqs[0:3], names, opfile = 'presentation/outputBRCA1_3_approx_exact.fasta')
    
    print("Seqs 1-3 (approx):")
    seqs, names = parse_fasta_multiple('testdata/brca1-testseqs.fasta')
    align_and_write_to_file(seqs[0:3], names, opfile = 'presentation/outputBRCA1_3_approx.fasta')
    
    print("Seqs 1-4:")
    seqs, names = parse_fasta_multiple('testdata/brca1-testseqs.fasta')
    align_and_write_to_file(seqs[0:3], names, opfile = 'presentation/outputBRCA1_4.fasta')
    
    print("Seqs 1-5:")
    seqs, names = parse_fasta_multiple('testdata/brca1-testseqs.fasta')
    align_and_write_to_file(seqs[0:3], names, opfile = 'presentation/outputBRCA1_5.fasta')
    
    print("Seqs 1-6:")
    seqs, names = parse_fasta_multiple('testdata/brca1-testseqs.fasta')
    align_and_write_to_file(seqs[0:3], names, opfile = 'presentation/outputBRCA1_6.fasta')

    print("Okay, bye!")
