from sp_exact_3 import sp_exact_3, backtrack_3seq
from sp_approx import sp_approx
from helpers.utils import parse_fasta_multiple, create_score_matrix, get_cost_3
from os import listdir
from matplotlib import pyplot as plt
from tqdm import tqdm  # loading bar

gap = 5
score_matrix = create_score_matrix("testdata/score_matrix.txt")

def run_first_exp():
    print("--- First experiment --- ")
    seqs, names = parse_fasta_multiple("testdata/brca1-testseqs.fasta")
    s1, n1 = seqs[0], names[0]
    s2, n2 = seqs[1], names[1]
    s3, n3 = seqs[2], names[2]
    print('The cost of alignment of the three sequences, brca1_bos_taurus, brca1_canis_lupus and brca1_gallus_gallus:')
    D = sp_exact_3(s1, s2, s3, score_matrix, gap)
    print(get_cost_3(D))
    alignments = backtrack_3seq(s1,s2,s3, D, score_matrix, gap)
    print("Alignments:\n> " + n1+"\n", alignments[0], "\n> "+n2+"\n", alignments[1], "\n> "+n3+"\n", alignments[2])
    print()

def run_second_exp():
    print("--- Second experiment ---")
    seqs, names = parse_fasta_multiple("testdata/brca1-testseqs.fasta")
    print('The cost of alignment of the three sequences, brca1_bos_taurus, brca1_canis_lupus and brca1_gallus_gallus:')
    cost, _, center_string = sp_approx(seqs[0:5], score_matrix, gap, return_center_string=True)
    print("The approximate cost of the first 5 seqs in brca1-testseqs.fasta is:")
    print(str(cost) + ".")
    print("The center string was:")
    print(names[center_string])
    print()

def run_third_exp():
    print("--- Third experiment ---")
    score_exact = []
    score_approx = []
    lengths = []
    for file in tqdm(listdir('testdata/testseqs')):
        # print('file is \t' + file)
        seqs, _ = parse_fasta_multiple('testdata/testseqs/' + file)
        score_exact.append(get_cost_3(sp_exact_3(seqs[0], seqs[1], seqs[2], score_matrix, gap)))
        score_approx.append(sp_approx(seqs, score_matrix, gap)[0])
        lengths.append(sum(map(len, seqs))/3) # log average length of sequence

    ratios = [a/b for a,b in zip(score_approx, score_exact)]
    average = sum(ratios)/len(ratios)

    fig, ax = plt.subplots()
    ax.set_title('Approximation ratio between sp_exact and sp_approx per sequence length')
    ax.set_xlabel('Sequence lengths')
    ax.set_ylabel('Approximation ratio')
    ax.bar(lengths, ratios, alpha=0.5, color="#D81B60")
    ax.axhline(y=4/3, color='#004D40', linestyle='-', label='upperbound')
    ax.axhline(y=average, color='#FFC107', linestyle='dashed', label='average')
    ax.axhline(y=1, color='#1E88E5', linestyle='-', label='lowerbound')
    ax.legend(loc="lower right")
    fig.savefig('figures/barplot', bbox_inches='tight')
    plt.show()

    fig, ax = plt.subplots()
    ax.set_title('Approximation ratio between sp_exact and sp_approx per sequence length')
    ax.set_xlabel('Sequence lengths')
    ax.set_ylabel('Approximation ratio')
    ax.scatter(lengths, ratios, c='#D81B60', marker='*')
    ax.axhline(y=4/3, color='#004D40', linestyle='-',label='upperbound')
    ax.axhline(y=average, color='#FFC107', linestyle='dashed', label='average')
    ax.axhline(y=1, color='#1E88E5', linestyle='-',label='lowerbound')
    ax.legend(loc="lower right")
    fig.savefig('figures/scatterplot', bbox_inches='tight')
    plt.show()

    print("On average, the approximation ratio was: ")
    print(average)
    print("Generated are plots in the figures-folder. :)")
    print()

if __name__ == "__main__":
    print("Beep boop!")
    print("Which experiment would you like to run, 1, 2, or 3?")
    while True:        
        userInput = input()
        if (userInput == "1"):
            run_first_exp()
        elif (userInput == "2"):
            run_second_exp()
        elif (userInput == "3"):
            run_third_exp()
        elif(userInput.lower() == "x"):
            break
        print("Experiment complete!")
        print("Which experiment would you like to run now, 1, 2, or 3? Press X to exit.")

    print("Okay, bye!")