from sys import setrecursionlimit
from pwa import PWA, AlignmentType, SubstitutionMatrix
if __name__ == '__main__':
    sequences = input("FASTA:\n").lower()
    pwa_test = PWA(AlignmentType.NUCLEOTIDE, sequences)
    pwa_test.do_alignment()

    # pwa_test.print_pretty()
    pwa_test.print_aligned_seq()

    # check score

    print("score: ", pwa_test.get_best_score())

    score = 0
    for i in range(len(pwa_test.aligned_sequences[0])):
        top = pwa_test.aligned_sequences[0][i]
        bot = pwa_test.aligned_sequences[1][i]
        if top is '-' or bot is '-':
            score += -2
        elif top is bot:
            score += 1
        elif top is not bot:
            score += -1
    print("calc'd score: ", score)


    # print("test substitution matrix Enum")
    # d = SubstitutionMatrix.BLOSUM62.get_dictionary
    # print(d['-','a'])
    # print(d['w', 'w'])