from pwa import PWA, AlignmentType, SubstitutionMatrix, AlignmentMethod

if __name__ == '__main__':
    sequences = input("FASTA:\n").lower()

    pwa_test = PWA(AlignmentType.PROTEIN, AlignmentMethod.LOCAL, sequences)
    # pwa_test.do_alignment()
    match = 1
    gap = -2
    pwa_test.do_alignment(match, -match, gap)

    print("\nShowing aligned sequences:")
    # pwa_test.print_pretty()
    pwa_test.print_aligned_seq()

    # check score
    print()
    print("score: ", pwa_test.get_final_score())
    print("best location: ", pwa_test.best_location)
    print("best score: ", pwa_test.best_score)

"""
    if pwa_test.alignment_type.name is AlignmentType.NUCLEOTIDE.name and \
        pwa_test.alignment_method.name is AlignmentMethod.LOCAL.name:
        score = 0
        for i in range(len(pwa_test.aligned_sequences[0])):
            top = pwa_test.aligned_sequences[0][i]
            bot = pwa_test.aligned_sequences[1][i]
            if top is '-' or bot is '-':
                score += gap
            elif top is bot:
                score += match
            elif top is not bot:
                score += -match
        print("calc'd score: ", score)

    print("\nTesting substitution matrix Enum:")
    d = SubstitutionMatrix.BLOSUM62.get_dictionary

    print(d['-', 'a']) if ('-','a') in d else print("d['-', 'a'] is not in the dictionary")
    print("d['w', 'w'] = ", d['w', 'w'])
"""