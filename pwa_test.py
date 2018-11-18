from pwa import PWA, AlignmentType, SubstitutionMatrix, AlignmentMethod


def get_test_sequences(get_aa: bool):
    fixed_sequence = list()

    if not get_aa:
        fixed_sequence.append(r">first\ngggggaccgtgggg\n>second\nacttgt")

        fixed_sequence.append(""">sequence A
ggtaagtcctctagtacaaacacccccaatattgtgatataattaaaattatattcatat
tctgttgccagaaaaaacacttttaggctatattagagccatcttctttgaagcgttgtc
>sequence B
ggtaagtgctctagtacaaacacccccaatattgtgatataattaaaattatattcatat
tctgttgccagattttacacttttaggctatattagagccatcttctttgaagcgttgtc
tatgcatcgatcgacgactg""")

        fixed_sequence.append("""
        > A
        1 GACTTTTTTT TTTTTTCCTT TGGGAAAGGT AGGGAGGTGT TCGTACGGGA GCAGCCTCGG
       61 GGACCCCTGC ACTGGGTCAG GGCTTATGAA GCTAGAAGCG TCCCTCTGTT CCCTTTGTGA
      121 GTTGGTGGGT TGTTGGTACA TTTGGTTGGA AGCTGTGTTG CTGGTTAGGG AGACTCGGTT
      181 TTGCTCCTTG GGTTCGAGGA AAGCTGGAGA ATAGAAGCCA TTGTTTGCCG TCTGTCGGCT
      241 TTGTCGACCA CGCTCACCCC CTCCTGTTCG TACTTTTTAA AGCAGTGAGG CGAGGTAGAC
        > B
      301 AGGGTGTGTC ACAGTACAGT TAAAGGGGTG AAGATCTAAA CGCCAAAAGA GAAGTTAATC
      361 ACAATAAGTG AGGTTTGGGA TAAAAAGTTG GGCTTGCCCC TTTCAAAGTC CCAGAAAGCT
      421 GGGAGGTAGA TGGAGAGGGG GCCATTGGGA AGTTTTTTTG GTGTAGGGAG AGGAGTAGAA
      481 GATAAAGGGT AAGCAGAGTG TTGGGTTCTG GGGGTCTTGT GAAGTTCCTT AAGGAAGGAG
      541 GGAGTGTGGC CCTGCAGCCC TCCCAAACTG CTCCAGCCTA TGCTCTCCGG CACCAGGAAG
      > C
      601 TTCCAAGGTT CCCTTCCCCT GGTCTCCAAA CTTCAGGTAT TCCTCTCCCC TCACACCCCT
      661 TCAACCTCAG CTCTTGGCCT CTACTCCTTA CTCCACTGTT CCTCCTGTTT CCCCCTTCCC
      721 CTTTTCCTGG TTCTTTATAT TTTTGCAAAG TGGGATCCGA ACTTGCTAGA TTTTCCAATT
      781 CTCCCAAGCC AGACCAGAGC AGCCTCTTTT AAAGGATGGA GACTTCTGTG GCAGATGCCG
      841 CTGAAAATGT GGGTGTAATG CTGGGACTTA GAGTTTGATG ACAGTTTGAC TGAGCCCTAG
      901 ATGCATGTGT TTTTCCTGAG AGTGAGGCTC AGAGAGCCCA TGGACGTATG CTGTTGAACC
      > D
      961 ACAGCTTGAT ATACCTTTTT CTCCTTCTGT TTTGTCTTAG GGGGAAGACT TTAACTAGGG
     1021 GCGCGCAGAT GTGTGAGGCC TTTTATTGTG AGAGTGGACA GACATCCGAG ATTTCAGGCA
     1081 AGTTCTGTGG TGGCTGCTTT GGGCT
     """)

    else:
        fixed_sequence.append(""">P29600|SUBS_BACLE Subtilisin Savinase - Bacillus lentus
        AQSVPWGISRVQAPAAHNRGLTGSGVKVAVLDTGISTHPDLNIRGGASFVPGEPSTQDGN
        GHGTHVAGTIAALNNSIGVLGVAPSAELYAVKVLGASGSGSVSSIAQGLEWAGNNGMHVA
        NLSLGSPSPSATLEQAVNSATSRGVLVVAASGNSGAGSISYPARYANAMAVGATDQNNNR
        ASFSQYGAGLDIVAPGVNVQSTYPGSTYASLNGTSMATPHVAGAAALVKQKNPSWSNVQI
        RNHLKNTATSLGSTNLYGSGLVNAEAATR

        >P41363|ELYA_BACHD Thermostable alkaline protease precursor - Bacillus halodurans
        MRQSLKVMVLSTVALLFMANPAAASEEKKEYLIVVEPEEVSAQSVEESYDVDVIHEFEEI
        PVIHAELTKKELKKLKKDPNVKAIEKNAEVTISQTVPWGISFINTQQAHNRGIFGNGARV
        AVLDTGIASHPDLRIAGGASFISSEPSYHDNNGHGTHVAGTIAALNNSIGVLGVAPSADL
        YAVKVLDRNGSGSLASVAQGIEWAINNNMHIINMSLGSTSGSSTLELAVNRANNAGILLV
        GAAGNTGRQGVNYPARYSGVMAVAAVDQNGQRASFSTYGPEIEISAPGVNVNSTYTGNRY
        VSLSGTSMATPHVAGVAALVKSRYPSYTNNQIRQRINQTATYLGSPSLYGNGLVHAGRAT
        Q
        """)

        fixed_sequence.append("""> dfdfd
        ettrgrksaqpesaalpdapastaptrsktpaqglarklhfstappnpdapwtprvagfnkrvfcaavgrlaamharmaavqlwdmsrprtdedlnellgittirvtvcegknllqranelvnpdvvqdvdaatatrgrsaasrpterpraparsasr
        > dfsafsdf
        dstgefcwichqpegplkrfcgckgscavshqdclrgwletsrrqtcalcgtpysmkwktkplrewtwgeeevlaameaclplvliplavlmivmgtwllvnhngflsprmqvvlvvivllamivfsasasyvmvegpgcldtctaknstvtvnsideaiatqqptktdlglaretlstrfrrgkcrsccr""")

        fixed_sequence.append(r">a\ntrvmk\n>b\ndvck")

        fixed_sequence.append(r">a\ntrvmpppp\n>b\ntrvmpppp")

    return fixed_sequence

if __name__ == '__main__':
    test_peptides = False

    fixed_sequences = get_test_sequences(test_peptides)
    for seq in fixed_sequences:
        align_type = AlignmentType.PROTEIN if test_peptides else AlignmentType.NUCLEOTIDE
        pwa_test = PWA(align_type, AlignmentMethod.LOCAL, seq)
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
        print("_____________________________________________________________________")

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