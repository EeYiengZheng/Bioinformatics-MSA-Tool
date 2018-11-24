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
      301 AGGGTGTGTC ACAGTACAGT TAAAGGGGTG AAGATCTAAA CGCCAAAAGA GAAGTTAATC
      361 ACAATAAGTG AGGTTTGGGA TAAAAAGTTG GGCTTGCCCC TTTCAAAGTC CCAGAAAGCT
      421 GGGAGGTAGA TGGAGAGGGG GCCATTGGGA AGTTTTTTTG GTGTAGGGAG AGGAGTAGAA
      481 GATAAAGGGT AAGCAGAGTG TTGGGTTCTG GGGGTCTTGT GAAGTTCCTT AAGGAAGGAG
      541 GGAGTGTGGC CCTGCAGCCC TCCCAAACTG CTCCAGCCTA TGCTCTCCGG CACCAGGAAG
      > B
      601 TTCCAAGGTT CCCTTCCCCT GGTCTCCAAA CTTCAGGTAT TCCTCTCCCC TCACACCCCT
      661 TCAACCTCAG CTCTTGGCCT CTACTCCTTA CTCCACTGTT CCTCCTGTTT CCCCCTTCCC
      721 CTTTTCCTGG TTCTTTATAT TTTTGCAAAG TGGGATCCGA ACTTGCTAGA TTTTCCAATT
      781 CTCCCAAGCC AGACCAGAGC AGCCTCTTTT AAAGGATGGA GACTTCTGTG GCAGATGCCG
      841 CTGAAAATGT GGGTGTAATG CTGGGACTTA GAGTTTGATG ACAGTTTGAC TGAGCCCTAG
      901 ATGCATGTGT TTTTCCTGAG AGTGAGGCTC AGAGAGCCCA TGGACGTATG CTGTTGAACC
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

        fixed_sequence.append("""
    >Pongo_pygmaeus
    MDLSAVRVEEVQNVINAMQKILECPICLELIKEPVSTKCDHIFCKFCMLKLLNQKKGPSQCPLCKNDITK
    RSLQESTRFSQLVEELLKIICAFQLDTGLQYANSYNFAKKENNSPEHLKDEVSIIQSMGYRNRAKRLLQS
    EPENPSLQETSPSVQLSNLGTVRTLRTKQRIQPQKKSVYIELGSDSSEDTVNKATYCSVGDQELLQITPQ
    GTSDEISLDSAKKAACEFSETDVTNTEHHQPSNNDLNTTEKRATERHPEKYQGSSVSNLHVEPCGTNTHA
    SSLQHENSSLLLTKDRMNVEKAEFCNKSKQPGLARSQHNRWAGSKETCNDRQTPSTEKKVDLNADPLCER
    KEWNKQKLPCSENPRDTEDVPWITLNSSIQKVNEWFSRSDELLGSDDSHDGRSESNAKVADVLDVLNEVD
    EYSGSSEKIDLLASDPHEALICKSERVHSKSVESNIEDKIFGKTYRRKASLPNLSHVTENLIIGAFVTEP
    QIIQERPLTNKLKRKRRATSGLHPEDFIKKADLAVQKTPEMINQGTNQMEQNGQVMNITNSGHENKTKGD
    SIQNEKNPNPIESLEKESAFKTKAEPISSSISNMELELNIHNSKAPKKNRLRRKSSTRHIHALELVVSRN
    LSPPNCTELQIDSCSSSEEIKKKKYNQMPVRHSRNLQLMEDKEPATGAKKSNKPNEQTSKRHDSDTFPEL
    KLTNAPGSFTNCSNTSELKEFVNPSLPREEKEEKLGTVKVSNNAKDPKDLMLSGERVLQTERSVESSSIS
    LVPGTDYGTQESISLLEVSTLGKAKTEPNKCVSQCAAFENPKELIHGCFKDTRNDTEGFKYPLGHEVNHS
    QETSIEMEESELDTQYLQNTFKVSKRQSFALFSNPGNPEEECATFSAHSRSLKKQSPKVTFECEQKEENQ
    GKNESNIKPVQTANITAGFPVVCQKDKPVDYAKCSIKGGSRFCLSSQFRGNETGLITPNKHGLSQNPYHI
    PPLFPIKSFVKTKCKKNLLEENSEEHSMSPEREMGNENIPSTVSIISRNNIRENVFKEASSSNINEVGSS
    TNEVGSSINEVGSSDENIQAELGRSRGPKLNAMLRLGVLQPEVYKQSFPGSNGKHPEIKKQEYEEVLQTV
    NTDFSPCLISDNLEQPMRSSHASQVCSETPNDLLDDGEIKEDTSFAENDIKESSAVFSKSVQRGELSRSP
    SPFTHTHLAQGYRRGAKKLESSEENLSSEDEELPCFQHLLFGKVSNIPSQSTRHSTVATECLSKNTEENL
    LSLKNSLNDYSNQVILVKASQEHHLSEETKCSASLFSSQCSELEDLTANTNTQDRFFIGSSKQMRHQSES
    QGVGLSDKELVSDDEERGTDLEENNQEEQGVDSNLGEAASGYESETSVSEDCSGLSSQSDILTTQQRDTM
    QDNLIKLQQEMAELEAVLEQHGSQPSNSYPSIISDSSALEDLRNPEQSTSEKAVLTSQKSSEYPISQNPE
    GLSADKFEVSADSSTNKNKEPGVERSSPSKCPSLDDRWYMHSCSGSLQNGNYPSQEELIKVVDVEKQQLE
    ESGPHDLTEPSYLPRQDLEGTPYLESGISLFSDDPESDASEDRAPESAHVGSIPSSTSALKVPQLKVAES
    AQSPAAAQTTNTAGYNAMEESVSREKPELTASTERVNKRMSMVVSGLTPEEFMLVYKFARKHHITLTNLI
    TEETTHVVMKTDAEFVCERTLKYFLGIAGGKWVVSYFWVTQSIKERKMLNEHDFEVRGDVVNGRNHQGPK
    RARESQDRKIFRGLEICCYGPFTNMPTDQLEWIVQLCGASVVKELSSFTLGTGVHPIVVVQPDAWTEDNG
    FHAIGQMCEAPVVTREWVLDSVALYQCQELDTYLIPQIPHSHY
    >Mus musculus
    MDLSAVQIQEVQNVLHAMQKILECPICLELIKEPVSTKCDHIFCKFCMLKLLNQKKGPSQCPLCKNEITK
    RSLQGSTRFSQLAEELLRIMAAFELDTGMQLTNGFSFSKKRNNSCERLNEEASIIQSVGYRNRVRRLPQV
    EPGNATLKDSLGVQLSNLGIVRSVKKNRQTQPRKKSVYIELDSDSSEETVTKPGDCSVRDQELLQTAPQE
    AGDEGKLHSAEEAACEFSEGIRNIEHHQCSDDLNPTENHATERHPEKCQSISISNVCVEPCGTDAHASSL
    QPETSSLLLIEDRMNAEKAEFCNKSKQPGIAVSQQSRWAASKGTCNDRQVPSTGEKVGPNADSLSDREKW
    THPQSLCPENSGATTDVPWITLNSSVEKVNEWFSRTGEMLTSDSASARRHESNAEAAVVLEVSNEVDGGF
    SSSRKTDLVTPDPHHTLMCKSGRDFSKPVEDNISDKIFGKSYQRKGSRPHLNHVTEIIGTFITEPQITQE
    QPFTNKLKRKRSTSLQPEDFIKKADSAGVQRTPDNINQGTDLMEPNEQAVSTTSNCQENQIAGSNLQKEK
    SAHPTESLRKEPASTAGAKSISNSVSDLEVELNVHSSKAPKKNRLRRKSSIRCALPLEPISRNPSPPTCA
    ELQIDSCGSSEETKKNHSNQQPAGHLREPQLIEDTEPAADAKKNEPNEHIRKRRASDAFPEEKLMNKAGL
    LTSCSSPRKSQGPVNPSPQRTGTEQLETRQMSDSAKELGDRVLGGEPSGKTTDRSEESTSVSLVPDTDYD
    TQNSVSVLDAHTVRYARTGSAQCMTQFVASENPKELVHGSNNAGSGTEGLKPPLRHALNLSQEKVEMEDS
    ELDTQYLQNTFQVSKRQSFALFSKPRSPQKDCAHSVPSKELSPKVTAKGKQKERQGQEEFEISHVQAVAA
    TVGLPVLCQEGKLAADTMCDRGSRLCPSSHYRSGENGLSATGKSGISQNSHFKQSVSPIRSSIKTDNRKP
    LTEGRFERHTSSTEMAVGNENILQSTVHTVSLNNRGNACQEAGSGSIHEVCSTGDSFPGQLGRNRGPKVN
    TVPPLDSMQPGVCQQSVPVSDKYLEIKKQEGEAVCADFSPCLFSDHLEQSMSGKVFQVCSETPDDLLDDV
    EIQGHTSFGEGDIMERSAVFNGSILRRESSRSPSPVTHASKSQSLHRASRKLESSEESDSTEDEDLPCFQ
    HLLSRISNTPELTRCRSAVTQGIPEKAEGTQAPWKGSSSDCNNEVIMIEASQEHQFSEDPRCSGRMFSSQ
    NSAAQGSTANANSQDSNFIPPSKQRSHQCGNEEAFLSDKELISDNEEMATCLEEDNDQEEDSIIPDSEAS
    GYESETNLSEDCSQSDILTTQQRATMKYNLIKLQQEMAHLEAVLEQRGNQPSGHSPSLLADPCALEDLPD
    LEPNMSGAAILTSKNINENPVSQNLKSACDDKFQLQHLEGPTSGDDESGMGRPSPFKSPLAGSRGSAHGC
    SRHLQKRNSPSQEELLQPAGSEASSEPHNSTGQSCLPRRELEGTPYLGSGISLFSSRDPESESPKEPAHI
    GTTPASTSALKIPQGQVAFRSAAAAGADKAVVGIVSKIKPELTSSEERADRDISMVVSGLTPKEVMTVQK
    FAEKYRLTLTDAITEETTHVIIKTDAEFVCERTLKYFLGIAGGKWIVSYSWVVRSIQERRLLNVHEFEVK
    GDVVTGRNHQGPRRSRESREKLFKGLQVYCCEPFTNMPKDELERMLQLCGASVVKELPSLTHDTGAHLVV
    IVQPSAWTEDSNCPDIGQLCKARLVMWDWVLDSLSSYRCRDLDAYLVQNITCDSSEPQDSND
    """)

    return fixed_sequence


def get_dq_sequences(A_or_B: str):
    if A_or_B is 'a':
        return ["""
        >DQA1 [human : Homo Sapiens]
            1 milnkalmlg alalttvmsp cggedivadh vasygvnlyq sygpsgqfth efdgdeefyv
           61 dlerketvwk lplfhrlrfd pqfaltniav lkhnlnilik rsnstaatne vpevtvfsks
          121 pvtlgqpntl iclvdnifpp vvnitwlsng hsvtegvset sflsksdhsf fkisyltflp
          181 sadeiydckv ehwgldepll khwepeipap mseltetvvc alglsvglvg ivvgtvliir
          241 glrsvgasrh qgpl
        >DQA1 [macaque : Macaca mulatta]
            1 milnkalllg alalttvmsp crgedivadh iasygvnlyq tyglsgqfth efdgdeqfyv
           61 dlerketvwr lplfskfggf dpqgalrnla vgkhnlnili khsnstaatn evpevtvfsk
          121 spvtlghpnt liclvdnifp pvvnitwlsn ghsvtegvse tsflsksdhs ffkisyltfl
          181 psadeiydcr vehwgldepl lkhwepeipa pmseltetvv calglsmglv givvgtvlii
          241 rglrsvgasr hqgpl
        >DQA1 [marmoset : Callithrix jacchus]
            1 milnkalmlg alvlttvmsp cggeniaadh vaacginlyq sygptgqyth efdgdeqfyv
           61 dlgrketlwr wpmlskfggf dpqgaltnia tmkhnldivi krsnctaatn evpeatvfsk
          121 spvalgqpnt liclvdnifp pvvnitwlsn ghsvtegvse tsflsksdhs ffkmsyltfl
          181 psadeiydck vehwgldepl lkhwepeipa ptseltetvv calglsvglv givvgtvfii
          241 qglrsvgasr hqgpl
        """]
    if A_or_B is 'b':
        return ["""
          >DQB1 [human : Homo Sapiens]
            1 mswkkalrip gglrvatvtl mlamlstsva egrdspedfv yqfkgmcyft ngtervrlvt
           61 ryiynreeya rfdsdvgvyr avtplgppda eywnsqkevl egtraeldtv crhnyqlelr
          121 ttlqrrvept vtispsrtea lnhhnllvcs vtdfypaqik vrwfrndqee ttgvvstpli
          181 rngdwtfqil vmlemtpqrg dvytchvehp slqtpitvew raqsesaqsk mlsgiggfvl
          241 gliflglgli ihhrsqkgll h
          >DQB1 [macaque : Macaca mulatta]
            1 mswkkalrip gglrvatvtl mlamlstpva egrdspedfv yqfkglcyft ngtervrsvt
           61 rhiynreeyv rfdsdwgeyr avtpqgrrsa eyfnsqkdil estraeldtv cmhnyevayr
          121 gilqrrvept vtislsrtea lnhhnllvcs vtdfypgqik vrwfrndqee ttgivstpli
          181 rngdwtfqil vmlemtpqrg dvytchvehp slqspitvew raqsesaqsk mlsgiggfvl
          241 gmiflglgli ihhrsqkglp h
          >DQB1 [marmoset : Callithrix jacchus]
            1 mswkkalwip gglraaavtl mlvmlsspva egrdspedfv yqfkflcyft ngtervrlvt
           61 eyvynreehv rfdsdvgeyw avtplglpda kywnsqkdil ertraeldtv crhnyevafr
          121 gilqwrvept viispsktea lnhhnllvcs vtdfypgqik vrwfrndqeq tagivstpli
          181 rngdwtfqil vmlemtpqrg dvytchvehp slqspitvew raqsesaqsk mlsgiggfvl
          241 gliflglgli irqrsrkgpq gpppagllh
      """]

def test_pwa():
    test_peptides = True

    seqs = get_test_sequences(test_peptides)
    for seq in seqs:
        align_type = AlignmentType.PROTEIN if test_peptides else AlignmentType.NUCLEOTIDE
        pwa_test = PWA(align_type, AlignmentMethod.LOCAL, seq)
        match = 5
        mismatch = -5
        gap = -3
        pwa_test.do_alignment(match, mismatch, gap)

        print("\nShowing aligned sequences:")
        # pwa_test.print_pretty()
        # pwa_test.print_aligned_seq()
        pwa_test.print_aligned_seq_formatted(show_conserved=True)

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


def test_dq():
    test_peptides = True

    dq = ["""
> DQB1 [pango]
        1 mswkkalrip gdlrvatvtl mlamlsslla egrdspedfv fqfkgmcyft ngtervrlvt
       61 ryiynreefv rfysdvgvyr avtpqgrpva eywnsqkevl ertraeldtv crhnyevayr
      121 gilqrrvept vtispsrtea lnhhnllvcs vtdfypgqik vrwfrndqee tagvvstpli
      181 rngdwtfqil vmlemtpqrg evytchvehp slqspitvew raqsesaqsk mlsgvggfvl
      241 gliflglgli irqrsqkgll h
> DQB1*020101 [human : Homo Sapiens]
MSWKKALRIPGGLRAATVTLMLSMLSTPVAEGRDSPEDFVYQFKGMCYFTNGTERVRLVSRSIYNREEIVRFDSDVGEFR
AVTLLGLPAAEYWNSQKDILERKRAAVDRVCRHNYQLELRTTLQRRVEPTVTISPSRTEALNHHNLLVCSVTDFYPAQIK
VRWFRNDQEETAGVVSTPLIRNGDWTFQILVMLEMTPQRGDVYTCHVEHPSLQSPITVEWRAQSESAQSKMLSGIGGFVL
GLIFLGLGLIIHHRSQKGLLH
> DQA1*050101 [human : Homo Sapiens]
MILNKALMLGALALTTVMSPCGGEDIVADHVASYGVNLYQSYGPSGQYTHEFDGDEQFYVDLGRKETVWCLPVLRQFRFD
PQFALTNIAVLKHNLNSLIKRSNSTAATNEVPEVTVFSKSPVTLGQPNILICLVDNIFPPVVNITWLSNGHSVTEGVSET
SFLSKSDHSFFKISYLTLLPSAEESYDCKVEHWGLDKPLLKHWEPEIPAPMSELTETVVCALGLSVGLVGIVVGTVFIIR
GLRSVGASRHQGPL
"""]

    for seq in dq:
        align_type = AlignmentType.PROTEIN if test_peptides else AlignmentType.NUCLEOTIDE
        pwa_test = PWA(align_type, AlignmentMethod.GLOBAL, seq)
        match = 5
        mismatch = -5
        gap = -10
        pwa_test.do_alignment(match, mismatch, gap)

        print("\nShowing aligned sequences:")
        # pwa_test.print_pretty()
        # pwa_test.print_aligned_seq()
        pwa_test.print_aligned_seq_formatted(show_conserved=True)

        # check score
        print()
        print("score: ", pwa_test.get_final_score())
        print("best location: ", pwa_test.best_location)
        print("best score: ", pwa_test.best_score)
        print("_____________________________________________________________________")

def test_x(seq: str, match=1, mismatch=-1, gap=-2):
    pwa_test = PWA(AlignmentType.PROTEIN, AlignmentMethod.LOCAL, seq)
    pwa_test.do_alignment(match, mismatch, gap)
    pwa_test.print_aligned_seq_formatted(max_len=60,show_mismatch_only=True)

if __name__ == '__main__':
    # test_dq()
    # test_pwa()
    print("DQB1 02:01 vs 03:02")
    test_x(""">_DQA1_[cow]
MVLNRALILGALALTTMMSSSGGEDIVADHVGSYGTEIYQSHGPSGQYTQEFDGDEMFYVDLGKKETVWR
LPMFSQFAGFDPQAALSEIATAKHNLDVLTKRSNFTPVINEVPEVTVFSKSPVMLGQPNTLICHVDNIFP
PVINITWLKNGHAVTEGVSETSFLPKDDHSFLKIGYLTFLPSDNDIYDCKVEHWGLDEPLLKHWEPEVPA
PMSELTETVVCALGLTVGLVGIVVGTIFIIQGLRSGGASRHQGPL
>_DQA1_outgroup_[western_claw_frog]
MLLLCALLALALKPSDAVTVDYFDYSALFYQTRRPTGEYLFDYNGNELFHVDLDSKSVVWTLPGLSEHES
FDPQGALQDINVARYNLDIGIKRSNSTAATNIPPLVNLYSSKPVVLGEPNILICFVKNIFPPVMNTTWIK
NGEKIEEGFTETSYLPAQDHSFRKLHYLAFIPNEHDIYTCEIEHWGLERPTRRVWQHDVPTPTSEAYQNV
ICALGLAVGIIGIIAGVMLIIKGMKQSAAQGRSQR
""")