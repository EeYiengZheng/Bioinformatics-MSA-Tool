"""
cs 123a bioinformatics

ee yieng zheng, raymond hong
"""

from enum import Enum


class AlignmentType(Enum):
    NUCLEOTIDE = 1
    PROTEIN = 2

class SubstitutionMatrix(Enum):
    BLOSUM62 = [
        [' ','a','r','n','d','c','q','e','g','h','i','l','k','m','f','p','s','t','w','y','v','b','z','x','-'],
        ['a','4','-1','-2','-2','0','-1','-1','0','-2','-1','-1','-1','-1','-2','-1','1','0','-3','-2','0','-2','-1','0','-4'],
        ['r','-1','5','0','-2','-3','1','0','-2','0','-3','-2','2','-1','-3','-2','-1','-1','-3','-2','-3','-1','0','-1','-4'],
        ['n','-2','0','6','1','-3','0','0','0','1','-3','-3','0','-2','-3','-2','1','0','-4','-2','-3','3','0','-1','-4'],
        ['d','-2','-2','1','6','-3','0','2','-1','-1','-3','-4','-1','-3','-3','-1','0','-1','-4','-3','-3','4','1','-1','-4'],
        ['c','0','-3','-3','-3','9','-3','-4','-3','-3','-1','-1','-3','-1','-2','-3','-1','-1','-2','-2','-1','-3','-3','-2','-4'],
        ['q','-1','1','0','0','-3','5','2','-2','0','-3','-2','1','0','-3','-1','0','-1','-2','-1','-2','0','3','-1','-4'],
        ['e','-1','0','0','2','-4','2','5','-2','0','-3','-3','1','-2','-3','-1','0','-1','-3','-2','-2','1','4','-1','-4'],
        ['g','0','-2','0','-1','-3','-2','-2','6','-2','-4','-4','-2','-3','-3','-2','0','-2','-2','-3','-3','-1','-2','-1','-4'],
        ['h','-2','0','1','-1','-3','0','0','-2','8','-3','-3','-1','-2','-1','-2','-1','-2','-2','2','-3','0','0','-1','-4'],
        ['i','-1','-3','-3','-3','-1','-3','-3','-4','-3','4','2','-3','1','0','-3','-2','-1','-3','-1','3','-3','-3','-1','-4'],
        ['l','-1','-2','-3','-4','-1','-2','-3','-4','-3','2','4','-2','2','0','-3','-2','-1','-2','-1','1','-4','-3','-1','-4'],
        ['k','-1','2','0','-1','-3','1','1','-2','-1','-3','-2','5','-1','-3','-1','0','-1','-3','-2','-2','0','1','-1','-4'],
        ['m','-1','-1','-2','-3','-1','0','-2','-3','-2','1','2','-1','5','0','-2','-1','-1','-1','-1','1','-3','-1','-1','-4'],
        ['f','-2','-3','-3','-3','-2','-3','-3','-3','-1','0','0','-3','0','6','-4','-2','-2','1','3','-1','-3','-3','-1','-4'],
        ['p','-1','-2','-2','-1','-3','-1','-1','-2','-2','-3','-3','-1','-2','-4','7','-1','-1','-4','-3','-2','-2','-1','-2','-4'],
        ['s','1','-1','1','0','-1','0','0','0','-1','-2','-2','0','-1','-2','-1','4','1','-3','-2','-2','0','0','0','-4'],
        ['t','0','-1','0','-1','-1','-1','-1','-2','-2','-1','-1','-1','-1','-2','-1','1','5','-2','-2','0','-1','-1','0','-4'],
        ['w','-3','-3','-4','-4','-2','-2','-3','-2','-2','-3','-2','-3','-1','1','-4','-3','-2','11','2','-3','-4','-3','-2','-4'],
        ['y','-2','-2','-2','-3','-2','-1','-2','-3','2','-1','-1','-2','-1','3','-3','-2','-2','2','7','-1','-3','-2','-1','-4'],
        ['v','0','-3','-3','-3','-1','-2','-2','-3','-3','3','1','-2','1','-1','-2','-2','0','-3','-1','4','-3','-2','-1','-4'],
        ['b','-2','-1','3','4','-3','0','1','-1','0','-3','-4','0','-3','-3','-2','0','-1','-4','-3','-3','4','1','-1','-4'],
        ['z','-1','0','0','1','-3','3','4','-2','0','-3','-3','1','-1','-3','-1','0','-1','-3','-2','-2','1','4','-1','-4'],
        ['x','0','-1','-1','-1','-2','-1','-1','-1','-1','-1','-1','-1','-1','-1','-2','0','0','-2','-1','-1','-1','-1','-1','-4'],
        ['-','-4','-4','-4','-4','-4','-4','-4','-4','-4','-4','-4','-4','-4','-4','-4','-4','-4','-4','-4','-4','-4','-4','-4','1']
    ]
    # other sub-matrices here

    def __init__(self, matrix_list: list):
        self.__matrix_dict = dict()
        for i in range(1, len(matrix_list)):
            for j in range(1, len(matrix_list[0])):
                self.__matrix_dict[(matrix_list[i][0], matrix_list[0][j])] = matrix_list[i][j]

    @property
    def get_dictionary(self):
        return self.__matrix_dict


class PWA(object):
    """
    Pairwise Alignment of FASTA sequences using
    the Needleman Wunsch algorithm

    supports nucleotide or peptide comparison
    """

    def __init__(self, alignment_type: AlignmentType, fasta: str):
        """
        Assuming all inputs are valid FASTA-string
        and using linux \n EOL
        :param alignment_type: the Type (Enum) of alignment
        :param seq1: FASTA sequences
        """
        self.alignment_type = alignment_type
        self.fasta_string = fasta.replace(" ", "")
        self.fasta_list = self.__to_list(fasta)
        self.traced_path = None
        self.aligned_sequences = None
        self.__matrix = [[" ", " "],
                       [" "]]

        if len(self.fasta_list) < 2:
            raise ValueError("Need at least 2 FASTA sequences")

    def do_alignment(self, match=1, mismatch=-1, gap=-2, sub_matrix=None):
        """
        performs pairwise alignment of two sequences and the traced path string at self.traced_path
        :param match: positive score for matching bit
        :param mismatch: negative score for mismatching bit
        :param gap: negative score for gap-penalty
        :param sub_matrix: substitution matrix to use. Default is BLOSUM62
        :return: void
        """
        if self.alignment_type == AlignmentType.NUCLEOTIDE:
            self.traced_path = self.__align_nuc(match, mismatch, gap, sub_matrix)
            self.aligned_sequences = self.__insert_indel()

        elif self.alignment_type == AlignmentType.PROTEIN:
            self.traced_path = self.__align_pro(match, mismatch, gap, sub_matrix)
        else:
            raise AttributeError("non-supported alignment Type. Supports: ", [t for t in AlignmentType])

    def get_best_score(self):
        end = self.__matrix[-1][-1]
        return end[0]

    def print_aligned_seq(self):
        if self.aligned_sequences:
            print(self.aligned_sequences[0])
            print(self.aligned_sequences[1])

    def print_pretty(self, separator=' '):
        print("rows: ", len(self.__matrix), "\ncolumns: ", len(self.__matrix[0]), "\n")
        for row in self.__matrix:
            if isinstance(row[0], str):
                print(separator.join(str(c) for c in row))

    def __align_nuc(self, match, mismatch, gap, sub_matrix):
        # copy the first 2 sequences to the matrix table
        seq1 = self.__strip_fasta_comment(self.fasta_list[0])
        seq2 = [list(c) for c in list(self.__strip_fasta_comment(self.fasta_list[1]))]
        self.__matrix[0].extend(list(seq1))
        self.__matrix.extend(seq2)

        # initialize gap penalty
        self.__matrix[1].append((0, 'd')) # 0 at row 2 column 2
        for i in range(1, len(seq1)+1):
            self.__matrix[1].append((i * gap, 'w'))
        for i in range(len(seq2)):
            self.__matrix[i+2].append(((i + 1) * gap, 'n'))

        # the algorithm part
        for row in range(len(seq2)):
            for col in range(len(seq1)):
                # diagonal (northwest)
                diag = self.__matrix[row+1][col+1]
                d = (diag[0] + (match if (self.__matrix[0][col+2] == self.__matrix[row+2][0]) else mismatch), 'd')

                # north and west
                north = self.__matrix[row+1][col+2]
                n = (north[0] + gap, 'n')

                west = (self.__matrix[row+2])[col+1]
                w = (west[0] + gap, 'w')

                # comparison
                m = d[0]
                if n[0] > m:
                    m = n[0]
                if w[0] > m:
                    m = w[0]

                best = filter(lambda tup: tup[0] == m, [d, n, w])
                arrow = ""
                for a in best:
                    arrow += a[1]

                self.__matrix[row+2].append((m, arrow))
        # trace back
        path = self.__trace_back()
        return path

    def __align_pro(self, match, mismatch, gap, matrix):

        return ''

    def __trace_back(self):
        row = len(self.__matrix) - 1
        col = len(self.__matrix[0]) - 1

        ret = ""

        while row > 1 or col > 1:
            score, path = self.__matrix[row][col]
            p = path[0]
            if p is 'd':
                row -= 1
                col -= 1
            elif p is 'n':
                row -= 1
            elif p is 'w':
                col -= 1
            else:
                raise RuntimeError("No path found. CHECK YOUR PWA ALGORITHM!!")
            ret = p + ret

        return ret

    def __insert_indel(self):
        """
        generate aligned sequences with indel insertion.
        by default uses the object's generated path
        :param path: a path different from the one generated by the algorithm
        :return: a list of aligned sequences by descending FASTA text
        """
        seq1 = ""
        seq2 = ""

        i = len(self.__matrix[0]) - 1
        j = len(self.__matrix) - 1
        for p in list(self.traced_path[::-1]):
            top = self.__matrix[0][i]
            bot = self.__matrix[j][0]
            if p is 'd':
                seq1 = top + seq1
                seq2 = bot + seq2
                i -= 1
                j -= 1
            elif p is 'n':
                seq1 = '-' + seq1
                seq2 = bot + seq2
                j -= 1
            elif p is 'w':
                seq1 = top + seq1
                seq2 = '-' + seq2
                i -= 1
            else:
                raise RuntimeError("traced path string contains invalid character. CHECK YOUR PWA ALGORITHM!!")

        return [seq1, seq2]

    def __to_list(self, seq: str):
        fastas = list()
        buffer = ""
        for line in seq.replace('\\n', '\n').splitlines():
            if line.startswith('>'):
                buffer = line
            else:
                fastas.append(buffer + "\n" + line)
        return fastas

    def __get_fasta_comment(self, seq: str):
        """
        find and return the FASTA comment from a single FASTA sequence
        :param seq: sequence
        :return: FASTA comment with a leading > and without the sequence
        """
        if seq.startswith('>'):
            return seq[:seq.find("\n")]
        else:
            return '>'

    def __strip_fasta_comment(self, seq: str):
        """
        remove the leading FASTA comment from a single FASTA sequence string
        :param seq: input sequence
        :return: sequence with FASTA comment removed
        """
        return seq.lstrip(self.__get_fasta_comment(seq)).strip()


if __name__ == '__main__':
    # sequences = input("FASTA:\n")
    # pwa = PWA(AlignmentType.NUCLEOTIDE, sequences)
    # pwa.do_alignment()

    # pwa.print_pretty()
    # pwa.print_aligned_seq()

    """ 
    # check score 
    
    print("score: ", pwa.get_best_score())
    
    score = 0
    for i in range(len(pwa.aligned_sequences[0])):
        top = pwa.aligned_sequences[0][i]
        bot = pwa.aligned_sequences[1][i]
        if top is '-' or bot is '-':
            score += -2
        elif top is bot:
            score += 1
        elif top is not bot:
            score += -1
    print(score)
    """

    # test substitution matrix Enum
    # d = SubstitutionMatrix.BLOSUM62.get_dictionary
    # print(d['-','a'])
    # print(d['w', 'w'])