"""
cs 123a bioinformatics

ee yieng zheng, raymond hong
"""

from enum import Enum
from math import inf
from re import compile
from warnings import warn

class AlignmentMethod(Enum):
    GLOBAL = 1
    LOCAL = 2


class AlignmentType(Enum):
    NUCLEOTIDE = 1
    PROTEIN = 2


class SubstitutionMatrix(Enum):
    BLOSUM62 = [
        ['','A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V'],
        ['A',4,-1,-2,-2,0,-1,-1,0,-2,-1,-1,-1,-1,-2,-1,1,0,-3,-2,0],
        ['R',-1,5,0,-2,-3,1,0,-2,0,-3,-2,2,-1,-3,-2,-1,-1,-3,-2,-3],
        ['N',-2,0,6,1,-3,0,0,0,1,-3,-3,0,-2,-3,-2,1,0,-4,-2,-3],
        ['D',-2,-2,1,6,-3,0,2,-1,-1,-3,-4,-1,-3,-3,-1,0,-1,-4,-3,-3],
        ['C',0,-3,-3,-3,9,-3,-4,-3,-3,-1,-1,-3,-1,-2,-3,-1,-1,-2,-2,-1],
        ['Q',-1,1,0,0,-3,5,2,-2,0,-3,-2,1,0,-3,-1,0,-1,-2,-1,-2],
        ['E',-1,0,0,2,-4,2,5,-2,0,-3,-3,1,-2,-3,-1,0,-1,-3,-2,-2],
        ['G',0,-2,0,-1,-3,-2,-2,6,-2,-4,-4,-2,-3,-3,-2,0,-2,-2,-3,-3],
        ['H',-2,0,1,-1,-3,0,0,-2,8,-3,-3,-1,-2,-1,-2,-1,-2,-2,2,-3],
        ['I',-1,-3,-3,-3,-1,-3,-3,-4,-3,4,2,-3,1,0,-3,-2,-1,-3,-1,3],
        ['L',-1,-2,-3,-4,-1,-2,-3,-4,-3,2,4,-2,2,0,-3,-2,-1,-2,-1,1],
        ['K',-1,2,0,-1,-3,1,1,-2,-1,-3,-2,5,-1,-3,-1,0,-1,-3,-2,-2],
        ['M',-1,-1,-2,-3,-1,0,-2,-3,-2,1,2,-1,5,0,-2,-1,-1,-1,-1,1],
        ['F',-2,-3,-3,-3,-2,-3,-3,-3,-1,0,0,-3,0,6,-4,-2,-2,1,3,-1],
        ['P',-1,-2,-2,-1,-3,-1,-1,-2,-2,-3,-3,-1,-2,-4,7,-1,-1,-4,-3,-2],
        ['S',1,-1,1,0,-1,0,0,0,-1,-2,-2,0,-1,-2,-1,4,1,-3,-2,-2],
        ['T',0,-1,0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1,1,5,-2,-2,0],
        ['W',-3,-3,-4,-4,-2,-2,-3,-2,-2,-3,-2,-3,-1,1,-4,-3,-2,11,2,-3],
        ['Y',-2,-2,-2,-3,-2,-1,-2,-3,2,-1,-1,-2,-1,3,-3,-2,-2,2,7,-1],
        ['V',0,-3,-3,-3,-1,-2,-2,-3,-3,3,1,-2,1,-1,-2,-2,0,-3,-1,4]
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


class UnsupportedCodeError(Exception):
    def __init__(self, option='', value="Sequence must contain only Amino acids: -ARNDCQEGHILKMFPSTWYV and Nucleotides: ATUCG. "):
        self.value = value + option

    def __str__(self):
        return repr(self.value)


class PWA(object):
    """
    Pairwise Alignment of FASTA sequences using
    the Needleman Wunsch algorithm

    supports nucleotide or peptide comparison
    """

    def __init__(self, alignment_type: AlignmentType, alignment_method: AlignmentMethod, fasta: str):
        """
        Assuming all inputs are valid FASTA-string
        and using linux \n EOL
        :param alignment_type: the AlignmentType (Enum) of alignment
        :param alignment_method: the AlignmentMethod (Enum); Local or global
        :param fasta: FASTA sequences
        """
        self.alignment_type = alignment_type
        self.alignment_method = alignment_method
        self.fasta_string = self.__conform_fasta_string(fasta)
        self.fasta_list = self.__to_list()
        self.traced_path = None
        self.aligned_sequences = None
        self.__matrix = [[" ", " "],
                       [" "]]

        self.best_location = None
        self.best_score = -inf

        if len(self.fasta_list) < 2:
            raise ValueError("Need at least 2 FASTA sequences")
        try:
            a = self.alignment_type
            b = self.alignment_method.name
        except AttributeError:
            raise AttributeError("non-supported alignment type or method.\nSupports only\nType: {}\nMethod: {}".format([t for t in AlignmentType], [m for m in AlignmentMethod]))

    def do_alignment(self, match=1, mismatch=-1, gap=-2, sub_matrix=SubstitutionMatrix.BLOSUM62):
        """
        performs pairwise alignment of two sequences and the traced path string at self.traced_path
        :param match: positive score for matching bit
        :param mismatch: negative score for mismatching bit
        :param gap: negative score for gap-penalty
        :param sub_matrix: substitution matrix to use. Default is BLOSUM62
        :return: void
        """
        self.traced_path = self.__align(match, mismatch, gap, sub_matrix)
        self.aligned_sequences = self.__insert_indel()

    def get_conserved_regions(self):
        if len(self.aligned_sequences) > 2:
            return self.aligned_sequences[2]
        else:
            s = self.aligned_sequences[0]
            t = self.aligned_sequences[1]
            v = ''.join(s[i] if s[i] == t[i] else ' ' for i in range(len(s)))

            self.aligned_sequences.append(v)
            return v

    def print_pretty(self, separator=' '):
        print("rows: ", len(self.__matrix), "\ncolumns: ", len(self.__matrix[0]), "\n")
        for row in self.__matrix:
            if isinstance(row[0], str):
                print(separator.join(str(c) for c in row))
        print()

    def print_aligned_seq(self):
        if self.aligned_sequences:
            print(self.aligned_sequences[0], " ", self.__get_fasta_comment(self.fasta_list[0]))
            print(self.aligned_sequences[1], " ", self.__get_fasta_comment(self.fasta_list[1]))

    def print_aligned_seq_formatted(self, max_len=60,show_conserved=True,show_mismatch_only=False):
        formatted_list = self.get_aligned_seq_formatted(max_len, show_conserved, show_mismatch_only)
        for row in formatted_list:
            print(row)

    def get_aligned_seq_formatted(self, max_len=60,show_conserved=True,show_mismatch_only=False):
        if len(self.aligned_sequences[0]) < max_len:
            max_len = len(self.aligned_sequences[0])

        seq1 = [self.aligned_sequences[0][i:i + max_len] for i in range(0, len(self.aligned_sequences[0]), max_len)]
        seq2 = [self.aligned_sequences[1][i:i + max_len] for i in range(0, len(self.aligned_sequences[1]), max_len)]

        c_region = self.get_conserved_regions()
        if show_mismatch_only:
            c_region = c_region.replace(' ', 'x')
            c_region = (compile('[A-Z]').sub(' ', c_region))
        seqC = [c_region[i:i + max_len] for i in range(0, len(c_region), max_len)]

        width = str(len(str(len(c_region))))
        s_len_total = 0
        t_len_total = 0

        formatted_list = list()

        row = 0
        for s, t, v in list(zip(seq1, seq2, seqC)):
            s_len = len(s) - s.count('-')
            t_len = len(t) - t.count('-')
            formatted = "query" + format(s_len_total + 1, ' <' + width + 'd') + s + str(s_len) + str(s_len_total) + '\n'
            formatted += "sbjct" + format(t_len_total + 1, ' <' + width + 'd') + t + str(t_len) + str(t_len_total) + '\n'
            if show_conserved:
                formatted += "     " + format('', ' <' + width + 's') + v + '\n'

            row += 1
            s_len_total += s_len
            t_len_total += t_len

            formatted_list.append(formatted)

        return formatted_list

    def get_final_score(self):
        end = self.__matrix[-1][-1]
        return end[0]

    def __init_matrix(self, gap, seq1, seq2):
        if self.alignment_method.name is AlignmentMethod.LOCAL.name:
            gap = 0

        # copy the first 2 sequences to the matrix table
        self.__matrix[0].extend(list(seq1))
        self.__matrix.extend(seq2)

        # initialize gap penalty
        self.__matrix[1].append((0, '')) # 0 at row 2 column 2
        for i in range(1, len(seq1)+1):
            self.__matrix[1].append((i * gap, 'w'))
        for i in range(len(seq2)):
            self.__matrix[i+2].append(((i + 1) * gap, 'n'))

    def __align(self, match, mismatch, gap, sub_matrix):
        seq1 = self.__strip_fasta_comment(self.fasta_list[0])
        seq2 = [list(c) for c in list(self.__strip_fasta_comment(self.fasta_list[1]))]

        # copy the first 2 sequences to the matrix table
        self.__init_matrix(gap, seq1, seq2)

        # the algorithm part
        self.__fill_matrix(match, mismatch, gap, sub_matrix, len(seq1), len(seq2))

        # trace back
        path = self.__trace_back()
        return path

    def __fill_matrix(self, match, mismatch, gap, sub_matrix, len_seq1, len_seq2):
        sub_matrix_dict = sub_matrix.get_dictionary if sub_matrix else None

        for row in range(len_seq2):
            for col in range(len_seq1):
                # diagonal (northwest)
                diag = self.__matrix[row+1][col+1]
                # north and west
                north = self.__matrix[row+1][col+2]
                west = (self.__matrix[row+2])[col+1]

                d = n = w = (0, '')
                if self.alignment_type.name is AlignmentType.NUCLEOTIDE.name:
                    d = (diag[0] + (match if (self.__matrix[0][col+2] == self.__matrix[row+2][0]) else mismatch), 'd')
                    n = (north[0] + gap, 'n')
                    w = (west[0] + gap, 'w')
                elif self.alignment_type.name is AlignmentType.PROTEIN.name:
                    key = (self.__matrix[0][col + 2], self.__matrix[row + 2][0])
                    sub = sub_matrix_dict[key] if key in sub_matrix_dict else gap
                    d = (diag[0] + sub, 'd')
                    n = (north[0] + gap, 'n')
                    w = (west[0] + gap, 'w')

                if not d[1] or not n[1] or not w[1]:
                    raise ValueError("the dynamic array is missing direction values.")

                # comparison
                m = d[0]
                if n[0] > m:
                    m = n[0]
                if w[0] > m:
                    m = w[0]

                best_direction = filter(lambda tup: tup[0] == m, [d, n, w])
                arrow = ""
                for a in best_direction:
                    arrow += a[1]

                self.__matrix[row+2].append((m, arrow))

                if m >= self.best_score:
                    self.best_score = m
                    self.best_location = (row + 2, col + 2)

    def __trace_back(self):
        # global
        row = len(self.__matrix) - 1
        col = len(self.__matrix[0]) - 1

        local = False
        if self.alignment_method.name is AlignmentMethod.LOCAL.name:
            row, col = self.best_location[0], self.best_location[1]
            local = True

        ret = ""
        while row > 1 or col > 1:
            score, path = self.__matrix[row][col]
            if score <= 0 and local:
                return ret

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
        seq1 = self.__strip_fasta_comment(self.fasta_list[0])
        seq2 = self.__strip_fasta_comment(self.fasta_list[1])
        new_seq1 = ""
        new_seq2 = ""

        # global
        i = len(self.__matrix[0]) - 1  # col, top
        j = len(self.__matrix) - 1     # row, bottom

        local = False
        if self.alignment_method.name is AlignmentMethod.LOCAL.name:
            local = True
            offset = 2
            i, j = i-offset, j-offset

            # copy suffix of first/second sequence
            fro, too = self.best_location[1]+1, i+1
            top_suffix = seq1[fro:too]

            fro, too = self.best_location[0]+1, j+1
            bot_suffix = seq2[fro:too]

            new_seq1 = top_suffix
            new_seq2 = bot_suffix

            # adding indels to end of sequence
            if len(new_seq1) is not len(new_seq2):
                diff = abs(len(new_seq1) - len(new_seq2))
                if len(new_seq2) < len(new_seq1):
                    new_seq2 = new_seq2 + '-' * diff
                elif len(new_seq1) < len(new_seq2):
                    new_seq1 = new_seq1 + '-' * diff

            i, j = self.best_location[1], self.best_location[0]

        for p in list(self.traced_path[::-1]):
            top = self.__matrix[0][i]
            bot = self.__matrix[j][0]

            if p is 'd':
                new_seq1 = top + new_seq1
                new_seq2 = bot + new_seq2
                i -= 1
                j -= 1
            elif p is 'n':
                new_seq1 = '-' + new_seq1
                new_seq2 = bot + new_seq2
                j -= 1
            elif p is 'w':
                new_seq1 = top + new_seq1
                new_seq2 = '-' + new_seq2
                i -= 1
            else:
                raise RuntimeError("traced path string contains invalid character. CHECK YOUR PWA ALGORITHM!!")

        if local:
            top_prefix = seq1[:i-1]
            bot_prefix = seq2[:j-1]

            new_seq1 = top_prefix + new_seq1
            new_seq2 = bot_prefix + new_seq2

            # adding indels to beginning of sequence
            if len(new_seq1) is not len(new_seq2):
                diff = abs(len(new_seq1) - len(new_seq2))
                if len(new_seq2) < len(new_seq1):
                    new_seq2 = '-' * diff + new_seq2
                elif len(new_seq1) < len(new_seq2):
                    new_seq1 = '-' * diff + new_seq1

        return [new_seq1, new_seq2]

    def __conform_fasta_string(self, raw_fasta: str):
        nucl_pat = compile("^[ATUCG]+$")
        prot_pat = compile("^[-ARNDCQEGHILKMFPSTWYV]+$")

        prot = False
        if self.alignment_type.name is AlignmentType.PROTEIN.name:
            prot = True

        buffer = ""
        lines = raw_fasta.replace("\\n", "\n").splitlines()
        for line in lines:
            line = line.strip()
            if line.startswith('>'):
                if buffer:
                    buffer += "\n"
                buffer += line + "\n"
            elif line:
                line = compile("[0-9\s]+").sub('', line).upper()
                illegal = prot_pat.sub('', line) if prot else nucl_pat.sub('', line)
                if illegal:
                    raise UnsupportedCodeError("illegal codes: " + illegal)
                else:
                    buffer += line
        return buffer

    def __to_list(self):
        fastas = self.fasta_string.split('>')
        func = lambda f: '>' + f
        fastas = list(map(func, filter(len, fastas)))
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
    alignment_method = ""
    alignment_type = ""
    sequence = ""
    while True:
        if not alignment_type:
            alignment_type = input("protein or nucleotides? p/n: ").lower()
        if not alignment_type or not alignment_type == "p" and not alignment_type == "n":
            print("incorrect type:", alignment_type)
            alignment_type = ""
            continue
        if not alignment_method:
            alignment_method = input("local or global? l/g: ").lower()
        if not alignment_method or not alignment_method == 'l' and not alignment_method == 'g':
            print("incorrect method:", alignment_method)
            alignment_method = ""
            continue

        alignment_t = AlignmentType.PROTEIN if alignment_type == 'p' else AlignmentType.NUCLEOTIDE
        alignment_m = AlignmentMethod.LOCAL if alignment_method == 'l' else AlignmentMethod.GLOBAL

        sequence = ""
        inp = input("FASTA sequences with >comment:\n")
        while inp:
            sequence += inp + '\n'
            inp = input()

        if not sequence.strip():
            print("please input correct sequences")
            continue

        print("Doing alignment on", alignment_t.name, "using method", alignment_m.name)

        pwa = PWA(alignment_t, alignment_m, sequence.strip())
        pwa.do_alignment()
        pwa.print_aligned_seq_formatted()