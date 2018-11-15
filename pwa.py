"""
cs 123a bioinformatics

ee yieng zheng, raymond hong
"""

class pwa(object):
    """
    Pairwise Alignment of FASTA sequences using
    the Needleman Wunsch algorithm

    supports nucleotide or peptide comparison
    """

    def __init__(self, seq1: str, seq2: str):
        """
        Parameters:
        seq1 (str): the first sequence. The 'row' -- left side.
        seq2 (str): The second sequence. The 'column' -- top.
        """

        self._seq1 = seq1.strip()
        self._seq2 = seq2.strip()

        self.matrix = [[" ", " "],
                       [" "    ]]

        self.matrix[0].extend(list(seq1))
        self.matrix.extend(list(list(c) for c in list(seq2)))

    def execute(self, match=1, mismatch=-1, gap=-2):
        """
        performs pairwise alignment of two sequences
        """
        pass

    def pretty(self, separator=' '):
        for row in self.matrix:
            print(separator.join(row))

if __name__ == '__main__':
    seqA = input("Seq1: ")
    seqB = input("Seq2: ")

    pwa = pwa(seqA, seqB)
    pwa.execute()
    pwa.pretty()
