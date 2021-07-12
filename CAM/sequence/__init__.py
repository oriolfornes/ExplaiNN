from Bio.Seq import Seq
import numpy as np

def one_hot_encode(seq):
    """One hot encodes a sequence."""

    # Initialize
    seq = seq.replace("A", "0")
    seq = seq.replace("C", "1")
    seq = seq.replace("G", "2")
    seq = seq.replace("T", "3")

    encoded_seq = np.zeros((4, len(seq)), dtype=float)

    for i in range(len(seq)):
        if seq[i].isdigit():
            encoded_seq[int(seq[i]), i] = 1
        else:
            # i.e. Ns
            encoded_seq[:, i] = .0

    return(encoded_seq)

def one_hot_decode(encoded_seq):
    """Reverts a sequence's one hot encoding."""

    seq = []
    code = list("ACGT")
 
    for i in encoded_seq.transpose(1, 0):
        try:
            seq.append(code[int(np.where(i == 1)[0])])
        except:
            # i.e. N?
            seq.append("N")

    return("".join(seq))

def rc_one_hot_encoding(encoded_seq):
    """Reverse complements one hot encoding for one sequence."""
    return(encoded_seq[::-1, ::-1])

def rc(seq):
    """Reverse complements one sequence."""
    return(str(Seq(seq).reverse_complement()))