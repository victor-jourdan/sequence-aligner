from readfa import readfq
from xopen import xopen
import matplotlib.pyplot as plt
def fileconcat(filename):
    with xopen(filename) as fasta:
        sum = []
        for _,seq,_ in readfq(fasta):
            sum += seq.upper().split("[^ATCGN]")
            print(sum)
        return sum


RC_TABLE = str.maketrans("ACTGNactgn", "TGACNtgacn")

def reverse_complement(seq: str) -> str:
    complement = seq.translate(RC_TABLE)
    return complement[::-1]


def canonical(kmer):
    return min(kmer, reverse_complement(kmer))



    
