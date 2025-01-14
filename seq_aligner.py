#! /usr/bin/env python3

from util import *
from readfa import readfq
from numpy import zeros
from numpy import argmin
import numpy as np
import sys


def interesting_positions(k, sequence, reference): #selects all positions that are the beginning a kmer (a contiguous seed) of the sequence
    indexes = []
    dic = {}
    #print("ref", len(sequence))
    for i in range (len(sequence) - k):
        dic[canonical(sequence[i:i+k]).upper()] = 1
    #print(dic)
    for i in range (len(reference) - k):
        errors = 0
        for j in range(k):
            if canonical(reference[i:i + k]).upper() in dic:
                indexes.append(i)

    return indexes




def cost(xc, yc): #every mistake is of cost 1, we have a credit of 4 mistakes
    if xc == yc: return 0
    else: return 1

def alignment(seq, ref): #where seq is the sequence

    distance = zeros((len(seq)+1, len(ref)+1), dtype=int)
    for j in range(1, len(ref)+1):
        distance[0, j] = 0 # we want semi global alignment, seq being the query and ref the ref, bigger, so the first line is 0s and not  distance[0, j-1] + cost('-', ref[j-1])
    for i in range(1, len(seq)+1):
        distance[i, 0] = distance[i-1, 0] + cost(seq[i-1], '-')
    for i in range(1, len(seq)+1):
        for j in range(1, len(ref)+1):
            distance[i, j] = min(distance[i-1, j-1] + cost(seq[i-1], ref[j-1]), distance[i-1, j ] + cost(seq[i-1], '-'), distance[i , j-1] + cost('-',ref[j-1]))


    return distance, min(distance[len(seq)]) # to allow 0 at the end of the query, we select the max on the last line.


def traceback(distance, x, y):

    # all possible directions, top = deletion, left = insertion, top-left = substitution
    directions = [(-1, 0, 'D'), (0, -1, 'I'), (-1, -1, 'S')]

    # nitialization of the path
    path = []

    while x != 0 or y != 0:

        min_cost = distance[x][y]

        next_move = 0 # will be the next move out of 3 possibilities, if there are multiple path i take the first one to come up
        for dx, dy, direction in directions:
            nx, ny = x + dx, y + dy
            if 0 <= nx < len(distance) and 0 <= ny < len(distance[0]):
                if distance[nx][ny] <= min_cost:
                    min_cost = distance[nx][ny]
                    next_move = (nx, ny, direction)

        if next_move:
            x, y, move = next_move
            path.append(move)

    # to get the path in correct order, we reverse it.
    path.reverse()
    return path



def argmins(line):
    mini = min(line)
    arr = []
    for i in range(len(line)):
        if line[i] == mini:
            arr.append(i)
    return arr

def compare(pos, sequence, reference):

    possible_alignments = []
    distance_array1, lev_dist1 = alignment(sequence, reference[max([0, pos - 100]): min(len(reference), pos + 101)])

    distance_array2, lev_dist2 = alignment(reverse_complement(sequence), reference[max([0, pos - 100]): min(len(reference), pos + 101)])
    # we know that a kmer matched
    distance_array = [distance_array1, distance_array2]
    lev_dist = [lev_dist1, lev_dist2]
    if min(lev_dist) / len(sequence) > 4 / 100: #if it goes further than 4% error it gets discarded.
        return []
    for i in [0,1]:
        starts = argmins(distance_array[i][-1])
        for y in starts:
            trace = traceback(distance_array[i], len(distance_array[i]) - 1, y)
            if lev_dist[i] <= lev_dist[1 - i] :
                offset = 0
                while trace[0] == "D": #calculate the offset from pas - 100
                    offset += 1
                    trace = trace[1:len(trace)]
                possible_alignments.append([sequence, max(pos - 100,0) + offset, y, trace, i, lev_dist[i]])

    return possible_alignments


def generate_reference(filename):
    with xopen(filename) as fasta:
        sum = []
        for name,seq,_ in readfq(fasta):

            sum += seq.upper().split("[^ATCGN]")
    return name, sum



def main_reading(k, out_file, ref_file, seq_file):
    reference_name, reference = generate_reference(ref_file)
    reference = reference[0]
    matches_final = []
    with xopen(seq_file) as fasta:
        matches = []
        best_score = np.inf
        for name,seq,_ in readfq(fasta):

            interesting_positions_array = interesting_positions(k, seq, reference)
            print(len(interesting_positions_array))

            for pos in interesting_positions_array:

                possible_alignments_with_best_score_for_pos =  compare(pos, seq, reference)
                for poss_align in possible_alignments_with_best_score_for_pos:

                    if poss_align[5] <= best_score: #if the edit distance is smaller or equal to the best one seen so far
                        best_score = poss_align[5]
                        matches.append ([name,  poss_align])
            best_only = [m for m in matches if m[1][5] > best_score]
            matches_final += best_only


    with open(out_file, "w") as out_file_s:
        for match in matches_final:

            seq_name, strand, start_pos, end_pos, score, transcript =  match[0], "F" if match[1][4] == 0 else "R", str(match[1][1]), str(match[1][1] + match[1][2]), str(match[1][5]), "".join(match[1][3])
            #it was not well organized, sorry :)
            output = reference_name + "\t" +seq_name + "\t" + strand + "\t" + start_pos + "\t" + end_pos + "\t" + score + "\t" + transcript + "\n"
            out_file_s.write(output)



def main():
    args = sys.argv[1:]
    if len(args) < 3:
        print("Usage: python seq_aligner.py <reference_file> <reads_file> <k>")
        sys.exit(1)
    reference_file = args[0]
    reads_file = args[1]
    k = int(args[2])
    main_reading(k, "output.txt", reference_file, reads_file)

main()
