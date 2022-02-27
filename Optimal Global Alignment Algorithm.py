import pandas as pd
import numpy as np

def optimal_global_alignment(path_fasta_1, path_fasta_2, path_score_mat, gap_penalty = -200):
    score = None
    aligned_1 = ""
    aligned_2 = ""
    
    # Read in data from path files
    with open(path_fasta_1) as f1:
        fasta_1 = f1.readlines()
    
    with open(path_fasta_2) as f2:
        fasta_2 = f2.readlines()
    
    with open(path_score_mat) as sm:
        score_mat = sm.read().split()
    
    # Assign inputted sequences as strings
    sequence_1 = str(fasta_1[1]).strip()
    sequence_2 = str(fasta_2[1]).strip()
    
    # Create a dictionary to determine match/mismatch score
    scores_dict = {
        "AA": score_mat[5],
        "AC": score_mat[6],
        "AT": score_mat[7],
        "AG": score_mat[8],
        "CA": score_mat[10],
        "CC": score_mat[11],
        "CT": score_mat[12],
        "CG": score_mat[13],
        "TA": score_mat[15],
        "TC": score_mat[16],
        "TT": score_mat[17],
        "TG": score_mat[18],
        "GA": score_mat[20],
        "GC": score_mat[21],
        "GT": score_mat[22],
        "GG": score_mat[23],
    }

    # Create matrix for optimal scores
    length_1 = len(sequence_1)
    length_2 = len(sequence_2)
        
    optimal_scores = np.zeros([length_1 + 1, length_2 + 1], dtype = float)
    
    optimal_scores[0:(length_1 + 1), 0] = [i * gap_penalty for i in range(length_1 + 1)]
    optimal_scores[0, 0:(length_2 + 1)] = [i * gap_penalty for i in range(length_2 + 1)]
    
    # Assign optimal score matrix values based on match, mismatch, and gaps
    i = 1
    while i <= length_1:
        j = 1
        while j <= length_2:
            if sequence_1[i - 1] == sequence_2[j - 1]:
                match = float(scores_dict[str(str(sequence_1[i - 1]) + str(sequence_2[j - 1]))])
                optimal_scores[i][j] = optimal_scores[i - 1][j - 1] + match
            else:
                mismatch = float(scores_dict[str(str(sequence_1[i - 1]) + str(sequence_2[j - 1]))])
                optimal_scores[i][j] = max(optimal_scores[i - 1][j - 1] + mismatch,
                            optimal_scores[i - 1][j] + gap_penalty,
                            optimal_scores[i][j - 1] + gap_penalty)
            j += 1
        i += 1
    
    # Assign score value
    score = float(optimal_scores[length_1][length_2])

    # Create array to store char of unicodes and corresponding index variables
    sq1 = np.zeros(length_1 + length_2 + 1, dtype = int)
    sq2 = np.zeros(length_1 + length_2 + 1, dtype = int)
     
    sq1_position = length_1 + length_2
    sq2_position = length_1 + length_2
    
    i = length_1
    j = length_2
    
    # Determine gap positions and put in correct location
    while not (i == 0 or j == 0):
        if sequence_1[i - 1] == sequence_2[j - 1]:       
            sq1[sq1_position] = ord(sequence_1[i - 1])
            sq2[sq2_position] = ord(sequence_2[j - 1])
            sq1_position -= 1
            sq2_position -= 1
            i -= 1
            j -= 1
        
        elif (optimal_scores[i - 1][j - 1] + float(scores_dict[str(str(sequence_1[i - 1]) + 
                                                                 str(sequence_2[j - 1]))])) == optimal_scores[i][j]:
            sq1[sq1_position] = ord(sequence_1[i - 1])
            sq2[sq2_position] = ord(sequence_2[j - 1])
            sq1_position -= 1
            sq2_position -= 1
            i -= 1
            j -= 1
         
        elif (optimal_scores[i - 1][j] + gap_penalty) == optimal_scores[i][j]:
            sq1[sq1_position] = ord(sequence_1[i - 1])
            sq2[sq2_position] = ord('-')
            sq1_position -= 1
            sq2_position -= 1
            i -= 1
         
        elif (optimal_scores[i][j - 1] + gap_penalty) == optimal_scores[i][j]:       
            sq1[sq1_position] = ord('-')
            sq2[sq2_position] = ord(sequence_2[j - 1])
            sq1_position -= 1
            sq2_position -= 1
            j -= 1
         
    # Finish char arrays
    while sq1_position > 0:
        if i > 0:
            i -= 1
            sq1[sq1_position] = ord(sequence_1[i])
            sq1_position -= 1
        else:
            sq1[sq1_position] = ord('-')
            sq1_position -= 1
     
    while sq2_position > 0:
        if j > 0:
            j -= 1
            sq2[sq2_position] = ord(sequence_2[j])
            sq2_position -= 1
        else:
            sq2[sq2_position] = ord('-')
            sq2_position -= 1
 
    # Create aligned strings for output from char sequences
    temp = 1
    i = length_1 + length_2
    while i >= 1:
        if (chr(sq1[i]) == '-') and chr(sq2[i]) == '-':
            temp = i + 1
            break
        i -= 1

    i = temp
    while i <= length_1 + length_2:
        aligned_1 += chr(sq1[i])
        i += 1
         
    i = temp
    while i <= length_1 + length_2:
        aligned_2 += chr(sq2[i])
        i += 1
    
    # Print output
#     print(f"The optimal alignment score between given sequences has score {score}")
#     print(aligned_1)
#     print(aligned_2)
    
    return score, aligned_1, aligned_2

# Test function using path files
fasta_1 = r"prob4_data\data_example\seq1.fasta"
fasta_2 = r"prob4_data\data_example\seq2.fasta"
score_mat = r"prob4_data\data_example\substitution_matrix.txt"

optimal_global_alignment(fasta_1, fasta_2, score_mat)





def test_optimal():
    i = 1
    for i in range(10):
        fasta_1 = r"prob4_data\data_" + str(i+1) + "\seq1.fasta"
        fasta_2 = r"prob4_data\data_" + str(i+1) + "\seq2.fasta"
        score_mat = r"prob4_data\data_1\substitution_matrix.txt"
        
        f = open("align_data_" + str(i+1) + ".txt","w+")
        f.write("The optimal alignment score between given sequences has score " + 
                str(optimal_global_alignment(fasta_1, fasta_2, score_mat)[0]) + "\n" 
                + optimal_global_alignment(fasta_1, fasta_2, score_mat)[1] + "\n"
                + optimal_global_alignment(fasta_1, fasta_2, score_mat)[2])
        f.close()

test_optimal()