
import re
import sys

#function to read in the substitution matrix from file.
subs = {}
def subs_matrix(subs_path): #subs_path is the substituion matrix file.
    fr = re.split('\n', open(subs_path).read())
    header = fr[0].split('\t')[1:]
    for subsline in fr[1:-1]:
        cols = subsline.strip().split('\t')
        aa_row = cols[0]
        idx = 1
        for aa_col in header:
            subs.update({'%s_%s' %(aa_row,aa_col):int(cols[idx])})#saving values from the file in subs matrix.
            idx += 1
# return the appropriate value out of three values i.e. gap, match or mismatch.

def match_score(firstchar, secondchar,gap_penalty ):
    if firstchar == '-' or secondchar == '-':
        return gap_penalty
    else:
        return subs['%s_%s' %(firstchar,secondchar)]


def needleWunsch(seq1, seq2,subs_path,gap_penalty ):
    subs_matrix(subs_path)
    m, n = len(seq1), len(seq2)  # length of two sequences to create Dynamic Programming table
    
    # Dynamic programming table generation
    score = zeros((m+1, n+1))      # the Dynamic Programming table initialization
   
    # Dynamic programming table calculations
    for i in range(0, m + 1):
        score[i][0] = gap_penalty * i
    for j in range(0, n + 1):
        score[0][j] = gap_penalty * j
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            match = score[i - 1][j - 1] + match_score(seq1[i-1], seq2[j-1],gap_penalty )#call to match_score function to get the value from substituion matrix
            delete = score[i - 1][j] + gap_penalty #add gap penalty for a deleted char
            insert = score[i][j - 1] + gap_penalty # add gap penalty for a insertion.
            score[i][j] = max(match, delete, insert)# choose maximum of the three values.

    # Traceback and compute the alignment 
    firstsequence, secondsequence = '', ''
    i,j = m,n # traceback starts from the end of the dynamic programming table
    while i > 0 and j > 0: # loop runs until reaches right top corner of the table.
        recent_score = score[i][j]
        match_mismatch_score = score[i-1][j-1]
        vertical_score = score[i][j-1]
        horizontal_score = score[i-1][j]
        # first if checks for the match/mismatch
        if recent_score == match_mismatch_score + match_score(seq1[i-1], seq2[j-1],gap_penalty ):
            firstsequence += seq1[i-1]
            secondsequence += seq2[j-1]
            i -= 1
            j -= 1
        #this conditon is for the horizontal left gap
        elif recent_score == horizontal_score + gap_penalty:
            firstsequence += seq1[i-1]
            secondsequence += '-'
            i -= 1
        #this condtion if for the vertical gap.
        elif recent_score == vertical_score + gap_penalty:
            firstsequence += '-'
            secondsequence += seq2[j-1]
            j -= 1

    # Finish tracing up to the top left cell
    while i > 0:
        firstsequence += seq1[i-1]
        secondsequence += '-'
        i -= 1
    while j > 0:
        firstsequence += '-'
        secondsequence += seq2[j-1]
        j -= 1

    reverse_seq(firstsequence, secondsequence,gap_penalty )
def reverse_seq(firstsequence, secondsequence,gap_penalty ):
    fout = open('alignment.txt','w')
    firstsequence = firstsequence[::-1]    #reverse sequence 1 to print the sequence in original order
    secondsequence = secondsequence[::-1]    #reverse sequence 2 to print the seqence in original order.
    
    i,j = 0,0
    
    #score calculation and aligned sequences generation
    Nucleotide = '' #initialization and defination of variables used for calculation.
    found = 0
    score = 0
    
    for i in range(0,len(firstsequence)):
        # condition for a match
        if firstsequence[i] == secondsequence[i]:
            Nucleotide = Nucleotide + firstsequence[i]
            score += match_score(firstsequence[i], secondsequence[i],gap_penalty )#call to the function defined to choose appropriate value out of three possible scores.
        
        # if condition for a mismatch
        elif firstsequence[i] != secondsequence[i] and firstsequence[i] != '-' and secondsequence[i] != '-': 
            score += match_score(firstsequence[i], secondsequence[i],gap_penalty )
            Nucleotide += ' '
            found = 0
        
        #if condtion for a gap
        elif firstsequence[i] == '-' or secondsequence[i] == '-':
            Nucleotide += ' '
            score += gap_penalty
    
    fout.write( 'The optimal alignment between given sequences has score %s\n\n' % score)
    # to write the sequence in fasta format with 60 nucleotides in each line.
    start = 0
    end = 60
    noOfLines = 0
    if len(firstsequence)%60 > 0:
        noOfLines = len(firstsequence)/60 + 1
    else:
        noOfLines = len(firstsequence)/60
    
    idx_num = 1
    while idx_num <= noOfLines:
        fout.write(firstsequence[start:end] + '\n')
        fout.write(secondsequence[start:end] + '\n\n')
        start = end
        end += 60
        idx_num += 1
    fout.close()


#initialization of table
def zeros(shape):
    returnvalue = []
    for x in range(shape[0]):
        returnvalue.append([])
        for y in range(shape[1]):
            returnvalue[-1].append(0)
    return returnvalue
