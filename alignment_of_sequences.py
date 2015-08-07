"""
Following code illustrates measuring the similarity between two sequences of characters.
Given an alphabet A and a scoring matrix M defined over (A union "-") ,
the dynamic programming method computed a score that measured the similarity of two sequences X and Y based on the values of this scoring matrix.
In particular, this method involved computing an alignment matrix S between X and Y whose entry S[i,j]
scored the similarity of the substrings X[0,...,i-i], Y[0,...,j-1]

Problem: Find the Global and Local Pairwise Alignment of two sequences:

    Aim is to find optimal globa/local alignment of two sequences, X,Y. (the most similar alignment as possible)
    ASSUMPTIONS: we assume that both sequences are above same alphabet and for simlicity we assume that alphabet and sequences
                 contains upper letters and "-"   
    CONDITIONS: "-" character can be mapped only to letter ("-" to "-" is not allowed).
                Optimal alignment is alignment with maximum score.
                Both sequences for global and local alignment have to be same length.

To determine score of alignment we use scoring matrix.
Because we want to find the most similar alignment thus same letters mapping has the highest score and
letter mapped to "-" ( ot "-" to letter) has the lowest score.
    
illustration of global allignment:
Let's have alphabet: {A,C,T,G}
and two sequences X= ACCT, Y = TACGGT
example of global alignment of X and Y are sequences X', Y' :
X' = -AC-CT, Y' = TACGGT


illustration of local allignment:
Let's have alphabet: {A,C,T,G}
and two sequences X= ACCT, Y = GA__TGATACCTG_ATGGAA__T
example of local alignment of X and Y are sequences X', Y' :
X' = ACCT, Y' = ACCT

I wrote code below as part of assignment in course Algorithmic Thinking 2 (available on coursera.org)

"""
def build_scoring_matrix(alphabet,diag_score, off_diag_score, dash_score):
    """
    function builds scoring matrix as dictionary of dictionaries
    INPUT: alphabet - set of letters of alphabet (doesn't contain "-")
            diag_score - score for same letters mapping (on diagonal in matrix)
            off_diag_score - score for different letters mapping (above and below main diagonal in matrix)
            dash_score - score for "-" and letter mapping
    OUTPUT: scoring matrix which is represent as dictionary of dictionaries
            contains record for "-" to "-" mapping despite of fact that it wont be used
    """
    scoring_matrix = dict()
    #add "-" to alphabet
    alphabet_ = alphabet.union("-")

    for letter1 in alphabet_:
        #for every letter in alphabet build appropriate dictionary of scores
        col_dict = dict()
        
        for letter2 in alphabet_:
            if letter1 == "-" or letter2 == "-":
                score = dash_score
            elif letter1 == letter2:
                score = diag_score
            else: 
                score = off_diag_score

            col_dict[letter2] = score
            
        scoring_matrix[letter1] = col_dict

    return scoring_matrix

def compute_alignment_matrix(seq_x,seq_y,scoring_matrix,global_flag):
    """
    function computes alignment matrix S whose entries S[i,j] are the maximum scores over all possible alignments
            for the pair of sequences  seq_x[0,...,i-1] and seq_y[0,...,j-1].
            Matrix is represented as list of list. 
    INPUT: seq_x, seq_y - strings whose elements share a common alphabet with the scoring matrix 
            scoring_matrix - scoring matrix for alignments represented as dictionary of dictionaries
            global_flag - boolean, if True ten global alignment is computed, otherwise local alignment is computed
    OUTPUT: list of list represents alignment matrix of global/local alignment.
            Let |seq_x| = m and |seq_y| = n, then its dimension is (m+1)x(n+1), where rows have indices 0,...,m and columns
            have indices 0,...,n
            If global alignment is computed entry S[m,n] contains score of optimal alignment.
            If local alignment is computed entry with maximal value represent score of local alignment.
    """
    
    rows = len(seq_x)
    cols = len(seq_y)
    #if sequences are empty return [[0]]
    if rows == 0 and cols == 0:
        return [[0]]
    
    #initialize of alignment matrix and other variables
    alignment_matrix = [[ 0 for col in range(cols+1)] for row in range(rows+1)]
    value = 0
    
    for row in range(rows+1):
        for col in range(cols+1):
            #for every entry its value is computed 
            if row == 0 and col == 0:
                #entry [0,0]
                alignment_matrix[row][col] = 0
            elif row == 0:
                #entry [0,j] is computed based on values [0,j-1] and score of ("-" and seq_y[j]) 
                value = alignment_matrix[row][col-1] + scoring_matrix["-"][seq_y[col-1]]
            elif col == 0:
                #entry [i,0] is computed based on values [i-1,0] and score of (seq_x[i] and "-")
                value = alignment_matrix[row-1][col] + scoring_matrix[seq_x[row-1]]["-"]
            else:
                #entry [i,j] is computed based of [i-1,j-1],[i,j-1],[i-1,j] as maximum of values
                val1 = alignment_matrix[row-1][col-1] + scoring_matrix[seq_x[row-1]][seq_y[col-1]]
                val2 = alignment_matrix[row-1][col] + scoring_matrix[seq_x[row-1]]["-"]
                val3 = alignment_matrix[row][col-1] + scoring_matrix["-"][seq_y[col-1]]

                value = max(val1,val2,val3)
           
            if not global_flag:
                #for local alignment negative score is replaced with 0
                value = max(value,0)
                
            alignment_matrix[row][col] = value    

    return alignment_matrix           

def compute_global_alignment(seq_x,seq_y,scoring_matrix,alignment_matrix):
    """
    function computes a global alignment of seq_x and seq_y using the global alignment matrix alignment_matrix.
    INPUT: seq_x,seq_y - string sequences whose elements share a common alphabet with the scoring matrix 
           scoring_matrix - scoring matrix for alignments represented as dictionary of dictionaries
           alignment_matrix - computed by function compute_alignment_matrix()
    OUTPUT: tuple of the form (score, result_seq_x, result_seq_y)
            where score is the score of the global alignment result_seq_x and result_seq_y.
            Note that result_seq_x and result_seq_y
            should have the same length and may include the padding character '-'
    """
    #initialization of start position as bottom-right corner of matrix
    x_pos = len(seq_x)
    y_pos = len(seq_y)

    #initialization of variables
    result_seq_x = ''
    result_seq_y = ''
    score = alignment_matrix[x_pos][y_pos]

    #start in bottom right corner of matrix and go upwards till we reach left or upper edge
    #in every iteration we reconstruct alignments based on value in alignment_matrix and scoring_matrix
    while x_pos != 0 or y_pos !=0:
        current_value = alignment_matrix[x_pos][y_pos]
        
        if current_value == alignment_matrix[x_pos-1][y_pos-1] + scoring_matrix[seq_x[x_pos-1]][seq_y[y_pos-1]] and x_pos > 0 and y_pos > 0:
            result_seq_x = seq_x[x_pos-1] + result_seq_x
            result_seq_y = seq_y[y_pos-1] + result_seq_y
            x_pos -= 1
            y_pos -= 1
        elif current_value == alignment_matrix[x_pos-1][y_pos] + scoring_matrix[seq_x[x_pos-1]]["-"]:
            result_seq_x = seq_x[x_pos-1] + result_seq_x
            result_seq_y = "-" + result_seq_y
            x_pos -= 1
        else:    
            result_seq_x = "-" + result_seq_x
            result_seq_y = seq_y[y_pos-1] + result_seq_y
            y_pos -= 1

    return (score,result_seq_x,result_seq_y)

def compute_local_alignment(seq_x,seq_y,scoring_matrix,alignment_matrix):
    """
    function computes a local alignment of seq_x and seq_y using the local alignment matrix alignment_matrix.
    INPUT: seq_x,seq_y - string sequences whose elements share a common alphabet with the scoring matrix 
           scoring_matrix - scoring matrix for alignments represented as dictionary of dictionaries
           alignment_matrix - computed by function compute_alignment_matrix()
    OUTPUT: tuple of the form (score, result_seq_x, result_seq_y)
            where score is the score of the local alignment result_seq_x and result_seq_y.
            Note that result_seq_x and result_seq_y should have the same length and
            may include the padding character '-'
    """
    #initialization of variables
    x_pos = -1
    y_pos = -1
    result_seq_x = ''
    result_seq_y = ''
    score = 0

    #determine start position in alignment_matrix as position with maximum value 
    for row in range(len(seq_x) + 1):
        for col in range(len(seq_y) + 1):
            if alignment_matrix[row][col] > score:
                score = alignment_matrix[row][col]
                x_pos = row
                y_pos = col

    #start in start position and go upwards till we reach first entry with value 0
    #in every iteration we reconstruct alignments based on value in alignment_matrix and scoring_matrix
    while x_pos != 0 and y_pos !=0:
        current_value = alignment_matrix[x_pos][y_pos]
        if current_value == 0:
            break
        
        if current_value == alignment_matrix[x_pos-1][y_pos-1] + scoring_matrix[seq_x[x_pos-1]][seq_y[y_pos-1]]:
            result_seq_x = seq_x[x_pos-1] + result_seq_x
            result_seq_y = seq_y[y_pos-1] + result_seq_y
            x_pos -= 1
            y_pos -= 1
        elif current_value == alignment_matrix[x_pos-1][y_pos] + scoring_matrix[seq_x[x_pos-1]]["-"]:
            result_seq_x = seq_x[x_pos-1] + result_seq_x
            result_seq_y = "-" + result_seq_y
            x_pos -= 1
        else:    
            result_seq_x = "-" + result_seq_x
            result_seq_y = seq_y[y_pos-1] + result_seq_y
            y_pos -= 1

    return (score,result_seq_x,result_seq_y)



