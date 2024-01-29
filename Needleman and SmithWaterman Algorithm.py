def initialize_matrix(rows, cols, gap_penalty):
    return [[0] * cols for _ in range(rows)]

def needleman_wunsch(seq1, seq2, match=1, mismatch=0, gap=-1):
    rows, cols = len(seq1) + 1, len(seq2) + 1
    matrix = initialize_matrix(rows, cols, gap)

    for i in range(1, rows):
        matrix[i][0] = matrix[i - 1][0] + gap

    for j in range(1, cols):
        matrix[0][j] = matrix[0][j - 1] + gap

    for i in range(1, rows):
        for j in range(1, cols):
            match_score = matrix[i - 1][j - 1] + (match if seq1[i - 1] == seq2[j - 1] else mismatch)
            delete_score = matrix[i - 1][j] + gap
            insert_score = matrix[i][j - 1] + gap
            matrix[i][j] = max(match_score, delete_score, insert_score)

    return matrix

def smith_waterman(seq1, seq2, match=1, mismatch=0, gap=-1):
    rows, cols = len(seq1) + 1, len(seq2) + 1
    matrix = initialize_matrix(rows, cols, gap)

    max_score = 0
    max_position = (0, 0)

    for i in range(1, rows):
        for j in range(1, cols):
            match_score = matrix[i - 1][j - 1] + (match if seq1[i - 1] == seq2[j - 1] else mismatch)
            delete_score = matrix[i - 1][j] + gap
            insert_score = matrix[i][j - 1] + gap
            matrix[i][j] = max(0, match_score, delete_score, insert_score)

            if matrix[i][j] > max_score:
                max_score = matrix[i][j]
                max_position = (i, j)

    return matrix, max_position

def print_alignment(matrix, seq1, seq2, alignment_type):
    print(f"\nAlignment Type: {alignment_type}")
    print("Score Matrix:")
    for row in matrix:
        print(row)

if __name__ == "__main__":
    sequence1 = input("Enter the first DNA sequence: ").upper()
    sequence2 = input("Enter the second DNA sequence: ").upper()

    alignment_type = input("Choose alignment type (global/local): ").lower()

    if alignment_type == "global":
        alignment_matrix = needleman_wunsch(sequence1, sequence2)
        print_alignment(alignment_matrix, sequence1, sequence2, "Needleman-Wunsch (Global)")
    elif alignment_type == "local":
        alignment_matrix, max_position = smith_waterman(sequence1, sequence2)
        print_alignment(alignment_matrix, sequence1, sequence2, "Smith-Waterman (Local)")
    else:
        print("Invalid alignment type. Please choose 'global' or 'local'.")