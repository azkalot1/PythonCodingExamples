def revcomp(dna_string):
    comp_d = {
        'A': 'T',
        'T': 'A',
        'C': 'G',
        'G': 'C',
    }
    dna_string = [comp_d[x]for x in dna_string[::-1]]
    dna_string = ''.join(dna_string)
    return(dna_string)


def rna(dna_string):
    rna_string = [x if x != 'T' else 'U' for x in dna_string]
    rna_string = ''.join(rna_string)
    return(rna_string)


def subs(dna_string, pattern):
    return([idx+1 for idx in range(len(dna_string)-len(pattern)) if dna_string[idx:idx+len(pattern)] == pattern])
    # 1-based


def findpalindrome(dna_string):
    # not the best implementation
    # sizes 4 to 12
    palindromes = []

    for i in range(len(dna_string)):
        for j in range(4, 13):
            if i + j > len(dna_string):
                continue
            s1 = dna_string[i:i+j]
            s2 = revcomp(s1)
            if s1 == s2:
                palindromes.append((i + 1, j))
    return(palindromes)


def hamm(s1, s2):
    return sum([x != y for x, y, in zip(s1, s2)])


def prot(rnastring):
    import numpy as np
    d = {
        'UUU': 'F',
        'UUC': 'F',
        'UUA': 'L',
        'UUG': 'L',
        'UCU': 'S',
        'UCC': 'S',
        'UCA': 'S',
        'UCG': 'S',
        'UAU': 'Y',
        'UAC': 'Y',
        'UAA': 'Stop',
        'UAG': 'Stop',
        'UGU': 'C',
        'UGC': 'C',
        'UGA': 'Stop',
        'UGG': 'W',
        'CUU': 'L',
        'CUC': 'L',
        'CUA': 'L',
        'CUG': 'L',
        'CCU': 'P',
        'CCC': 'P',
        'CCA': 'P',
        'CCG': 'P',
        'CAU': 'H',
        'CAC': 'H',
        'CAA': 'Q',
        'CAG': 'Q',
        'CGU': 'R',
        'CGC': 'R',
        'CGA': 'R',
        'CGG': 'R',
        'AUU': 'I',
        'AUC': 'I',
        'AUA': 'I',
        'AUG': 'M',
        'ACU': 'T',
        'ACC': 'T',
        'ACA': 'T',
        'ACG': 'T',
        'AAU': 'N',
        'AAC': 'N',
        'AAA': 'K',
        'AAG': 'K',
        'AGU': 'S',
        'AGC': 'S',
        'AGA': 'R',
        'AGG': 'R',
        'GUU': 'V',
        'GUC': 'V',
        'GUA': 'V',
        'GUG': 'V',
        'GCU': 'A',
        'GCC': 'A',
        'GCA': 'A',
        'GCG': 'A',
        'GAU': 'D',
        'GAC': 'D',
        'GAA': 'E',
        'GAG': 'E',
        'GGU': 'G',
        'GGC': 'G',
        'GGA': 'G',
        'GGG': 'G'
    }
    return ''.join([d[rnastring[idx:idx+3]] for idx in np.arange(start=0, stop=len(rnastring)-3, step=3)])


def orf(dnastring):
    # find start-end codons, extract RNA substring, pass them thru rna2prot
    d = {
        'UUU': 'F',
        'UUC': 'F',
        'UUA': 'L',
        'UUG': 'L',
        'UCU': 'S',
        'UCC': 'S',
        'UCA': 'S',
        'UCG': 'S',
        'UAU': 'Y',
        'UAC': 'Y',
        'UAA': 'Stop',
        'UAG': 'Stop',
        'UGU': 'C',
        'UGC': 'C',
        'UGA': 'Stop',
        'UGG': 'W',
        'CUU': 'L',
        'CUC': 'L',
        'CUA': 'L',
        'CUG': 'L',
        'CCU': 'P',
        'CCC': 'P',
        'CCA': 'P',
        'CCG': 'P',
        'CAU': 'H',
        'CAC': 'H',
        'CAA': 'Q',
        'CAG': 'Q',
        'CGU': 'R',
        'CGC': 'R',
        'CGA': 'R',
        'CGG': 'R',
        'AUU': 'I',
        'AUC': 'I',
        'AUA': 'I',
        'AUG': 'M',
        'ACU': 'T',
        'ACC': 'T',
        'ACA': 'T',
        'ACG': 'T',
        'AAU': 'N',
        'AAC': 'N',
        'AAA': 'K',
        'AAG': 'K',
        'AGU': 'S',
        'AGC': 'S',
        'AGA': 'R',
        'AGG': 'R',
        'GUU': 'V',
        'GUC': 'V',
        'GUA': 'V',
        'GUG': 'V',
        'GCU': 'A',
        'GCC': 'A',
        'GCA': 'A',
        'GCG': 'A',
        'GAU': 'D',
        'GAC': 'D',
        'GAA': 'E',
        'GAG': 'E',
        'GGU': 'G',
        'GGC': 'G',
        'GGA': 'G',
        'GGG': 'G'
    }

    prot_orf = []
    rnastring = rna(dnastring)

    for idx_orf in [0, 1, 2]:
        protein_str = ''
        for i in range(idx_orf, len(rnastring)-2, 3):
            protein_str = protein_str + str(d[rnastring[i:i+3]])
        for i in range(len(protein_str)):
            if protein_str[i] == 'M' and protein_str.find('Stop', i) != -1:
                prot_orf.append(protein_str[i:protein_str.find('Stop', i)])

    rnastring = rna(revcomp(dnastring))
    for idx_orf in [0, 1, 2]:
        protein_str = ''
        for i in range(idx_orf, len(rnastring)-2, 3):
            protein_str = protein_str + str(d[rnastring[i:i+3]])
        for i in range(len(protein_str)):
            if protein_str[i] == 'M' and protein_str.find('Stop', i) != -1:
                prot_orf.append(protein_str[i:protein_str.find('Stop', i)])

    return(set(prot_orf))


def recur_factorial(n):
    if n == 1:
        return 1
    if n == 0:
        return 1
    else:
        return n*recur_factorial(n-1)


def perm(n):
    print(recur_factorial(n))
    from itertools import permutations
    for permutation in permutations(list(range(n))):
        print(' '.join([str(x+1) for x in permutation]))


def comb(n, k):
    return recur_factorial(n) / (recur_factorial(n-k) * recur_factorial(k))


def sign(n):
    # it will be \sum_k=0^n choice(n, k) * n!
    from itertools import combinations
    from itertools import permutations
    import numpy as np
    combs = [comb(n, k) * recur_factorial(n) for k in range(n+1)]
    print(int(sum(combs)))
    for k in range(n+1):  # how many with -1
        for x in combinations(range(n), k):
            list_comb = np.array(range(n)) + 1
            list_comb[list(x)] *= -1
            for permutation in permutations(list_comb):
                print(' '.join([str(perm) for perm in permutation]))


def lcsm_naive(dna_strings):
    # we generate all possible substring from
    # the first string
    # and find if they are present in the others
    from Bio import SeqIO
    records = list(SeqIO.parse(dna_strings, "fasta"))
    print([record.seq for record in records])
    lcsm = ''
    for idx_lo in range(len(records[0].seq)-1):
        idx_hi = idx_lo + 1
        commong_substr_is_growing = True
        while idx_hi < len(records[0].seq) and commong_substr_is_growing:
            # there is no need to check for patterns which 
            # are shorter than existing pattern
            # also if we already know that we don't have 
            # a match we can move idx_lo up
            pattern = records[0].seq[idx_lo:idx_hi]
            if len(pattern) > len(lcsm):
                matching_n = sum([len(subs(record.seq, pattern)) > 0 for record in records[1:]])
                if matching_n == len(records[1:]):
                    lcsm = pattern
                else:
                    commong_substr_is_growing = False
            idx_hi += 1
    print(lcsm)
    return lcsm

def lcsm(dna_strings):
    # a faster verstion
    # 1. take the shortest seq and generate motifs using that
    # 2. generate all possible motifs from that shortes sequence
    # 3. iterate over sequences
    # 4 for m in motifs, if m not in seq - remove m
    # 5. we end up with the all shared motifs, pick longes
    def get_shortest_seq(records):
        shortest_len = 10000 # we know they are up to 1kbp
        for record in records:
            if len(record.seq) < shortest_len:
                shortest_len = len(record)
                shortest_seq = record.seq
        return shortest_seq
    
    from Bio import SeqIO
    records = list(SeqIO.parse(dna_strings, "fasta"))   
    # get shortest seq
    shortest_seq = get_shortest_seq(records)
    # now generate all possible motifs
    possible_motifs = set()
    for idx_lo in range(len(shortest_seq)-1):
        for idx_hi in range(idx_lo+1, len(shortest_seq)):
            subseq = shortest_seq[idx_lo:idx_hi]
            possible_motifs.add(subseq)
    # now iterate over sequences
    for record in records:
        update_possible_motifs = list(possible_motifs)
        for motif in update_possible_motifs:
            if motif not in record.seq:
                possible_motifs.remove(motif)

    # get the longest
    n = 0
    longest_motif = ''
    for motif in possible_motifs:
        if len(motif) > n:
            longest_motif = motif
            n = len(motif)
    print(longest_motif)
    return longest_motif


def splc(dna_strings):
    import numpy as np
    # splice out all introns
    exons_only = ''
    from Bio import SeqIO
    records = list(SeqIO.parse(dna_strings, "fasta"))
    template = records[0].seq
    introns = records[1:]
    introns_locations = []
    introns_lengths = []
    for intron in introns:
        intron_pos = subs(template, intron.seq)
        for pos in intron_pos:
            pos = pos - 1
            introns_locations.append(pos)
            introns_lengths.append(len(intron.seq))

    introns_lengths = np.array(introns_lengths)
    introns_locations = np.array(introns_locations)

    introns_lengths = introns_lengths[np.argsort(introns_locations)]
    introns_locations = np.sort(introns_locations)
    print(introns_locations)
    print(introns_lengths)
    pointer = 0
    for idx in range(len(introns_locations)):
        exons_only += template[pointer:introns_locations[idx]]
        pointer = introns_locations[idx] + introns_lengths[idx]
    exons_only += template[pointer:len(template)]
    return(prot(rna(exons_only)))


def sseq(str1, str2):
    # find if str2 is a subseq of str1
    m = len(str2)
    n = len(str1)

    j = 0
    i = 0
    pos = []
    while j < m and i < n:
        if str2[j] == str1[i]:
            j = j+1  # look for the next
            pos.append(i)
        i = i + 1

    # if j==m mean we used all characters in the str2
    # so yes, we found it
    print(''.join([str1[x] for x in pos]))
    print(str2)
    return((j == m, ' '.join([str(x+1) for x in pos])))


def lgis(n, sequence):
    def lis(n, sequence):
        lis = [1]*n
        lis_seq = [[sequence[x]] for x in range(n)]
        for i in range(1, n):
            for j in range(0, i):
                if sequence[i] > sequence[j] and lis[i] < lis[j] + 1:
                    lis[i] = lis[j] + 1
                    lis_seq[i] =  lis_seq[j] + [sequence[i]]

        maximum = 0
        max_lis_seq = []
        for i in range(n):
            if maximum < lis[i]:
                maximum = lis[i]
                max_lis_seq = lis_seq[i]

        return maximum, max_lis_seq

    def lds(n, sequence):
        lds = [1]*n
        lds_seq = [[sequence[x]] for x in range(n)] 
        for i in range(1, n):
            for j in range(0, i):
                if sequence[i] < sequence[j] and lds[i] < lds[j] + 1:
                    lds[i] = lds[j] + 1
                    lds_seq[i] =  lds_seq[j] + [sequence[i]]

        maximum = 0
        max_lds_seq = []
        for i in range(n):
            if maximum < lds[i]:
                maximum = lds[i]
                max_lds_seq = lds_seq[i]

        return maximum, max_lds_seq

    sequence = sequence.split(' ')
    sequence = [int(x) for x in sequence]
    maximum_lis, max_lis_seq = lis(n, sequence)
    maximum_lds, max_lds_seq = lds(n, sequence)
    print(' '.join([str(x) for x in max_lis_seq]))
    print(' '.join([str(x) for x in max_lds_seq]))


def lexf(alphabet, n, output_file='lexf.txt'):
    from itertools import product
    with open(output_file, 'w') as f:
        for x in list(product(alphabet, repeat=n)):
            f.write(''.join(x)+'\n')


def lcsq(X, Y):
    m = len(X)
    n = len(Y)
    L = [[0 for x in range(n+1)] for x in range(m+1)] 

    # Following steps build L[m+1][n+1] in bottom up fashion. Note 
    # that L[i][j] contains length of LCS of X[0..i-1] and Y[0..j-1]  
    for i in range(m+1):
        for j in range(n+1):
            if i == 0 or j == 0:
                L[i][j] = 0
            elif X[i-1] == Y[j-1]:
                L[i][j] = L[i-1][j-1] + 1  # go diagonal
            else:
                L[i][j] = max(L[i-1][j], L[i][j-1])  # add max from left \ top

    # Following code is used to print LCS
    index = L[m][n]

    # Create a character array to store the lcs string
    lcs = [""] * (index+1) 
    lcs[index] = ""

    # Start from the right-most-bottom-most corner and 
    # one by one store characters in lcs[] 
    i = m
    j = n
    # Now print
    while i > 0 and j > 0:
        if X[i-1] == Y[j-1]:
            lcs[index-1] = X[i-1]
            i -= 1
            j -= 1
            index -= 1
            # Goind up diagonal

        elif L[i-1][j] > L[i][j-1]:
            i -= 1
        else:
            j -= 1
        # Followind the max

    print("".join(lcs))


if __name__ == '__main__':
    lcsq(
        'GGGACTCGCCCCAGATAATATATCGTAGCGCTCGACCTAGCAAAGATCAAGTGATTTCCGGCAAGCTGCACTTTGTACAGGGCTTCAGATCCCCCCCAGTAGTTTTCAACGAGCTCTCGCATTCAACCCCGTCCGGTAAATCCGGAGGTCATCATTTTTATAACTGAAGAAGCGGCACGATACCGCTATACATACGGATGAGTTCGCAGTTGTCATCGGTAACAGACCTTCAATTTCTAACAGGGCAACAGCAACCTGCGTGCTTGCGCTGATGTCCATGCTCGCTTGTAGGATCACTGACGCCAAAGGAGGCGTTCCGAGGCGTCGCTGACCGCTAGGTACAATTGATGCCACAGACTGTTTGGAAGTGAAGGCAGGTGGTGCGATCCGTTTCGCGGGATGTCTCAGAATCGGCTTTAATCCTAGGTAGTAGAGAATCCAATGCTACGGACTTTATCTTAGGATTAACCGCAATCTATCCAACATGATTGCAGGAGCTCGGCACTGTGACTGTTCTGTAAGGCTACTGCACGCCAAGCTAGGCGGATACTGGCTGCCCAAACAAGGTACTAGATCCCGGATTGAGGAGCACTAGCAAGCTTAAACCCGAGACAGGTAAAGCCTGCGCTGTGGCCCCATTGGCAGCGTCCTGCTCCTGTAGGCGGAGCAAAAGCACCCAGGCAGAGCACAACTCGTGGATGGCTTCACACCTTAGATTTTCGGAAGTCCCCCGCACATTCATGCGAACGCAGCCCCGGCGATGAATAGGTTGTTAATATTTGTGCAAGACACGGGGGTTCCGTCCAAAATGTTCTGTGTAAGTTCATCTTCAAGCAGAATACAACTCCAGTGGTTTTGTGGCTGCGACCAAGCCGGTAGCGCCACGATTTTTTATTTTAGTCCGGGCCAGGATAA', 
        'CCAGATAAGATGTGTGACCGGTCTAACCTGTATAGATTCACTCCTGCCTCGTTCCTTCCCCGTAGTTGGGGAGAAGACTCTGTTTGAGTATAGCGCGTTTTAGAACCGGCGGCCGGTGATATGGTTTAATTGACACTAAAACTCACCAGGAGGGTACTCTGAGCCGAGTCACCCCTTTTCCCTTTGCTCACGAGAGCGCACCCCAGTCCCTTTTAAGTGGACGTCTTCAGTGCGCACCATAAAGATCTCGGGCCGGCGTGTAGTTCGTGACTACCGTGGCAGAAATGCTTCAACTGATCAGTACACACGGGGGACCGCTGTGTATATGTTTGAGTGCAACGGCATCCGAAAACCCAGCCCATGATGTCAAAGTAACACAGCGGCAACAATGCATTCAAGCAGAGGTTGTAACCGGCCGCCATTCGACCCGGCTAAACACCTAGATTAAGCATGTAATCAGAACGAATTCCACTAGATAAGGGGTGTCCCATTCTAACCCGTGATTTACGGCCGGACCTACTGAATGCTTTTTGGAGTGGCATAACACGAGTTTGTTCCCGTAGTACAGTCACGTATAGGTCGCCAGAAGAAAAAGAGCGGCGCTTAAGCATCTAGCGTTCGGTCTATAGATTCCAAAATCAAGATTTAATGCGTCAGTTGTGATACGTCTCTCGGTTTGTCTGGGAAAACTTGACAGCCGCCTAAAGCCTGGATTTTATGACGCAGGCTGGTCGAGGATCCATTAACAGTGGATACGATAATCGGACGGAAAGGTGCTACTGCCCTCCCCTCTAAGCATAGGATCACGTGGATTTCTAATTCGCGAGAGATCACCGATCGGGAGACCGCTTTTTGATCTCTCATATCTTAAAATGACGCTGCTGCAACGCTACGGAGTGGCTTGTTCTGTTAATTTTGTACTTTATGTGAGACTG'
        )

