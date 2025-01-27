from hmm import HiddenMarkovModel

# Counts non-gap residues ('-') in a given sequence
def count_residues(seq):
    count = 0
    for char in seq:
        if char != "-":
            count += 1
    return count

# Returns a list of residue counts for each column in the multiple sequence alignment (MSA)
def columns_length(msa):
    count_list = []
    for i in range(len(msa)):
        count_list.append(count_residues(msa.column(i)))  # Counts residues in each column
    return count_list

# Assigns labels (Match 'M' or Insert 'I') to each column based on the number of residues
def column_labels(msa):
    thres = int(len(msa.listseqs) / 2)  # Threshold: half the number of sequences
    col_len = columns_length(msa)
    state_list = []
    count = 0  # Counter for match states
    for i in col_len:
        if i <= thres:
            state_list.append(("I", count))  # Columns with fewer residues are Inserts
        else:
            count += 1
            state_list.append(("M", count))  # Columns with more residues are Matches
    return state_list

# Returns the length of the HMM (number of Match states)
def length_HMM(col_labels):
    return col_labels[-1][1]  # Last match state index

# Builds the list of all possible states in the HMM (Match, Insert, Delete)
def states_HMM(msa):
    labels = column_labels(msa)
    length = length_HMM(labels)
    states = []

    # Add Match states
    M_state = 0
    while M_state <= length + 1:
        states.append(("M", M_state))
        M_state += 1

    # Add Insert states
    I_state = 0
    while I_state <= length:
        states.append(("I", I_state))
        I_state += 1

    # Add Delete states
    D_state = 1
    while D_state <= length:
        states.append(("D", D_state))
        D_state += 1

    return states

# Modifies emission probabilities for Insert states by setting gaps ('-') to zero
def emission_insert(prob):
    res = prob
    res["-"] = 0
    return res

# Modifies emission probabilities for Delete states by setting all probabilities to zero except gaps ('-')
def emission_delete(prob):
    res = {i: 0 for i in prob}  # Zero out all probabilities
    res["-"] = 1  # Gaps are fully probable in a delete state
    return res

# Initializes emission probabilities for the Start state
def emission_start(alphabet):
    res = {}
    for i in alphabet:
        if i != "-":
            res[i] = 0  # No emission probability for residues
        else:
            res[i] = 1  # Only emits gaps ('-')
    return res

# Initializes emission probabilities for the End state (all probabilities are zero)
def emission_end(alphabet):
    res = {i: 0 for i in alphabet}
    return res

# Counts the occurrences of each residue in a column of the MSA
def emission_column(alphabet, column): 
    res = {i: 0 for i in alphabet}
    for char in column:
        if char != "-":
            res[char] += 1
    return res

# Calculates emission probabilities for a Match state using pseudocounts
def emission_match(alphabet, column_ind, msa):
    count = emission_column(alphabet, msa.column(column_ind))  # Residue counts in column
    residues = count_residues(msa.column(column_ind))  # Number of non-gap residues
    count = count.items()
    dic = {}

    for i in count:
        emission_prob = (i[1] + 1) / (residues + len(alphabet))  # Add pseudocounts
        dic[i[0]] = emission_prob
    dic["-"] = 0  # No emissions of gaps in a Match state
    return dic

# Identifies the final Match state in the HMM
def end_state(col_labels):
    return ("M", length_HMM(col_labels) + 1)

# Calculates all emission probabilities for the HMM
def emission_probabilities(msa, background_probs):
    probs = {}
    alphabet = list(background_probs.keys())  # Amino acid alphabet
    alphabet.append("-")  # Add gap symbol
    col_labels = column_labels(msa)
    states = states_HMM(msa)
    end = end_state(col_labels)

    for state in states:
        if state[0] != "D":  # No emissions for Delete states
            if state[0] == "M":  # Match state emissions
                if state == ("M", 0):
                    probs[state] = emission_start(alphabet)  # Start state
                elif state == end:
                    probs[state] = emission_end(alphabet)  # End state
                else:
                    col_index = col_labels.index(state)
                    probs[state] = emission_match(alphabet, col_index, msa)
            else:
                probs[state] = emission_insert(background_probs)  # Insert state
        else:
            probs[state] = emission_delete(background_probs)  # Delete state
    return probs

# Generates the state transitions for a single sequence based on the MSA and column labels
def sequence_transitions(sequence, labels, length):
    transitions = [("M", 0)]  # Start state
    for i, label in enumerate(labels):
        if label[0] == "M":
            if sequence[i] == "-":
                transitions.append(("D", label[1]))  # Transition to Delete state
            else:
                transitions.append(label)  # Match state
        else:
            if sequence[i] != "-":
                transitions.append(label)  # Insert state
    transitions.append(("M", length + 1))  # End state
    return transitions

# Initializes transition counts for all HMM states
def state_transitions(length):
    dic = {}
    # Match states
    for i in range(length + 1):
        dic[("M", i)] = {("M", i+1): 0, ("I", i): 0, ("D", i+1): 0}
    dic[("M", length + 1)] = {}

    # Insert states
    for i in range(length + 1):
        dic[("I", i)] = {("M", i+1): 0, ("I", i): 0, ("D", i+1): 0}

    # Delete states
    for i in range(1, length + 1):
        dic[("D", i)] = {("M", i+1): 0, ("I", i): 0, ("D", i+1): 0}
    return dic

# Aligns sequences in the MSA and counts state transitions
def align_transitions(msa, labels):
    length = length_HMM(labels)
    transitions = state_transitions(length)
    for seq in msa:
        seq_trans = sequence_transitions(seq, labels, length)
        for i in range(len(seq_trans) - 1):
            transitions[seq_trans[i]][seq_trans[i+1]] += 1
    return transitions

# Converts transition counts to probabilities using pseudocounts
def align_probabilities(trans):
    probs = {}
    for start_state in trans:
        probs[start_state] = {}
        total_transitions = 3  # Pseudocount
        for end_state in trans[start_state]:
            total_transitions += trans[start_state][end_state]
        for end_state in trans[start_state]:
            probs[start_state][end_state] = (trans[start_state][end_state] + 1) / total_transitions
    return probs

# Calculates transition probabilities between all HMM states using a multiple sequence alignment (MSA).
def transition_probabilities(msa):
    labels = column_labels(msa)  # Label columns as Match (M) or Insert (I)
    states = states_HMM(msa)  # Generate all possible states (M, I, D)
    trans = align_transitions(msa, labels)  # Count transitions between states
    probs = align_probabilities(trans)  # Convert counts into probabilities with pseudocounts

    dic = {}  # Final transition probabilities dictionary
    for start_state in probs:  
        dic[start_state] = {}
        for end_state in states:  
            dic[start_state][end_state] = probs[start_state].get(end_state, 0)
    return dic

# Calculates initial state probabilities
def initial_probabilities(msa):
    states = states_HMM(msa)
    return {i: 1 if i == ("M", 0) else 0 for i in states}

# Builds the final HMM profile from the MSA
def build_HMM_profile(msa,background_probs):
    initial_probs = initial_probabilities(msa)
    emission_probs = emission_probabilities(msa,background_probs)
    transition_probs = transition_probabilities(msa)
    HMMprofile = HiddenMarkovModel(initial_probs, emission_probs, transition_probs, msa.name)
    return HMMprofile