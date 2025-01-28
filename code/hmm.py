import math
import numpy as np
import warnings
import numpy as np

with warnings.catch_warnings():
    # Suppress all RuntimeWarnings.
    warnings.filterwarnings("ignore", category=RuntimeWarning)

def save_dict_to_txt_file(data, file_path):
    """
    Saves a dictionary to a .txt file.

    Args:
        data (dict): The dictionary to save.
        file_path (str): The path of the .txt file to save to.
    """
    try:
        with open(file_path, 'w') as file:
            for key, value in data.items():
                file.write(f"{key}: {value}\n")
    except Exception as e:
        print(f"An error occurred: {e}")

class HiddenMarkovModel:

    def __init__(self, init_probs, emission_probs, trans_probs, name):
        """Create a constructor based on three different attributes: probability of start states; emission probabilities matrix; transition probabilities matrix.
        Both emission and transition probability matrices can be implemented as dictionaries of dictionaries. States and symbols are represented as lists and can be infered from the probabilities.
        """
        self.initstate_probs = init_probs
        self.emission_probs = emission_probs
        self.transition_probs = trans_probs
        self.name=name
        self.states = emission_probs.keys()
        self.symbols = list(emission_probs[list(self.states)[0]].keys())

    def get_init_prob(self, state):
        '''Get initial probability of a given state'''
        if state in self.states:
            return (self.initstate_probs[state])
        else:
            return 0

    def get_emission_prob(self, state, symbol):
        '''Get probability of a given state to emit a symbol'''
        if state in self.states and symbol in self.symbols:
            return (self.emission_probs[state][symbol])
        else:
            return 0

    def get_transition_prob(self, state_orig, state_dest):
        '''Get probability of transition from a origin state to destination state'''
        if state_orig in self.states and state_dest in self.states:
            return (self.transition_probs[state_orig][state_dest])
        else:
            return 0

    def set_init_prob(self, state, p):
        '''Set initial probability of a given state'''
        if state in self.states:
            self.initstate_probs[state] = p

    def set_emission_prob(self, state, symbol, p):
        '''Set probability of a given state emit a symbol'''
        if state in self.states and symbol in self.symbols:
            self.emission_probs[state][symbol] = p

    def set_transition_prob(self, state_orig, state_dest, p):
        '''Set probability of transition from a origin state to destination state'''
        if state_orig in self.states and state_dest in self.states:
            self.transition_probs[state_orig][state_dest] = p
    '''
    def forward(self, sequence):
        """Calculate the forward probabilities for a given sequence using optimized matrix operations."""
        seq_len = len(sequence)
        if seq_len == 0:
            print("SEQUÊNCIA VAZIA!")
            return 0.0

        # Convert transition probabilities into a matrix
        transition_matrix = np.array([
            [self.get_transition_prob(state_orig, state_dest) for state_dest in self.states]
            for state_orig in self.states
        ])

        # Initialize probabilities for the first sequence symbol
        initial_probs = np.array([self.get_init_prob(state) for state in self.states])
        first_emissions = np.array([self.get_emission_prob(state, sequence[0]) for state in self.states])
        current_probs = initial_probs * first_emissions

        # Iterate through the sequence using matrix operations
        for i in range(1, seq_len):
            # Get emission probabilities for the current symbol
            current_symbol = sequence[i]
            emission_probs = np.array([self.get_emission_prob(state, current_symbol) for state in self.states])

            # Update probabilities: matrix multiplication + element-wise multiplication for emissions
            current_probs = np.dot(current_probs, transition_matrix) * emission_probs
        # print(sequence)
        #for prob in current_probs:
            # print(prob)

        # Return the total probability
        return np.sum(current_probs)

    def forward_log(self, sequence):
        """Calculate the log forward probabilities for a given sequence using log-space operations."""
        seq_len = len(sequence)
        if seq_len == 0:
            print("SEQUÊNCIA VAZIA!")
            return float('-inf')  # Log of 0 is -inf

        # Create mappings for state indices
        # state_index = {state: idx for idx, state in enumerate(self.states)}
        num_states = len(self.states)

        # Convert transition probabilities into a log matrix
        log_transition_matrix = np.log([
            [self.get_transition_prob(state_orig, state_dest) for state_dest in self.states]
            for state_orig in self.states
        ])

        # Initialize log probabilities for the first sequence symbol
        log_initial_probs = np.log([self.get_init_prob(state) for state in self.states])
        log_first_emissions = np.log([self.get_emission_prob(state, sequence[0]) for state in self.states])
        log_current_probs = log_initial_probs + log_first_emissions

        # Iterate through the sequence in log-space
        for i in range(1, seq_len):
            current_symbol = sequence[i]
            log_emission_probs = np.log([self.get_emission_prob(state, current_symbol) for state in self.states])

            # Update probabilities in log-space
            log_next_probs = np.full(num_states, float('-inf'))  # Initialize with -inf for log-space summation
            for state_dest in range(num_states):
                log_sum = np.logaddexp.reduce(
                    log_current_probs + log_transition_matrix[:, state_dest]
                )  # Log-sum-exp trick
                log_next_probs[state_dest] = log_sum + log_emission_probs[state_dest]

            log_current_probs = log_next_probs

        # Log-sum-exp for the final probabilities
        return np.logaddexp.reduce(log_current_probs)
    

    
    def forward_log_vectorized(self, sequence):
        seq_len = len(sequence)
        if seq_len == 0:
            return float('-inf')

        num_states = len(self.states)

        # Precompute log transition matrix, initial, and emission probabilities
        log_transition_matrix = np.log([
            [self.get_transition_prob(s_from, s_to) for s_to in self.states]
            for s_from in self.states
        ])
        log_initial_probs = np.log([self.get_init_prob(s) for s in self.states])
        log_first_emissions = np.log([self.get_emission_prob(s, sequence[0]) for s in self.states])

        # Initialize the DP array of shape (num_states,)
        log_current_probs = log_initial_probs + log_first_emissions

        # Forward pass
        for i in range(1, seq_len):
            symbol = sequence[i]
            # Emission log-probs for all states
            log_emission_probs = np.log([self.get_emission_prob(s, symbol) for s in self.states])

            # Vectorized log-sum-exp:
            # 1) shape (num_states, 1) + shape (num_states, num_states) => (num_states, num_states)
            tmp = log_current_probs[:, np.newaxis] + log_transition_matrix
            # 2) reduce across the "origin" dimension => shape (num_states,)
            log_next_probs = np.logaddexp.reduce(tmp, axis=0)
            # 3) add emission probabilities => shape (num_states,)
            log_next_probs += log_emission_probs

            log_current_probs = log_next_probs

        # Final log-sum-exp
        return np.logaddexp.reduce(log_current_probs)

    '''
    def forward(self, sequence):
        """
        Given an observed sequence, calculate the forward probability of the sequence
        using scaling at each time step (to prevent underflow), and NumPy for speed.
        """
        import math
        import numpy as np

        T = len(sequence)
        if T == 0:
            return 0.0

        states = self.states
        n_states = len(states)

        # 1. Build NumPy arrays from your existing accessors
        #    (No changes to the class; we just query it here.)
        init_prob = np.array([self.get_init_prob(s) for s in states], dtype=float)
        trans_prob = np.zeros((n_states, n_states), dtype=float)
        for i, s_orig in enumerate(states):
            for j, s_dest in enumerate(states):
                trans_prob[i, j] = self.get_transition_prob(s_orig, s_dest)

        # alpha[t,i] = scaled forward prob of state i at time t
        alpha = np.zeros((T, n_states), dtype=float)
        # scale_factors[t] = sum of unscaled alpha[t]
        scale_factors = np.zeros(T, dtype=float)

        # 2. Initialization at t = 0
        obs0 = sequence[0]
        for i, s in enumerate(states):
            alpha[0, i] = init_prob[i] * self.get_emission_prob(s, obs0)
        c0 = alpha[0].sum()
        if c0 == 0.0:
            return 0.0
        alpha[0] /= c0
        scale_factors[0] = c0

        # 3. Recursion: for t = 1..T-1
        for t in range(1, T):
            obs_t = sequence[t]
            # Vector/matrix multiply: alpha[t-1] · trans_prob => shape (n_states,)
            unscaled = alpha[t - 1].dot(trans_prob)
            # Multiply each state-prob by that state's emission probability
            for j, s in enumerate(states):
                unscaled[j] *= self.get_emission_prob(s, obs_t)

            ct = unscaled.sum()
            if ct == 0.0:
                return 0.0
            alpha[t] = unscaled / ct  # scale
            scale_factors[t] = ct

        # 4. Final probability = product of scale_factors
        #    We accumulate in log-space to avoid underflow:
        log_prob = np.sum(np.log(scale_factors))
        return log_prob

    def log_odds(self, sequence, log_background_prob):
        log_hmm_prob = self.forward(sequence)
        # print(f"Log HMM Probability: {log_hmm_prob}")
        log_odds_score = log_hmm_prob - log_background_prob

        return log_odds_score


    def save_model(self):
        d={"initial_probs":self.initstate_probs,
           "emission_probs":self.emission_probs,
           "transition_probs":self.transition_probs}
        save_dict_to_txt_file(d, f"{self.name}.txt")


