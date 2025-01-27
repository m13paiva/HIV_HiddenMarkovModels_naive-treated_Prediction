import os
from hmm_profile import build_HMM_profile
from prepare_hmm_data import prepare_hmm_data,prepare_hmm_data_test
from background_probs import BACKGORUND_PROBS
from classify import *
from auxiliary import *
from navigate import *
import random

def build_profiles(folder_path):
    files=os.listdir(folder_path)
    profiles = {}
    total = len(files)
    print_progress_bar(0, total, "Building HMM profiles:")
    for i, file in enumerate(files):
        msa = read_msa(os.path.join(folder_path, file))
        hmm_profile = build_HMM_profile(msa, BACKGORUND_PROBS)
        profiles[msa.name] = hmm_profile
        print_progress_bar(i + 1, total, "Building HMM profiles:")
    return profiles

def score_sequences(seqs,profiles):
    seq_scores = {}
    total = len(seqs)
    print_progress_bar(0, total, "Scoring sequences:    ")
    for i, seq in enumerate(seqs):
        seq_label = seq[0]
        sequence = seq[1]
        seq_scores[seq_label] = {}
        log_bg = log_background_prob(sequence, BACKGORUND_PROBS)
        for profile in profiles:
            score = profiles[profile].log_odds(sequence, log_bg)
            seq_scores[seq_label][profile] = score
        print_progress_bar(i + 1, total, "Scoring sequences:    ")
    return seq_scores

def score_sequences_test(seqs,profiles):
    seq_scores = {}
    total = len(seqs)
    print_progress_bar(0, total, "Scoring sequences:    ")
    for i, seq in enumerate(seqs):
        seq_label = (seq[0],seq[1])
        sequence = seq[2]
        seq_scores[seq_label] = {}
        log_bg = log_background_prob(sequence, BACKGORUND_PROBS)
        for profile in profiles:
            score = profiles[profile].log_odds(sequence, log_bg)
            seq_scores[seq_label][profile] = score
        print_progress_bar(i + 1, total, "Scoring sequences:    ")
    return seq_scores

def execute_model(fasta_file, subtype, pred_subtype, dis_drm, output_file,PI_codons, RT_codons ):
    anim = DotAnimation("Preparing profile data", interval=0.7)
    anim.start()
    if dis_drm:
        codons = compile_codons(PI_codons, RT_codons)
    else:
        codons = None

    split_subtype = True
    if (subtype == "unknown") and (not pred_subtype):
        split_subtype = False

    seqs = read_fasta(fasta_file, dis_drm, codons)
    rawdata_path = get_abs_path("model_data", "raw_data", "global_AA_table_wts.csv")
    proccessdata_path = get_abs_path("model_data", "data")
    prepare_hmm_data(
        rawdata_path,
        proccessdata_path,
        codons,
        subtype,
        dis_drm,
        split_subtype
    )
    anim.stop()
    profiles = build_profiles(proccessdata_path)
    seq_scores = score_sequences(seqs,profiles)
    seq_classified = classify(seq_scores)
    save_output(seq_classified, output_file)
    print_output(seq_classified)

def test_model(subtype, n_iter, test_frac, pred_subtype, dis_drm,PI_codons, RT_codons,
               test_samples=None, train_samples=None, eq=True):

    rawdata_path = get_abs_path("model_data", "raw_data", "global_AA_table_wts.csv")
    train_path=get_abs_path("model_data", "train")
    test_path = get_abs_path("model_data", "test")
    output_file= get_abs_path("output_data", "test_output.csv")

    if dis_drm:
        codons = compile_codons(PI_codons, RT_codons)
    else:
        codons = None

    split_subtype = True
    if (subtype == "unknown") and (not pred_subtype):
        split_subtype = False

    random_states = [random.randint(1, 1000) for _ in range(n_iter)]
    seq_count=1
    total_seq_scores={}
    for i,random_state in enumerate(random_states):
        seq_scores={}
        print(f"\n»»»»»»»»»»»»»»»»»»»»»»» ITERATION {i+1}/{n_iter} »»»»»»»»»»»»»»»»»»»»»»\n")
        #print(f"test {test_samples} train {train_samples}")
        anim = DotAnimation("Preparing profile data", interval=0.7)
        anim.start()
        prepare_hmm_data_test(rawdata_path,
                              train_path,
                              test_path,
                              codons,
                              subtype,
                              test_frac,
                              test_samples,
                              train_samples,
                              dis_drm,
                              split_subtype,
                              eq,
                              random_state
        )
        anim.stop()
        profiles=build_profiles(train_path)
        seqs,seq_count=read_test(test_path,seq_count)
        seq_scores=score_sequences_test(seqs,profiles)
        write_output(profiles, seq_scores, output_file)
        read_output(output_file)
        total_seq_scores.update(seq_scores)
    #print(seq_scores)
    write_output(profiles,total_seq_scores,output_file)
    read_output(output_file)

def test():
    test_model(subtype="B",
               n_iter=50,
               test_frac=0.2,
               pred_subtype=False,
               dis_drm=True,
               PI_codons=[30, 32, 33, 46, 47, 48, 50, 54, 76, 82, 84, 88, 90],
               RT_codons=[184, 65, 70, 74, 115, 41, 67, 70, 210, 215, 219, 100, 101, 103, 106, 138, 181, 188, 190, 230],
               test_samples=None,
               train_samples=None,
               eq=True)

#test()







