import argparse
import getpass
from navigate import *
from run_model import *

# =========== Constants and Global Variables ===========
SUBTYPES = ["A", "B", "02_AG", "UNKNOWN"]
PI_codons = [30, 32, 33, 46, 47, 48, 50, 54, 76, 82, 84, 88, 90]
RT_codons = [184, 65, 70, 74, 115, 41, 67, 70, 210, 215, 219, 100, 101, 103, 106, 138, 181, 188, 190, 230]
SEP_BAR = "\n==============================================================\n"
PROMPT_ARROW = "»»» "
CODON_FORMAT="List of codon positions must be given in the format 'pos1, pos2, ..., pos N'.\n"
WRONG_FORMAT="\nWrong format! Try again...\n"
INVALID_ANS="\nInvalid answer! Try again...\n"

# =========== Custom Multi-line Help ===========
HELP = r"""
==============================================================
     HIV-1 ART-naïve/treated prediction with HMM - v1.0
==============================================================

 USAGE:
    python hiv_hmm_pred.py [OPTIONS]

 DESCRIPTION:
    This program takes HIV-1 "pol" gene sequences (in FASTA format) 
    and predicts if they belong to a patient who has experienced 
    antiretroviral treatment (ART-treated) or not (ART-naïve),
    using Hidden Markov Models.

 OPTIONS:
    -h, --help
        Display this help message and exit.

    -c, <IN_FILE>, --classify <IN_FILE>
        Classify sequences in the provided FASTA file.

        (Sub-option)
        -o, <OUT_FILE>, --output <OUT_FILE>
            Specify an output CSV file path for classification.
            (Default: output_data/output.csv)

    -t, --test
        Test performance of HMM algorythm.

 ARGUMENTS:
    IN_FILE    Fasta file with sequences to classify.

==============================================================
 Developed by : Manuel Almeida
 Email address: mpa.almeida@campus.fct.unl.pt
 Available in : 
https://github.com/m13paiva/IHMT-HIV-naive-experienced-prediction
==============================================================
"""


# =========== Existing Functions (Unmodified) ===========
def test_cmd():
    """
    Interactive flow to test sequences as ART-naïve or treated.
    """
    global PI_codons
    global RT_codons
    test_samples=None
    train_samples=None
    print(SEP_BAR, "\nThis section trains and then tests the model based is your")
    print("input and returns performance info (for each iteration and for\noverall performance).", SEP_BAR)
    proceed = False
    while not proceed:
        pred_subtype = False
        print("Please write which HIV subtype you intend to model.")
        print('("A", "B", "02_AG" or "unknown")\n')
        subtype = input(PROMPT_ARROW).upper()
        if subtype not in SUBTYPES:
            print("\nInvalid subtype!\n")
        elif subtype == "UNKNOWN":
            subtype=subtype.lower()
            print("\nI see that you want to model unknown subtypes.\n")
            proceed2 = False
            while not proceed2:
                print("Would you like the model to predict subtype as well?")
                print("(y/n)\n")
                ans = input(PROMPT_ARROW).upper()
                if ans == "Y":
                    pred_subtype = True
                    proceed2 = True
                    proceed = True
                    print(SEP_BAR)
                elif ans == "N":
                    pred_subtype = False
                    proceed2 = True
                    proceed = True
                    print(SEP_BAR)
                else:
                    print(INVALID_ANS)
        else:
            proceed = True
            print(SEP_BAR)

    proceed = False
    while not proceed:
        print("\nDo you want the model to discard the major DRM codons?")
        print("(y/n)\n")
        ans = input(PROMPT_ARROW).upper()
        if ans == "Y":
            print("\nThe major DRM codons considered by the model are the following:")
            print(f"» PI codons: {PI_codons}")
            print(f"» RT codons: {RT_codons}\n")
            proceed = True
            dis_drm = True
            proceed2 = False
            while not proceed2:
                # spcecify which DRM codons user wants to discard
                print("\nDo you want to specify a different list of DRM codons?\n(y/n)\n")
                ans = input(PROMPT_ARROW).upper()
                if ans == "Y":
                    proceed3 = False
                    while not proceed3:
                        print("\nSpecify Protease Inhibitor (PI) codons.\n",
                              CODON_FORMAT)
                        codons = input("PI codons " + PROMPT_ARROW)
                        codons_ls = codons_prompt(codons)
                        if not verify_codon_format(codons_ls):
                            print(WRONG_FORMAT)
                        else:
                            PI_codons = codons_ls
                            proceed3 = True
                    proceed3 = False
                    while not proceed3:
                        print("\nSpecify Reverse Transcriptase (RT) codons.\n",
                              CODON_FORMAT)
                        codons = input("RT codons " + PROMPT_ARROW)
                        codons_ls = codons_prompt(codons)
                        if not verify_codon_format(codons_ls):
                            print(WRONG_FORMAT)
                        else:
                            RT_codons = codons_ls
                            proceed3 = True
                    proceed2 = True
                elif ans == "N":
                    proceed2 = True
                else:
                    print(INVALID_ANS)

            print(SEP_BAR)
        elif ans == "N":
            proceed = True
            dis_drm = False
            print(SEP_BAR)
        else:
            print(INVALID_ANS)
    # specify test fraction
    proceed=False
    while not proceed:
        print("\nWhich fraction of the data do you want to use for testing?\n",
              "The answer must be a number between 0 and 1 and the decimal\n separator must be a dot '.'\n")
        test_frac=input(PROMPT_ARROW)
        try:
            test_frac=float(test_frac)
            if not (0<test_frac<1):
                print("\nThe answer must be a number between 0 and 1! Try again...\n")
            else:
                proceed=True
                print(SEP_BAR)
        except:
            print(WRONG_FORMAT)

    # ask if user wants to downsize data set
    proceed = False
    while not proceed:
        print("\nDo you want to use a downsized versison of the dataset?\n(y/n)\n")
        ans = input(PROMPT_ARROW).upper()
        if ans=="Y":
            proceed2=False
            while not proceed2:
                print("\nIndicate de maximum number of samples for training:")
                print("(if you don't want to limit training input 'none'\n")
                train_samples=input(PROMPT_ARROW).upper()
                if train_samples=="NONE":
                    train_samples=None
                    proceed2=True
                else:
                    try:
                        train_samples=int(train_samples)
                        proceed2=True
                    except:
                        print(WRONG_FORMAT)
            proceed2 = False
            while not proceed2:
                print("\nIndicate de maximum number of samples for testing:")
                print("(if you don't want to limit testing input 'none'\n")
                test_samples = input(PROMPT_ARROW).upper()
                if test_samples=="NONE":
                    test_samples=None
                    proceed2 = True
                else:
                    try:
                        test_samples = int(test_samples)
                        proceed2 = True
                    except:
                        print(WRONG_FORMAT)
            proceed=True
            print(SEP_BAR)
        elif ans=="N":
            proceed=True
            print(SEP_BAR)
        else:
            print(INVALID_ANS)
    proceed=False
    while not proceed:
        print("\nDo you want to use a balanced dataset?")
        print("(same number of naive and treated sequences but, consequently less data)\n")
        ans = input(PROMPT_ARROW).upper()
        if ans == "Y":
            eq=True
            proceed=True
            print(SEP_BAR)
        elif ans=="N":
            eq=False
            proceed=True
            print(SEP_BAR)
        else:
            print(INVALID_ANS)
    proceed=False
    while not proceed:
        print("\nIndicate the number of iterations:\n")
        n_iter=input(PROMPT_ARROW)
        try:
            n_iter=int(n_iter)
            proceed=True
            print(SEP_BAR)
        except:
            print(WRONG_FORMAT)
    #print(f"test {test_samples} train {train_samples}")
    test_model(subtype, n_iter, test_frac, pred_subtype, dis_drm, PI_codons, RT_codons, test_samples, train_samples, eq)








    """test_model(subtype=subtype,
               n_iter=50,
               test_frac=0.2,
               pred_subtype=pred_subtype,
               dis_drm=True,
               PI_codons=[30, 32, 33, 46, 47, 48, 50, 54, 76, 82, 84, 88, 90],
               RT_codons=[184, 65, 70, 74, 115, 41, 67, 70, 210, 215, 219, 100, 101, 103, 106, 138, 181, 188, 190,
                          230],
               test_samples=None,
               train_samples=None,
               eq=True)"""

def classify_cmd(fasta_file, output_file):
    """
    Interactive flow to classify sequences as ART-naïve or treated.
    """
    global PI_codons
    global RT_codons
    username = getpass.getuser()
    print(SEP_BAR, f"Hello Dr. {username}! Ready to find who's naïve and who's treated?", SEP_BAR)
    use_prev_choices = False

    pred_subtype = False
    choices_file = get_abs_path("model_data", "choices.txt")
    use_prev_choices = False
    # verify if the previous choices were saved and if so ask if user wants to use them
    if os.path.exists(choices_file):
        proceed = False
        while not proceed:
            print("\nDo you want to use the choices previously saved?\n(y/n)\n")
            ans = input(PROMPT_ARROW).upper()
            if ans == "Y":
                d = access_choices()
                pred_subtype = d["pred_subtype"]
                dis_drm = d["dis_drm"]
                PI_codons = d["PI_codons"]
                RT_codons = d["RT_codons"]
                use_prev_choices = True
                proceed = True
                print(SEP_BAR)
            elif ans == "N":
                proceed = True
                print(SEP_BAR)
            else:
                print(INVALID_ANS)
    proceed = False
    # ask subtype of samples
    while not proceed:
        print("Please write the HIV subtype of the sample(s).")
        print('("A", "B", "02_AG" or "unknown")\n')
        subtype = input(PROMPT_ARROW).upper()
        if subtype not in SUBTYPES:
            print("\nInvalid subtype!\n")
        elif subtype == "UNKNOWN":
            subtype = subtype.lower()
            if not use_prev_choices:
                print("\nI see that you don't know the subtype of your sample(s).")
                proceed2 = False
                while not proceed2:
                    print("Would you like to predict the subtype as well?")
                    print("(y/n)\n")
                    ans = input(PROMPT_ARROW).upper()
                    if ans == "Y":
                        pred_subtype = True
                        proceed2 = True
                        print(SEP_BAR)
                    elif ans == "N":
                        pred_subtype = False
                        proceed2 = True
                        print(SEP_BAR)
                    else:
                        print(INVALID_ANS)
            proceed = True
        else:
            proceed = True


    if not use_prev_choices:
        # ask if user wants to discard major DRM codons
        print(SEP_BAR)
        proceed = False
        while not proceed:
            print("\nDo you want the model to discard the major DRM codons?")
            print("(y/n)\n")
            ans = input(PROMPT_ARROW).upper()
            if ans == "Y":
                print("\nThe major DRM codons considered by the model are the following:")
                print(f"» PI codons: {PI_codons}")
                print(f"» RT codons: {RT_codons}\n")
                proceed = True
                dis_drm = True
                proceed2 = False
                while not proceed2:
                    # spcecify which DRM codons user wants to discard
                    print("\nDo you want to specify a different list of DRM codons?\n(y/n)\n")
                    ans = input(PROMPT_ARROW).upper()
                    if ans=="Y":
                        proceed3=False
                        while not proceed3:
                            print("\nSpecify Protease Inhibitor (PI) codons.\n",
                                  CODON_FORMAT)
                            codons=input("PI codons "+PROMPT_ARROW)
                            codons_ls=codons_prompt(codons)
                            if not verify_codon_format(codons_ls):
                                print(WRONG_FORMAT)
                            else:
                                PI_codons=codons_ls
                                proceed3=True
                        proceed3=False
                        while not proceed3:
                            print("\nSpecify Reverse Transcriptase (RT) codons.\n",
                                  CODON_FORMAT)
                            codons=input("RT codons "+PROMPT_ARROW)
                            codons_ls=codons_prompt(codons)
                            if not verify_codon_format(codons_ls):
                                print(WRONG_FORMAT)
                            else:
                                RT_codons=codons_ls
                                proceed3=True
                        proceed2 = True
                    elif ans=="N":
                        proceed2 = True
                    else:
                        print(INVALID_ANS)


                print(SEP_BAR)
            elif ans == "N":
                proceed = True
                dis_drm = False
                print(SEP_BAR)
            else:
                print(INVALID_ANS)
        # ask if user wants to save their choices
        proceed = False
        while not proceed:
            print("Would you like the program to remember your choices?\n(y/n)\n")
            ans = input(PROMPT_ARROW).upper()
            if ans == "Y":
                choices = {"pred_subtype": pred_subtype, "dis_drm": dis_drm, "PI_codons": PI_codons,
                           "RT_codons": RT_codons}
                save_choices(choices)
                print("\nYour choices have been registered!",
                      "\nThe next time the program is ran, you'll be asked if you want \nto execute your previous choices.")
                proceed = True
            elif ans == "N":
                delete_choices()
                proceed = True
            else:
                print(INVALID_ANS)
    print(SEP_BAR)
    execute_model(fasta_file, subtype, pred_subtype, dis_drm, output_file, PI_codons, RT_codons)

def main():
    # Create a parser, but don't rely on argparse's default -h/--help message
    parser = argparse.ArgumentParser(add_help=False)

    # Custom -h/--help (store_true) to print your multi-line HELP
    parser.add_argument(
        "-h", "--help",
        action="store_true",
        help="Display this help message and exit."
    )

    # -c/--classify expects the FASTA file path
    parser.add_argument(
        "-c", "--classify",
        metavar="IN_FILE",
        type=str,
        help="Classifies the sequences in the given F_FILE (FASTA)."
    )

    # -o/--output is relevant only when classifying
    parser.add_argument(
        "-o", "--output",
        metavar="OUT_FILE",
        type=str,
        default="output.csv",
        help="(Used with -c) Specify the output CSV file path."
    )

    # -t/--test runs the test_cmd function
    parser.add_argument(
        "-t", "--test",
        action="store_true",
        help="Run the test suite."
    )

    args, unknown = parser.parse_known_args()

    # 1) Handle unknown commands
    if unknown:
        print("\nNo valid options provided. Run with -h or --help for usage.\n")
        sys.exit(1)

    # 2) If user requests help, print your custom HELP
    if args.help:
        print(HELP)
        sys.exit(0)

    # 3) If user wants to run tests
    if args.test:
        test_cmd()
        sys.exit(0)

    # 4) If user wants to classify, check file validity
    if args.classify:
        file_path = get_abs_path("input_data",args.classify)

        # Check file extension
        if not str(file_path).endswith(".fasta"):
            print("ERROR: Wrong file format!")
            print("Please indicate a FASTA file (must end with .fasta).")
            sys.exit(1)

        if not os.path.exists(file_path):
            print(f"The file '{file_path}' does not exist!")
            print("Please indicate a valid FASTA file.")
            sys.exit(1)

        # Use the custom or default output file
        classify_cmd(fasta_file=file_path, output_file=get_abs_path("output_data",args.output))

    else:
        # If no args are provided (or no valid flags given), guide the user:
        print("\nNo valid options provided. Run with -h or --help for usage.\n")


if __name__ == "__main__":
    main()
