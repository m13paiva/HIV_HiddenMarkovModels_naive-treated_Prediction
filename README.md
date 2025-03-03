# HIV-1 ART-naïve/treated prediction with HMM - v1.0

## SETUP
1. Navigate to the directory `code`:
   ```bash
   cd your_directory\HIV_HiddenMarkovModels_naive-treated_Prediction\code
   ```
2. Create a virtual environment (Python 3.12 required):
   ```bash
   python3.12 -m venv hiv1_hmm
   ```
3. Activate the virtual environment:
   - On Windows:
     ```bash
     .\hiv1_hmm\Scripts\activate
     ```
   - On macOS/Linux:
     ```bash
     source hiv1_hmm/bin/activate
     ```
4. Install required packages within the virtual environment:
   ```bash
   pip install -r requirements.txt
   ```
5. Verify the installation (optional):
   ```bash
   python --version   # Should display Python 3.12.x
   pip list           # To see installed packages
   ```

## USAGE
```bash
python hiv_hmm_pred.py [OPTIONS]
```

## DESCRIPTION
This program takes HIV-1 **pol** gene sequences (in FASTA format) and predicts if they belong to a patient who has experienced antiretroviral treatment (**ART-treated**) or not (**ART-naïve**), using Hidden Markov Models.

## OPTIONS
- `-h, --help`  
  Display this help message and exit.

- `-c, <IN_FILE>, --classify <IN_FILE>`  
  Classify sequences in the provided FASTA file.  
  **Sub-option**  
  - `-o, <OUT_FILE>, --output <OUT_FILE>`  
    Specify an output CSV file path for classification.  
    *(Default: `output_data/output.csv`)*

- `-t, --test`  
  Test performance of HMM algorithm.

## ARGUMENTS
- **IN_FILE**  
  Fasta file with sequences to classify.

- **OUT_FILE**  
  Personalized CSV for output data.

---

**Developed by**: Manuel Almeida  
**Email address**: [mpa.almeida@campus.fct.unl.pt](mailto:mpa.almeida@campus.fct.unl.pt)
