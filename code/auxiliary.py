import math
import time
import sys
import threading
import itertools


def log_background_prob(sequence, background_probs):
    log_background = 0.0
    for symbol in sequence:
        if (symbol in background_probs) and (symbol != "-"):
            log_background += math.log(background_probs[symbol])
    return log_background

def compile_codons(PI_codons, RT_codons):
    RT_codons_acc_PI = [x + 99 for x in RT_codons]
    final_codons = PI_codons + RT_codons_acc_PI
    return final_codons

def print_progress_bar(current, total, message, bar_length=20):
    """
    Prints and updates a progress bar in the console on the same line.

    Args:
        current (int): The number of processes completed.
        total (int): The total number of processes.
        bar_length (int): The character width of the progress bar (default 20).
    """
    # Protect against division by zero or invalid values
    if total <= 0:
        total = 1

    fraction = current / total
    # Ensure the fraction never goes beyond 1
    if fraction > 1:
        fraction = 1.0

    # Calculate the number of '#' characters
    filled = int(bar_length * fraction)
    # Build the bar
    bar = '#' * filled + ' ' * (bar_length - filled)
    # Calculate percentage
    percent = fraction * 100

    # Write the progress bar to stdout, using \r to return to the start of the line
    sys.stdout.write(f"\r{message} [{bar}] {percent:6.2f}%")
    sys.stdout.flush()
    # Optionally, when done, print a newline
    if current >= total:
        print()  # Moves to the next line once complete

def codons_prompt(prompt):
    return sorted([int(x) for x in prompt.replace(" ","").split(",")])

def verify_codon_format(ls):
    for x in ls:
        if not isinstance(x, int):  # Check if the element is not an integer
            return False
    return True  # All elements are integers


class DotAnimation:
    def __init__(self, message="Processing", interval=0.5):
        self.message = message
        self.interval = interval
        self._stop = threading.Event()
        self._thread = threading.Thread(target=self._run, daemon=True)

    def start(self):
        self._thread.start()

    def stop(self):
        self._stop.set()
        self._thread.join()
        # Optionally overwrite the line with spaces and a carriage return
        sys.stdout.write("\r" + " " * 80 + "\r")
        sys.stdout.flush()

    def _run(self):
        frames = [".", ". .", ". . ."]
        prev_length = 0
        for frame in itertools.cycle(frames):
            if self._stop.is_set():
                break

            # Build the new text
            new_text = f"{self.message} {frame}"

            # Overwrite old text fully
            sys.stdout.write("\r" + new_text + " " * (prev_length - len(new_text)))
            sys.stdout.flush()

            prev_length = len(new_text)
            time.sleep(self.interval)
