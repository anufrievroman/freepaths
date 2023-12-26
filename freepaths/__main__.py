"""FreePATHS - Free Phonon and THermal Simulator"""

import argparse
import colorama
from colorama import Fore, Style

import freepaths.main_tracing
import freepaths.main_mfp_sampling

__version__ = "1.8"

colorama.init()

# Parse user arguments:
parser = argparse.ArgumentParser(
                prog = 'FreePATHS',
                description = 'Phonon Monte Carlo simulator',
                epilog = 'For more information, examples, and tutorials, visit: https://anufrievroman.gitbook.io/freepaths'
                )
parser.add_argument('input_file', nargs='?', default=None, help='The input file')
parser.add_argument("-s", "--sampling", help="Run in MFP sampling mode", action="store_true")
args = parser.parse_args()


def run():
    """Run the program depending on the mode"""
    print(f"\n{Fore.BLUE}FreePATHS v{__version__}{Style.RESET_ALL}")
    if args.sampling:
        freepaths.main_mfp_sampling.main(args.input_file)
    else:
        freepaths.main_tracing.main(args.input_file)


if __name__ == "__main__":
    run()
