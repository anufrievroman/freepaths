"""FreePATHS - Free Phonon and THermal Simulator"""

import argparse
import colorama
from colorama import Fore, Style


__version__ = "2.2"

colorama.init()

# Parse user arguments:
parser = argparse.ArgumentParser(
                prog = 'FreePATHS',
                description = 'Phonon Monte Carlo simulator',
                epilog = 'For more information, examples, and tutorials, visit: https://anufrievroman.gitbook.io/freepaths'
                )
parser.add_argument('input_file', nargs='?', default=None, help='The input file')
parser.add_argument("-s", "--sampling", help="Run in MFP sampling mode", action="store_true")
parser.add_argument("-e", "--electron", help="Run simulation for electrons", action="store_true")
args = parser.parse_args()


def run():
    """Run the program depending on the mode"""
    print(f"\n{Fore.BLUE}FreePATHS v{__version__}{Style.RESET_ALL}")
    if args.sampling:
        import freepaths.main_mfp_sampling
        freepaths.main_mfp_sampling.main(args.input_file)

    elif args.electron:
        import freepaths.main_tracing
        from freepaths.particle_types import ParticleType
        freepaths.main_tracing.main(args.input_file, ParticleType.ELECTRON)

    else:
        import freepaths.main_tracing
        from freepaths.particle_types import ParticleType
        freepaths.main_tracing.main(args.input_file, ParticleType.PHONON)


if __name__ == "__main__":
    run()
