"""FreePATHS - Free Phonon and THermal Simulator"""

import sys
import argparse
import colorama
from colorama import Fore, Style

from freepaths.particle_types import ParticleType

__version__ = "2.3.3"

colorama.init()

# Parse user arguments:
parser = argparse.ArgumentParser(
                prog = 'FreePATHS',
                description = 'Phonon Monte Carlo simulator',
                epilog = 'For more information, examples, and tutorials, visit: https://anufrievroman.gitbook.io/freepaths'
                )
parser.add_argument('input_file', nargs='?', default=None, help='The input file')
parser.add_argument("-v", "--version", action="version", version=f"FreePATHS {__version__}")
parser.add_argument("-s", "--sampling", help="Run in phonon MFP sampling mode", action="store_true")
parser.add_argument("-e", "--electron", help="Run simulation for electrons", action="store_true")
parser.add_argument("--demo", help="Run a demo simulation with default parameters", action="store_true")
args = parser.parse_args()


def run():
    """Run the program depending on the mode"""
    print(f"\n{Fore.BLUE}FreePATHS v{__version__}{Style.RESET_ALL}")

    if not args.input_file and not args.demo:
        print("\nUsage:")
        print(f"  freepaths <input_file>              run phonon simulation")
        print(f"  freepaths <input_file> -e           run electron simulation")
        print(f"  freepaths <input_file> -s           run in MFP sampling mode")
        print(f"  freepaths --demo                    run a demo with default parameters")
        print(f"  freepaths --version                 show version")
        print(f"\nExample input files:")
        print(f"  https://github.com/anufrievroman/freepaths/tree/master/examples")
        print(f"\nDocumentation:")
        print(f"  https://anufrievroman.gitbook.io/freepaths")
        sys.exit(0)

    try:
        if args.sampling:
            import freepaths.main_mfp_sampling
            freepaths.main_mfp_sampling.main(args.input_file, ParticleType.PHONON)

        elif args.electron:
            import freepaths.main_tracing
            freepaths.main_tracing.main(args.input_file, ParticleType.ELECTRON)

        else:
            import freepaths.main_tracing
            freepaths.main_tracing.main(args.input_file, ParticleType.PHONON)

    except KeyboardInterrupt:
        print(f"\n{Fore.RED}Simulation interrupted by user.{Style.RESET_ALL}")
        sys.exit(0)


if __name__ == "__main__":
    run()
