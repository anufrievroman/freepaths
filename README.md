# FreePATHS - Free Phonon And THermal Simulator

This Monte Carlo algorithm simulates trajectories of phonons in 3D models of nanostructures, which consists of a box with holes or pillars of various shapes. The algorithm outputs phonon trajectories, heat fluxes, temperature maps and profiles, the thermal conductivity, scattering maps and statistics and other information. See [the wiki pages](https://github.com/anufrievroman/Monte-Carlo/wiki/General-algorithm-flow) for the details of the simulation.

![Screenshot](https://github.com/anufrievroman/Monte-Carlo/blob/master/screenshot.png)


## Installation

FreePATHS requires python 3. On Linux and MacOS, it is probably already installed. On Windows, you may choose to install [Anaconda](https://www.anaconda.com) package, which will install everything for you.

Install the package from PyPi repository by entering this command into a terminal or a python console:

`pip install --upgrade freepaths`

For more details about installation, see [this wiki page](https://github.com/anufrievroman/freepaths/wiki/Installation).

## Usage

FreePATHS is a command line application, so it runs inside Linux, MacOS, or Windows terminal. It takes an input file from the user, which contains all the settings, and outputs the results in a new folder.

There are two modes of using the program. Main mode traces a large number of phonons through a structure and collects statistics about their paths. The MFP sampling mode measures phonon mean free paths using a small number of phonons and calculates the thermal conductivity by integrating phonon dispersion. 

### Demo

If you simply run `freepaths` without specifying an input file, the program will run a demo simulation.

### Main mode

In the main mode, the program traces large number of phonons through a structure and calculates various statistical distributions and maps. In this mode, the thermal conductivity will be calculated via Fourier law. Note that the thermal conductivity will be correct only in the absence of holes.

Run the program as:

`freepaths your_input_file.py`

See [wiki page](https://github.com/anufrievroman/freepaths/wiki/Creating-input-files) for explanations about creating your own input files. In the [examples](https://github.com/anufrievroman/freepaths/tree/master/examples) folder, you will find examples of various input files. Try using one of them, for instance as:

`freepaths simple_nanowire.py`

After the simulation, see the results in a newly created **Results** folder.


### MFP sampling mode

Alternatively, you can run FreePATHS in the mean free path sampling mode, which is designed to calculate the thermal conductivity by integrating phonon dispersion. To run the program in this mode, it is advised to reduce the number of phonons to about 30 and add `-s` flag in the command:

`freepaths -s simple_nanowire.py`

The calculated thermal conductivity will be output in the terminal. However, other statistical quantities and plots will still be calculated and output in the `Results` folder.


## Troubleshooting

- [Troubles with installation](https://github.com/anufrievroman/freepaths/wiki/Installation)
- [Troubles with usage](https://github.com/anufrievroman/freepaths/wiki/Usage)


## Disclaimer

The code is still in development and provided as is. It likely contains bugs or might be inappropriate for your research. It is your responsibility to understand the underlying physics, test the code, and verify that all the equations and the code are correct. See [the wiki pages](https://github.com/anufrievroman/Monte-Carlo/wiki/General-algorithm-flow) and the references below for more details on the code.


## References and acknowledgments

The code has been developed by [Roman Anufriev](https://anufrievroman.com) in [Nomura lab](https://www.nlab.iis.u-tokyo.ac.jp/index-e.html) at the University of Tokyo in 2018-2022.
If you would like to use this code for your research, please see the disclaimer above and consider citing the papers below, if it is appropriate.
Details of the code and examples of the output can be found in the following papers:

1. Anufriev et al. [Materials Today Physics 15, 100272 (2021)](https://www.sciencedirect.com/science/article/pii/S2542529320300961)
2. Anufriev et al. [Nanoscale, 11, 13407-13414 (2019)](https://pubs.rsc.org/en/content/articlehtml/2019/nr/c9nr03863a)
3. Anufriev et al. [ACS Nano 12, 11928 (2018)](https://pubs.acs.org/doi/abs/10.1021/acsnano.8b07597)
4. Huang et al. [ACS Applied Materials & Interfaces 12, 25478 (2020)](https://pubs.acs.org/doi/10.1021/acsami.0c06030)

Development of this code was funded by the following grants:

- PRESTO JST (No. JPMJPR19I1)
- CREST JST (No. JPMJCR19Q3)
- Kakenhi (15H05869, 15K13270, and 18K14078)
- Postdoctoral Fellowship program of Japan Society for the Promotion of Science.
