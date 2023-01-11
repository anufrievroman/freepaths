# FreePATHS - Free Phonon And THermal Simulator

This Monte Carlo algorithm simulates trajectories of phonons in 3D models of nanostructures, which consists of a box with holes or pillars of various shapes. The algorithm outputs phonon trajectories, heat fluxes, temperature maps and profiles, the thermal conductivity, scattering maps and statistics and other information. See [the wiki pages](https://github.com/anufrievroman/Monte-Carlo/wiki/General-algorithm-flow) for the details of the simulation.

![Screenshot](https://github.com/anufrievroman/Monte-Carlo/blob/master/screenshot.png)

## Installation

- Install python 3. On Linux and MacOS it is probably already installed. On Windows, you may choose to install [Anaconda](https://www.anaconda.com).
- Install `numpy matplotlib scipy` libraries. On linux and MacOS `pip install numpy matplotlib scipy`, on Windows it is already installed with [Anaconda](https://www.anaconda.com).
- Download this repository. Unzip it. You only need the `freepaths` folder.

## Usage

- Adjust parameters in `parameters.py` file.
- See `Examples` folder for typical structures and how you may build arrays of holes.
- Start the simulation by running `main.py` file. On Linux and MacOS `python main.py`. On Windows, open Spider that was installed with [Anaconda](https://www.anaconda.com), open `main.py` and run it.
- See the results in a newly created `Results` folder.

## Disclaimer

The code is provided as is. Note that it likely contains bugs or might be inappropriate for your use case. It is your responsibility to understand the underlying physics, test the code, and verify that all the equations and the code are correct. See [the wiki pages](https://github.com/anufrievroman/Monte-Carlo/wiki/General-algorithm-flow) and the references below for more details on the code.

## References

Details of the code and examples of the output can be found in the following papers:

1. Anufriev et al. [Materials Today Physics 15, 100272 (2021)](https://www.sciencedirect.com/science/article/pii/S2542529320300961)
2. Anufriev et al. [Nanoscale, 11, 13407-13414 (2019)](https://pubs.rsc.org/en/content/articlehtml/2019/nr/c9nr03863a)
3. Anufriev et al. [ACS Nano 12, 11928 (2018)](https://pubs.acs.org/doi/abs/10.1021/acsnano.8b07597)

## Credits

The code has been developed by [Roman Anufriev](https://anufrievroman.com) in [Nomura lab](https://www.nlab.iis.u-tokyo.ac.jp/index-e.html) at the University of Tokyo in 2018-2022. If you would like to use this code for your research, please cite the papers above, if it is appropriate.

## Acknowledgments

Development of this code was funded by the following grants:

- PRESTO JST (No. JPMJPR19I1)
- CREST JST (No. JPMJCR19Q3)
- Kakenhi (15H05869, 15K13270, and 18K14078)
- Postdoctoral Fellowship program of Japan Society for the Promotion of Science.
