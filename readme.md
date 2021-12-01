# Phonon Monte Carlo

The algorithm simulates trajectories of phonons in 3D models of nanostructures, which consists of a box (membrane) with holes and/or pillars of various shapes. The algorithm can calculate the heat fluxes, temperature profiles, the thermal conductivity of the structure, and lots of other things. See [the wiki pages](https://github.com/anufrievroman/Monte-Carlo/wiki/General-algorithm-flow) for the details of the simulation.

![Screenshot](https://github.com/anufrievroman/Monte-Carlo/blob/master/screenshot.png)

## Installation and usage

- Install python 3
- Download this repository
- Adjust parameters in `parameters.py` file and run the `monte-carlo.py` file. You will find the results in a newly created folder.

## Disclaimer

The code its provided as is. It may contain bugs or be inappropriate for your use case. It is your responsibility to understand the underlying physics, test the code for your use case, and verify that all the equations and the code are correct.

## Credits

The code is created by [Roman Anufriev](https://anufrievroman.com) in the University of Tokyo in 2018-2021. If you use this code for your research, please cite the papers below, if it is appropriate. If you need a consultation, I'll be happy to help you, but then you'll need to consider adding me as a co-author.

## References

- Anufriev et al. [Materials Today Physics 15, 100272 (2021)](https://www.sciencedirect.com/science/article/pii/S2542529320300961)
- Anufriev et al. [Nanoscale, 11, 13407-13414 (2019)](https://pubs.rsc.org/en/content/articlehtml/2019/nr/c9nr03863a)
- Anufriev et al. [ACS Nano 12, 11928 (2018)](https://pubs.acs.org/doi/abs/10.1021/acsnano.8b07597)
