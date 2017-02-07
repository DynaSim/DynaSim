# DynaSim

[![Join the chat at https://gitter.im/DynaSim/DynaSim](https://badges.gitter.im/DynaSim/DynaSim.svg)](https://gitter.im/DynaSim/DynaSim?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)

DynaSim toolbox for modeling and simulating dynamical systems in Matlab or Octave

Installation:

1. Download the DynaSim toolbox: `git clone
   https://github.com/dynasim/dynasim.git`
2. Create a file named `startup.m` located in
    - If Mac/Linux `<home folder>/Documents/MATLAB`
    - If Windows `<home folder>\Documents\MATLAB`
3. Put the below in the file:
    - If Mac/Linux `addpath(genpath('/path/to/dynasim'))`
    - If Windows `addpath(genpath('\path\to\dynasim'))`

Documentation:
- Get started with the demos: [demos/demos.m](https://github.com/DynaSim/DynaSim/blob/master/demos/demos.m)
- For more details and examples walk through the tutorial: [demos/tutorial.m](https://github.com/DynaSim/DynaSim/blob/master/demos/tutorial.m)
- Example modeling projects: [PFC networks](https://github.com/jsherfey/PFC_models), [Thalamus](https://github.com/asoplata/ching2010_tcre_dynasim_mechanisms)

Mailing lists:
- Join the [user mailing list](https://groups.google.com/forum/#!forum/dynasim-users) to ask questions, request features, report bugs, and discuss DynaSim related issues
- Join the [developer mailing list](https://groups.google.com/forum/#!forum/dynasim-developers) or email [Jason Sherfey](http://jasonsherfey.com/) if you are interested in becoming a collaborator on DynaSim

Branches:
+ master - a direct fork of Jasons code
+ dave_mods - the modded branch

This fork contains Dave Stanley's customizations to DynaSim. Changes incliude:

Core simulator
+ **SimulateModel.m**: 'parallel_flag' now works (it will produce some warnings due to lock files on server, but these can be ignored)

Plotting
+ **PlotFR2.m**: New plotting command, an expansion of Jason's PlotFR command with better functionality for handling parameter sweeps (only works with plotting single variables for now).

Utilites (functions that operate on DynaSim Data Structures)
+ **ThevEquiv.m**: Calculates the Thévenin equivalent voltage and conductance for a given set of M specified ionic channels.
+ **CalcAverages.m**: which operates on a DynaSim data structure and averages over all cells. This structure can then be passed to PlotData as a means to force PlotData to produce average plots.
+ **DownsampleData**: Downsamples all data in a DynaSim data structre
+ **CalcSumOverFields.m**: Creates a new field that is the sum of a bunch of other fields (specified by the "fields" cell array). Useful for adding multiple ionic currents together.

