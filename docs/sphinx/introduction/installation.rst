Installation
=================================

Just follow these steps:

0. Install Matlab. If you are at BU, you can get it for free `from
   here <http://www.bu.edu/tech/services/cccs/desktop/distribution/mathsci/matlab/>`__
1. Download the DynaSim toolbox using either:
      1. `Zip
         download <https://github.com/DynaSim/DynaSim/archive/master.zip>`__
      2. `Git clone <https://github.com/DynaSim/DynaSim.git>`__ (so you
         have access to all the latest updates)
   -  If using Git on the Mac/Linux command line (aka Terminal), run the
      command ``git clone https://github.com/dynasim/dynasim.git``
   -  For updates using git, run the command ``git pull`` from inside
      the DynaSim directory.
   -  Alternatively, you can also use `GitHub
      Desktop <https://desktop.github.com/>`__

2. Create a file named ``startup.m`` in your ``Documents/MATLAB`` folder
3. Add this line to the ``startup.m`` file:
   ``addpath(genpath(fullfile('your', 'custom', 'path', 'to', 'DynaSim')))``.

   -  Example: If I downloaded the DynaSim code to
      ``/home/austin/Dropbox/DynaSim``, then I would add the following
      to the ``startup.m`` file:
      ``addpath(genpath(fullfile('/home/austin/Dropbox/DynaSim')``
   -  IMPORTANT: Make sure that your startup file DOES NOT change the
      current directory (i.e., no uses of ``cd`` to change directory).
      This prevents problems during submission of batch jobs on a
      cluster.

Optional: If you want to use your own models/mechanisms, you can put
them inside the ``models/personal`` folder inside the DynaSim code
folder, or inside any folder on the MATLAB path.

That's it! This should work identically on your local computer AND the
SCC cluster.
