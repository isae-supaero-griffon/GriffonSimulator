The Simulator is available at https://github.com/isae-supaero-griffon/GriffonSimulator

You should use the most recent branch named PIR-1D-Felix-Dev
We started a new branch named AEther 0.1 for the purpose of improving the code for use in project AEther

To understand this project better, you can start by reading the 2-page Summary. Then you can read the first and second paragraph of Final Version 2 Semester Project Report. And then the third paragraph of S3_Prject_Report.

# GriffonSimulator
GriffonSimulator aims to provide a simulation of the Griffon hybrid rocket engine developped in the Supaero Space Section.

# Installation
You need to install Python 3.7, RocketCEA, GCC, G+, GFORTRAN and OpenCV

## On Linux
The standard way to do it would be through a pycharm professional python distribution, for which you can get a free student license.

### Install conda
Download the installer at https://www.anaconda.com/products/individual#linux
Execute it : From the terminal type "./Anaconda3-2020.07-Linux-x86_64.sh" (with the name of the installer)
    After scrolling through the license, type "yes", then choose where to install (your home directory is nice)

### Create a new conda environment
In a terminal, type "conda create -n name_of_your_environment python=3.8" (or later version of python)
You will need to type "conda activate name_of_your_environment" each time you wish to work on this project

### Install the needed packages
In a terminal, type "conda install matplotlib numpy scipy pandas pip"
Then "pip install rocketcea"
  You can check the correct installation of rocketcea by running in command prompt :
  python -c "from rocketcea.cea_obj import CEA_Obj; C=CEA_Obj(oxName='LOX', fuelName='LH2'); print(C.get_Isp())"
  It Should print 374.30361765576265
Type "conda install -c conda-forge opencv" (cf https://anaconda.org/conda-forge/opencv)

### Configure pycharm
In pycharm, close any opened project. In the "Project" tab, click "Get from VCS" (Version Control System), then synchronize pycharm with your github account and select the Griffon Simulator.
Open it, then click File->Setting->Project:"..."->Python interpreter->click on the gearwheel->add->Conda environment->Existing environment->Select the one you just created

Then in the Python interpreter tab of the Settings pannel, you should see every package that you installed (including numpy, opencv, but maybe not rocketcea ?)

## On Windows
Please add information if you can !

# Get started
In order to run your first simulation, simply run the testIntegration.py. The simulator will calculate the parameters at different timesteps and then display the results in several graphs.
# Any docs?
Several reports from Jose Felix Zapata are available on the S³ nextcloud (5.4 Commande Moteurs) explaining how he conceived this simulator.
# What next?

* Provide a user-friendly way to interact with the simulator.
* Writing an extensive documentation of the inner working of the program.
* Extend its simulation capabilities.
