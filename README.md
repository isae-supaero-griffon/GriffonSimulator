The Simulator is available at https://github.com/isae-supaero-griffon/GriffonSimulator

You should use the most recent branch named PIR-1D-Felix-Dev
We started a new branch named AEther 0.1 for the purpose of improving the code for use in project AEther

To understand this project better, you can start by reading the 2-page Summary. Then you can read the first and second paragraph of Final Version 2 Semester Project Report. And then the third paragraph of S3_Prject_Report.

# GriffonSimulator
GriffonSimulator aims to provide a simulation of the Griffon hybrid rocket engine developped in the Supaero Space Section.

# Installation
## On Windows
While it's technically possible to run the simulator on Windows natively, this configuration hasn't been tested. The recommended way to run GriffonSimulator on Windows is therefore to either install WSL or a Linux virtual machine and follow the Linux instructions below.

## On Linux
The following instructions might work on macOS too. As a package manager, use [brew](https://brew.sh).
This guide is written for Debian/Ubuntu, if you are running a different distribution please change the package installation commands to use the distro's own package manager.

Make sure to have Python3 installed, by running
```bash
python --version
```
If that returns 2.X.X, try replacing `python` with `python3`.

### Install conda
Go to the [installers page](https://www.anaconda.com/products/individual#linux) and scroll to the list at the bottom.
Copy the link to the Linux installer for your architecture (usually x86, if the installation fails try ARM64 instead), and download it
```bash
sudo apt install libgl1-mesa-glx libegl1-mesa libxrandr2 libxrandr2 libxss1 libxcursor1 libxcomposite1 libasound2 libxi6 libxtst6
wget "https://repo.anaconda.com/archive/Anaconda3-2021.05-Linux-x86_64.sh"
```
Install conda by running
```bash
sudo chmod +x Anaconda3-2021.05-Linux-x86_64.sh
./Anaconda3-2021.05-Linux-x86_64.sh
```
After scrolling through the license, type "yes", then choose where to install it (your home directory is alright, just press enter).
When asked whether to initialize Anaconda3, type "yes".

### Create a new conda environment
Optional step.  
Run (you can replace the python version with the output of `python3 --version`):
```bash
conda create -n <insert name of the new environment here> python=3.8
```
You will need to execute the command
```bash
conda activate <insert name of the environment previously created here>
```
each time you wish to work on this project.

### Install needed packages
Run
```bash
conda update --all
conda install matplotlib numpy scipy pandas pip
pip install rocketcea
python -c "from rocketcea.cea_obj import CEA_Obj; C=CEA_Obj(oxName='LOX', fuelName='LH2'); print(C.get_Isp())"
```
The last command should print `374.30361765576265` as its output.
Then, run
```bash
conda install -c conda-forge opencv
```
Now go to the folder where you want to save the simulator:
* Linux:
```bash
cd /home/<your linux username>/
```
* WSL (if you want the project folder to show up in Windows, else follow the Linux instructions)
```bash
cd /mnt/c/Users/<your windows username>/
```
And run (on WSL you may need to prepend `sudo` to `git` commands in case of errors)
```bash
git clone "https://github.com/isae-supaero-griffon/GriffonSimulator.git"
cd GriffonSimulator
git checkout AEther-0.1
export PYTHONPATH="${PYTHONPATH}:$(pwd)"
```
To avoid having to type the last command every time you open a new terminal, run
```bash
pwd
```
and take note of the output. Then run
```bash
nano ~/.bashrc
```
and append `export PYTHONPATH="${PYTHONPATH}:<output of pwd>"` to it (make sure to replace the last part).

### Configure pycharm
Optional step.  
The professional edition is recommended (check with your school on whether they offer free licenses).  
In pycharm, close any open project. In the "Project" tab, click "Get from VCS" (Version Control System), then synchronize pycharm with your github account and select the Griffon Simulator.
Open it, then click File->Setting->Project:"..."->Python interpreter->click on the gearwheel->add->Conda environment->Existing environment->Select the one you just created

Then in the Python interpreter tab of the Settings pannel, you should see every package that you installed (including numpy, opencv, but maybe not rocketcea ?)

# Get started
To test the simulator, run (use `python3` if the output of `python --version` isn't version 3 or later)
```bash
cd tests
python testIntegration.py
```
The simulator will calculate the parameters at different timesteps and then display the results in several graphs.
# Any docs?
Several reports from Jose Felix Zapata are available on the S³ nextcloud (5.4 Commande Moteurs) explaining how he conceived this simulator.
# What next?

* Provide a user-friendly way to interact with the simulator.
* Writing an extensive documentation of the inner working of the program.
* Extend its simulation capabilities.
