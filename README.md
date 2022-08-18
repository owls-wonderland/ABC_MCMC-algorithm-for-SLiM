# About this repository

This repository contains self-contained script created for my research project for the Undergraduate Research Support Scheme at the University of Warwick in the summer 2021. The presentation poster is available at [URSS Showcase 2021 website](https://urss.warwick.ac.uk/items/show/160). Please familiarise yourself with the poster to find out about the project.

#### Contents

- `background.slim` - script for Background selection model in SLiM. SLiM is an evolutionary simulation package that provides facilities for very easily and quickly constructing genetically explicit individual-based evolutionary models.
- `findingsfs.py` - python script that executed ABC-MCMC to find parameters for SLiM Background model.
- `requirements.txt` - file contains libraries required for python script
- `run_mcmc.sh` - bash script to run `findingsfs.py`

## Before executing the script, SLiM needs to be installed

Installation of SLiM on Ubuntu is described in [SLiM manual](http://benhaller.com/slim/SLiM_Manual.pdf) in section 2.2

There are 2 methods for that:

#### 1. Using installation script

There is an installation script created by Bryce Carson. The script is hosted in the [SLiM-Extras](https://raw.githubusercontent.com/MesserLab/SLiM-Extras/master/installation/DebianUbuntuInstall.sh) repository on GitHub. Note
that it requires cmake, qmake, Qt, and either curl or wget to be installed on your system first; if they
are not, it will print instructions on how to install them. Follow those instructions to install the
needed packages, and then try running the install script again. 
Note that because this script installs the built components into /usr/bin, the sudo command is
used to run the script with root privileges, so you will need to enter your system’s root password. 

If you have wget installed (which appears to be present by default on Ubuntu), you can run the
install script directly from the web with the following single command line (which should be
entered as a single line in Terminal):

'''
wget --output-document /dev/stdout - --quiet https://
raw.githubusercontent.com/MesserLab/SLiM-Extras/master/installation/
DebianUbuntuInstall.sh | sudo bash -s
'''

or if you have curl installed and prefer to use it, you can execute this single command line to
install directly from the web:
'''
curl --silent https://raw.githubusercontent.com/MesserLab/SLiM-Extras/
master/installation/DebianUbuntuInstall.sh | sudo bash -s
'''

These two options should be identical for all practical purposes. They download the install
script from its URL and pipe it directly into bash to be executed. 

#### 2.Building SLiM from sources on Linux 

Refer to section 2.2.2 in [SLiM manual](http://benhaller.com/slim/SLiM_Manual.pdf)

## Installing pyenv(in case it is not installed)

Installation is performed on a per-user basis, using the [pyenv-installer](https://github.com/pyenv/pyenv-installer):

`curl https://pyenv.run | bash`

We then need to configure our bash profile for the pyenv installation. Replace `<user>` below with your username.

```bash
echo '
export PATH="/home/<user>/.pyenv/bin:$PATH"
eval "$(pyenv init -)"
eval "$(pyenv virtualenv-init -)"' >> ~/.bashrc
```

Finally, activate the updated bash profile:

`source ~/.bashrc`

You should now be able to execute pyenv instructions:

`pyenv`


## Usage

All files from this branch should be saved in one directory. 

Navigate to the directory, where files are saved. Once in that directory, give permision to bash file:

`chmod 755 run_mcmc.sh`

Then all that have to be done is to run bash file:

`./run_mcmc.sh`

Bash file should take care of everything. It will create virtual environment and download all requared packages for python code to run.
Bash file also gives directory's path such that all future files would be created in the same directory.

Note that after the first usage of bash file the virtual environment has been created and all required packages has been installed. Therefore, to avoid repeated installation lines 2 & 6 can be commented out:

'2 #pyenv virtualenv 3.8.2 mcmc-slim'

'6 #pip install -r requirements.txt'


