# cod-model



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

`./run_msms.sh`

Bash file should take care of everything. It will create virtual environment and download all requared pachaged for python code to run.
Bash file also gives directory's path such that all future files would be creted in the same directory.
