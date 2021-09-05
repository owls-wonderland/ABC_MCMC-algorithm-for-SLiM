#!/usr/bin/env bash
pyenv virtualenv 3.8.2 mcmc-slim
eval "$(pyenv init -)"
eval "$(pyenv virtualenv-init -)"
pyenv activate mcmc-slim
pip install -r requirements.txt
python findingsfs.py - $(pwd)

