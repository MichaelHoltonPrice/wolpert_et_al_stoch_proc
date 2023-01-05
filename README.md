# Overview
This repository contains code to generate a number of the figures for the
following article:

Wolpert et al. -- The past as a stochastic process

This README is for Windows Powershell, but will work with minor
modifications for other platforms or approaches.

# Setup
Clone this repository and change directory into the new folder:

```bash
git clone https://github.com/MichaelHoltonPrice/wolpert_et_al_stoch_proc
cd wolpert_et_al_stoch_proc
```

Create and activate a virtual environment called stoch:

```bash
python -m venv stoch
Set-ExecutionPolicy Unrestricted -Scope Process
.\stoch\Scripts\activate
```

Install requirements:

```bash
pip install -r requirements.txt
```

In addition, the Python package bighist must be installed; its github
repository is here:

https://github.com/MichaelHoltonPrice/bighist

Change directory one level up, clone the bighist repository, change directory
into bighist, install the requirements for bighist, and install bighist:

```bash
cd ..
git clone https://github.com/MichaelHoltonPrice/bighist
cd bighist
pip install -r requirements.txt
python setup.py install
```

Change directory back to wolpert_et_al_stoch_proc, then proceed to the next
section to run the analysis code.

```bash
cd ..
cd wolpert_et_al_stoch_proc
```

# Running the code
The working directory must be wolpert_et_al_stoch_proc. This will be the case
if the steps in the preceding section have been followed. If necessary,
activate the virtual environment

```bash
Set-ExecutionPolicy Unrestricted -Scope Process
.\stoch\Scripts\activate
```