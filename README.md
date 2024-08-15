# LOCNESS-pH

LOCNESS-pH contains code to take polled readings from an MBARI mFET pH sensor/logger, store them in a file and database, and convert raw readings to pH.

mFET hardware and internal software by Yui Takeshita et al.

Calculation code is from https://github.com/SUPScientist/pH-with-the-Honeywell-Durafet which is in turn from https://pubs.acs.org/doi/abs/10.1021/acs.analchem.5b04653

# Installation

After cloning the github repo:

* Create a python venv for the project: 'python -m venv .venv`. 
* Activate the venv: `source .venv/bin/activate`. 
* Install dependencies. `pip install -r requirements.txt`.
* Copy `default_config.yaml` to `config.yaml` and edit as needed, 
especially serial port and file config. 

# Use

Run the python program. Data and communications should be logged to the terminal and logfile, and data should be added to the files specified.

Quit using `Ctrl-C`.
