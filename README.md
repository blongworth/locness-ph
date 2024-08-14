# LOCNESS-pH

LOCNESS-pH contains code to take polled readings from an MBARI mFET pH sensor/logger, store them in a file and database, and convert raw readings to pH.

mFET hardware and internal software by Yui Takeshita et al.

Calculation code is from https://github.com/SUPScientist/pH-with-the-Honeywell-Durafet which is in turn from https://pubs.acs.org/doi/abs/10.1021/acs.analchem.5b04653

# Installation

After cloning the github repo, create a python venv for the project and install dependencies.

# Use

Determine the correct serial port for the mFET using a terminal program (terraterm, CoolTerm).
Set serial port and other parameters near the top of `read_ph_sensor.py`
Run the python program. Data and communications should be logged to the terminal, and data should be added to the files specified.
Quit using `Ctrl-C`.
