#!/usr/bin/python

import os

dir0=os.getcwd()
os.system("echo export PATH=" + dir0 + "/bin:\$PATH >> ~/.bashrc")
os.system("echo export PYTHONPATH=" + dir0 + "/bin:\$PYTHONPATH >> ~/.bashrc")
os.system("source ~/.bashrc")
