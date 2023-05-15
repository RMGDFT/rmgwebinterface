import streamlit as st
import os
import sys
import string
import copy
import math
from datetime import datetime
from optparse import OptionParser, OptionGroup
import warnings
import tarfile

class config_part():
  def __init__(self, dirname, a,b,c, nx, ny, nz, num_atom, num_orb):
    self.dirname = dirname
    self.a = a
    self.b = b
    self.c = c
    self.nx = nx
    self.ny = ny
    self.nz = nz
    self.num_atoms = num_atom
    self.num_orb = num_orb
