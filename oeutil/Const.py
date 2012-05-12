#!/usr/bin/env python
# encoding: utf-8
"""
Module level constants

File: Const.py
Author: MKayala
Created 2012-05-12
Copyright 2012. 



"""

import logging
import Env

"""Application name","for example to identify a common logger object"""
APPLICATION_NAME = "chemdbfields.oeutil"

"""Default level for application logging.  Modify these for different scenarios.  See Python logging package documentation for more information"""
LOGGER_LEVEL = Env.LOGGER_LEVEL

"""Default format of logger output"""
LOGGER_FORMAT = "[%(asctime)s %(levelname)s] %(message)s"

"""SMILES separator to indicate distinct molecules"""
SMILES_MOL_DELIM = "."
