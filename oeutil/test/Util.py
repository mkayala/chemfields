#!/usr/bin/env python
# encoding: utf-8
"""
Logging for the test module.

File: Util.py
Author: MKayala
Created 2012-05-12
Copyright 2012. 

"""
import Const
import sys, os
import logging
import unittest

log = logging.getLogger(Const.APPLICATION_NAME)
log.setLevel(Const.LOGGER_LEVEL)

handler = logging.StreamHandler(sys.stderr)
formatter = logging.Formatter(Const.LOGGER_FORMAT)

handler.setFormatter(formatter)
log.addHandler(handler)
