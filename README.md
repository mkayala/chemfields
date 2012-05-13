chemfields: Django interface module for OpenEye OEChem objects
===========================================================

This repository is a Django module to provide a database interface 
for chemical datastructures using the [OpenEye OEChem 
Toolkit](http://eyesopen.com/oechem-tk).  Essentially, provides
an interface mapping OEChem structures to/from canonical 
[SMILES](http://www.daylight.com/dayhtml/doc/theory/theory.smiles.html)
strings in a database.  Handles special canonicalization
issues like kekulization and atom mapping.                                        

Provides: 

* Django model fields for Molecule and Atom Mapped Molecules
  with options for having specified kekulization and aromatization.
* Example Django models for Molecule, AtomMappedMolecule
* Utilities to canonicalize kekulizations, atom mappings.

Written by Matt Kayala at University of California, Irvine.  Code 
licensed under GPL.

Requirements
------------

* Python 2.7
* Django 1.3+
* OpenEye OEChem Toolkit 1.7.6+ (Note non-free software!)

Usage
----- 

Look in the chemfields/tests.py and oeutil/test/ module for some examples.  To run the tests

  python manage.py test
  python oeutil/test/TestCanonicalAtomMapSmiles.py

Acknowledgements
----------------

Some of the oeutil module code is based off of low-level utilities developed by
Jonathan Chen and Josh Swamidass at the University of California, Irvine.