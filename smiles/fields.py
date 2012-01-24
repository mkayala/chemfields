#!/usr/bin/env python
# encoding: utf-8
"""
Fields for different types of smiles.

File: fields.py
Created 2012-01-24
Copyright 2012. Institute of Genomics and Bioinformatics

"""
from django.db import models;
from django import forms;
from django.core import exceptions

from openeye.oechem import OEGraphMol, OEParseSmiles, OECreateIsoSmiString
from openeye.oechem import OEAssignAromaticFlags, OEClearAromaticFlags, OEKekulize
from openeye.oechem import OECanonicalOrderAtoms, OECanonicalOrderBonds
from CHEM.Common.Util import standardizeSmiles, createAtomMapSmiString
from CHEM.Common.Util import createStdAtomMapSmiString
from CHEM.Common.CanonicalAtomMapSmiles import (canonicalizeAtomMapSmiString,
                                                createCanonicalAtomMapSmiString
                                                )
from CHEM.Common.MolExt import clearAtomMaps, removeNonsenseStereo

from widgets import SmilesWidget

class InvalidSmilesException(Exception):
    """Basic class for invalid smiles"""
    pass;

def canonicalKekule(smi):
    """Function to make sure that the kekule structure is canonical"""
    mol = value_to_oemol(smi)
    removeNonsenseStereo(mol)
    OEAssignAromaticFlags(mol)
    for bond in mol.GetBonds():
        if bond.IsAromatic():
            bond.SetIntType(5)
    OECanonicalOrderAtoms(mol)
    OECanonicalOrderBonds(mol)
    OEClearAromaticFlags(mol)
    OEKekulize(mol)
    return createCanonicalAtomMapSmiString(mol)


def value_to_oemol(value):
    """Convenience to turn a value into an oemol
    Raises proper validation errors on parsing problems

    Note:  Have to coerce the value to non-unicode string for OEChem
    compatability.

    Note the canon value.  Set this to TRUE to have the OEParseSmiles
    NOT Kekulize.  
    """
    if value == '':
        return None;
    mol = OEGraphMol()
    try:
        success = OEParseSmiles(mol, str(value), False, True)
        if not success:    
            raise InvalidSmilesException(value);
    except (InvalidSmilesException, ):
        raise exceptions.ValidationError('Value %r is an invalid smiles' % value)

    return mol;


def clean_iso_smiles(value):
    """Basic function to clean a smiles to a canonical iso_smiles
    """
    mol = value_to_oemol(value);
    if mol is None:
        return ''
    clearAtomMaps(mol)
    removeNonsenseStereo(mol)
    OEAssignAromaticFlags(mol)
    return OECreateIsoSmiString(mol)

def clean_atom_map_smiles(value, doKekule=True, doArom=False):
    """Basic function to clean a smiles to a atom_map_smiles
    """
    mol = value_to_oemol(value)
    if mol is None:
        return '';
    removeNonsenseStereo(mol)
    if doArom:
        OEAssignAromaticFlags(mol)
    if doKekule:
        OEAssignAromaticFlags(mol)
        smi = createStdAtomMapSmiString(mol)
        smi = canonicalKekule(smi)
        mol = value_to_oemol(smi)
    return createCanonicalAtomMapSmiString(mol)


### Form fields        
class SmilesFormField(forms.CharField):
    """Basic setup for a Form field for smiles with validation"""
    def __init__(self, *args, **kwargs):
        """Setup to make sure running with the same widget"""
        kwargs['widget'] = SmilesWidget(attrs={'size':'80'})
        super(SmilesFormField, self).__init__(*args, **kwargs)
        
    def clean(self, value):
        clean_value = super(SmilesFormField, self).clean(value)
        return self.to_python(value)
    
class IsoSmilesFormField(SmilesFormField):
    """Class for handling form field for iso_smiles"""
    def to_python(self, value):
        """Basic iso smiles cleaning """
        return clean_iso_smiles(value)

class AtomMapSmilesFormField(SmilesFormField):
    """Class for handling form fields for atom_map_smiles"""
    def to_python(self, value):
        """Basic atm map smiles cleaning"""
        return clean_atom_map_smiles(value)

class AtomMapSmilesNonKekuleFormField(SmilesFormField):
    """Class for handling form fields for atom_map_smiles where we do NOT want
    canonical Kekulization"""
    def to_python(self, value):
        """Basic atm map smiles cleaning"""
        return clean_atom_map_smiles(value, doKekule=False)

class AtomMapSmilesAromaticFormField(SmilesFormField):
    """Class for handling form fields for atom_map_smiles where we do NOT want
    canonical Kekulization"""
    def to_python(self, value):
        """Basic atm map smiles cleaning"""
        return clean_atom_map_smiles(value, doKekule=False, doArom=True)


### DB Fields
class SmilesField(models.CharField):
    """Validate smiles"""
    __metaclass__ = models.SubfieldBase

    formfield_class = SmilesFormField
    
    default_error_messages = {
        'invalid_smiles': 'Value %r is an invalid Smiles'
        }

    def get_prep_value(self, value):
        """Simply call to_python"""
        return self.to_python(value)

    def formfield(self, **kwargs):
        """Set to use the form fields from below"""
        kwargs['form_class'] = self.formfield_class
        return super(SmilesField, self).formfield(**kwargs)

    
class IsoSmilesField(SmilesField):
    """Canonicalize Iso Smiles"""
    formfield_class = IsoSmilesFormField
    
    def to_python(self, value):
        """Canonicalize a non atom mapped smiles"""
        return clean_iso_smiles(value)
        
class AtomMapSmilesField(SmilesField):
    """Canonicalize atom mapped smiles"""
    formfield_class = AtomMapSmilesFormField
    
    def to_python(self, value):
        """Ensure canonical atom mapped smiles"""
        return clean_atom_map_smiles(value)

class AtomMapSmilesNonKekuleField(SmilesField):
    """Canonicalize atom mapped smiles without canonical Kekulization"""
    formfield_class = AtomMapSmilesNonKekuleFormField
    
    def to_python(self, value):
        """Ensure canonical atom mapped smiles"""
        return clean_atom_map_smiles(value, doKekule=False)


class AtomMapSmilesAromaticField(SmilesField):
    """Canonicalize atom mapped smiles without canonical Kekulization"""
    formfield_class = AtomMapSmilesAromaticFormField
    
    def to_python(self, value):
        """Ensure canonical atom mapped smiles"""
        return clean_atom_map_smiles(value, doKekule=False, doArom=True)
