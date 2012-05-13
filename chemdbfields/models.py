#!/usr/bin/env python
# encoding: utf-8
"""
Some basic super simple models to capture common things with chemicals.

A Model for a molecule, and another for an atom.

File: models.py
Created 2012-01-24

"""
from fields import  AtomMapSmilesField, IsoSmilesField
from django.db import models

class BaseCDBFieldsMeta:
    """Simple meta class to be subclassed in spresidb db models"""
    app_label='chemdbfields'


class Molecule(models.Model):
    """Simple representation of a molecule"""
    iso_smiles = IsoSmilesField(max_length=400, unique=True)


    class Meta(BaseCDBFieldsMeta):
        pass;

    def __unicode__(self):
        return self.iso_smiles



class AtomMappedMolecule(models.Model):
    """Simple rep to capture an atom mapped smiles"""
    atm_mapped_smiles = AtomMapSmilesField(max_length=400, unique=True)
    molecule = models.ForeignKey(Molecule)

    class Meta(BaseCDBFieldsMeta):
        pass;
    
    def __unicode__(self):
        return self.atm_mapped_smiles

    @classmethod
    def fromAtmMappedSmiles(cls, atm_mapped_smiles):
        """Convenience that handles the look up of the Molecule."""
        mol, isNew = Molecule.objects.get_or_create(iso_smiles=atm_mapped_smiles)

        aMol, isNew = cls.objects.get_or_create(atm_mapped_smiles=atm_mapped_smiles,
                                                molecule=mol)
        return aMol



from django.db.models.signals import pre_save 

def validate_model(sender, **kwargs): 
    if 'raw' in kwargs and not kwargs['raw']: 
        kwargs['instance'].full_clean() 
pre_save.connect(validate_model, sender=AtomMappedMolecule, dispatch_uid='validate_models')
pre_save.connect(validate_model, sender=Molecule, dispatch_uid='validate_models') 

