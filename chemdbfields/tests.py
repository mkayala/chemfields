#!/usr/bin/env python
# encoding: utf-8
"""
Testing the general purpose chemdbfields.

In particular, we want to test that saving Kekule and non-Kekule version of the same atom mapped molecule should
give the same molecule back.

File: tests.py
Created 2012-03-01
Copyright 2012. Institute of Genomics and Bioinformatics


"""

import sys, os
from pprint import pformat
from django.test import TestCase
from django.core import exceptions

from Util import log;
from models import Molecule, AtomMappedMolecule
from fields import clean_atom_map_smiles

class TestMolecules(TestCase):
    """Testing all kinds of stuff about molecules.

    Diff ways of writing atommapped and iso smiles should give same mols back.
    """
    
    def setUp(self):
        """Basic Setup"""
        super(TestMolecules, self).setUp()
        
    
    def tearDown(self):
        """Basic Tear Down"""
        super(TestMolecules, self).tearDown()
        

    def test_isoSmilesCanonicalized(self):
        """Check that different iso smiles end up pointing to same mol"""
        smi1 = 'CC1CCC1'
        smi2 = 'C1C(C)CC1'
        smi3 = 'CC1C[CH2:1]C1'

        mol1 = Molecule(iso_smiles=smi1)
        mol2 = Molecule(iso_smiles=smi2)
        mol3 = Molecule(iso_smiles=smi3)

        self.assertEquals(mol1.iso_smiles, mol2.iso_smiles)
        self.assertEquals(mol1.iso_smiles, mol3.iso_smiles)

        log.info("Mols are %s" % pformat([mol1, mol2, mol3]))

        mol1.save()
        mol2.clean();
        mol3.clean();

        self.assertEquals(mol1.iso_smiles, mol2.iso_smiles)
        self.assertEquals(mol1.iso_smiles, mol3.iso_smiles)


    def test_isoSmilesKekuleCanonicalized(self):
        """Check diff iso smiles in diff kekule forms look same"""
        smi1 = 'C1=CC=CC=C1'
        smi2 = 'c1ccccc1'
        smi3 = 'c1cc[cH:1]cc1'

        mol1 = Molecule(iso_smiles=smi1)
        mol2 = Molecule(iso_smiles=smi2)
        mol3 = Molecule(iso_smiles=smi3)

        self.assertEquals(mol1.iso_smiles, mol2.iso_smiles)
        self.assertEquals(mol1.iso_smiles, mol3.iso_smiles)

        log.info("Mols are %s" % pformat([mol1, mol2, mol3]))

        mol1.save()
        mol2.clean();
        mol3.clean();

        self.assertEquals(mol1.iso_smiles, mol2.iso_smiles)
        self.assertEquals(mol1.iso_smiles, mol3.iso_smiles)


    def test_atommappedSmilesCanonicalized(self):
        """Check that atom mapped smiles get canonicalized """
        smi1 = 'CC1C[CH2:1]C1'
        smi2 = 'C1C(C)C[CH2:1]1'
        smi3 = 'CC1C[CH2:1]C1'

        amol1 = AtomMappedMolecule.fromAtmMappedSmiles(smi1);
        amol2 = AtomMappedMolecule.fromAtmMappedSmiles(smi2);
        amol3 = AtomMappedMolecule.fromAtmMappedSmiles(smi3);

        self.assertEquals(amol1.atm_mapped_smiles, amol2.atm_mapped_smiles)
        self.assertEquals(amol1.atm_mapped_smiles, amol3.atm_mapped_smiles)

        log.info("Mols are %s" % pformat([amol1, amol2, amol3]))

        amol1.save()
        amol2.clean();
        amol3.clean();

        self.assertEquals(amol1.molecule, amol2.molecule)
        self.assertEquals(amol1.molecule, amol3.molecule)


    def test_atmMapKekuleCanonicalize(self):
        """Check that atom maps from kekule, non kekule forms work out."""
        smi1 = 'C1=CC=[CH:1]C=C1'
        smi2 = 'c1[cH:1]cccc1'
        smi3 = 'c1cc[cH:1]cc1'

        amol1 = AtomMappedMolecule.fromAtmMappedSmiles(smi1);
        amol2 = AtomMappedMolecule.fromAtmMappedSmiles(smi2);
        amol3 = AtomMappedMolecule.fromAtmMappedSmiles(smi3);

        self.assertEquals(amol1.atm_mapped_smiles, amol2.atm_mapped_smiles)
        self.assertEquals(amol1.atm_mapped_smiles, amol3.atm_mapped_smiles)

        log.info("Mols are %s" % pformat([amol1, amol2, amol3]))

        amol1.save()
        amol2.clean();
        amol3.clean();

        self.assertEquals(amol1.molecule, amol2.molecule)
        self.assertEquals(amol1.molecule, amol3.molecule)

    def test_atmMapKekuleCanonicalizeLarge(self):
        """Check that atom maps from kekule, non kekule forms work out with a large setup."""

        #smi1 = '[CH3:23][c:22]1[cH:17][cH:14][c:13]([cH:15][cH:18]1)[S:16](=[O:19])(=[O:20])[OH:21].O'
        #smi2 = '[CH3:23][c:22]1[cH:17][cH:14][c:13]([cH:15][cH:18]1)[S:16](=[O:19])(=[O:20])[OH:21].O'
        #smi3 = '[CH3:23][C:22]1=[CH:18][CH:15]=[C:13]([CH:14]=[CH:17]1)[S:16](=[O:19])(=[O:20])[OH:21].O'

        smi1 = '[CH3:23][c:22]1[cH:17][cH:14][c:13]([cH:15][cH:18]1)[S:16](=[O:19])(=[O:20])[OH:21]'
        smi2 = '[CH3:23][c:22]1[cH:17][cH:14][c:13]([cH:15][cH:18]1)[S:16](=[O:19])(=[O:20])[OH:21]'
        smi3 = '[CH3:23][C:22]1=[CH:18][CH:15]=[C:13]([CH:14]=[CH:17]1)[S:16](=[O:19])(=[O:20])[OH:21]'

        amol1 = AtomMappedMolecule.fromAtmMappedSmiles(smi1);
        log.info('amol1 : %s, mol:%s' % (amol1, amol1.molecule))

        mol, isNew = Molecule.objects.get_or_create(iso_smiles = smi2)
        self.assertEquals(amol1.molecule, mol)

        cleanSmiles = clean_atom_map_smiles(smi2)
        self.assertEquals(amol1.atm_mapped_smiles, cleanSmiles)
        log.info('cleanSmiles : %s' % cleanSmiles)
        
        amol2, isNew = AtomMappedMolecule.objects.get_or_create(atm_mapped_smiles=cleanSmiles)
        
        self.assertEquals(amol2, amol1)


        amol2 = AtomMappedMolecule.fromAtmMappedSmiles(smi2);
        amol3 = AtomMappedMolecule.fromAtmMappedSmiles(smi3);

        self.assertEquals(amol1.atm_mapped_smiles, amol2.atm_mapped_smiles)
        self.assertEquals(amol1.atm_mapped_smiles, amol3.atm_mapped_smiles)

        log.info("Mols are %s" % pformat([amol1, amol2, amol3]))

        amol1.save()
        amol2.clean();
        amol3.clean();

        self.assertEquals(amol1.molecule, amol2.molecule)
        self.assertEquals(amol1.molecule, amol3.molecule)
        
        self.assertEquals(amol1, amol2)
        self.assertEquals(amol1, amol3)
