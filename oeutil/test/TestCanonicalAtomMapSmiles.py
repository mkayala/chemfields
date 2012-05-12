#!/usr/bin/env python
# encoding: utf-8
"""
TestCanonicalAtomMapSmiles.py

Created by Matt Kayala on 2011-05-10.
Copyright (c) 2011 Institute for Genomics and Bioinformatics. All rights reserved.

A Test Case in the CHEM Module.
"""
import sys, os;
import unittest;
import cStringIO;
import tempfile;
from pprint import pformat

from chemfields.oeutil import molBySmiles, createAtomMapSmiString, createStdAtomMapSmiString

from chemfields.oeutil.CanonicalAtomMapSmiles import createCanonicalAtomMapSmiString, createPossibleAtomMapSmiles
from chemfields.oeutil.CanonicalAtomMapSmiles import canonicalizeAtomMapSmiString
from chemfields.oeutil.CanonicalAtomMapSmiles import layBreadCrumbHydrogens, followBreadCrumbHydrogens
from openeye.oechem import OEGraphMol, OEParseSmiles, OECreateIsoSmiString
from openeye import oechem as oe;


import Const, Util;
from Util import log;

class TestCanonicalAtomMapSmiles(unittest.TestCase):
    def setUp(self):
        """Set up anything for the tests.., """
        super(TestCanonicalAtomMapSmiles, self).setUp();
                
    def tearDown(self):
        """Restore state"""
        super(TestCanonicalAtomMapSmiles, self).tearDown();
    
    
    def test_substructMatchesReached(self):
        """Test of the canonicalizeAtomMapSmiString stuff where we are going too deep recursively"""
        smiList = \
            [
                '[Al-](Cl)(Cl)(Cl)Cl.[Al-](Cl)(Cl)(Cl)Cl.[CH3:11]C(=O)C1[CH+]C=C(O1)C(F)(F)F.[CH3:21]C(=O)C1[CH+]C=C(O1)C(F)(F)F',
                '[Al-](Cl)(Cl)(Cl)Cl.[Al-](Cl)(Cl)(Cl)Cl.[CH3:11]C(=O)C1[CH+]C=C(O1)C(F)(F)F.[CH3:21][C:22](=[O:23])C1[CH+]C=C(O1)C(F)(F)F',
                '[Al-](Cl)(Cl)(Cl)Cl.[Al-](Cl)(Cl)(Cl)Cl.[CH3:11]C(=O)C1[CH+]C=C(O1)C(F)(F)F.[CH3:21][C:22](=O)[CH:23]1[CH+]C=C(O1)C(F)(F)F',
            ]
        for smi in smiList:
            #log.info('About to look at %s' % smi)
            expSmi = smi
            self.assertEquals(canonicalizeAtomMapSmiString(smi), expSmi)
    
    
    
    def test_PossibleSmilesEqNone(self):
        """Testing that the warning of possibleSMILES == None does not come up."""
        smiList = \
            [
                '[H-].[H-].[H:11][CH2:10]C(=O)Cl.[H:20][CH2:21]C(=O)Cl',
                '[H-].[H-].[H:11][CH2:10]C(=O)Cl.[H:20][CH2:21][C:22](=[O:23])Cl',
                '[H-].[H-].[H:11][CH2:10]C(=O)Cl.[H:20][CH2:21][C:22](=O)[Cl:23]',
                '[H-].[H-].[H:11][CH2:10]C(=O)Cl.[H:21][CH2:20]C(=O)Cl',
            ]
        
        for smi in smiList:
            #log.info('About to look at %s' % smi)
            expSmi = smi;
            self.assertEquals(canonicalizeAtomMapSmiString(smi), expSmi)
        
    
    def test_hydrogenAtomMaps(self):
        """Test with a few smiles where we have hydrogen atom maps"""
        smi = '[H-:10]'
        expSmi = smi;
        self.assertEquals(canonicalizeAtomMapSmiString(smi), expSmi)
        
        smi = '[H-:10].[H:20][C:21]1(C(=O)CC[C:22]1=[O:23])C'
        expSmi = smi 
        self.assertEquals(canonicalizeAtomMapSmiString(smi), expSmi)
        
        
    
    def _testCanonicalSmi(self, smi):
        """Convenience to help with testing"""
        mol1 = molBySmiles(smi)
        smi2 = createCanonicalAtomMapSmiString(mol1)
        mol2 = molBySmiles(smi2)
        smi3 = createCanonicalAtomMapSmiString(mol2)
        mol3 = molBySmiles(smi3)
        smi4 = createCanonicalAtomMapSmiString(mol3)
        
        #log.info('smi2, smi3, smi4 : %s' % str([smi2, smi3, smi4]))
        # smi2, smi3, smi4 should ALL be the same!!!
        self.assertEquals(smi2, smi3)
        self.assertEquals(smi3, smi4)
    
    
    def test_createCanonicalAtomMapSmiString(self):
        """Test that new code to do this mapping works"""
        smi1 = 'C1=C[CH:1]=CC=C1'
        self._testCanonicalSmi(smi1)
        smi2 = 'C1=[CH:1]C=CC=C1'
        self._testCanonicalSmi(smi2)
        
        smi3 = '[H][C]1=[CH:1]C=CC=C1'
        self._testCanonicalSmi(smi2)
        
        self._test_multipleSame([smi1, smi2, smi3])
        
        aroSmi1 = '[cH:1]1ccccc1'
        self._testCanonicalSmi(aroSmi1)
        aroSmi2 = 'c1ccc[cH:1]c1'
        self._testCanonicalSmi(aroSmi2)
        
        self._test_multipleSame([aroSmi1, aroSmi2])
        
        smi1 = 'C1=C[CH:10]=[CH:20]C=C1'
        self._testCanonicalSmi(smi1)
        smi2 = '[CH:20]1=[CH:10]C=CC=C1'
        self._testCanonicalSmi(smi2)
        
        self._test_multipleSame([smi1, smi2])
        
        ## And then a bunch of others
        smi1 = 'ClC1=C(Cl)[C:1](Cl)=C(Cl)C(Cl)=C1Cl'
        self._testCanonicalSmi(smi1)
        smi2 = 'ClC1=[C:1](Cl)C(Cl)=C(Cl)C(Cl)=C1Cl'
        self._testCanonicalSmi(smi2)
        
        self._test_multipleSame([smi1, smi2])
        
        smi1 = 'C1=CC2=[C:1](C=C1)C=CC=C2'
        self._testCanonicalSmi(smi1)
        smi1 = 'C1=CC2=C(C=C1)[CH:1]=CC=C2'
        self._testCanonicalSmi(smi1)
        smi1 = 'C1=CC2=C(C=C1)C=[CH:1]C=C2'
        self._testCanonicalSmi(smi1)
        smi1 = 'C1CC[CH2:1]CC1'
        self._testCanonicalSmi(smi1)
        
        
        smi1 = 'ClC1=C[CH:1]=CC=C1'
        self._testCanonicalSmi(smi1)
        smi1 = 'ClC1=C[CH:1]=C(C=C1)C1=CC=C(Cl)C=C1'
        self._testCanonicalSmi(smi1)
        smi1 = 'ClC1=C[C:1](Cl)=CC(Cl)=C1'
        self._testCanonicalSmi(smi1)
        smi1 = 'ClC1=CC(Cl)=[CH:1]C(Cl)=C1'
        self._testCanonicalSmi(smi1)
        smi1 = 'ClC1=[CH:1]C2=CC(=C1)C1=CC(=CC(Cl)=C1)C1=CC2=CC(Cl)=C1'
        self._testCanonicalSmi(smi1)
    
    
    def test_multipleSame(self):
        """Test that all of these atom mapped smiles are the same"""
        smiList = [
                'C1CCCCC1[Br:1]',
                'C1CC([Br:1])CCC1',
            ]
        self._test_multipleSame(smiList)
        
        smiList = [
                'ClC1=C(Cl)[C:1]([Cl:2])=C(Cl)C(Cl)=C1Cl',
                'ClC1=[C:1]([Cl:2])C(Cl)=C(Cl)C(Cl)=C1Cl'
            ]
        self._test_multipleSame(smiList)
        
        smiList = [
            '[H:10][CH:11](C(=O)OCC)[C:12](=[O:13])OCC',
            '[H:10][CH:11]([C:12](=[O:13])OCC)C(=O)OCC'
        ]
        self._test_multipleSame(smiList)
        
        smiList = [
            'CC(C)(C)[C:13]1=[CH:12][CH:11]=[CH:10]C=C1',
            'CC(C)(C)[C:13]1=[CH:12][CH:11]=[CH:10]C=C1',
        ]
        self._test_multipleSame(smiList)
        
        smiList = [
            '[H:10][C:11]1(C=C[C+:14]([CH:13]=[CH:12]1)C(C)(C)C)C(C)C',
            '[H:10][C:11]1(C=C[C+:14]([CH:13]=[CH:12]1)C(C)(C)C)C(C)C',
        ]
        self._test_multipleSame(smiList)
        
        smiList = [
            'CC1=C[CH:1]=CC=C1C',
            'CC1=CC=[CH:1]C=C1C'
        ]
        self._test_multipleSame(smiList)
    
    
    def test_breadCrumbCode(self):
        """Test that code to remove hydrogens, but retain a breadcrumb works"""
        smi = '[H:10][C:11]1(C=C[C+:14]([CH:13]=[CH:12]1)C(C)(C)C)C(C)C'
        expSmi = 'CC(C)[CH:11]1C=C[C+:14]([CH:13]=[CH:12]1)C(C)(C)C'
        expBreadCrumbTuple = ({11:set([10])}, set([]))
        
        mol = molBySmiles(smi)
        origSmi = createStdAtomMapSmiString(mol)
        nMol, breadCrumbTuple = layBreadCrumbHydrogens(mol)
        endSmi = createStdAtomMapSmiString(nMol)
        
        self.assertEquals(endSmi, expSmi);
        self.assertEquals(breadCrumbTuple[0], expBreadCrumbTuple[0])
        self.assertEquals(breadCrumbTuple[1], expBreadCrumbTuple[1])
        
        
        backMol = followBreadCrumbHydrogens(nMol, breadCrumbTuple)
        finalSmi = createStdAtomMapSmiString(backMol)
        
        self.assertEquals(finalSmi, origSmi)
    
    
    def test_breadCrumbOnlyHMap(self):
        """Test that code to handle bread crumb laying/following with hydrogens works"""
        smi = '[H:1]CC(=O)C1=CC=C2C=C[CH+]C(C2=C1)C(C)C'
        expSmi = 'CC(C)C1[CH+]C=CC2=CC=C(C=C12)C(=O)[CH3:2]'
        expBreadCrumbTuple = ({2:set([1])}, set([2]))
        
        mol = molBySmiles(smi)
        origSmi = createStdAtomMapSmiString(mol)
        nMol, breadCrumbTuple = layBreadCrumbHydrogens(mol)
        endSmi = createStdAtomMapSmiString(nMol)
        
        self.assertEquals(endSmi, expSmi);
        self.assertEquals(breadCrumbTuple[0], expBreadCrumbTuple[0])
        self.assertEquals(breadCrumbTuple[1], expBreadCrumbTuple[1])
        
        
        backMol = followBreadCrumbHydrogens(nMol, breadCrumbTuple)
        finalSmi = createStdAtomMapSmiString(backMol)
        #log.info('FinalSmi: %s' % finalSmi)
        self.assertEquals(finalSmi, origSmi)
    
    
    def _test_multipleSame(self, smiList):
        """Convenience to test that a number are all the same"""
        nSmiList = [canonicalizeAtomMapSmiString(smi) for smi in smiList]
        
        for aSmi in nSmiList[1:]:
            #log.info('Getting into test!')
            self.assertEquals(aSmi, nSmiList[0])
    
    
    def test_retainsAtomMaps(self):
        """Test that the code retains the correct atom mapping"""
        # A single test here.
        smi = '[H:10][CH:11](C(=O)OCC)[C:12](=[O:13])OCC'
        mol = molBySmiles(smi)
        numAtomMaps = 0;
        for a in mol.GetAtoms():
            if a.GetMapIdx() > 0:
                numAtomMaps += 1;
        
        #possibleSmiles = createPossibleAtomMapSmiles(smi)
        #log.info('possibleSmiles: %s' % pformat(possibleSmiles))
        
        nSmi = canonicalizeAtomMapSmiString(smi)
        nMol = molBySmiles(nSmi)
        nNumAtomMaps = 0;
        for a in nMol.GetAtoms():
            if a.GetMapIdx() > 0:
                nNumAtomMaps += 1;
        
        #log.info('smi: %s\nnSmi:%s' % (smi, nSmi))
        self.assertEquals(numAtomMaps, nNumAtomMaps)
    
    

    

def suite():
    """Returns the suite of tests to run for this test class / module.
    Use unittest.makeSuite methods which simply extracts all of the
    methods for the given class whose name starts with "test"

    Actually, since this is mostly all basic input / output functions,
    can do most of it with doctests and DocTestSuite
    """
    suite = unittest.TestSuite();
    suite.addTest(unittest.makeSuite(TestCanonicalAtomMapSmiles));
    return suite;

if __name__=="__main__":
    unittest.TextTestRunner(verbosity=Const.RUNNER_VERBOSITY).run(suite())
