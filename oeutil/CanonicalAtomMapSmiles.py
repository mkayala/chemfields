#!/usr/bin/env python
# encoding: utf-8
"""
CanonicalAtomMapSmiles.py

Created by Matt Kayala on 2011-05-15.
Copyright (c) 2011 Institute for Genomics and Bioinformatics. All rights reserved.
"""

import sys
import os

from openeye import oechem as oe

from chemfields.oeutil import molBySmiles, createStdAtomMapSmiString
from chemfields.oeutil import clearAtomMaps, hasAtomMaps
from chemfields.oeutil  import splitCompositeSmilesToList, joinSmilesListToCompositeSmiles

from Util import log

def symmetricAtomsExist(mol):
    """Test for whether or not all the atoms involved are actually symmetric.
    
    If not, then can skip canonicalization
    """
    oe.OEPerceiveSymmetry(mol)
    aMapList = []
    countBySymClass = {};
    
    for atom in mol.GetAtoms():
        if atom.GetMapIdx() > 0:
            aMapList.append(atom)
        symClass = atom.GetSymmetryClass();
        if symClass not in countBySymClass:
            countBySymClass[symClass] = 1;
        else:
            countBySymClass[symClass] += 1;
    
    for atom in aMapList:
        if countBySymClass[atom.GetSymmetryClass()] > 1:
            return True;
    return False;


def makeNonMappedHydrogensImplicit(mol):
    """To ensure common format, we need to make nonMapped Hydrogens to be implicit"""
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 1 and atom.GetMapIdx() == 0:
            bonds = [b for b in atom.GetBonds()]
            if len(bonds) == 1:
                nbr = bonds[0].GetNbr(atom)
                implHCount = nbr.GetImplicitHCount()
                mol.DeleteAtom(atom)
                nbr.SetImplicitHCount(implHCount + 1)

def createStdAtomMapSmiString(mol, semiImplicitH=True):
    """Just a wrapper to ensure proper flags are set to handle hydrogens"""
    nMol = oe.OEGraphMol(mol)
    if semiImplicitH:
        makeNonMappedHydrogensImplicit(nMol)
    flavor = oe.OESMILESFlag_DEFAULT|oe.OESMILESFlag_Hydrogens|oe.OESMILESFlag_AtomStereo|oe.OESMILESFlag_BondStereo
    return oe.OECreateSmiString(nMol, flavor)


def hasAtomMappedHydrogen(mol):
    """Convenience to check if we have atom mapped hydrogens.
    
    Good to know because if not, we can skip a large number of matches
    """
    hasH = False;
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 1 and atom.GetMapIdx() > 0:
            hasH = True;
            break;
    return hasH;


def layBreadCrumbHydrogens(mol):
    """If we have a smiles string with atom mapped hydrogens.
    
    We want to run substructure matching to canonicalize.  However, doing so can 
    cause a recursion limit if ALL hydrogens are explicit
    
    So, here, go through all hydrogens.  Unless bonded to another H, 
    make implicit, but save a breadcrumb dict object with map of adjacent heavy 
    atom with set of maps to recreate
    
    Assumes that adjacent heavy atom has ATOM MAP
    
    If the adjacent heavy atom DOES NOT have an atom map, we need to map it, 
    and save an extra breadCrumb to tell us to delete it later.
    
    """
    nMol = oe.OEGraphMol(mol)
    makeNonMappedHydrogensImplicit(nMol)
    breadCrumb = {};
    toDeleteBreadCrumb = set([]);
    usedMaps = set([0])
    
    # First, if any hydrogen is bonded to another H, or has degree != 1, 
    # add Explicit hydrogens return with no breadcrumbs
    for atom in nMol.GetAtoms():
        usedMaps.add(atom.GetMapIdx());
        if atom.IsHydrogen():
            bonds = [b for b in atom.GetBonds()]
            areNbrsH = any([b.GetNbr(atom).IsHydrogen() for b in bonds])
            if len(bonds) != 1 or areNbrsH:
                oe.OEAddExplicitHydrogens(nMol)
                return (nMol, ({}, {}))
    
    # Then we know can make all H's implicit.
    for atom in nMol.GetAtoms():
        if atom.IsHydrogen():
            # From code above, guaranteed to be only 1!
            nbr = [a for a in atom.GetAtoms()][0]
            
            bCrumbMap = nbr.GetMapIdx()
            if bCrumbMap == 0:
                # Have to add a map to neighbor 
                bCrumbMap = max(usedMaps) + 1;
                usedMaps.add(bCrumbMap)
                nbr.SetMapIdx(bCrumbMap)
                toDeleteBreadCrumb.add(bCrumbMap)
                #raise Exception('Laying breadcrumbs, with non-mapped adjacent atom!')
            if bCrumbMap not in breadCrumb:
                breadCrumb[bCrumbMap] = set([atom.GetMapIdx()])
            else:
                breadCrumb[bCrumbMap].add(atom.GetMapIdx())
            
            nMol.DeleteAtom(atom)
            nbr.SetImplicitHCount(nbr.GetImplicitHCount() + 1);
    return (nMol, (breadCrumb, toDeleteBreadCrumb))


def followBreadCrumbHydrogens(mol, breadCrumbTuple):
    """Function to return a mol back to its previous state with explicit hydrogens"""
    (breadCrumb, toDeleteBreadCrumb) = breadCrumbTuple
    nMol = oe.OEGraphMol(mol)
    for atom in nMol.GetAtoms():
        if atom.GetMapIdx() in breadCrumb:
            for nMapIdx in breadCrumb[atom.GetMapIdx()]:
                nAtom = nMol.NewAtom(1);
                nBond = nMol.NewBond(nAtom, atom, 1)
                nAtom.SetMapIdx(nMapIdx)
                atom.SetImplicitHCount(atom.GetImplicitHCount() - 1)
            if atom.GetMapIdx() in toDeleteBreadCrumb:
                atom.SetMapIdx(0)
    return nMol;


def createPossibleAtomMapSmiles(smi):
    """Use substructure search on self to return all ways to write out atom mapped smiles
    
    Some tricky stuff going on here.  If there are atom mapped hydrogens, then we drop these 
    but save a breadcrumb to be able to get back to them.
    
    The reason for this is that substructure search with explicit hydrogens can be HORRIBLY slow.
    """
    mol = molBySmiles(smi)
    
    #Check for explicit atom mapped Hydrogens:
    hAtomMaps = hasAtomMappedHydrogen(mol);
    breadCrumb = {}
    if hAtomMaps:
        nMol, breadCrumb = layBreadCrumbHydrogens(mol)
        nSmi = createStdAtomMapSmiString(nMol)
        #log.info('mol (with maps) after explicit hydrogens: %s' % nSmi)
    else:
        nSmi = createStdAtomMapSmiString(mol)
    
    #log.info('Canonical Srch: nSmi : %s' % nSmi)
    
    nMol = molBySmiles(nSmi)
    ss = oe.OESubSearch(nSmi)
    #ss.SetMaxMatches(0)
    clearAtomMaps(mol)    
    
    atomByIdxDict = {};
    for atom in mol.GetAtoms():
        atomByIdxDict[atom.GetIdx()] = atom;
    potSmiList = []
    for match in ss.Match(mol):
        for ma in match.GetAtoms():
            mapidx = ma.pattern.GetMapIdx()
            if mapidx > 0:
                atomByIdxDict[ma.target.GetIdx()].SetMapIdx(mapidx)
        potSmiList.append(createStdAtomMapSmiString(mol))
        clearAtomMaps(mol);
    
    # Finally add the hAtoms back in if needed
    if hAtomMaps:
        nPotSmiList = [];
        for smi in potSmiList:
            mol = molBySmiles(smi)
            nMol = followBreadCrumbHydrogens(mol, breadCrumb)
            nPotSmiList.append(createStdAtomMapSmiString(nMol))
        potSmiList = nPotSmiList;
    return potSmiList;

def canonicalizeAtomMapSmiString(atmMapSmi, arom=False):
    """Wrapper to call createCanonicalAtomMapSmiString starting with a smiles string"""
    mol = molBySmiles(atmMapSmi)
    #log.info('inCanon.  Smi: %s' %createStdAtomMapSmiString(mol))
    return createCanonicalAtomMapSmiString(mol, arom);


def createCanonicalAtomMapSmiString(mol, arom=False):
    """TRULY canonicalize an atom mapped mol
    
    First, check for symmetry.  If none, then simply write out
    If there is some, then split up the mol by components.
    
    For each atmMapped component, enumerate ALL possible matches 
    using substructure search on self.  
    
    Then recombine the components.
    """
    copyMol = oe.OEGraphMol(mol)
    if arom:
        oe.OEAssignAromaticFlags(copyMol)
    
    sym = symmetricAtomsExist(copyMol)
    if not sym:
        return createStdAtomMapSmiString(copyMol)
        
    smi = createStdAtomMapSmiString(copyMol)
    smiList = splitCompositeSmilesToList(smi)
    
    newSmiList = []
    for aSmi in smiList:
        aMol = molBySmiles(aSmi)
        newSmiList.append(createCanonicalAtomMapSmiString_singleComponent(aMol))
    # Finally, take the smiles from all and combine, relative ordering will have been set from 
    # original write out!
    ## Canonicalize by sorting here!
    newSmiList.sort()
    nSmi = joinSmilesListToCompositeSmiles(newSmiList)
    #nSmi = createStdAtomMapSmiString(molBySmiles(nSmi))
    return nSmi
    


def createCanonicalAtomMapSmiString_singleComponent(mol):
    """return the canon atm map string for a single component"""
    sym = symmetricAtomsExist(mol)
    if not sym:
        return createStdAtomMapSmiString(mol)
    isAtomMapped = hasAtomMaps(mol)
    if not isAtomMapped:
        return createStdAtomMapSmiString(mol)
    
    ## Otherwise, run through
    possibleSmiles = createPossibleAtomMapSmiles(createStdAtomMapSmiString(mol))
    possibleSmiles.sort();
    if len(possibleSmiles) > 0:
        return possibleSmiles[0]
    
    log.warning('Possible SMILES = NONE, for %s' % createStdAtomMapSmiString(mol))
    return ''


def _standardAtomMapSmi(smi):
    """Helper function to read into a mol, write out with atom maps"""
    flavor = oe.OESMILESFlag_DEFAULT|oe.OESMILESFlag_AtomMaps
    mol = molBySmiles(smi);
    return oe.OECreateSmiString(mol, flavor)

def smiToSmirks(smi):
    return "%s>>%s" % (smi, smi)

def selfMatches(atmMapSmi):
    """Make a SMIRKS string of self.  Return all atm mapped matches"""
    libgen = oe.OELibraryGen()
    libgen.Init(smiToSmirks(atmMapSmi))
    mol = molBySmiles(atmMapSmi)
    matches = libgen.SetStartingMaterial(mol, 0, False)
    
    prodSmiList = []
    for prod in libgen.GetProducts():
        hasAtomMap = False;
        for atom in prod.GetAtoms():
            if atom.GetMapIdx() > 0:
                hasAtomMap = True;
                break;
        if hasAtomMap:
            oe.OEClearAromaticFlags(prod)
            prodSmiList.append(oe.OECreateIsoSmiString(prod))
    return prodSmiList



def enumeratePotentialSmilesStringsForMol( mol):
    """Given a mol with atom maps, enumerate all the different types we get 
    by writing as smiles, reading back in.
    
    If we call standardize Smiles multiple times this will sometimes 
    return different smiles. Especially with highly symmetric mols
    and only one or two mapped atoms.  
    
    However, this renaming is circular, and eventually comes back on itself
    so, to canonicalize, enumerate all of these and then string sort, take the 
    first.
    """
    flavor = oe.OESMILESFlag_DEFAULT|oe.OESMILESFlag_AtomMaps;
    potSmiList = [oe.OECreateSmiString(mol, flavor)]
    #log.info('First potSmi: %s' % potSmiList[0])
    while True:
        nMol = molBySmiles(potSmiList[-1])
        nSmi = oe.OECreateSmiString(nMol, flavor)
        if nSmi not in potSmiList:
            #log.info('nSmi : %s' % nSmi)
            potSmiList.append(nSmi)
        else:
            break;
    return potSmiList
