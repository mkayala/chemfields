#!/usr/bin/env python
# encoding: utf-8
"""
Allow module level import

File: __init__.py
Author: MKayala
Created 2012-05-12
Copyright 2012. 

Module with utilities for oechem interface.
"""
from openeye import oechem as oe
from MolStdValue import METAL_ATOMIC_NUMS
from Const import SMILES_MOL_DELIM

def molBySmiles(smiles):
    """Convenience method for common construct of creating a
    new molecule by SMILES string.
    """
    mol = oe.OEGraphMol();
    oe.OEParseSmiles(mol, smiles);
    return mol;


def createStandardSmiString( mol ):
    """Strange.  Should be perfectly redundant with OECreateIsoSmiString, but in some cases,
    the direct call messes up, such as alkene stereochemistry giving C\\C=C\\C instead of C/C=C/C.
    Functionally the same, but distinct when doing unique SMILES checks.
    """
    if mol is None:
        return "";
    
    smi = oe.OECreateIsoSmiString(mol);
    return standardizeSmiles(smi);


def standardizeSmiles( smi, ensureAtomMaps=False, kekulize=False, aromatize=False ):
    """Given a SMILES string, generate a standardized isomeric SMILES version.
    Option to use createAtomMapSmiString to ensure any hydrogens with atom mappings
    are retained in the generated SMILES string.
    Also add option to kekulize.
    """
    mol = molBySmiles(smi);
    if kekulize:
        oe.OEClearAromaticFlags(mol);
    
    if aromatize:
        oe.OEAssignAromaticFlags(mol);
    
    stdSmi = None;
    if ensureAtomMaps:
        stdSmi = createAtomMapSmiString(mol);
    else:
        stdSmi = oe.OECreateIsoSmiString(mol);
    return stdSmi;

def createAtomMapSmiString( mol ):
    """Should be the same as createStandardSmiString, but beware of molecules 
    that should generate a SMILES string like [H:1][Br:2].  
    The standard canonization function will ignore
    "implicit hydrogens" and thus yield [BrH:2] as the result, losing the atom
    mapping information on the hydrogen.  Can hack around this, though beware
    that doing so will mess up any hydrogens specified with an isotope of 1
    like [1H]Br.
    """
    originalIsotopes = dict();
    for atom in mol.GetAtoms():
        originalIsotopes[atom.GetIdx()] = atom.GetIsotope();
        if atom.GetAtomicNum() == 1 and atom.GetMapIdx() > 0:   
            # Hydrogen atom with atom mapping.  Artificially set isotope, 
            #   otherwise will be lost as "implicit hydrogen" upon canonization
            atom.SetIsotope(1);

    smi = createStandardSmiString(mol);
    smi = smi.replace("[1H","[H"); # Get rid of the artificial 1H isotopes at this point
    
    # Revert the molecule object back to its original state
    for atom in mol.GetAtoms():
        atom.SetIsotope( originalIsotopes[atom.GetIdx()] );
    
    return smi;

def createStdAtomMapSmiString(mol, semiImplicitH=True):
    """Just a wrapper to ensure proper flags are set to handle hydrogens.
    
    This is a replacement for the createAtomMapSmiString function (that should work better)
    ie, it should NOT alter aromaticity or kekulization.
    """
    nMol = oe.OEGraphMol(mol)
    if semiImplicitH:
        makeNonMappedHydrogensImplicit(nMol)
    flavor = oe.OESMILESFlag_DEFAULT|oe.OESMILESFlag_Hydrogens|oe.OESMILESFlag_AtomStereo|oe.OESMILESFlag_BondStereo
    return oe.OECreateSmiString(nMol, flavor)


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

def removeNonsenseStereo(mol):
    """Method to remove stereo information about any non-chiral atoms or bonds"""
    oe.OEPerceiveChiral(mol)
    for atom in mol.GetAtoms():
        if not atomIsChiral(atom):
            clearAtomStereo(atom);
    for bond in mol.GetBonds():
        if not bondIsChiral(bond):
            clearBondStereo(bond);

def bondIsChiral( bond ):
    """Extension of OpenEye OEBondBase.IsChiral method.
    Manual corrections for
        - Vinylic carbocations should not have stereochemistry, should be linear, sp hybridized
    """
    chiral = bond.IsChiral();
    chiral = chiral and not (bond.GetBgn().GetDegree() < 3 and bond.GetBgn().GetFormalCharge() != 0 );
    chiral = chiral and not (bond.GetEnd().GetDegree() < 3 and bond.GetEnd().GetFormalCharge() != 0 );
    return chiral;


def atomIsChiral( atom ):
    """Extension of OpenEye OEAtomBase.IsChiral method.
    Manual corrections for
        - Organometallics do not have configurationally stable stereochemistry
        - Carbanions should be readily invertible (similar to nitrogen) so achiral
        - onium ions should be considered achiral as well. 
        
    """
    chiral = atom.IsChiral();
    for neighbor in atom.GetAtoms():
        if neighbor.GetAtomicNum() in METAL_ATOMIC_NUMS:
            chiral = False;
            break;  # No point in checking more
    if atom.IsCarbon() and atom.GetFormalCharge() == -1:
        chiral = False;
    if atom.IsOxygen() and atom.GetFormalCharge() == +1:
        chiral = False;
    return chiral;

def clearAtomStereo(atom):
    if atom.HasStereoSpecified():
        neighbors = [];
        for nbr in atom.GetAtoms():
            neighbors.append(nbr);
        atom.SetStereo(neighbors, oe.OEAtomStereo_Tetrahedral, oe.OEAtomStereo_Undefined);

def clearBondStereo(bond):
    if bond.HasStereoSpecified():
        neighbors = [];
        # Find any adjacent neighbor atom, other than the same bonded atom
        for atom in bond.GetBgn().GetAtoms():
            if atom.GetIdx() != bond.GetEnd().GetIdx():
                neighbors.append(atom);
                break;
        for atom in bond.GetEnd().GetAtoms():
            if atom.GetIdx() != bond.GetBgn().GetIdx():
                neighbors.append(atom);
                break;
                
        bond.SetStereo(neighbors, oe.OEBondStereo_CisTrans, oe.OEBondStereo_Undefined);    

def hasAtomMaps(mol):
    """Predicate to determine if there are atom maps in this mol"""
    for atom in mol.GetAtoms():
        if atom.GetMapIdx() != 0:
            return True;
    return False;

def splitCompositeSmilesToList(compositeSmi, retainCounterIons=False):
    """retainCounterIons option will try to not separate out components
    if they have opposing (and neutralizing) formal charges.
    """
    smilesList = compositeSmi.split(SMILES_MOL_DELIM);
    if retainCounterIons:
        # Sort out by net formal charges on each component
        chargeSmiList = [];
        mol = oe.OEGraphMol();
        for smi in smilesList:
            oe.OEParseSmiles(mol, smi);
            chargeSmiList.append( (OENetCharge(mol), smi) );
            mol.Clear();
        chargeSmiList.sort();
        
        # Constants to extract tuple components of the chargeSmiList items
        CHARGE = 0;
        SMILES = 1;
        
        # Now put these back to the smilesList, but combine any charged components
        #   from either end of the sorted list
        smilesList = [];
        iFirst = 0;
        iLast = len(chargeSmiList)-1;
        
        currComponent = [];
        currNetCharge = 0;
        while iFirst <= iLast:
            if len(currComponent) < 1:
                # Haven't added any components yet, just add the first one
                currComponent.append( chargeSmiList[iFirst][SMILES] );
                currNetCharge += chargeSmiList[iFirst][CHARGE];
                iFirst += 1;
            elif currNetCharge < 0:   
                # Negative charge, look for a positive addition at end of the list to neutralize
                if chargeSmiList[iLast][CHARGE] > 0:
                    currComponent.append( chargeSmiList[iLast][SMILES] );
                    currNetCharge += chargeSmiList[iLast][CHARGE];
                    iLast -= 1;
                else:
                    # No positive charged components to add for neutralization, have to just accept the lone ion then
                    smilesList.append( str.join(SMILES_MOL_DELIM, currComponent) );
                    currComponent = [];
                    currNetCharge = 0;
            elif currNetCharge > 0:
                # Positive charge, look for a negative addition at front of the list to neutralize
                if chargeSmiList[iFirst][CHARGE] < 0:
                    currComponent.append( chargeSmiList[iFirst][SMILES] );
                    currNetCharge += chargeSmiList[iFirst][CHARGE];
                    iFirst += 1;
                else:
                    # No negative charged components to add for neutralization, have to just accept the lone ion then
                    smilesList.append( str.join(SMILES_MOL_DELIM, currComponent) );
                    currComponent = [];
                    currNetCharge = 0;
            elif currNetCharge == 0:
                # Found enough components to get a neutral component
                smilesList.append( str.join(SMILES_MOL_DELIM, currComponent) );
                currComponent = [];
                currNetCharge = 0;
            else:
                raise Exception("Sanity check, this shouldn't happen");
        # Whatever's left as a component, have to accept it
        smilesList.append( str.join(SMILES_MOL_DELIM, currComponent) );
    return smilesList;

def joinSmilesListToCompositeSmiles(smiList):
    """Simple convenience to join together a set of smiles"""
    return SMILES_MOL_DELIM.join(smiList);

def clearAtomMaps(mol):
    """Convenience to clear the atom mapping on a whole mol  
    Removes isotope settings on the atom mapped hydrogens as well."""
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 1 and atom.GetMapIdx() > 0:   
            # Hydrogen atom with atom mapping.  Artificially set isotope, to zero 
            atom.SetIsotope(0);
    [atm.SetMapIdx(0) for atm in mol.GetAtoms()];
