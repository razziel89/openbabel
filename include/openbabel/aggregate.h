/**********************************************************************
aggregate.h - Handle aggregates. Declarations of OBAggregate.

Copyright (C) 2016 by Torsten Sachse

This file is part of the Open Babel project modified by Torsten Sachse.
For more information about Open Babel, see <http://openbabel.org/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/

#ifndef OB_AGGREGATE_H
#define OB_AGGREGATE_H

#include <openbabel/atom.h>
#include <openbabel/mol.h>
#include <vector>
#include <algorithm>

namespace OpenBabel
{
   bool IdxSort(OBAtom* i, OBAtom* j){
       return i->GetIdx() < j->GetIdx();
   }

  //class OBAggregate;

  // In the following, the word "molecule" means "a covalently bound entity".
  // class introduction in aggregate.cpp
  class OBAPI OBAggregate: public OBMol
  {
  protected:
    int                                  _nrMolecules; //!< number of molecules in this aggregate
    std::vector<std::vector<OBAtom*> >   _connections; //!< vector of molecules (each a vector of atoms)
    bool                                 _needRefresh; //!< whether or not _connections must be refreshed
    //! Tags can be used to manipulate multiple molecules in an
    //! aggregate at the same time after adding them to one tag
    std::vector<std::vector<OBAtom*> >   _tags;    //!< all tags in this aggregate
    std::vector<std::vector<int> >       _molTags; //!< which molecule belongs to which tag
    bool                                 _useTag;  //!< whether or not to use tags to determine which molecules to move
    //all the other private members are inherited from OBMol

  public:

    //! \name Initialization and data (re)size methods
    //@{
    //! Constructor
    OBAggregate();
    //! Copy constructor, copies atoms,bonds and OBGenericData
    OBAggregate(const OBAggregate &agg); // only needed because swig does not seem to know implicit casts
    OBAggregate(const OBMol &mol);
    //! Destructor
    ~OBAggregate();

    //! Assignment, copies atoms and bonds but no OBGenericData
    //! If given an OBAggregate, inplicit casting is used
    OBAggregate &operator=(const OBMol &mol);
    OBAggregate &operator=(const OBAggregate &agg);
    void Assign(const OBMol &mol){ *this = mol; }
    void Assign(const OBAggregate &agg){ *this = agg; }

	//! Append one molecule to this aggregate
	void AppendMolecule(const OBMol &source);
	//! Append one molecule to this aggregate and translate and rotate it prior to appending it without
    //! changing the original
	void AppendMolecule(const OBMol &source, const double v[3], const double axis[3], const double angle);
    //@}
    
    //! \name Aggregate analysis methods
    //@{
	//! A very simple steric filter without the need to use the steric filter class.
    //! This filter ignores steric clashes between atoms within each molecule.
    //! vdw_factor: an atoms vdw-radius is multiplied by this factor to get the radius witin which clashes are considered to happen
    //! vdw_added: add this much to each vdw-radius if it is >0.0
	bool IsGoodVDW(double vdw_factor=0.9, double vdw_added=-1.0);
    //! Get the minimum distance between the vdW spheres of all atoms in all
    //! non-covalently bound entities in the aggregate.
    //! This means that the result is <0 if there are vdW clashes and >0 otherwise.
    //! If break_on_clash==true, the function will return as soon as a vdW clash has been detected
    //! thus not returning the true minimum vdW distance.
    //! If vdw_factor>0.0, it will be used to scale all vdw_radii.
	float MinVDWDist(bool break_on_clash=false, double vdw_factor=-1.0);
	//! Find all atoms that are connected to the given atom
	OBBitVec FindConnectedChildren(std::vector<OBAtom*> &children,OBAtom *check);
	//! Find all connections within the aggregate and optimally sort the atoms by their Idx
	void FindAllConnections(bool sort_it=true);
    //@}

    //! \name Aggregate information printing methods
    //@{
    //! Print out a simple representation of which atoms belong to which molecule
	void PrintConnections();
    //! Get a simple representation of which atoms belong to which molecule.
    //!
    //! The string can easily be parsed. Newlines separate molecules and spaces the
    //! atoms in each molecule. The first number per line is the molecule count.
    std::string GetConnections();
    //! Print out a simple representation of which molecules belong to which tag
    //! Returns true on success and false otherwise
    bool PrintTags();
    //! Get a simple representation of which molecules belong to which tag.
    //!
    //! The string can easily be parsed. Newlines separate tags and spaces
    //! the molecules per tag. The special tag "-1" denotes untagged molecules.
    //! Returns true on success and false otherwise
    std::string GetTags();
    //@}
    
    //! \name Aggregate information getting methods
    //@{
    //!Return the number of molecules in the aggregate
    int GetNrMolecules();
    //!Return the number of tags. -1 is returned if aggregate has not been configures to use tags.
    int GetNrTags();
    //!Return the value of _useTag
    bool GetUseTag();
    //! Get the center of a part (i.e. molecule or tag) of the aggregate.
    //!
    //!Part numbers start at 0.
    vector3 GetCenterPart(int part);
    //@}

    //! \name Aggregate tagging methods
    //@{
    //! Set the member _useTag to true. If _useTag==true, all member function ending in Part will
    //! manipulate the respective tag and not the molecule.
    void EnableTags();
    //! Set the member _useTag to false. Will clear all tag data
    void DisableTags();
    //! Create a new tag. Returns the index of the new tag (can be used to add a molecule to the tag)
    int CreateTag(int nr_elements=0);
    //! Add a molecule to a tag. Return true on success and false on error
    bool AddToTag(int molnr, int tag);
    //@}

    //! \name Manipulation of parts of aggregate methods
    //@{
    //! A part can either be specified by the number of the molecule in the aggregate or the number
    //! of the tag if the boolean value "_useTag" is set to true by calling EnableTags().
    //! Rotate a molecule using only doubles
    void RotatePart(const int part, const double axis[3], const double angle);
    //! Rotate a molecule around one of its main axes
    void RotatePart(const int part, const int main_axis_nr, const double angle);
    //! Translate a molecule using only doubles
    void TranslatePart(const int part, const double vector[3]);
    void TranslatePart(const int part, const vector3 vec);
    //! Align a molecule with a plane as good as possible.
    //! The plane is given by the point p which will be the molecule's new centre
    //! and two axes v1 and v2 which are the third main axis and the second main axis respectively.
    void AlignPart(const int part, const double p[3], const double v1[3], const double v2[3]);
    void AlignPart(const int part, const vector3 p, const vector3 v1, const vector3 v2);
    //! Mirrors the molecule at a point (inversion, if normal=(0,0,0)) or plane in Hessian normal form
    void MirrorPart(const int part, const double normal[3], const double point[3], bool center_it=false);
    void MirrorPart(const int part, const vector3 normal, const vector3 point, bool center_it=false);
    //! Move two parts of a molecule closer together in the direction given by vector.
    //! This will stop once a vdw-clash has been detected or when the distance between
    //! the centers of the parts increases when following the vector.
    bool MovePartsCloser(const double vector[3], const int part1, const int part2,
            const double stepsize=0.1, const double vdw_factor=0.9, double vdw_added=-1.0);
    //! Move two parts of a molecule closer together in the direction given by the centers of both parts.
    //! This will stop once a vdw-clash has been detected or when the distance between
    //! the centers of the parts increases when following the vector.
    bool MovePartsCloser(const int part1, const int part2, 
            const double stepsize=0.1, const double vdw_factor=0.9, double vdw_added=-1.0);
    //@}
    
  };

} // end namespace OpenBabel

#endif // OB_AGGREGATE_H
