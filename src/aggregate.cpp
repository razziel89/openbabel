/**********************************************************************
aggregate.cpp - Handle aggregates.

Copyright (C) 2016 by Torsten Sachse

This file is part of the Open Babel project modified by Torsten Sachse.
For more information on Open Babel, see <http://openbabel.org/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/
#include <openbabel/babelconfig.h>

#include <openbabel/mol.h>
#include <openbabel/aggregate.h>
#include <openbabel/rotamer.h>
#include <openbabel/phmodel.h>
#include <openbabel/bondtyper.h>
#include <openbabel/builder.h>
#include <openbabel/math/matrix3x3.h>

#include <openbabel/stereo/tetrahedral.h>
#include <openbabel/stereo/cistrans.h>

#include <sstream>
#include <set>

#include <algorithm>
#include <limits>
#include <cmath>

using namespace std;

namespace OpenBabel
{
    
  /** \class OBAggregate aggregate.h <openbabel/aggregate.h>
      \brief Aggregate Class

      An aggregate is an object consisting of several molecules (i.e.,
      covalently bound entities) and, hence, is based on OBMol. The OBAggregate
      class is designed to store the additional information associated with an
      aggregate, to enable manipulations of single molecules in the aggregate,
      the creation of such aggregates, and determine which atoms belong to the
      same molecule.

  */

  //
  // OBAggregate member functions
  //

  //! \brief Append a new molecule.
  //!
  //! Append all atoms in a molecule to this aggregate as a new molecule.
  void OBAggregate::AppendMolecule(const OBMol &source)
  {
#define nullvector {0,0,0}
#define standardaxis {1,0,0}
    const double n[3] = nullvector;
    const double a[3] = standardaxis;
#undef nullvector
#undef standardaxis
    AppendMolecule(source,n,a,0);
  }

  //! \brief Append a new molecule and move it.
  //!
  //! Append all atoms in a molecule to this aggregate as a new molecule.
  //! The to-be-appended molecule can be moved (translation and rotation) before being appended.
  //! The original data is left unchanged.
  void OBAggregate::AppendMolecule(const OBMol &source, const double v[3], const double axis[3], const double angle)
  {
    OBMol src;
    if (this == &source){
      src = OBMol(source);
    }
    else{
      src = (OBMol &)source;
    }
    vector<OBAtom*>::iterator i;
    vector<OBBond*>::iterator j;
    OBAtom *atom, *a;
    OBBond *bond;
    unsigned int at_count=1;
    const unsigned int nr_atoms=NumAtoms();
    const unsigned int nr_bonds=NumBonds();
    unsigned int begin_atom_idx, end_atom_idx, order, flags;

    src.Rotate(axis,angle);
    src.Translate(v);

    BeginModify();
    _vatom.reserve(nr_atoms+src.NumAtoms());
    _atomIds.reserve(nr_atoms+src.NumAtoms());
    _vbond.reserve(nr_bonds+src.NumBonds());
    _bondIds.reserve(nr_bonds+src.NumBonds());

    for (atom = src.BeginAtom(i);atom;atom = src.NextAtom(i), ++at_count)
    {
      a = CreateAtom();
      a->Duplicate(atom);
      AddAtom(*a,true);
    }
  
    for (bond = src.BeginBond(j);bond;bond = src.NextBond(j))
    {
      begin_atom_idx = bond->GetBeginAtomIdx() + nr_atoms;
      end_atom_idx = bond->GetEndAtomIdx() + nr_atoms;
      order = bond->GetBondOrder();
      flags = bond->GetFlags();
      AddBond(begin_atom_idx, end_atom_idx, order, flags, -1);
    }
    EndModify();
    _needRefresh = true;
  }

  //! Determine the minimum van-der-Waals distance between all atoms in different molecules
  float OBAggregate::MinVDWDist(bool break_on_clash, double vdw_factor)
  {
    //floats are used for faster calculation of a square root since it needs fewer iterations
    float result = std::numeric_limits<float>::infinity(); 
    vector3 v;
    float sum_vdw, dist;
    vector<OBAtom*>::iterator i, j;
    OBAtom *a1, *a2;
    int first, second;
    vector<OBAtom*> first_molecule, second_molecule;
    int length;
    int at_nr_a1;
    vector3 at_vec_a1;

    vector<float> vdw_radii_temp;
    vdw_radii_temp.reserve(etab.GetNumberOfElements()+1);
    //to be able to directly use an atoms Idx (which start at 1), the
    //element with number 0 is filled with bogus data
    vdw_radii_temp.push_back(0.0f);
    if (vdw_factor<=0.0){
      for (int at_nr=1; at_nr<=etab.GetNumberOfElements();++at_nr){
        vdw_radii_temp.push_back((float)etab.GetVdwRad(at_nr));
      }
    }
    else{
      for (int at_nr=1; at_nr<=etab.GetNumberOfElements();++at_nr){
        vdw_radii_temp.push_back((float)vdw_factor * (float)etab.GetVdwRad(at_nr));
      }
    }

    if (_needRefresh){
      FindAllConnections();
    }
    length = _connections.size();

    if (break_on_clash){
      bool loop = true;
      for (first = 0; loop && first < length-1; ++first)
      {
        for (second = first+1; loop && second < length; ++second)
        {
          first_molecule  = *(_connections[first]);
          second_molecule = *(_connections[second]);
          for (i = first_molecule.begin(); loop && i != first_molecule.end();++i)
          {
            a1=*i;
            at_nr_a1 = a1->GetAtomicNum();
            at_vec_a1 = a1->GetVector();
            for (j = second_molecule.begin(); loop && j != second_molecule.end();++j)
            {
              a2=*j;
              v = at_vec_a1 - a2->GetVector();
              sum_vdw = vdw_radii_temp[at_nr_a1] + vdw_radii_temp[a2->GetAtomicNum()];
              dist = sqrtf((float)v.length_2()) - sum_vdw;
              if ( dist < result ){ 
                result=dist;
                if (result<0.0){loop = false;}
              }
            }
          }
        }
      }
    }
    else{
      for (first = 0; first < length-1; ++first)
      {
        for (second = first+1; second < length; ++second)
        {
          first_molecule  = *(_connections[first]);
          second_molecule = *(_connections[second]);
          for (i = first_molecule.begin();i != first_molecule.end();++i)
          {
            a1=*i;
            at_nr_a1 = a1->GetAtomicNum();
            at_vec_a1 = a1->GetVector();
            for (j = second_molecule.begin();j != second_molecule.end();++j)
            {
              a2=*j;
              v = at_vec_a1 - a2->GetVector();
              sum_vdw = vdw_radii_temp[at_nr_a1] + vdw_radii_temp[a2->GetAtomicNum()];
              dist = sqrtf((float)v.length_2()) - sum_vdw;
              if ( dist < result ){ 
                result=dist;
              }
            }
          }
        }
      }
    }
    return (double)result;
  }

  //! \brief Determine whether or not van-der-Waals surfaces of molecules overlap (i.e., clash).
  //!
  //! All van-der-Waals radii can be scaled by vdw_factor and vdw_added will be added
  //! to the result before determining clashes.
  bool OBAggregate::IsGoodVDW(double vdw_factor, double vdw_added)
  {
    vector3 v;
    bool result = true;
    bool loop = true;
    double sum_vdw;
    vector<OBAtom*>::iterator i, j;
    OBAtom *a1, *a2;
    int first, second;
    vector<OBAtom*> first_molecule, second_molecule;
    int length;
    vdw_factor *= vdw_factor;
    int at_nr_a1;
    vector3 at_vec_a1;

    vector<double> vdw_radii_temp;
    vdw_radii_temp.reserve(etab.GetNumberOfElements()+1);
    vdw_radii_temp.push_back(0.0);
    for (int at_nr=1; at_nr<=etab.GetNumberOfElements();++at_nr){
      vdw_radii_temp.push_back(etab.GetVdwRad(at_nr));
    }

    if (_needRefresh){
      FindAllConnections();
    }
    length = _connections.size();
    
    if (vdw_added > 0.0){
      vdw_added += vdw_added;
      for (first = 0; loop && first < length-1; ++first)
      {
        for (second = first+1; loop && second < length; ++second)
        {
          first_molecule  = *(_connections[first]);
          second_molecule = *(_connections[second]);
          for (i = first_molecule.begin(); loop && i != first_molecule.end();++i)
          {
            a1=*i;
            at_nr_a1 = a1->GetAtomicNum();
            at_vec_a1 = a1->GetVector();
            for (j = second_molecule.begin(); loop && j != second_molecule.end();++j)
            {
              a2=*j;
              v = at_vec_a1 - a2->GetVector();
              sum_vdw = vdw_radii_temp[at_nr_a1] + vdw_radii_temp[a2->GetAtomicNum()] + vdw_added;
              sum_vdw *= sum_vdw;
              if ( v.length_2() < vdw_factor*sum_vdw){
                result=false;
                loop = false;
              }
            }
          }
        }
      }
    }
    else{
      for (first = 0; loop && first < length-1; ++first)
      {
        for (second = first+1; loop && second < length; ++second)
        {
          first_molecule  = *(_connections[first]);
          second_molecule = *(_connections[second]);
          for (i = first_molecule.begin(); loop && i != first_molecule.end();++i)
          {
            a1=*i;
            at_nr_a1 = a1->GetAtomicNum();
            at_vec_a1 = a1->GetVector();
            for (j = second_molecule.begin(); loop && j != second_molecule.end();++j)
            {
              a2=*j;
              v = a1->GetVector() - a2->GetVector();
              sum_vdw = vdw_radii_temp[at_nr_a1] + vdw_radii_temp[a2->GetAtomicNum()];
              sum_vdw *= sum_vdw;
              if ( v.length_2() < vdw_factor*sum_vdw){
                result=false;
                loop = false;
              }
            }
          }
        }
      }
    }
    return result;
  }

  //! Find all molecules in an aggregate.
  void OBAggregate::FindAllConnections(bool sort_it)
  {
    _needRefresh = false;
    OBAtom *atom, *a;
    int nr_at = NumAtoms();
    int atoms_checked = 0;
    vector<OBAtom*>::iterator i;
    OBBitVec bitvec, tempvec;
    bitvec.Resize((unsigned)nr_at+1);
    tempvec.Resize((unsigned)nr_at+1);

    bitvec.Clear();
    bitvec.Negate();
    bitvec.SetBitOff(0);
    for (vector<vector<OBAtom*>* >::iterator it = _connections.begin(); it != _connections.end(); ++it){
      (*it)->clear();
    }
    _connections.clear();
    for (vector<vector<OBAtom*>* >::iterator it = _tags.begin(); it != _tags.end(); ++it){
      (*it)->clear();
    }
    _tags.clear();

    for (;;)
    {
      atom=GetAtom(bitvec.NextBit(-1));
      
      vector<OBAtom*> *tempchildren = new vector<OBAtom*>();
      tempvec = FindConnectedChildren(*tempchildren,atom);
      tempchildren->push_back(atom);
      tempvec.Negate();
      bitvec&= tempvec;
      atoms_checked +=tempchildren->size();
      if (sort_it){
        sort(tempchildren->begin(), tempchildren->end(), IdxSort);
      }
      _connections.push_back(tempchildren);

      if (nr_at <= atoms_checked){break;}; 
      //find first atom in aggregate that does not belong to any of the children found so far
    }
    _nrMolecules = _connections.size();
    if (_useTag){
      vector<vector<int> >::iterator outer_it;
      vector<int>::iterator inner_it;
      //refresh the data in the vector
      _tags.reserve(_molTags.size());
      for (outer_it = _molTags.begin(); outer_it != _molTags.end(); ++outer_it){
        if ( outer_it->size() > 0 ){
          vector<OBAtom*> *adjustvec = new vector<OBAtom*>();
          for (inner_it = outer_it->begin(); inner_it != outer_it->end(); ++inner_it){
            adjustvec->reserve(adjustvec->size()+(_connections[*inner_it])->size());
            adjustvec->insert(adjustvec->end(), (_connections[*inner_it])->begin(), (_connections[*inner_it])->end());
          }
          _tags.push_back(adjustvec);
        }
      }
    }
  }

  //! Find all atoms that can be reached from a given one by traversing bonds.
  OBBitVec OBAggregate::FindConnectedChildren(vector<OBAtom*> &children,OBAtom *check)
  {
    int nr_at = NumAtoms();
    OBBitVec used,curr,next;
    used.Resize((unsigned)nr_at+1);
    curr.Resize((unsigned)nr_at+1);
    next.Resize((unsigned)nr_at+1);

    used |= check->GetIdx();
    curr |= check->GetIdx();
    children.clear();

    int i;
    OBAtom *atom,*nbr;
    vector<OBBond*>::iterator j;

    for (;;)
    {
      next.Clear();
      for (i = curr.NextBit(-1);i != curr.EndBit();i = curr.NextBit(i))
      {
        atom = GetAtom(i);
        for (nbr = atom->BeginNbrAtom(j);nbr;nbr = atom->NextNbrAtom(j))
          if (!used[nbr->GetIdx()])
          {
            children.push_back(nbr);
            next |= nbr->GetIdx();
            used |= nbr->GetIdx();
          }
      }
      if (next.Empty())
        break;
      curr = next;
    }
    return used;
  }

  OBAggregate &OBAggregate::operator=(const OBMol &source)
  {
    if (this == &source)
      return *this;

    OBMol &src = (OBMol &)source;
    vector<OBAtom*>::iterator i;
    vector<OBBond*>::iterator j;
    OBAtom *atom;
    OBBond *bond;

    Clear();
    BeginModify();

    _vatom.reserve(src.NumAtoms());
    _atomIds.reserve(src.NumAtoms());
    _vbond.reserve(src.NumBonds());
    _bondIds.reserve(src.NumBonds());

    for (atom = src.BeginAtom(i);atom;atom = src.NextAtom(i))
      AddAtom(*atom);
    for (bond = src.BeginBond(j);bond;bond = src.NextBond(j))
      AddBond(*bond);

    this->_title  = src.GetTitle();
    this->_energy = src.GetEnergy();
    this->_dimension = src.GetDimension();
    this->SetTotalCharge(src.GetTotalCharge()); //also sets a flag
    this->SetTotalSpinMultiplicity(src.GetTotalSpinMultiplicity()); //also sets a flag

    EndModify(); //zeros flags!

    //!!!FLAGS ARE NOT COPIED OVER!!!

    //Copy Residue information
    unsigned int NumRes = src.NumResidues();
    if (NumRes)
      {
        unsigned int k;
        OBResidue *src_res=NULL;
        OBResidue *res=NULL;
        OBAtom *src_atom=NULL;
        OBAtom *atom=NULL;
        vector<OBAtom*>::iterator ii;
        for (k=0 ; k<NumRes ; ++k)
          {
            res = NewResidue();
            src_res = src.GetResidue(k);
            res->SetName(src_res->GetName());
            res->SetNum(src_res->GetNumString());
            res->SetChain(src_res->GetChain());
            res->SetChainNum(src_res->GetChainNum());
            for (src_atom=src_res->BeginAtom(ii) ; src_atom ; src_atom=src_res->NextAtom(ii))
              {
                atom = GetAtom(src_atom->GetIdx());
                res->AddAtom(atom);
                res->SetAtomID(atom,src_res->GetAtomID(src_atom));
                res->SetHetAtom(atom,src_res->IsHetAtom(src_atom));
                res->SetSerialNum(atom,src_res->GetSerialNum(src_atom));
              }
          }
      }

    //Copy conformer information
    if (src.NumConformers() > 1) {
      int k;//,l;
      vector<double*> conf;
      int currConf = -1;
      double* xyz = NULL;
      for (k=0 ; k<src.NumConformers() ; ++k) {
        xyz = new double [3*src.NumAtoms()];
        memcpy( xyz, src.GetConformer(k), sizeof( double )*3*src.NumAtoms() );
        conf.push_back(xyz);

        if( src.GetConformer(k) == src.GetCoordinates() ) {
          currConf = k;
        }
      }

      SetConformers(conf);
      if( currConf >= 0 && _vconf.size() ) {
        _c = _vconf[currConf];
      }
    }

    //Copy all the OBGenericData, providing the new molecule, this,
    //for those classes like OBRotameterList which contain Atom pointers
    //OBGenericData classes can choose not to be cloned by returning NULL
    vector<OBGenericData*>::iterator itr;
    for(itr=src.BeginData();itr!=src.EndData();++itr)
      {
        OBGenericData* pCopiedData = (*itr)->Clone(this);
        SetData(pCopiedData);
      }

    // copy chiral data for all atoms
    FOR_ATOMS_OF_MOL (atom, src) {
      if (atom->HasData(OBGenericDataType::ChiralData)) {
        OBChiralData* cd = (OBChiralData*) atom->GetData(OBGenericDataType::ChiralData);
        OBGenericData* pCopiedData = cd->Clone(NULL); // parent not used in OBChiralData::Clone()
        GetAtom(atom->GetIdx())->SetData(pCopiedData);
      }
    }

    if (src.HasChiralityPerceived())
      SetChiralityPerceived();

     _needRefresh = true;
     FindAllConnections();
     _nrMolecules = GetNrMolecules();

    return(*this);
  }


  //! Get the number of molecules in the aggregate.
  int OBAggregate::GetNrMolecules(){
      if (_needRefresh){
          FindAllConnections();
      }
      return _connections.size();
  }

  //! Print which atoms belong to which molecules.
  void OBAggregate::PrintConnections(){
      if (_needRefresh){
          FindAllConnections();
      }
      vector<vector<OBAtom*>* >::iterator outer_it;
      vector<OBAtom*>::iterator inner_it;
      int molcount = 0;
      for (outer_it = _connections.begin(); outer_it != _connections.end(); ++outer_it){
        std::cout<< "MolNr. " << molcount << ": ";
        int id = -1;
        int elements = 0;
        for (inner_it = (*outer_it)->begin(); inner_it != (*outer_it)->end(); ++inner_it){
            OBAtom* atom = *inner_it;
            if ( atom->GetIdx() != id+1 ){
                if ( id != -1 ){
                    std::cout<< " - " << id << ", " << std::endl;
                    std::cout<< atom->GetIdx();
                }
                id = atom->GetIdx();
                std::cout<< id;
            }
            else{
                ++id;
            }
            ++elements;
        }
        std::cout<< " - " << id;
        std::cout<< " = " << elements  << " atoms" << std::endl;
        ++molcount;
      }
      std::cout<< std::endl;
  }

  //! \brief Enable the use of Tags.
  //!
  //! This causes all methods that treat "parts" to use all part-indices
  //! as indices for tagged molecule assemblies.
  void OBAggregate::EnableTags(){
      _useTag = true;
  }

  //! \brief Disable the use of Tags.
  //!
  //! This causes all methods that treat "parts" to use all part-indices
  //! as indices for molecules in the aggregate.
  void OBAggregate::DisableTags(){
      _useTag = false;
      for (vector<vector<OBAtom*>* >::iterator it = _tags.begin(); it != _tags.end(); ++it){
          (*it)->clear();
      }
      _tags.clear();
      _molTags.clear();
  }

  //! \brief Get the number of tags
  //!
  //! Makes it easy to apply something to all tags. The value "-1"
  //! is returned if tags are not being used.
  int OBAggregate::GetNrTags(){
      if (_useTag){
          return _molTags.size();
      }
      else{
          return -1;
      }
  }

  //! Determine whether or not this aggregate uses tagging.
  bool OBAggregate::GetUseTag(){
      return _useTag;
  }

  //! Print which molecules belong to which tags.
  bool OBAggregate::PrintTags(){
      if (not(_useTag)){
          std::cerr << "ERROR: Aggregate not configured to use tags.";
          return false;
      }
      if (_needRefresh){
          FindAllConnections();
      }
      
      std::vector<bool> checked;
      vector<bool>::iterator bool_it;
      checked.reserve(_nrMolecules);
      for (int i=0; i<_nrMolecules; ++i){
          checked.push_back(false);
      }

      vector<vector<int> >::iterator outer_it;
      vector<int>::iterator inner_it;

      int tagged = 0;
      int tagcount = 0;
      for (outer_it = _molTags.begin(); outer_it != _molTags.end(); ++outer_it){
          if ( outer_it->size() > 0 ){
              std::cout << "TagNr. " << tagcount << ": ";
              for (inner_it = outer_it->begin(); inner_it != outer_it->end(); ++inner_it){
                  std::cout << *inner_it << " ";
                  checked[*inner_it] = true;
                  ++tagged;
              }
              ++tagcount;
              std::cout << std::endl;
          }
      }

      if ( tagged < _nrMolecules ){
          std::cout << "Untagged: ";
          int molcount = 0;
          for (bool_it = checked.begin(); bool_it != checked.end(); ++bool_it){
              if ( not(*bool_it) ){
                  std::cout << molcount << " ";
              }
              ++molcount;
          }
          std::cout << std::endl;
      }

      return true;
  }

  //! Create a new tag (index is returned) and optionally reserve space
  //! for a certain number of molecules.
  int OBAggregate::CreateTag(int nr_elements){
      if (not(_useTag)){
          return -1;
      }

      int newsize = _molTags.size()+1;
      _molTags.reserve(newsize);

      std::vector<int> tempvec;
      tempvec.reserve(nr_elements);

      _molTags.push_back(tempvec);

      _needRefresh = true;
      return newsize-1;
  }

  //! Add a molecule to a tag. Both are declared using their indices.
  bool OBAggregate::AddToTag(int molnr, int tag){
      if (not(_useTag)){
          return false;
      }
      if (tag >= _molTags.size()){
          return false;
      }

      _molTags[tag].push_back(molnr);

      _needRefresh = true;
      return true;
  }

  void OBAggregate::RotatePart(const int part, const double axis[3], const double angle){
      matrix3x3 matrix;
      matrix.RotAboutAxisByAngle(vector3(axis[0],axis[1],axis[2]),angle);
      if (_needRefresh){
          FindAllConnections();
      }
      if ( _useTag ){
          Rotate(*(_tags.at(part)), matrix);
      }
      else{
          Rotate(*(_connections.at(part)), matrix);
      }
  }

  void OBAggregate::RotatePart(const int part, const int main_axis_nr, const double angle){
      if (_needRefresh){
          FindAllConnections();
      }
      if ( _useTag ){
          Rotate(*(_tags.at(part)), main_axis_nr, angle);
      }
      else{
          Rotate(*(_connections.at(part)), main_axis_nr, angle);
      }
  }

  void OBAggregate::TranslatePart(const int part, const double vector[3]){
      if (_needRefresh){
          FindAllConnections();
      }
      if ( _useTag ){
          Translate(*(_tags.at(part)), vector);
      }
      else{
          Translate(*(_connections.at(part)), vector);
      }
  }

  void OBAggregate::TranslatePart(const int part, const vector3 vec){
      double v[3];
      v[0]=vec[0]; v[1]=vec[1]; v[2]=vec[2];
      TranslatePart(part, v);
  }

  void OBAggregate::AlignPart(const int part, const double p[3], const double v1[3], const double v2[3]){
      AlignPart(part, vector3(p[0],p[1],p[2]), vector3(v1[0],v1[1],v1[2]), vector3(v2[0],v2[1],v2[2]));
  }

  void OBAggregate::AlignPart(const int part, const vector3 p, const vector3 v1, const vector3 v2){
      if (_needRefresh){
          FindAllConnections();
      }
      if ( _useTag ){
          Align(*(_tags.at(part)), p, v1, v2);
      }
      else{
          Align(*(_connections.at(part)), p, v1, v2);
      }
  }

  void OBAggregate::MirrorPart(const int part, const double normal[3], const double point[3], bool center_it){
      MirrorPart(part, vector3(normal[0],normal[1],normal[2]), vector3(point[0],point[1],point[2]), center_it);
  }

  void OBAggregate::MirrorPart(const int part, const vector3 normal, const vector3 point, bool center_it){
      if (_needRefresh){
          FindAllConnections();
      }
      if ( _useTag ){
          Mirror(*(_tags.at(part)), normal, point, center_it);
      }
      else{
          Mirror(*(_connections.at(part)), normal, point, center_it);
      }
  }

  vector3 OBAggregate::GetCenterPart(int part){
    std::vector<OBAtom*>::iterator it;
    std::vector<OBAtom*> tempvec;
    vector3 center;
    OBAtom *atom;
    int elements;

    if (_needRefresh){
      FindAllConnections();
    }
    
    if ( _useTag ){
      tempvec = *(_tags.at(part));
    }
    else{
      tempvec = *(_connections.at(part));
    }

    elements = 0;
    for (it = tempvec.begin(); *it != NULL && it != tempvec.end(); ++it)
    {
      ++elements;
      atom = *it;
      center += atom->GetVector();
    }
    center /= elements;
    return center;
  }

  //! Move two parts of an aggregate closer together until their vdW surfaces overlap or the distante
  //! between the two increases.
  bool OBAggregate::MovePartsCloser(const double vector[3], const int part1, const int part2,
      const double stepsize, const double vdw_factor, double vdw_added){

    if (_needRefresh){
      FindAllConnections();
    }

    bool result = true;
    vector3 center1       = GetCenterPart(part1);
    vector3 center2       = GetCenterPart(part2);
    double distance_2     = (center1-center2).length_2();
    double new_distance_2;
    
    vector3 vec1 = vector3( vector[0],  vector[1],  vector[2]);
    vector3 vec2 = vector3(-vector[0], -vector[1], -vector[2]);

    if (vec1.CanBeNormalized() && vec2.CanBeNormalized()){

      vec1.normalize();
      vec2.normalize();

      vec1 *= stepsize/2.0;
      vec2 *= stepsize/2.0;
      
      bool loop = true;
      while ( loop ){
        center1        += vec1;
        center2        += vec2;
        new_distance_2  = (center1-center2).length_2();
        if ( new_distance_2 <= distance_2 ){ 
          TranslatePart(part1, vec1);
          TranslatePart(part2, vec2);
          if (IsGoodVDW(vdw_factor, vdw_added)){
            distance_2 = new_distance_2;
          }
          else{
            TranslatePart(part1, -vec1);
            TranslatePart(part2, -vec2);
            loop = false;
          }
        }
        else{ loop = false; }
      }
    }
    else{
      result=false;
    }
    return result;

  }

  bool OBAggregate::MovePartsCloser(const int part1, const int part2, 
      const double stepsize, const double vdw_factor, double vdw_added){

    vector3 center1 = GetCenterPart(part1);
    vector3 center2 = GetCenterPart(part2);
    double vector[3] = {0.0,0.0,0.0};

    vector[0] = center2[0] - center1[0];
    vector[1] = center2[1] - center1[1];
    vector[2] = center2[2] - center1[2];

    return MovePartsCloser(vector, part1, part2, stepsize, vdw_factor, vdw_added);
  }

  //parent constructor is called automaticaly for default constructor
  OBAggregate::OBAggregate()
  {
    _useTag=false;
    _nrMolecules = 0;
    _connections.clear();
    _needRefresh = false;
    _tags.clear();
    _molTags.clear();
  }

  OBAggregate::OBAggregate(const OBAggregate &agg) : OBMol(agg)
  {
    _useTag = agg._useTag;
    _needRefresh = true;
    _tags.clear();
    _molTags.clear();
    _molTags.reserve(agg._molTags.size());
    for (int i = 0; i<agg._molTags.size(); ++i){
        std::vector<int> tempvec;
        tempvec.reserve((agg._molTags[i]).size());
        for (int j=0; j<(agg._molTags[i]).size(); ++j){
            tempvec.push_back((agg._molTags[i])[j]);
        }
        _molTags.push_back(tempvec);
    }
    _connections.clear();
    _natoms = _nbonds = 0;
    _mod = 0;
    _totalCharge = 0;
    _dimension = 3;
    _vatom.clear();
    _atomIds.clear();
    _vbond.clear();
    _bondIds.clear();
    _vdata.clear();
    _title = "";
    _c = (double*)NULL;
    _flags = 0;
    _vconf.clear();
    _autoPartialCharge = true;
    _autoFormalCharge = true;
    _energy = 0.0;
    *this = agg;
    FindAllConnections();
    _nrMolecules = _connections.size();
  }

  OBAggregate::OBAggregate(const OBMol &mol) : OBMol(mol)
  {
    _useTag=false;
    _needRefresh = true;
    _tags.clear();
    _molTags.clear();
    _connections.clear();
    _natoms = _nbonds = 0;
    _mod = 0;
    _totalCharge = 0;
    _dimension = 3;
    _vatom.clear();
    _atomIds.clear();
    _vbond.clear();
    _bondIds.clear();
    _vdata.clear();
    _title = "";
    _c = (double*)NULL;
    _flags = 0;
    _vconf.clear();
    _autoPartialCharge = true;
    _autoFormalCharge = true;
    _energy = 0.0;
    *this = mol;
    FindAllConnections();
    _nrMolecules = _connections.size();
  }

  //parent destructor is called automaticaly for default destructor
  OBAggregate::~OBAggregate()
  {
      vector<vector<OBAtom*>* >::iterator it;
      for (it = _connections.begin(); it != _connections.end(); ++it){
          delete *it;
      }
      for (it = _tags.begin(); it != _tags.end(); ++it){
          delete *it;
      }
  }

} // end namespace OpenBabel
