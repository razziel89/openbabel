/**********************************************************************
confsearch.cpp - Conformer Searching Routines (see also forcefield.cpp)

Copyright (C) 2010 Noel O'Boyle <baoilleach@gmail.com>
Some portions Copyright (C) 2016 by Torsten Sachse

This file is part of the Open Babel project.
For more information, see <http://openbabel.sourceforge.net/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

***********************************************************************/

#include <openbabel/babelconfig.h>

#include <openbabel/forcefield.h>
#include <openbabel/rotamer.h>
#include <openbabel/rotor.h>
#include <openbabel/math/align.h>
#include <openbabel/tree/tree.hh>
#include <openbabel/tree/tree_util.hh>
#include <openbabel/math/vector3.h>

#include <float.h> // For DBL_MAX
#include <algorithm> // For min and std::sort
#include <limits.h> // For UINTS_MAX with certain old GCC4
#include <map> //for std::map
#include <utility> //for std::pair
#include <functional> //for std::hash
#include <deque> //for std::deque

#include <iomanip> // For setprecision

namespace OpenBabel
{

  class LFSR
  {
    /**
     * Usage:
     *     LFSR lfsr(N); // where N < 2^31
     *     unsigned int d;
     *     do {
     *       d = lfsr.GetNext();
     *       // .... do something with d ...
     *     } while (d != 1);
     **/
  public:
    LFSR(unsigned int range, unsigned int start);
    unsigned int GetNext(); // Return 1 when finished
  private:
    unsigned int _range, _lfsr, _poly;
    static const unsigned int _polynomials[31];
  };

  const unsigned int LFSR::_polynomials[31] = {0x3, 0x6, 0xc, 0x14, 0x30, 0x60, 0xb8, 0x110, 0x240, 0x500,
       0x829, 0x100d, 0x2015, 0x6000, 0xd008, 0x12000, 0x20400, 0x40023, 0x90000, 0x140000,
       0x300000, 0x420000, 0xe10000, 0x1200000, 0x2000023, 0x4000013, 0x9000000, 0x14000000,
       0x20000029, 0x48000000, 0x80200003};

  LFSR::LFSR(unsigned int range, unsigned int start = 1): _range(range), _lfsr(start)
  {
    assert ( _range < (1 << 31) ); // We can only handle up to 2^31 - 1
    // (Currently don't use a start value)
    // assert ( start < _range ); // Otherwise the _start value will never be returned

    int i = 0;
    unsigned int tot = 4;
    while (tot <= _range) {
      i++;
      tot <<= 1;
    }
    _poly = _polynomials[i];
  }

  inline unsigned int LFSR::GetNext()
  {
    do {
      _lfsr = (_lfsr >> 1) ^ (unsigned int)(0 - (_lfsr & 1u) & _poly);
    } while (_lfsr > _range);

    return _lfsr;
  }

  class OBDiversePoses {
    public:
      OBDiversePoses(const OBMol &ref, double RMSD, bool percise=false);
      ~OBDiversePoses() {
        delete palign; // Allocated with 'new'
      }
      bool AddPose(double* coords, double energy, int conf_nr);
      bool AddPose(vector<vector3> coords, double energy, int conf_nr);
      // The PosePair is used to remember a conformer's geometry and its energy.
      // However, for some purposes, it is beneficial to also remember
      // additional information (here: a number).
      typedef triple<vector<vector3>, double, int> PoseData;
      //typedef pair<vector<vector3>, double> PoseData;
      typedef tree<PoseData> Tree;
      Tree* GetTree() { return &poses; }
      typedef tree<PoseData>::iterator Tree_it;
      typedef tree<PoseData>::sibling_iterator Tree_sit;
      size_t GetSize();
      inline int GetNRMSD() {
        return n_rmsd;
      }
      inline double GetCutoff() {
        return cutoff;
      }
      bool KeepHydrogens(OBBitVec mask);

    private:
      bool _percise;
      vector<vector3> GetHeavyAtomCoords(const vector<vector3> &all_coords);
      int natoms;
      Tree poses;
      std::vector<double> levels;
      OBAlign* palign;
      const double cutoff;
      int n_rmsd;
      OBBitVec hydrogens;
    };

  bool OBDiversePoses::KeepHydrogens(OBBitVec mask){
      OBBitVec tmp = hydrogens;
      tmp.Negate();
      if (!(tmp & mask).IsEmpty())
          return false;
      hydrogens ^= mask;
      return true;
  }

  OBDiversePoses::OBDiversePoses(const OBMol &ref, double RMSD, bool percise):
          palign(new OBAlign(false, percise)), cutoff(RMSD), _percise(percise)
  {
    natoms = ref.NumAtoms();
    palign->SetRefMol(ref);
    palign->SetMethod(OBAlign::QCP);
    n_rmsd = 0;

    const double arr[] = {3.0, 2.0, 1.5, 1.0, 0.5, 0.25};
    std::vector<double> vec (arr, arr + sizeof(arr) / sizeof(arr[0]) );
    vec.erase(std::remove_if(vec.begin(), vec.end(), std::bind2nd(std::less<double>(), (cutoff + 0.1) )), vec.end());
    vec.push_back(cutoff);

    levels = vec;
    vector<vector3> pdummy;
    poses.insert(poses.begin(), PoseData(pdummy, 0.0, 0)); // Add a dummy top node

    // Remember the hydrogens
    hydrogens.Resize(natoms);
    for (int i=1; i<=natoms; i++)
      if (ref.GetAtom(i)->IsHydrogen())
        hydrogens.SetBitOn(i - 1);
  }

  bool OBDiversePoses::AddPose(double* coords, double energy, int conf_nr=0) {
    vector<vector3> vcoords;
    vcoords.reserve(natoms);
    for (unsigned int a = 0; a < natoms; ++a)
      vcoords.push_back(vector3(coords[a*3], coords[a*3+1], coords[a*3+2]));
    return AddPose(vcoords, energy, conf_nr);
  }

  bool OBDiversePoses::AddPose(vector<vector3> vcoords, double energy, int conf_nr=0) {
    Tree_it node = poses.begin();
    int level = 0;
    bool first_time = true;

    // Convert coords to vector<vector3>

    // Only use the heavy-atom coords for the alignment, but store
    // the full set of coordinates in the tree
    vector<vector3> vcoords_hvy = GetHeavyAtomCoords(vcoords);
    palign->SetRef(vcoords_hvy);

    vector<Tree_it> nodes, min_nodes;
    vector<double> min_nodes_rmsds;
    nodes.push_back(node);
    vector<int> stack_levels;
    stack_levels.push_back(level);

    vector<Tree_it> insert_pt;
    vector<PoseData> insert_data;
    vector<int> insert_level;

    while(nodes.size() > 0) { // Using stack-based recursion
      node = nodes.back();
      nodes.pop_back();
      level = stack_levels.back();
      stack_levels.pop_back();

      // Find whether the molecule is similar to any of the children of this node.
      // - min_node will hold the result of this search

      min_nodes.clear();
      min_nodes_rmsds.clear();
      double rmsd;
      double min_rmsd = DBL_MAX;

      Tree_sit sib = poses.begin(node);
      // Skip the first child after the first time through this loop
      // - it will already have been tested against at the end of the previous loop
      if (!first_time)
        ++sib;
      for (; sib != poses.end(node); ++sib) { // Iterate over children of node
        vector<vector3> tcoords = (*sib).first;
        vector<vector3> tcoords_hvy = GetHeavyAtomCoords(tcoords);
        palign->SetTarget(tcoords_hvy);

        palign->Align();
        rmsd = palign->GetRMSD();
        n_rmsd++;
        if (rmsd < levels.at(level)) {
          if (rmsd < cutoff)
            return false;

          min_nodes.push_back(sib);
          min_nodes_rmsds.push_back(rmsd);

          if (!_percise) // Exit as soon as one is found
            break;

        }
      } // end of for loop

      if (min_nodes.size() == 0) {
        // No similar molecule found, so remember it for later so that we can
        // append it the children and add it as the first child all the way down
        // through the levels. The reason we don't add it now is that the molecule
        // could still be rejected for addition to the tree.
        insert_pt.push_back(node);
        insert_level.push_back(level);
        insert_data.push_back(PoseData(vcoords, energy, conf_nr));
        continue;
      }

      // If we reach here, then a similar molecule was found
      int startlevel = level;
      vector<double>::const_iterator n_rmsd = min_nodes_rmsds.begin();
      for (vector<Tree_it>::iterator n = min_nodes.begin(); n != min_nodes.end(); ++n, ++n_rmsd) {
        node = *n;
        rmsd = *n_rmsd;
        level = startlevel + 1;
        while (rmsd < levels.at(level)) {
          node = poses.child(node, 0); // Get the first child
          level++;
        }
        nodes.push_back(node);
        stack_levels.push_back(level);
      }

      first_time = false;

    } // end of while loop


    // If we get here, then the molecule has been accepted for addition to the tree
    vector<PoseData>::iterator b = insert_data.begin();
    vector<int>::iterator c = insert_level.begin();
    for (vector<Tree_it>::iterator a = insert_pt.begin(); a != insert_pt.end(); ++a, ++b, ++c) {
      node = *a;
      for (int k = *c; k < levels.size(); ++k) {
        node = poses.append_child(node, *b);
      }
    }

    return true;
  }

  size_t OBDiversePoses::GetSize() {
    return poses.size() - 1; // Remove the dummy
  }

  vector<vector3> OBDiversePoses::GetHeavyAtomCoords(const vector<vector3> &all_coords) {
    vector<vector3> v_hvyatoms;
    for (unsigned int a = 0; a < natoms; ++a)
      if (!hydrogens.BitIsSet(a))
        v_hvyatoms.push_back(all_coords[a]);
    return v_hvyatoms;
  }

  template <typename T>
  bool sort_by_second(const T& a, const T& b) {
    return (a.second < b.second);
  }
  template <typename T>
  bool sort_by_first(const T& a, const T& b) {
    return (a.first < b.first);
  }

vector<vector3> GetHeavyAtomCoords(const OBMol* mol, const vector<vector3> &all_coords) {
  vector<vector3> v_hvyatoms;
  for (unsigned int a = 1; a <= mol->NumAtoms(); ++a)
    if (!mol->GetAtom(a)->IsHydrogen())
      v_hvyatoms.push_back(all_coords[a]);
  return v_hvyatoms;
}

void UpdateConformersFromTree(OBMol* mol, vector<double> &energies, OBDiversePoses* divposes, bool verbose,
        bool precise=true, bool sort=true, std::vector<int>* new_confs=NULL)
  {

  OBDiversePoses::Tree* poses = divposes->GetTree();
  double cutoff = divposes->GetCutoff();

  vector <OBDiversePoses::PoseData> confs, newconfs;

  // The leaf iterator will (in effect) iterate over the nodes just at the loweset level
  for (OBDiversePoses::Tree::leaf_iterator node = poses->begin(); node != poses->end(); ++node)
    if (node->first.size() > 0) // Don't include the dummy head node
      confs.push_back(*node);

  // Sort the confs by energy (lowest first)
  if (sort)
    std::sort(confs.begin(), confs.end(), sort_by_second<OBDiversePoses::PoseData>);

  if(verbose)
    cout << "....tree size = " << divposes->GetSize() <<  " confs = " << confs.size() << "\n";

  typedef vector<OBDiversePoses::PoseData> vpp;

  // Loop through the confs and filter using a tree
  newconfs.clear();
  OBDiversePoses newtree(*mol, cutoff, precise);
  for (vpp::iterator conf = confs.begin(); conf!=confs.end(); ++conf) {
    if (newtree.AddPose(conf->first, conf->second)) {
      newconfs.push_back(*conf);
    }
  }
  if (verbose)
    cout << "....new tree size = " << newtree.GetSize() <<  " confs = " << newconfs.size() << "\n";

  if (new_confs!=NULL){
    new_confs->reserve(newconfs.size());
    // Only add the indices of the conformers that were kept and do not modify
    // the molecule.
    for (vpp::iterator chosen = newconfs.begin(); chosen!=newconfs.end(); ++chosen) {
        new_confs->push_back(chosen->third);
        //energies.push_back(chosen->second);

        //// To avoid making copies of vectors or vector3s, I am using pointers throughout
        //vector<vector3> *tmp = &(chosen->first);
        //double *confCoord = new double [mol->NumAtoms() * 3];
        //for(unsigned int a = 0; a<mol->NumAtoms(); ++a) {
        //  vector3* pv3 = &(*tmp)[a];
        //  confCoord[a*3] = pv3->x();
        //  confCoord[a*3 + 1] = pv3->y();
        //  confCoord[a*3 + 2] = pv3->z();
        //}
        //mol->AddConformer(confCoord);
    }
  }else{
    // Add confs to the molecule's conformer data and add the energies to molecules's energies
    energies.reserve(energies.size()+newconfs.size());
    for (vpp::iterator chosen = newconfs.begin(); chosen!=newconfs.end(); ++chosen) {
      energies.push_back(chosen->second);

      // To avoid making copies of vectors or vector3s, I am using pointers throughout
      vector<vector3> *tmp = &(chosen->first);
      double *confCoord = new double [mol->NumAtoms() * 3];
      for(unsigned int a = 0; a<mol->NumAtoms(); ++a) {
        vector3* pv3 = &(*tmp)[a];
        confCoord[a*3] = pv3->x();
        confCoord[a*3 + 1] = pv3->y();
        confCoord[a*3 + 2] = pv3->z();
      }
      mol->AddConformer(confCoord);
    }
  }
}

int OBForceField::DiverseConfGen(double rmsd, unsigned int nconfs, double energy_gap, bool verbose)
  {
    _energies.clear(); // Wipe any energies from previous conformer generators

    // Remove all conformers (e.g. from previous conformer generators) even the current conformer
    double *initialCoord = new double [_mol.NumAtoms() * 3]; // initial state
    double *store_initial = new double [_mol.NumAtoms() * 3]; // store the initial state
    memcpy((char*)initialCoord,(char*)_mol.GetCoordinates(),sizeof(double)*3*_mol.NumAtoms());
    memcpy((char*)store_initial,(char*)_mol.GetCoordinates(),sizeof(double)*3*_mol.NumAtoms());
    std::vector<double *> newConfs(1, initialCoord);
    _mol.SetConformers(newConfs);

    if (_mol.NumRotors() == 0) {
      SetupPointers();
      _energies.push_back(Energy(false));
      delete [] store_initial;
      return 0;
    }

    // Get estimate of lowest energy conf using FastRotorSearch
    FastRotorSearch(true);
    double lowest_energy = Energy(false);

    int origLogLevel = _loglvl;

    OBRotorList rl;
    OBBitVec fixed = _constraints.GetFixedBitVec();
    rl.SetFixAtoms(fixed);
    if (_loglvl == 0)
      rl.SetQuiet(); // Don't print info on symmetry removal
    rl.Setup(_mol);

    OBRotorIterator ri;
    OBRotamerList rotamerlist;
    rotamerlist.SetBaseCoordinateSets(_mol);
    rotamerlist.Setup(_mol, rl);

    // Can take shortcut later, as 4 components of the energy will be constant
    SetupPointers();
    double energy_offset = E_Bond(false) + E_Angle(false) + E_StrBnd(false) + E_OOP(false);
    lowest_energy -= energy_offset;
    _energies.push_back(lowest_energy);

    OBRotorKeys rotorKeys;
    OBRotor* rotor = rl.BeginRotor(ri);
    unsigned int combinations = 1;
    vector<size_t> rotor_sizes;
    for (int i = 1; i < rl.Size() + 1; ++i, rotor = rl.NextRotor(ri)) { // foreach rotor
      size_t size = rotor->GetResolution().size();
      rotorKeys.AddRotor(size);
      combinations *= size;
      rotor_sizes.push_back(size);
      if(verbose) {
        cout << "....rotor " << i << " from " << rotor->GetBond()->GetBeginAtomIdx() << " to "
             << rotor->GetBond()->GetEndAtomIdx() << " has " << size << " values" << endl;
      }
    }
    if (rotor_sizes.size() > 0 && combinations == 0) { // Overflow!
      combinations = UINT_MAX;
    }
    cout << "..tot conformations = " << combinations << "\n";

    if (nconfs == 0)
      nconfs = 1 << 20;
    unsigned int max_combinations = min<unsigned int>(nconfs , combinations);
    LFSR lfsr(max_combinations); // Systematic random number generator
    if (verbose && combinations > max_combinations) {
      cout << "....Using a cutoff of "
           << nconfs << " we will only explore " << std::fixed << setprecision(1)
           << static_cast<float>(nconfs * 100)/static_cast<float>(combinations) << "% of these\n";
    }

    unsigned int combination;
    OBDiversePoses divposes(_mol, rmsd, false);
    vector<int> my_rotorkey(rotor_sizes.size() + 1, 0);
    int counter = 0;

    // Main loop over rotamers
    unsigned int N_low_energy = 0;
    do {
      _mol.SetCoordinates(store_initial);

      combination = lfsr.GetNext();
      unsigned int t = combination;
      // Convert the combination number into a rotorkey
      for (int i = 0 ; i < rotor_sizes.size(); ++i) {
        my_rotorkey[i + 1] = t % rotor_sizes[i];
        t /= rotor_sizes[i];
      }

      rotamerlist.SetCurrentCoordinates(_mol, my_rotorkey);
      SetupPointers();
      double currentE = E_VDW(false) + E_Torsion(false) + E_Electrostatic(false);
      if (currentE < lowest_energy + energy_gap) { // Don't retain high energy poses
        divposes.AddPose(_mol.GetCoordinates(), currentE);
        N_low_energy++;
        if (currentE < lowest_energy)
          lowest_energy = currentE;
      }
      counter++;
    } while (combination != 1 && counter < nconfs); // The LFSR always terminates with a 1
    cout << "..tot confs tested = " << counter << "\n..below energy threshold = " << N_low_energy << "\n";

    // Reset the coordinates to those of the initial structure
    _mol.SetCoordinates(store_initial);

    // Get results from the tree
    UpdateConformersFromTree(&_mol, _energies, &divposes, verbose);

    // Add back the energy offset
    transform(_energies.begin(), _energies.end(), _energies.begin(), bind2nd(plus<double>(), energy_offset));

    // Clean up
    delete [] store_initial;

    return 0;
 }


// This class allows to easily round coordinates. Here, this has to be used
// since conformers, that are actually identical, do not have to be bit-
// identical. Since identity is checked by computing a hash of a string-
// representation of the molecule, rouding has to be employed.
class iatom {
    public:
        int el, x ,y, z, fac; 
        iatom(int _el, double _x, double _y, double _z, int _fac): el(_el), fac(_fac)
        {
            x = (int)(_x*fac*10);
            y = (int)(_y*fac*10);
            z = (int)(_z*fac*10);
        };
        // Return a string representation of this atom's rounded coordinate
        // either by actually rounding (rounding==true) or by truncating
        // (rounding==false).
        std::string str(bool rounding){
            std::stringstream ss;
            if (rounding){
                ss << el << " " << int((x+5)/10) << " " << int((y+5)/10) << " " << int((z+5)/10);
            }else{
                ss << el << " " << int(x/10) << " " << int(y/10) << " " << int(z/10);
            }
            return ss.str();
        }
  };
  
inline bool operator< (const iatom& lhs, const iatom& rhs){
      if (lhs.el<rhs.el)
          return true;
      if (lhs.el>rhs.el)
          return false;
      if (lhs.x<rhs.x)
          return true;
      if (lhs.x>rhs.x)
          return false;
      if (lhs.y<rhs.y)
          return true;
      if (lhs.y>rhs.y)
          return false;
      if (lhs.z<rhs.z)
          return true;
      if (lhs.z>rhs.z)
          return false;
      return false;
  }
inline bool operator> (const iatom& lhs, const iatom& rhs){ return rhs < lhs; }

// Variables of type AtIdx can easily be sorted with respect to the operator< defined
// above, thus enforcing the strictly weak ordering required for the identity screening.
typedef std::pair<iatom,int> AtIdx;
typedef std::size_t Hash;
typedef std::pair<Hash,Hash> HashPair;

// Return two hashes of string representations of this molecule (one: truncating, two: rounding).
// The parameter "prec" declares to how many decimal places the conversion shall be performed.
HashPair GetSortedConformerHash(OBMol& mol, int conf_nr, int prec, std::vector<unsigned short>& atom_order, OBBitVec mask){
    int fac = 1;
    for (int i=prec; i>=1; --i)
        fac *= 10;
    int nr_atoms = mol.NumAtoms();
    std::hash<std::string> hash_fn;
    atom_order.reserve(nr_atoms);
    std::vector<AtIdx> symsort;
    symsort.reserve(nr_atoms);
    double* coords = mol.GetConformer(conf_nr);
    OBAtomIterator it=mol.BeginAtoms();
    for (int count=0; it!=mol.EndAtoms(); ++it, coords+=3, ++count){
        //if (not((*it)->IsHydrogen())){
        if (mask.BitIsSet(count)){
            iatom at = iatom((*it)->GetAtomicNum(), *(coords+0), *(coords+1), *(coords+2), fac);
            symsort.push_back(std::make_pair(at, (*it)->GetIdx()-1));
        }
    }

    // in order to check whether 2 conformers are symmetry equivalent,
    // a definitive order of atoms has to be enforced to avoid pairwise
    // checks of identity of atoms
    std::sort(symsort.begin(), symsort.end(), sort_by_first<AtIdx>);

    // two conformers are likely identical if the hashes of their
    // hydrogen-less string representations are equivalent
    std::stringstream ss_trunc;
    std::stringstream ss_round;
    for (std::vector<AtIdx>::iterator it = symsort.begin(); it!=symsort.end(); ++it){
        ss_trunc << it->first.str(false) << "\n";
        ss_round << it->first.str(true) << "\n";
        atom_order.push_back(it->second);
    }
    HashPair geom_hash = std::make_pair(hash_fn(ss_trunc.str()) , hash_fn(ss_round.str()));
    return geom_hash;
  }

// In order to perform identity screening (see below, search for "symscreening"), 
// the atoms have to be sorted in a strictly weak ordering. This type can hold
// a conformer number as well as a vector that allows for sorting the atoms in
// said ordering.
typedef std::pair< int,std::vector<unsigned short> > HashInfo;

// Compute the RMSD between two conformers of the same molecule. The coordinates
// of both conformers are taken in a sorted order specifyed by the second element
// of the HashInfo type.
double DirectRMSDComparison(OBMol& mol, const HashInfo& cur, const HashInfo& ref)
  {
    // cur_it will iterate over the sorted atom indices in the conformer that might be added to the order_map (cur)
    // ref_it will iterate over the sorted atom indices in the conformers already in the order_map (ref)
    // ref_it_end is just a helper variable to check for the end of the iteration
    // No more checks for the end of the iteration have to be performed since
    // all conformers must have the same number of non-hydrogen atoms.
    // This process DOES ignore hydrogen atoms.
    std::vector<unsigned short>::const_iterator cur_it, ref_it, ref_it_end;
    cur_it     = cur.second.begin();
    ref_it     = ref.second.begin();
    ref_it_end = ref.second.end();
    // This allows for comparing the geometries of conformers without
    // explicitly setting all atom coordinates using SetConformer.
    double* geom_cur = mol.GetConformer((int)cur.first);
    double* geom_ref = mol.GetConformer((int)ref.first);
    // Will contain the squared rmsd
    double rmsd_2 = 0.0;
    // The conformers will be considered "identical" if their RMSD is
    // below a certain threshold.
    for (; ref_it!=ref_it_end; ++cur_it, ++ref_it){
        double dist;
        // these now point to the x-coordinates of the two atoms
        double* coord_cur = geom_cur+(3*(*cur_it));
        double* coord_ref = geom_ref+(3*(*ref_it));
        // add contribution to RMSD
        dist = *(coord_cur+0) - *(coord_ref+0);
        rmsd_2 += dist*dist;
        dist = *(coord_cur+1) - *(coord_ref+1);
        rmsd_2 += dist*dist;
        dist = *(coord_cur+2) - *(coord_ref+2);
        rmsd_2 += dist*dist;
    }
    rmsd_2 /= mol.NumAtoms();
    return rmsd_2;
  }

// The main similarity screening function
int OBForceField::ScreenByRMSD(double rmsd, double egap, short prec, OBBitVec hydrogen_mask, bool verbose)
  {
    
    if (_mol.NumConformers()<=0){
        std::cerr << "ERROR: no conformers in given molecule, cannot perform similarity screening." << std::endl;
        return 1;
    }

    _energies.clear(); // Wipe any energies from previous algorithms
    _energies.reserve(_mol.NumConformers());
    
    // determine which kinds of screenings shall be performed
    bool escreening = egap > 0.0 ? true : false;
    bool rmsdscreening = rmsd > 0.0 ? true : false;
    if (!rmsdscreening){rmsd = 0.0;}
    bool symscreening = prec >= 0 ? true : false;

    if (!escreening && !rmsdscreening && !symscreening){
        if (verbose){
            cout << "None of the 3 screening procedures activated." << endl;
        }
        return 0;
    }

    // energy screening: do not keep conformers whose associated energy
    // is higher than that of the global minimum plus a given gap (egap)
    double lowest_energy = 0.0;
    if (escreening){
        _mol.SetConformer(0);
        SetupPointers();
        lowest_energy = Energy(false); // Energy(false) means "do not evaluate gradients"
        _energies.push_back(lowest_energy);
        //Determine energies of all conformers if screening by energies is desired.
        //Add all zeros otherwise.
        for (int conf_count = 1; conf_count < _mol.NumConformers(); ++conf_count){
            _mol.SetConformer(conf_count);
            SetupPointers();
            double temp_energy = Energy(false);
            if (temp_energy < lowest_energy){
                lowest_energy = temp_energy;
            }
            _energies.push_back(temp_energy);
        }
    }
    else{
        // _energies has been resized to hold at least _mol.NumConformers() numbers
        // but the energy is not to be considered for screening (hence fill _energies
        // will all zeroes to be able to use the same algorithm for everything else)
        std::fill(_energies.begin(), _energies.begin()+_mol.NumConformers(), 0.0);
    }

    if (!rmsdscreening && symscreening){
        std::cerr << "WARNING: no RMSD set but one is required for symmetry screening. Will use 0.10."
            << std::endl;
        rmsd = 0.10;
    }

    // symscreening: do not keep conformers that are identical This is done by
    // sorting all non-hydrogen atoms in a strict- weak order whilst keeping
    // the hydrogens at their places.  Here, the algorithm determines which
    // atoms can be sorted.  First, hashes are computed that are associated
    // with the coordinates of each conformer. Then, only such conformers whose
    // hashes are identical are considered to be possibly identical. Then, an
    // explicit check for identity (with respect to an RMSD) is performed. To
    // avoid rounding errors, two hashes are created per conformer (one by
    // rounding and one by truncating).
    // All this has to be done since the tree (divposes) cannot cope with very
    // large data sets that contain a large number of almost identical
    // molecules (runs out of RAM very quickly).  Hence, the identity screening
    // (also called symmetry screening) is performed beforehand to keep the
    // required RAM to a minimum.
    std::vector<unsigned short> sort_important;
    OBBitVec sort_important_mask = OBBitVec();
    sort_important_mask.Resize(_mol.NumAtoms());
    sort_important_mask.SetRangeOff(0,_mol.NumAtoms()-1);
    {
        unsigned short count = 0;
        sort_important.reserve(_mol.NumAtoms());
        for (OBAtomIterator it=_mol.BeginAtoms(); it!=_mol.EndAtoms(); ++it, ++count){
            if (not((*it)->IsHydrogen()) || hydrogen_mask.BitIsSet(count)){
                sort_important.push_back(count);
                sort_important_mask.SetBitOn(count);
            }
        }
    }

    // This tree will be used to screen by RMSD
    OBDiversePoses divposes(_mol, rmsd, false);

    bool success = divposes.KeepHydrogens(hydrogen_mask);
    if (!success){
        std::cerr << "ERROR: given hydrogen_mask also masks non-hydrogen atoms." << std::endl;
        return 1;
    }

    if (verbose){
        if (escreening){
            puts("Using energy screening.");
        }else{
            puts("Not using energy screening.");
        }
        if (rmsdscreening){
            puts("Using rmsd screening.");
        }else{
            puts("Not using rmsd screening.");
        }
        if (symscreening){
            puts("Using symmetry screening.");
        }else{
            puts("Not using symmetry screening.");
        }
    }

    // Save the number of conformers currently present in the molecule
    // since all of them have to be removed later on.
    int nr_confs = _mol.NumConformers();

    // This is the type of the map that will save the sorted atom order
    // (allowing for easy checks of identity of conformers).  Information about
    // every conformer that has the same hash has to be saved hence the
    // std::deque.
    typedef std::map<Hash,std::deque<HashInfo> > MapType;
    MapType order_map;

    struct sScreenCount;
    typedef struct sScreenCount {
        int e/*nergy*/, si /*strict_identity*/, wi /*weak_identity*/;
        sScreenCount() : e(0), si(0), wi(0){};
        int sum(){
            return e+si+wi;
        }
    } ScreenCount;
    ScreenCount sc = ScreenCount();

    // Main screening loop over conformers
    for (int conf_count = 0; conf_count < _mol.NumConformers(); ++conf_count){

        bool screened = false;
        _mol.SetConformer(conf_count);
        SetupPointers();

        double currentE = 0.0;
        // screen such conformers whose energy is too high
        if (escreening){
            currentE = _energies.at(conf_count);
            screened = (currentE >= lowest_energy + egap);
            if (screened)
                sc.e += 1;
        }

        // symmetry screening (see above for explanation)
        double* coords = _mol.GetCoordinates();
        // True if coords has to be free'd after having been added to divposes.
        // Coords can also contain atom positions in a rearranged order.
        bool free_me = false;
        // sort the atoms in the conformer to screen for symmetry
        if (!screened && symscreening){
            // contains the strictly weak ordering
            std::vector<unsigned short> atom_order;
            HashPair geom_hash = GetSortedConformerHash(_mol, conf_count, prec, atom_order, sort_important_mask);

            HashInfo order_element = std::make_pair(conf_count,atom_order);

            MapType::const_iterator it1 = order_map.find(geom_hash.first);
            MapType::const_iterator it2 = order_map.find(geom_hash.second);
            if (it1->first==it2->first){
                it2 = order_map.end();
            }
            if (it1 != order_map.end()){
                bool found = false;
                // If they share a hash, they still have to be checked for strict identity
                // with respect to the given RMSD (since the hashes can be equivalent just by accident).
                // This check has to be performed for every conformer with the determined hash.
                for (std::deque<HashInfo>::const_iterator data_it = it1->second.begin(); !found && data_it!=it1->second.end(); ++data_it){
                    double rmsd_2 = DirectRMSDComparison(_mol, order_element, *data_it);
                    found = (rmsd_2 <= rmsd*rmsd);
                }
                if (!found){
                    order_map[it1->first].push_back(order_element);
                }
                screened = (screened || found);
            }else{
                order_map[geom_hash.first].push_back(order_element);
            }
            if (it2 != order_map.end()){
                bool found = false;
                for (std::deque<HashInfo>::const_iterator data_it = it2->second.begin(); !found && data_it!=it2->second.end(); ++data_it){
                    double rmsd_2 = DirectRMSDComparison(_mol, order_element, *data_it);
                    found = (rmsd_2 <= rmsd*rmsd);
                }
                if (!found){
                    order_map[it2->first].push_back(order_element);
                }
                screened = (screened || found);
            }else{
                order_map[geom_hash.second].push_back(order_element);
            }

            // If this conformer is not strictly identical to any of the previous ones, it can
            // still be almost strictly identical (i.e., hashes are not equivalent). Hence, use
            // divposes to check for that (but use the correctly sorted atomic coordinates).
            // This reordering of coordinates is performed here.
            if (!screened){
                double* sorted_coords = (double*)malloc(sizeof(double)*3*_mol.NumAtoms());
                memcpy(sorted_coords, coords, sizeof(double)*3*_mol.NumAtoms());
                int count=0;
                // Only overwrite those which are to be resorted (i.e., non-hydrogen or important hydrogens)
                for (std::vector<unsigned short>::iterator it=atom_order.begin(); it!=atom_order.end(); ++it, ++count){
                    memcpy(sorted_coords+(3*(sort_important[count])), coords+(3*(*it)), 3*sizeof(double));
                }
                coords = sorted_coords;
                free_me = true;
            }
            if (screened)
                sc.si += 1;
        }

        // don't retain poses if at least one of:
        //  - energy too high
        //  - pose is symmetry equivalent to one that has already been processed
        if (!screened) {
            divposes.AddPose(coords, currentE, conf_count);
            if (free_me)
                free(coords);
        }

    }
    order_map.clear();

    std::vector<int> new_confs;
    _mol.SetConformer(0);
    SetupPointers();
    // This only saves the indices of the conformers that were kept to new_confs.
    UpdateConformersFromTree(&_mol, _energies, &divposes, verbose, false, escreening, &new_confs);
    sc.wi = nr_confs - sc.sum() - new_confs.size();
    // Hence they still have to be manually appended to the molecule
    for (std::vector<int>::iterator it = new_confs.begin(); it!=new_confs.end(); ++it){
        _mol.AddConformer(_mol.GetConformer(*it),true);
    }

    //Delete all old conformers from molecule and their associated energies
    _mol.DeleteConformers(0,nr_confs-1);
    _energies.clear();
    //_energies.erase(_energies.begin(),_energies.begin()+nr_confs);

    _mol.SetConformer(0);
    SetupPointers();

    std::cout << "..screened conformers due to:\n";
    if (escreening){
        std::cout << "    energy: " << sc.e << "\n";
    }
    if (symscreening){
        std::cout << "    strict identity: " << sc.si << "\n";
    }
    if (rmsdscreening){
        std::cout << "    rmsd: " << sc.wi << "\n";
    }
    if (!rmsdscreening && symscreening){
        std::cout << "    weak identity: " << sc.wi << "\n";
    }
    if (int(escreening) + int(symscreening) + int(rmsdscreening) + int(!rmsdscreening && symscreening) > 1){
        std::cout << "    total: " << sc.sum() << "\n";
        std::cout << "..screening was performed in the above order\n";
    }

    return 0;
 }

} // end of namespace OpenBabel

//! \file confsearch.cpp
//! \brief Conformer searching and screening routines
