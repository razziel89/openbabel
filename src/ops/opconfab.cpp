/**********************************************************************
opconfab.cpp - Confab, the conformer generator described in
               Journal of Cheminformatics, 3, 1, 8
               http://www.jcheminf.com/content/3/1/8

Copyright (C) 2013 by Noel O'Boyle

This file also contains an OBOp called SimScreen to screen molecules depending
on how similar they are, RMSD-wise. The code has been added below that
concerning Confab.

Copyright for SimScreen (C) 2015 by Torsten Sachse (torsten.sachse@uni-jena.de)

This file is part of the Open Babel project.
For more information, see <http://openbabel.org/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/

#include <openbabel/babelconfig.h>
#include <iostream>
#include <sstream>
#include <openbabel/op.h>
#include <openbabel/mol.h>
#include <openbabel/forcefield.h>
#include <openbabel/obconversion.h>
#include <openbabel/generic.h>
#include <openbabel/bitvec.h>

#define CONFAB_VER "1.1.0"

namespace OpenBabel
{
  using namespace std;

  //////////////////////////////////////////////////////////
  //
  //  Confab
  //
  //////////////////////////////////////////////////////////

  class Confab
  {
  public:
    Confab() {};

  };

  class OpConfab : public OBOp
  {
    public:
      OpConfab(const char* ID) : OBOp(ID, false) {
      }

      const char* Description()
      {
        return "Confab, the diverse conformer generator\n"
          "Typical usage: obabel infile.xxx -O outfile.yyy --confab --conf 1000000\n"
          "  options:\n"
          "    --conf #     Max number of conformers to test (default is 1 million)\n"
          "    --rcutoff #  RMSD cutoff (default 0.5 Angstrom)\n"
          "    --ecutoff #  Energy cutoff (default 50.0 kcal/mol)\n"
          "    --original   Include the input conformation as the first conformer\n"
          "    --verbose    Verbose output\n"
          ;
      }

      virtual bool WorksWith(OBBase* pOb) const
      {
        return dynamic_cast<OBMol*>(pOb) != NULL;
      }
      virtual bool Do(OBBase* pOb, const char* OptionText, OpMap* pmap, OBConversion*);
      
      void DisplayConfig(OBConversion* pConv);
      void Run(OBConversion* pConv, OBMol* pmol);
      double rmsd_cutoff;
      double energy_cutoff;
      unsigned int conf_cutoff;
      bool verbose;
      bool include_original;
      unsigned int N;
      OBForceField *pff;
  };

  //////////////////////////////////////////////////////////
  OpConfab theConfab("confab"); //Global instance



  //////////////////////////////////////////////////////////
  bool OpConfab::Do(OBBase* pOb, const char* OptionText, OpMap* pmap, OBConversion* pConv=NULL)
  {
    OBMol* pmol = dynamic_cast<OBMol*>(pOb);
    if(!pmol)
      return false;

    if(pConv->IsFirstInput())
    {
      pConv->AddOption("writeconformers", OBConversion::GENOPTIONS);
      rmsd_cutoff = 0.5;
      energy_cutoff = 50.0;
      conf_cutoff = 1000000; // 1 Million
      verbose = false;
      include_original = false;

      OpMap::const_iterator iter;
      iter = pmap->find("rcutoff");
      if(iter!=pmap->end())
        rmsd_cutoff = atof(iter->second.c_str());
      iter = pmap->find("ecutoff");
      if(iter!=pmap->end())
        energy_cutoff = atof(iter->second.c_str());
      iter = pmap->find("conf");
      if(iter!=pmap->end())
        conf_cutoff = atoi(iter->second.c_str());
      iter = pmap->find("verbose");
      if(iter!=pmap->end())
        verbose = true;
      iter = pmap->find("original");
      if(iter!=pmap->end())
        include_original = true;

      cout << "**Starting Confab " << CONFAB_VER << "\n";
      cout << "**To support, cite Journal of Cheminformatics, 2011, 3, 8.\n";
      pff = OpenBabel::OBForceField::FindType("mmff94");
      if (!pff) {
        cout << "!!Cannot find forcefield!" << endl;
        exit(-1);
      }
      DisplayConfig(pConv);
    }

    Run(pConv, pmol);

    return false;
  }

  void OpConfab::Run(OBConversion* pConv, OBMol* pmol)
  {
    OBMol mol = *pmol;
    
    N++;
    cout << "**Molecule " << N << endl << "..title = " << mol.GetTitle() << endl;
    cout << "..number of rotatable bonds = " << mol.NumRotors() << endl;
    mol.AddHydrogens();
    bool success = pff->Setup(mol);
    if (!success) {
      cout << "!!Cannot set up forcefield for this molecule\n"
           << "!!Skipping\n" << endl;
      return;
    }
    pff->DiverseConfGen(rmsd_cutoff, conf_cutoff, energy_cutoff, verbose);

    pff->GetConformers(mol);
    int nconfs = include_original ? mol.NumConformers() : mol.NumConformers() - 1;

    cout << "..generated " << nconfs << " conformers" << endl;

    unsigned int c = include_original ? 0 : 1;

    for (; c < mol.NumConformers(); ++c) {
      mol.SetConformer(c);
      if(!pConv->GetOutFormat()->WriteMolecule(&mol, pConv))
        break;
    }
    cout << endl;

  }

  void OpConfab::DisplayConfig(OBConversion* pConv)
  {
    cout << "..Input format = " << pConv->GetInFormat()->GetID() << endl;
    cout << "..Output format = " << pConv->GetOutFormat()->GetID() << endl;
    cout << "..RMSD cutoff = " << rmsd_cutoff << endl;
    cout << "..Energy cutoff = " << energy_cutoff << endl;
    cout << "..Conformer cutoff = " << conf_cutoff << endl;
    cout << "..Write input conformation? " << (include_original ? "True" : "False") << endl;
    cout << "..Verbose? " << (verbose ? "True" : "False") << endl;
    cout << endl;
  }

}//namespace


//
namespace OpenBabel
{
  using namespace std;

  //////////////////////////////////////////////////////////
  //
  //  SimScreen
  //
  //////////////////////////////////////////////////////////

  class SimScreen
  {
  public:
    SimScreen() {};
  };

  class OpSimScreen : public OBOp
  {
    public:

      enum ALIGNING {
          never = 0,
          before = 1,
          after = 2,
          always = 3,
          nr_elements = 4
      };

      OpSimScreen(const char* ID) : OBOp(ID, false) {
      }

      const char* Description()
      {
        return "SimScreen, screen given molecules with respect to their RMSD and energy.\n"
          "Typical usage: obabel infile.xxx -O outfile.yyy --simscreen --rcutoff 5 --ecutoff 50\n"
          "  options:\n"
          "    --rcutoff #  RMSD cutoff (default: do not use)\n"
          "    --ecutoff #  Energy cutoff (default: do not use)\n"
          "    --ffname #   Name of the force field to use (default mmff94)\n"
          "    --prec #     Number of decimal places used for determination\n"
          "                 of conformer equivalence due to symmetry. A high value\n"
          "                 reduces computational time but increases risk of getting\n"
          "                 duplicates (default: do not use symmetry screening,\n"
          "                 highly recommended: 2)\n"
          "    --ssalign #  Automatically align all conformers prior to or after the\n"
          "                 similarity screening procedure. The third and second main axes\n"
          "                 would be aligned to the x and y axes, respectively and the\n"
          "                 molecules' centers would be aligned to the origin.\n"
          "                 Accepted values are arbitrary combinations of 'b' and 'a'\n"
          "                 for 'before' and 'after' respectively. Note that the\n"
          "                 symmetry screening relies on the conformers being aligned\n"
          "                 in some uniform fashion.\n"
          "    --imp-H #    Declare a comma-separated list of indices (counting starts\n"
          "                 at 1) that indicate \"important\" hydrogens, i.e., such H-\n"
          "                 atoms that should not be discarded when computing an RMSD or\n"
          "                 when judging equivalence due to symmetry. This is useful in\n"
          "                 cases where, e.g., different tautomers of a molecule exist.\n"
          ;
      }

      virtual bool WorksWith(OBBase* pOb) const
      {
        return dynamic_cast<OBMol*>(pOb) != NULL;
      }
      virtual bool Do(OBBase* pOb, const char* OptionText, OpMap* pmap, OBConversion*);
      
      void DisplayConfig(OBConversion* pConv);
      bool Run(OBConversion* pConv, OBMol* pmol);
      double rmsd_cutoff;
      double energy_cutoff;
      bool verbose;
      ALIGNING ssalign;
      short int prec;
      string ffname;
      OBForceField *pff;
      OBBitVec hydrogen_mask;
      string str_hydrogen_mask;
  };

  //////////////////////////////////////////////////////////
  OpSimScreen theSimScreen("simscreen"); //Global instance
  //////////////////////////////////////////////////////////

  void OpSimScreen::DisplayConfig(OBConversion* pConv)
  {
    cout << "..Input format  = " << pConv->GetInFormat()->GetID() << endl;
    cout << "..Output format = " << pConv->GetOutFormat()->GetID() << endl;
    if (rmsd_cutoff >= 0.0){
        cout << "..RMSD cutoff   = " << rmsd_cutoff << endl;
    }else{
        cout << "..RMSD will be ignored for screening" << endl;
    }
    if (energy_cutoff >= 0.0){
        cout << "..Energy cutoff = " << energy_cutoff << endl;
        cout << "..Lowest energy will be determined" << endl;
    }
    else{
        cout << "..Energy will be ignored for screening" << endl;
    }
    if (prec>=0){
        cout << "..Check for symmetry will use " << prec << " decimal places" << endl;
    }else{
        cout << "..Not considering symmetry for screening" << endl;
    }
    switch (ssalign){
        case never:
            cout << "..Not aligning molecules prior to or after screening" << endl;
            break;
        case before:
            cout << "..Aligning molecules prior to screening" << endl;
            break;
        case after:
            cout << "..Aligning molecules after screening" << endl;
            break;
        case always:
            cout << "..Aligning molecules prior to and after screening" << endl;
            break;
    }
    if (!str_hydrogen_mask.empty())
        cout << "..Considering important hydrogens for screening: " << str_hydrogen_mask << endl;
    else
        cout << "..Ignoring all hydrogens for screening." << endl;
    cout << "..forcefield    = " << ffname << endl;
    cout << "..Verbose? " << (verbose ? "True" : "False") << endl;
  }

  bool OpSimScreen::Do(OBBase* pOb, const char* OptionText, OpMap* pmap, OBConversion* pConv=NULL)
  {
    OBMol* pmol = dynamic_cast<OBMol*>(pOb);

    if(pConv && pConv->IsFirstInput())
    {
      rmsd_cutoff = -1.0;
      energy_cutoff = -1.0;
      ffname = "mmff94";
      prec = -1;
      verbose = false;
      ssalign = never;
      str_hydrogen_mask = "";
      hydrogen_mask = OBBitVec();
      hydrogen_mask.Resize(pmol->NumAtoms());
      hydrogen_mask.SetRangeOff(0,pmol->NumAtoms()-1);

      OpMap::const_iterator iter;
      iter = pmap->find("rcutoff");
      if(iter!=pmap->end())
        rmsd_cutoff = atof(iter->second.c_str());
      iter = pmap->find("ecutoff");
      if(iter!=pmap->end())
        energy_cutoff = atof(iter->second.c_str());
      iter = pmap->find("ffname");
      if(iter!=pmap->end())
        ffname = iter->second;
      iter = pmap->find("prec");
      if(iter!=pmap->end())
        prec = atoi(iter->second.c_str());
      iter = pmap->find("ssalign");
      if(iter!=pmap->end()){
        for (std::string::const_iterator it = iter->second.begin(); it!=iter->second.end(); ++it){
            switch (*it){
                case 'b':
                    if       (ssalign == never){
                        ssalign = before;
                    }else{if (ssalign == after){
                        ssalign = always;
                    }}
                    break;
                case 'a':
                    if       (ssalign == never){
                        ssalign = after;
                    }else{if (ssalign == before){
                        ssalign = always;
                    }}
                    break;
                default:
                    cerr << "Accepted values for '--ssalign' are 'b' and 'a' and arbitrary combinations. Aborting." << endl;
                    return false;
                    break;
            }
        }
      }
      iter = pmap->find("verbose");
      if(iter!=pmap->end())
        verbose = true;
      iter = pmap->find("imp-H");
      if(iter!=pmap->end()){
          str_hydrogen_mask = iter->second;
          stringstream ss(str_hydrogen_mask);
          int i;
          while(ss >> i){
              if (pmol->GetAtom(i)->GetAtomicNum() != 1){
                cerr << "ERROR: atom with index " << i << " is no hydrogen." << endl;
                return false;
              }
              hydrogen_mask.SetBitOn(i-1);
              if (ss.peek() == ',' || ss.peek() == EOF)
                  ss.ignore();
              else{
                cerr << "ERROR: argument to --imp-H must be comma-separated integers, but found " << ss.peek() << endl;
                return false;
              }
          }
      }

      cout << "**Starting SimScreen" << endl;
      pff = OpenBabel::OBForceField::FindType(ffname.c_str());
      if (!pff) {
        cout << "!Cannot find forcefield " << ffname << "!" << endl;
        exit(-1);
      }
      DisplayConfig(pConv);
    }

    return Run(pConv, pmol);
  }

  bool OpSimScreen::Run(OBConversion* pConv, OBMol* pmol)
  {
    
    cout << endl << "..conformers before screening: " << pmol->NumConformers() << endl;

    if (ssalign == before || ssalign == always){
      if (verbose)
        cout << "..aligning conformers prior to screening" << endl;
      for (int i=pmol->NumConformers()-1; i>=0; --i){
          pmol->SetConformer(i);
          pmol->Align({0.0,0.0,0.0},{1.0,0.0,0.0},{0.0,1.0,0.0});
      }
    }
    //force reinitializing the ff's internal OBMol and take over this molecule's conformers
    bool success = pff->Setup(*pmol,true,true);
    if (!success) {
      cout << "!!Cannot set up forcefield for this molecule" << endl;
      return false;
    }

    success = (pff->ScreenByRMSD(rmsd_cutoff, energy_cutoff, prec, hydrogen_mask, verbose) == 0);
    if (!success) {
      cout << "!!Error in screening procedure for this molecule" << endl;
      return false;
    }

    success = pff->GetConformers(*pmol);
    if (!success) {
      cout << "!!Error updating conformers for this molecule after screening" << endl;
      return false;
    }

    cout << "..conformers after screening: " << pmol->NumConformers() << endl;

    if (ssalign == after || ssalign == always){
      if (verbose)
        cout << "..aligning conformers after screening" << endl;
      for (int i=pmol->NumConformers()-1; i>=0; --i){
          pmol->SetConformer(i);
          pmol->Align({0.0,0.0,0.0},{1.0,0.0,0.0},{0.0,1.0,0.0});
      }
    }

    OBFormat* multixyz = pConv->FindFormat("multixyz");
    if (multixyz!=NULL && multixyz==pConv->GetOutFormat()){
        pConv->GetOutFormat()->WriteMolecule(pmol, pConv);
    }else{
        for (unsigned int c = 0; c < pmol->NumConformers(); ++c) {
          pmol->SetConformer(c);
          if(!pConv->GetOutFormat()->WriteMolecule(pmol, pConv))
            break;
        }
    }
    if (pmol->NumConformers()>0) {
        pmol->SetConformer(0);
    }
    cout << endl;

    //Never have OpenBabel output anything after this. However, the boolean
    //that is returned can still be used to check whether the screening
    //finished successfully or not.
    pConv->SetOutFormat("nul");

    return true;

  }

} //namespace
