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
#include<openbabel/op.h>
#include<openbabel/mol.h>
#include<openbabel/forcefield.h>
#include <openbabel/obconversion.h>
#include<openbabel/generic.h>

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
      OpSimScreen(const char* ID) : OBOp(ID, false) {
      }

      const char* Description()
      {
        return "SimScreen, screen given molecules with respect to their RMSD and energy.\n"
          "Typical usage: obabel infile.xxx -O outfile.yyy --simscreen --rcutoff 5 --ecutoff 50\n"
          "  options:\n"
          "    --rcutoff #  RMSD cutoff (default 0.5 Angstrom)\n"
          "    --ecutoff #  Energy cutoff (default -100 kcal/mol)\n"
          "                 Negative values switch off energy screening\n"
          "    --emin #     Manually declare the global minimum energy \n"
          "    --ffname #   Name of the force field to use (default mmff94)\n"
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
      double lowest_energy;
      bool lowest_energy_given;
      bool verbose;
      string ffname;
      OBForceField *pff;
  };

  //////////////////////////////////////////////////////////
  OpSimScreen theSimScreen("simscreen"); //Global instance
  //////////////////////////////////////////////////////////

  void OpSimScreen::DisplayConfig(OBConversion* pConv)
  {
    cout << "..Input format  = " << pConv->GetInFormat()->GetID() << endl;
    cout << "..Output format = " << pConv->GetOutFormat()->GetID() << endl;
    cout << "..RMSD cutoff   = " << rmsd_cutoff << endl;
    if (energy_cutoff >= 0.0){
        cout << "..Energy cutoff = " << energy_cutoff << endl;
        if (lowest_energy_given){
            cout << "..Lowest energy = " << lowest_energy << endl;
        }
        else{
            cout << "..Lowest energy will be determined" << endl;
        }
    }
    else{
        cout << "..Energy will be ignored for screening" << endl;
    }
    cout << "..forcefield    = " << ffname << endl;
    cout << "..Verbose? " << (verbose ? "True" : "False") << endl;
  }

  bool OpSimScreen::Do(OBBase* pOb, const char* OptionText, OpMap* pmap, OBConversion* pConv=NULL)
  {
    OBMol* pmol = dynamic_cast<OBMol*>(pOb);

    if(pConv && pConv->IsFirstInput())
    {
      pConv->AddOption("writeconformers", OBConversion::GENOPTIONS);
      rmsd_cutoff = 0.5;
      energy_cutoff = -100.0;
      ffname = "mmff94";
      lowest_energy_given = false;
      lowest_energy = 0.0;
      verbose = false;

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
      iter = pmap->find("emin");
      if(iter!=pmap->end()){
        lowest_energy_given = true;
        lowest_energy = atof(iter->second.c_str());
      }
      iter = pmap->find("verbose");
      if(iter!=pmap->end())
        verbose = true;

      cout << "**Starting SimScreen" << endl;
      pff = OpenBabel::OBForceField::FindType(ffname.c_str());
      if (!pff) {
        cout << "!Cannot find forcefield " << ffname << "!" << endl;
        exit(-1);
      }
      DisplayConfig(pConv);
    }

    Run(pConv, pmol);

    return false;
  }

  void OpSimScreen::Run(OBConversion* pConv, OBMol* pmol)
  {
    //OBMol mol = *pmol;
    
    //mol.AddHydrogens(); //this should not be necessary for the screening
    bool success = pff->Setup(*pmol);
    if (!success) {
      cout << "!!Cannot set up forcefield for this molecule\n"
           << "!!Skipping\n" << endl;
      return;
    }
    cout << endl << "..conformers before screening: " << pmol->NumConformers() << endl;

    pff->ScreenByRMSD(rmsd_cutoff, energy_cutoff, lowest_energy, lowest_energy_given, verbose);

    pff->GetConformers(*pmol);

    cout << "..conformers after screening: " << pmol->NumConformers() << endl;

    for (unsigned int c = 0; c < pmol->NumConformers(); ++c) {
      pmol->SetConformer(c);
      if(!pConv->GetOutFormat()->WriteMolecule(pmol, pConv))
        break;
    }
    if (pmol->NumConformers()>0) {
        pmol->SetConformer(0);
    }
    cout << endl;

  }

} //namespace
