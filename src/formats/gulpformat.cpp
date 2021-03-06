/**********************************************************************
Copyright (C) 2004 by Chris Morley for template
Copyright (C) 2009 by David C. Lonie for GULP

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/

#include <openbabel/babelconfig.h>
#include <openbabel/obmolecformat.h>

#define EV_TO_KCAL_PER_MOL 23.060538

using namespace std;
namespace OpenBabel {
  class GULPFormat : public OBMoleculeFormat
  {
  public:

    GULPFormat()
    {
      OBConversion::RegisterFormat("got",this);
    }

    virtual const char* Description()
    {
      return "GULP format\n\n";
    };

    virtual const char* SpecificationURL(){return "http://cms.mpi.univie.ac.at/vasp/vasp/vasp.html";};

    /* Flags() can return be any of the following combined by |
       or be omitted if none apply
       NOTREADABLE  READONEONLY  NOTWRITABLE  WRITEONEONLY  DEFAULTFORMAT
       READBINARY  WRITEBINARY  READXML  ZEROATOMSOK */
    virtual unsigned int Flags()
    {
      return READONEONLY | NOTWRITABLE;
    };

    virtual int SkipObjects(int n, OBConversion* pConv)
    {
      return 0;
    };

    ////////////////////////////////////////////////////
    /// Declarations for the "API" interface functions. Definitions are below
    virtual bool ReadMolecule(OBBase* pOb, OBConversion* pConv);
    //    virtual bool WriteMolecule(OBBase* pOb, OBConversion* pConv);

  private:
    /* Add declarations for any local function or member variables used.
       Generally only a single instance of a format class is used. Keep this in
       mind if you employ member variables. */
  };
  ////////////////////////////////////////////////////

  //Make an instance of the format class
  GULPFormat theGULPFormat;

  /////////////////////////////////////////////////////////////////

  bool GULPFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv)
  {
    OBMol* pmol = pOb->CastAndClear<OBMol>();
    if(pmol==NULL)
      return false;

    //Define some references so we can use the old parameter names
    istream &ifs = *pConv->GetInStream();
    OBMol &mol = *pmol;
    const char* title = pConv->GetTitle();

    char buffer[BUFF_SIZE], tag[BUFF_SIZE];
    double x,y,z,a,b,c,alpha,beta,gamma;
    vector<string> vs;
    matrix3x3 ortho;
    int atomicNum;
    OBUnitCell *cell = new OBUnitCell();
    pmol->BeginModify();

    while (ifs.getline(buffer,BUFF_SIZE)) {
      // Unit cell info
      if (strstr(buffer, "Final cell parameters and derivatives :")) {
        ifs.getline(buffer,BUFF_SIZE); // Blank
        ifs.getline(buffer,BUFF_SIZE); // -----

        ifs.getline(buffer,BUFF_SIZE); // a
        tokenize(vs, buffer);
        a = atof(vs.at(1).c_str());

        ifs.getline(buffer,BUFF_SIZE); // b
        tokenize(vs, buffer);
        b = atof(vs.at(1).c_str());

        ifs.getline(buffer,BUFF_SIZE); // c
        tokenize(vs, buffer);
        c = atof(vs.at(1).c_str());

        ifs.getline(buffer,BUFF_SIZE); // alpha
        tokenize(vs, buffer);
        alpha = atof(vs.at(1).c_str());

        ifs.getline(buffer,BUFF_SIZE); // beta
        tokenize(vs, buffer);
        beta = atof(vs.at(1).c_str());

        ifs.getline(buffer,BUFF_SIZE); // gamma
        tokenize(vs, buffer);
        gamma = atof(vs.at(1).c_str());

        // Build unit cell
        cell->SetData(a, b, c, alpha, beta, gamma);
      }

      // Unit cell info
      if (strstr(buffer, "Cell parameters (Angstroms/Degrees):")) {
        ifs.getline(buffer,BUFF_SIZE); // Blank

        ifs.getline(buffer,BUFF_SIZE); // a
        tokenize(vs, buffer);
        a = atof(vs.at(2).c_str());
        alpha = atof(vs.at(5).c_str());

        ifs.getline(buffer,BUFF_SIZE); // b
        tokenize(vs, buffer);
        b = atof(vs.at(2).c_str());
        beta = atof(vs.at(5).c_str());

        ifs.getline(buffer,BUFF_SIZE); // c
        tokenize(vs, buffer);
        c = atof(vs.at(2).c_str());
        gamma = atof(vs.at(5).c_str());

        // Build unit cell
        cell->SetData(a, b, c, alpha, beta, gamma);
      }

      // Atoms info
      if (strstr(buffer, "Final fractional coordinates of atoms :")) {
        // Clear old atoms from pmol
        vector<OBAtom*> toDelete;
        FOR_ATOMS_OF_MOL(a, *pmol)
          toDelete.push_back(&*a);
        for (int i = 0; i < toDelete.size(); i++)
          pmol->DeleteAtom(toDelete.at(i));

        // Load new atoms from molecule
        ifs.getline(buffer,BUFF_SIZE); // Blank
        ifs.getline(buffer,BUFF_SIZE); // Header
        ifs.getline(buffer,BUFF_SIZE); // Header
        ifs.getline(buffer,BUFF_SIZE); // Header
        ifs.getline(buffer,BUFF_SIZE); // Header

        ifs.getline(buffer,BUFF_SIZE); // First entry
        tokenize(vs, buffer);
        int size = vs.size();
        while (size >= 7 && size <= 10) {
          atomicNum = etab.GetAtomicNum(vs[1].c_str());

          // Gulp sometimes places extra chars between the coords, so
          // it's not so straight-forward to parse them...
          x = y = z = 0;
          int set = 0;
          for (unsigned i = 3; i < size; i++) {
            if (strstr(vs[i].c_str(), "*")) continue; // Skip "*" in output
            // Else assign x,y,z based on how many coords have been
            // set already. These are currently fractional, we'll
            // convert all at the end of the run.
            switch (set) {
            case 0:
              x = atof((char*)vs[i].c_str());
              set++;
              break;
            case 1:
              y = atof((char*)vs[i].c_str());
              set++;
              break;
            case 2:
              z = atof((char*)vs[i].c_str());
              set++;
              break;
            default:
              break;
            }
          }
          // Add atom
          OBAtom *atom = pmol->NewAtom();
          atom->SetAtomicNum(atomicNum);
          vector3 coords (x,y,z);
          atom->SetVector(coords);

          // Reset vars
          ifs.getline(buffer,BUFF_SIZE); // First entry
          tokenize(vs, buffer);
          size = vs.size();
        }
      }

      // Free energy
      if (strstr(buffer, "Final energy")) {
        tokenize(vs, buffer);
        pmol->SetEnergy(atof(vs[3].c_str()) * EV_TO_KCAL_PER_MOL);
      }

      // Enthalphy
      if (strstr(buffer, "Components of enthalpy :")) {
        bool hasPV = false;
        float en, en_eV, pv, pv_eV;
        OBPairData *enthalpy = new OBPairData();
        OBPairData *enthalpy_pv = new OBPairData();
        OBPairData *enthalpy_eV = new OBPairData();
        OBPairData *enthalpy_pv_eV = new OBPairData();
        enthalpy->SetAttribute("Enthalpy (kcal/mol)");
        enthalpy_pv->SetAttribute("Enthalpy PV term (kcal/mol)");
        enthalpy_eV->SetAttribute("Enthalpy (eV)");
        enthalpy_pv_eV->SetAttribute("Enthalpy PV term (eV)");

        ifs.getline(buffer,BUFF_SIZE);

        while (strstr(buffer, "kJ/(mole unit cells)") == 0) {
          if (strstr(buffer, "Pressure*volume")) {
            tokenize(vs, buffer);
            pv_eV = static_cast<float> (atof(vs[2].c_str()));
            pv = static_cast<float> (pv_eV * EV_TO_KCAL_PER_MOL);
            snprintf(tag, BUFF_SIZE, "%f", pv);
            enthalpy_pv->SetValue(tag);
            snprintf(tag, BUFF_SIZE, "%f", pv_eV);
            enthalpy_pv_eV->SetValue(tag);
            pmol->SetData(enthalpy_pv);
            pmol->SetData(enthalpy_pv_eV);
            hasPV = true;
          }

          if (strstr(buffer, "Total lattice enthalpy")) {
            tokenize(vs, buffer);
            en_eV = static_cast<float> (atof(vs[4].c_str()));
            en = static_cast<float> (en_eV * EV_TO_KCAL_PER_MOL);
            snprintf(tag, BUFF_SIZE, "%f", en);
            enthalpy->SetValue(tag);
            snprintf(tag, BUFF_SIZE, "%f", en_eV);
            enthalpy_eV->SetValue(tag);
            pmol->SetData(enthalpy);
            pmol->SetData(enthalpy_eV);
          }
          
          ifs.getline(buffer,BUFF_SIZE);
        }
        if (hasPV)
          pmol->SetEnergy(en - pv);
      }
    }

    // Convert coords to cartesian
    FOR_ATOMS_OF_MOL(atom, pmol) {
      atom->SetVector(cell->FractionalToCartesian(atom->GetVector()));
    }

    // set final unit cell
    pmol->SetData(cell);

    pmol->EndModify();

    return true;
  }

} //namespace OpenBabel
