#include <openbabel/babelconfig.h>
#include <openbabel/obmolecformat.h>

#include "cansmiles.h"

using namespace std;
namespace OpenBabel
{

class Can2Format : public OBMoleculeFormat
{
  public:
    Can2Format()
    {
      OBConversion::RegisterFormat("Can2",this);
      OBConversion::RegisterOptionParam("f", this, 1);
      OBConversion::RegisterOptionParam("n", this);
      OBConversion::RegisterOptionParam("s", this, 0, OBConversion::INOPTIONS);
    }

    virtual const char* Description() //required
    {
      return "";
    }

    virtual const char* SpecificationURL()
    {
      return "";
    }

    virtual const char* GetMIMEType() 
    { return "chemical/x-xxx"; }


    virtual unsigned int Flags()
    {
      return NOTREADABLE;
    };

    virtual int SkipObjects(int n, OBConversion* pConv)
    {
      return 0;
    }

    virtual bool WriteMolecule(OBBase* pOb, OBConversion* pConv);
  private:
};  

Can2Format theCan2Format;

////////////////////////////////////////////////////////////////

bool Can2Format::WriteMolecule(OBBase* pOb, OBConversion* pConv)
{
  OBMol* pmol = dynamic_cast<OBMol*>(pOb);
  if(pmol==NULL)
      return false;

  ostream& ofs = *pConv->GetOutStream();

  return true; //or false to stop converting
}

} //namespace OpenBabel

