#include <openbabel/babelconfig.h>
#include <openbabel/op.h>
#include <openbabel/obconversion.h>
#include <vector>

namespace OpenBabel
{
/**
DeferredFormat class is intended to assist writing ops that influence the
conversion of multiple molecules with the OBConversion Convert interface.
See, for instance, OpSort. Although it is a format, it does not registered
itself and is used in a different way from normal formats. An op makes an
instance of DeferredFormat, probably when it is first called in its Do()
function.
\code
  if(pConv && pConv->IsFirstInput())
    new DeferredFormat(pConv, this); //it will delete itself
\endcode
Output objects (probably molecules) will then all be diverted to the
DeferredFormat instance and pointers to them stored there. After the last
object, it calls the op's ProcessVec(std::vector<OBBase*>& vec) function.
The objects can be manipulated or deleted (call delete with their pointer).
When the function returns, the remaining molecules in the vector will be
output to the normal output format. No conversion options are applied at
this stage, since they already have been earlier.
**/
class DeferredFormat : public OBFormat
{
public:
  DeferredFormat(OBConversion* pConv, OBOp* pOp=NULL)
  {
    _pRealOutFormat = pConv->GetOutFormat();
    pConv->SetOutFormat(this);
    _pOp = pOp;
  }
  virtual const char* Description() { return "Read and write an OBBase* array"; }
  virtual bool ReadChemObject(OBConversion* pConv)
  {
    if(_obvec.empty())
    {
      delete this;//self destruction; was made in new in an OBOp
      return false;
    }
    //returned in reverse order, because its easier with a vector
    pConv->AddChemObject(_obvec.back());
    _obvec.pop_back();
    return true;
  }

  virtual bool WriteChemObject(OBConversion* pConv)
  {
    //Store the object pointer.
    //Unlike most formats, no deletion of object here or object constuction in ReadChemObject.
    _obvec.push_back(pConv->GetChemObject());
    if(pConv->IsLast())
    {
      //At the end, sort, or whatever, the vector
      if(_pOp)
      {
        if(_pOp->ProcessVec(_obvec))
        {
          //Now output the processed vector
          std::reverse(_obvec.begin(),_obvec.end()); //because DeferredFormat outputs in reverse order
          pConv->SetInAndOutFormats(this, _pRealOutFormat);

          std::ifstream ifs; // get rid of gcc warning
          pConv->SetInStream(&ifs);//Not used, but Convert checks it is ok
          pConv->GetInStream()->clear();

          //clear the options - they have already been applied
          pConv->SetOptions("",OBConversion::GENOPTIONS);
          pConv->SetOutputIndex(0);
          pConv->Convert();
        }
      }
    }
    return true;
  }
private:
  OBFormat* _pRealOutFormat;
  std::vector<OBBase*> _obvec;
  public:
  OBOp* _pOp;
};

} //namespace

