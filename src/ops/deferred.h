#include <openbabel/babelconfig.h>
#include <openbabel/op.h>
#include <openbabel/obconversion.h>
#include <vector>
#include <algorithm>

namespace OpenBabel
{
/**
DeferredFormat class is intended to assist writing ops that influence the
conversion of multiple molecules with the OBConversion Convert interface.
See, for instance, OpSort. Although it is a format, it does not registered
itself, an object is not constructed in ReadChemObject() or deleted in
WriteChemObject(). It is used in a different way from normal formats. An 
op makes an instance of DeferredFormat, probably when it is first called
in its Do() function.
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

Constructing with a third boolean parameter set true allows an alternative mode
of operation. Before storing the pointer to an object, DeferredFormat calls
the op's Do() function and stores the object pointer only if this returns true.
This has the effect of allowing the op to act after all the other options,
rather than before most of them. See OpLargest for an example.
**/
class DeferredFormat : public OBFormat
{
public:
  DeferredFormat(OBConversion* pConv, OBOp* pOp=NULL, bool CallDo=false)
  {
    _pRealOutFormat = pConv->GetOutFormat();
    pConv->SetOutFormat(this);
    _pOp = pOp;
    _callDo = CallDo;
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
    OBBase* pOb = pConv->GetChemObject();
    if(!_callDo || _pOp->Do(pOb, "", pConv->GetOptions(OBConversion::GENOPTIONS), pConv))  
      _obvec.push_back(pOb); // Store the object pointer.

    if(pConv->IsLast())
    {
      //At the end, sort, or whatever, the vector
      if(_pOp)
      {
        //clear the options if return is true - they have already been applied
        if(_pOp->ProcessVec(_obvec))
          pConv->SetOptions("",OBConversion::GENOPTIONS);

        //Now output the processed vector, unless it is empty
        if(!_obvec.empty())
        {
          std::reverse(_obvec.begin(),_obvec.end()); //because DeferredFormat outputs in reverse order
          pConv->SetInAndOutFormats(this, _pRealOutFormat);

          std::ifstream ifs; // get rid of gcc warning
          pConv->SetInStream(&ifs);//Not used, but Convert checks it is ok
          pConv->GetInStream()->clear();

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
  bool _callDo;
};

} //namespace

