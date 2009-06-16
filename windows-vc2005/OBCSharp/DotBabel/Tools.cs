using System;
using System.Collections;
using System.Collections.Generic;
using System.Collections.Specialized;
using System.Linq;
using System.Text;
using OpenBabel;

namespace DotBabel
{
    static class Tools
    {
        public static BitArray AsBitArray(this OBFingerprint fp)
        {
            BitArray b = new BitArray(3);
            BitVector32 v = new BitVector32(3);
            
            return null;
        }
    }
}
