#ifndef Analysis_Core_HLT_h
#define Analysis_Core_HLT_h

#include "Object.h"

namespace analysis
{
    namespace core
    {
        class HLT : public Object
        {
            public:
                HLT() : Object()
                    {this->reset();}
                    virtual ~HLT() {}
                    virtual void reset()
                    {
                        _IsoMu24 = 0;
                    }

                    bool _IsoMu24;
        };
        
		typedef std::vector<analysis::core::HLT> HLTs;
    }
}

#endif
