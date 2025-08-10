#ifndef Analysis_Core_L1_h
#define Analysis_Core_L1_h

#include "Object.h"

namespace analysis
{
    namespace core
    {
        class L1 : public Object
        {
            public:
                L1() : Object()
                    {this->reset();}
                    virtual ~L1() {}
                    virtual void reset()
                    {
                        _DoubleMu_12_5 = 0;
                        _SingleMu12er1p5 = 0;
                    }

                    bool _DoubleMu_12_5;
                    bool _SingleMu12er1p5;
        };
        
		typedef std::vector<analysis::core::L1> L1s;
    }
}

#endif
