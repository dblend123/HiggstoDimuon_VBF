#ifndef Analysis_Core_TrigObj_h
#define Analysis_Core_TrigObj_h

#include <cstdint>
#include "Object.h"
#include "defs/Track.h"

namespace analysis
{
    namespace core
    {
        class TrigObj : public Object
        {
            public:
                TrigObj() : Object()
                    {this->reset();}
                    virtual ~TrigObj() {}
                    virtual void reset()
                    {
                        _eta = 0.f;
                        _pt = 0.f;
                        _phi = 0.f;
                        _id = 0;
                        _filterBits = 0;
                    }

                    float       _eta;
                    float       _pt;
                    float       _phi;
                    uint16_t    _id;
                    uint16_t    _filterBits;
        };
        
		typedef std::vector<analysis::core::TrigObj> TrigObjs;
    }
}

#endif
