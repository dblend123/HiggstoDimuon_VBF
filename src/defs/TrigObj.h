#ifndef Analysis_Core_TrigObj_h
#define Analysis_Core_TrigObj_h

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
                        _eta = 0;
                        _pt = 0;
                        _phi = 0;
                        _id = 0;
                        _filterBits = 0;
                    }

                    float _eta;
                    float _pt;
                    float _phi;
                    short _id;
                    int _filterBits;
        };
        
		typedef std::vector<analysis::core::TrigObj> TrigObjs;
    }
}

#endif
