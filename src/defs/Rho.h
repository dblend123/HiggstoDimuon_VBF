#ifndef Analysis_Core_Rho_h
#define Analysis_Core_Rho_h

#include "Object.h"

namespace analysis
{
    namespace core
    {
        class Rho : public Object
        {
            public:
                Rho() : Object() 
                {this->reset();}
                virtual ~Rho() {}
                virtual void reset()
                {
                    _fixedGridRhoFastjetAll = 0;
                }

                float _fixedGridRhoFastjetAll;

        };

        typedef std::vector<analysis::core::Rho> Rhos;
    }
}

#endif
