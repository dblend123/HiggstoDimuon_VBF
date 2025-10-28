#ifndef Analysis_Core_Flag_h
#define Analysis_Core_Flag_h

#include "Object.h"

namespace analysis
{
	namespace core
	{
		class Flag : public Object
		{
			public:
				Flag() : Object() 
                {this->reset();}
				virtual ~Flag() {}
				virtual void reset()
				{
				    _goodVertices                       = 0;
                    _globalSuperTightHalo2016Filter     = 0;
                    _EcalDeadCellTriggerPrimitiveFilter = 0;
                    _BadPFMuonFilter                    = 0;
                    _BadPFMuonDzFilter                  = 0;
                    _eeBadScFilter                      = 0;
                    _ecalBadCalibFilter                 = 0;
                    _hfNoisyHitsFilter                  = 0;
                }

                bool _goodVertices;
                bool _globalSuperTightHalo2016Filter;
                bool _EcalDeadCellTriggerPrimitiveFilter;
                bool _BadPFMuonFilter;
                bool _BadPFMuonDzFilter;
                bool _eeBadScFilter;
                bool _ecalBadCalibFilter;
                bool _hfNoisyHitsFilter;
		};

		typedef std::vector<analysis::core::Flag> Flags;
	}
}

#endif
