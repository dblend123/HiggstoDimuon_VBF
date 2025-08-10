#ifndef Analysis_Core_FsrPhoton_h
#define Analysis_Core_FsrPhoton_h

#include "Object.h"

namespace analysis
{
	namespace core
	{
		class FsrPhoton : public Object 
		{
			public:
				FsrPhoton() : Object() 
                {this->reset();}
				virtual ~FsrPhoton() {}
				virtual void reset()
				{
					_pt = 0;
					_eta = 0;
					_phi = 0;
					_muonIdx = 0;
					_relIso03 = 0;
                    _dROverEt2 = 0;
                    _electronIdx = 0;
                
				}
				
				int _muonIdx;
				float _relIso03;
				float _pt;
				float _eta;
				float _phi;
                float _dROverEt2;
                float _electronIdx;

		};

		typedef std::vector<analysis::core::FsrPhoton> FsrPhotons;
	}
}

#endif
