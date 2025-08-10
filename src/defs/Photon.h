#ifndef Analysis_Core_Photon_h
#define Analysis_Core_Photon_h

#include "Object.h"

namespace analysis
{
	namespace core
	{
		class Photon : public Object 
		{
			public:
				Photon() : Object() 
                {this->reset();}
                virtual ~Photon() {}
				virtual void reset()
				{
                    _pt = 0;
                    _eta = 0;
                    _phi = 0;
					_hoe = 0;
                    _isScEtaEB = 0;
                    _isScEtaEE = 0;
                    _mvaID = 0;
                    _pixelSeed = 0;
				}

				float _pt;
                float _eta;
                float _phi;
                float _hoe;
                float _isScEtaEB;
                bool _isScEtaEE;
                float _mvaID;
                bool _pixelSeed; 


		};

		typedef std::vector<analysis::core::Photon> Photons;
	}
}

#endif
