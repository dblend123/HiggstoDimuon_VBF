#ifndef Analysis_Core_GenJet_h
#define Analysis_Core_GenJet_h

#include "Object.h"

namespace analysis
{
	namespace core
	{
		class GenJet : public Object
		{
			public:
				GenJet() : Object()
				{this->reset();}
				virtual ~GenJet() {}
				virtual void reset()
				{
					_pt = 0;
					_eta = 0;
					_phi = 0;
					_px = 0;
					_py = 0;
					_pz = 0;
					_mass = 0;
				}

				float _px;
				float _py;
				float _pz;
				float _pt;
				float _eta;
				float _phi;
				float _mass;

		};

		typedef std::vector<analysis::core::GenJet> GenJets;
	}
}

#endif
