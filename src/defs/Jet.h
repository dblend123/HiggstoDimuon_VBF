#ifndef Analysis_Core_Jet_h
#define Analysis_Core_Jet_h

#include "GenJet.h"

namespace analysis
{
	namespace core
	{
		class Jet : public Object
		{
			public:
				Jet() : Object() 
                {this->reset();}
				virtual ~Jet() {}
				virtual void reset()
				{
					_px = 0;
					_py = 0;
					_pz = 0;
					_pt = 0;
					_eta = 0;
					_phi = 0;
					_mass = 0;
					_charge = 0;
					_partonFlavour = 0;
					_chf = 0;
					_nhf = 0;
					_cef = 0;
					_nef = 0;
					_muf = 0;
					_hfhf = 0;
					_hfef = 0;
					_cm = 0;
					_chm = 0;
					_nhm = 0;
					_cem = 0;
					_nem = 0;
					_mum = 0;
					_hfhm = 0;
					_hfem = 0;
					_jecf = 0;
					_jecu = 0;
					
                    _btag.clear();
					
                    _puid = 0;

                    _uncAK5 = 0;
                    _uncAK4 = 0;
                    _pt_upAK5 = 0;
                    _pt_upAK4 = 0;
                    _pt_downAK5 = 0;
                    _pt_downAK4 = 0;

					_genjet.reset();
					
                    _genMatched = 0;
					_genemf = 0;
					_genhadf = 0;
					_geninvf = 0;
					_genauxf = 0;
				    
                    _jetId = 0;
                    _btagDeepFlavB = 0;
                    _area = 0;
                    
                    _chHEF = 0;
                    _neHEF = 0;
                    _chEmEF = 0;
                    _neEmEF = 0;
                    _muEF = 0;
                }

				float _px;
				float _py;
				float _pz;
				float _pt;
				float _eta;
				float _phi;
				float _mass;
				float _charge;
				float _partonFlavour;

				float _chf;
				float _nhf;
				float _cef;
				float _nef;
				float _muf;
				float _hfhf;
				float _hfef;
				float _cm;
				float _chm;
				float _nhm;
				float _cem;
				float _nem;
				float _mum;
				float _hfhm;
				float _hfem;
				float _jecf;
				float _jecu;
                std::vector<float> _btag;
				float _puid;

                double _uncAK5;
                double _uncAK4;
                double _pt_upAK5;
                double _pt_upAK4;
                double _pt_downAK5;
                double _pt_downAK4;

				//	
				GenJet _genjet;
				bool _genMatched;
				float _genemf;
				float _genhadf;
				float _geninvf;
				float _genauxf;
                
                unsigned char  _jetId;
                float _btagDeepFlavB;
                float _area;

                float _chHEF;
                float _neHEF;
                float _chEmEF;
                float _neEmEF;
                float _muEF;
		};

		typedef std::vector<analysis::core::Jet> Jets;
	}
}

#endif
