#ifndef Analysis_Core_Electron_h
#define Analysis_Core_Electron_h

#include "Track.h"

namespace analysis
{
    namespace core
    {
        class Electron : public Track
        {
            public:
                Electron() : Track() 
                {this->reset();}
                virtual ~Electron() {}
                virtual void reset()
                {
                    Track::reset();
                    _ids.clear();
                    _sumChargedHadronPt = 0;
                    _sumNeutralHadronEt = 0;
                    _sumPhotonEt = 0;
                    _sumPUPt = 0;
                    _sumChargedParticlePt = 0;
                    _dz = 0;
                    _isPF = false;
                    _convVeto = false;

                    _mvagp_value = 0;
                    _mvagp_category = 0;
                    _mvagpid_medium = false;
                    _mvagpid_tight = false;

                    _mvahzz_value = 0;
                    _mvahzz_category = 0;
                    _mvahzzid_loose = false;
                    
                    _cutBased = 0;
                    _mvaIso_WP90 = false;
                }

                std::vector<bool> _ids;
                float _sumChargedHadronPt;
                float _sumNeutralHadronEt;
                float _sumPhotonEt;
                float _sumPUPt;
                float _sumChargedParticlePt;
                double _dz;
                bool _isPF;
                bool _convVeto;

                float _mvagp_value;
                int _mvagp_category;
                bool _mvagpid_medium;
                bool _mvagpid_tight;

                float _mvahzz_value;
                int _mvahzz_category;
                bool _mvahzzid_loose;

                char _cutBased;
                bool _mvaIso_WP90;

        };

        typedef std::vector<analysis::core::Electron> Electrons;
    }
}

#endif
