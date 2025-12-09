#pragma once

#include <vector>

#include "defs/Track.h"

namespace analysis
{
	namespace core
	{
		class Muon : public Track
		{
			public:
				Muon() : Track() 
                {this->reset();}
				virtual void reset()
				{
					_dxy = 0;
                    _dz = 0;

                    _isTracker = 0;
					_isStandAlone = 0;
					_isGlobal = 0;
					_tightId = 0;
					_mediumId = 0;
					_looseId = 0;

					_isPFcand = 0;

					_normChi2 = 0; 
					_d0BS = 0;
					_dzBS = 0;
					_d0PV = 0;
					_dzPV = 0;

					_nPLs = 0;
					_nTLs = 0;
					_nSLs = 0;

					_vfrTrk = 0;
					_nvMHits = 0;
					_nvPHits = 0;
					_nvTHits = 0;
					_nvSHits = 0;
					_nSegMts =0 ;
					_nMtsStations = 0;

					_trackIsoSumPt = 0;
					_trackIsoSumPtCorr = 0;

					_hIso = 0;
					_eIso = 0;
					_relCombIso = 0;
					Track::reset();
					_track.reset();
					_isHLTMatched.clear();

					_segmentCompatibility = 0;
					_combinedQChi2LocalPosition = 0;
					_combinedQTrkKink = 0;

					_sumChargedHadronPtR03 = 0;
					_sumChargedParticlePtR03 = 0;
					_sumNeutralHadronEtR03 = 0;
					_sumPhotonEtR03 = 0;
					_sumPUPtR03 = 0;
					_sumChargedHadronPtR04 = 0;
					_sumChargedParticlePtR04 = 0;
					_sumNeutralHadronEtR04 = 0;
					_sumPhotonEtR04 = 0;
					_sumPUPtR04 = 0;
				    
                    _pfRelIso04_all = 0;
                    _nTrackerLayers = 0;
                    _bsConstrainedPt = 0; 
                }
				virtual ~Muon() {}

                float _dxy;
                float _dz;

				bool _isTracker;
				Track _track;
				bool _isStandAlone;
				bool _isGlobal;
				bool _tightId;
				bool _mediumId;
				bool _looseId;

				bool _isPFcand;

				float _normChi2;
				float _d0BS;
				float _dzBS;
				float _d0PV;
				float _dzPV;

				int _nPLs;
				int _nTLs;
				int _nSLs;

				float _vfrTrk;
				int _nvMHits;
				int _nvPHits;
				int _nvTHits;
				int _nvSHits;
				int _nSegMts;
				int _nMtsStations;

				float _trackIsoSumPt;
				float _trackIsoSumPtCorr;

				float _hIso;
				float _eIso;
				float _relCombIso;

				float _segmentCompatibility;
				float _combinedQChi2LocalPosition;
				float _combinedQTrkKink;
				std::vector<bool> _isHLTMatched;

				float _sumChargedHadronPtR03;
				float _sumChargedParticlePtR03;
				float _sumNeutralHadronEtR03;
				float _sumPhotonEtR03;
				float _sumPUPtR03;
				float _sumChargedHadronPtR04;
				float _sumChargedParticlePtR04;
				float _sumNeutralHadronEtR04;
				float _sumPhotonEtR04;
				float _sumPUPtR04;

                float _pfRelIso04_all;
                char _nTrackerLayers;
                float _bsConstrainedPt;
		};

		typedef std::vector<analysis::core::Muon> Muons;
	}
}
