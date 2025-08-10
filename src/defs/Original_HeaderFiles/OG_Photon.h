#pragma once

#include <vector>

#include "defs/Track.h"

namespace analysis
{
	namespace core
	{
		class Photon : public Track
		{
			public:
				Photon() : Track() {this->reset();}

				virtual void reset()
				{
					_hoe = 0;
                    _isScEtaEB = 0;
                    _isScEtaEE = 0;
                    _mvaID = 0;
                    _pixelSeed = 0;
                    Track::reset();
					_track.reset();

				}
				virtual ~Photon() {}

				Track _track;
				float _hoe;
                float _isScEtaEB;
                bool _isScEtaEE;
                float _mvaID;
                bool _pixelSeed; 


		};

		typedef std::vector<analysis::core::Photon> Photons;
	}
}
