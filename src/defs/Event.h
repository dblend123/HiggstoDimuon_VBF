#ifndef Analysis_Core_Event_h
#define Analysis_Core_Event_h

#include "Object.h"
#include <vector>

namespace analysis
{
	namespace core
	{
		class Event : public Object
		{
			public:
				Event() : Object() {this->reset();}

        virtual void reset()
        {
            _run              = 0;
            _lumi             = 0;
            _event            = 0;
            _bx               = 0;
            _orbit            = 0;

            _genWeight        = 1.0f;
            _LHEWeight        = 1.0f;

            _nLHEScaleWeight  = 0;
            _nLHEPdfWeight    = 0;
            _nPSWeight        = 0;

            _LHEScaleWeight.clear();
            _LHEPdfWeight.clear();
            _PSWeight.clear();
        }

        virtual ~Event() {}

        int _run;
        int _lumi;
        long long int _event;
        int _bx;
        int _orbit;

        float _genWeight;
        float _LHEWeight;

        int _nLHEScaleWeight;
        int _nLHEPdfWeight;
        int _nPSWeight;

        std::vector<float> _LHEScaleWeight;
        std::vector<float> _LHEPdfWeight;
        std::vector<float> _PSWeight;
};
		class EventAuxiliary : public Object
		{
			public:
				EventAuxiliary() :  Object() {this->reset();}
				virtual ~EventAuxiliary() {}

				virtual void reset()
				{
					_nPU = 0;
					_genWeight = 0;
					_hasHLTFired.clear();
				}

				int _nPU;
				int _genWeight;
				std::vector<bool> _hasHLTFired;

		};

		typedef std::vector<analysis::core::EventAuxiliary> EventAuxiliaries;
		typedef std::vector<analysis::core::Event> Events;
	}
}

#endif
