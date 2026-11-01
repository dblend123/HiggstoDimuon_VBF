#pragma once

#include "TString.h"
#include "TChain.h"

#include <fstream>
#include <string>
#include <iostream>

#include "TError.h"
#include <cstdlib>

namespace analysis
{
	namespace processing
	{
		class Streamer
		{
			public:
				Streamer() : _inputname(""), _tpathname(""), _chain(NULL)
				{}
				Streamer(std::string const& inputname, 
					std::string const& tpathname):
					_inputname(inputname), _tpathname(tpathname)
				{
					_chain = new TChain(_tpathname.c_str());
				}
				virtual ~Streamer() {}

				void chainup();
				void chainup(std::string const& inputname, 
					std::string const& tpathname);

				std::string			_inputname;
				std::string			_tpathname;
				TChain				*_chain;
                Long64_t _expectedEntries = -1;
		};
	}
}

namespace analysis
{
	namespace processing
	{
		inline void Streamer::chainup(std::string const& inputname, 
			std::string const& tpathname)
		{
			_inputname = inputname;
			_tpathname = tpathname;
			_chain = new TChain(_tpathname.c_str());
			chainup();
		}

		inline void Streamer::chainup()
		{
			gErrorIgnoreLevel = kWarning;
            std::cout << "Chaining All The Data Files" << std::endl;
			std::cout << "\n";
            std::ifstream input(_inputname);
			for (std::string line; std::getline(input, line); )
			{
			
            const int added = _chain->Add(TString(line.c_str()), 0);

            if (added <= 0) {
                std::cerr
                    << "[FATAL INPUT FILE] Failed to add: "
                    << line << std::endl;
                std::exit(EXIT_FAILURE);
            }
            
            std::cout << "[CHAIN CHECK] Files registered = "
                << _chain->GetListOfFiles()->GetEntries()
                << std::endl; 
            
            }
            Long64_t TotalEntries = _chain->GetEntries();
            if (TotalEntries <= 0) {
                std::cerr << "[FATAL] Invalid input entry count: "
                          << TotalEntries << std::endl;
                std::exit(EXIT_FAILURE);
            }
            _expectedEntries = TotalEntries;
            std::cout << "Total Entries = " << TotalEntries 
				<< std::endl;
			std::cout << "\n";
		}
	}
}
