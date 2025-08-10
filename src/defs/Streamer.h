#pragma once

#include "TString.h"
#include "TChain.h"

#include <fstream>
#include <string>
#include <iostream>

#include "TError.h"

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
			gErrorIgnoreLevel = kFatal;
            std::cout << "Chaining All The Data Files" << std::endl;
			std::cout << "\n";
            std::ifstream input(_inputname);
			for (std::string line; std::getline(input, line); )
			{
				_chain->Add(TString(line.c_str()));
			}
            int TotalEntries = _chain->GetEntries();
			std::cout << "Total Entries = " << TotalEntries 
				<< std::endl;
			std::cout << "\n";
		}
	}
}
