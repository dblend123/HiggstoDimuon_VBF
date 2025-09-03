#include <boost/program_options/options_description.hpp>
#include <boost/program_options/variables_map.hpp>
#include <iostream>
#include <string>

#include "TChain.h"

#include "defs/Macros.h"

#include "boost/program_options.hpp"

namespace bpo = boost::program_options;

using namespace std::string_literals;

#define PRINT_VAR(NAME) std::cout << #NAME "=" << NAME << std::endl

int main(int argc, char** argv) {
    std::string root_file_path;

    bpo::options_description desc{"options"};
    desc.add_options()
        ("help", "print help message")
        ("root-file", bpo::value<std::string>(&root_file_path), "input root file")
        ;
    bpo::variables_map vm;
    bpo::store(bpo::parse_command_line(argc, argv, desc), vm);
    bpo::notify(vm);

    if (vm.count("help")) {
        std::cout << desc << std::endl;
        return 0;
    }

    auto root_chain = new TChain{"Events"};
    root_chain->Add(root_file_path.c_str());
    std::cout << "nevents = "<< root_chain->GetEntries() << std::endl;

#define ROOT_CHAIN root_chain

    SET_BRANCH_UINT(nMuon);
    SET_BRANCH_FLOAT_ARRAY(Muon_pt);
    SET_BRANCH_FLOAT_ARRAY(Muon_phi);
    SET_BRANCH_FLOAT_ARRAY(Muon_eta);

    for (size_t i=0; i<root_chain->GetEntries(); i++) {
        root_chain->GetEntry(i);
        if (nMuon > 1) {
            for (size_t imu=0; imu<nMuon; imu++) {
                std::cout << "Muon_pt[" << imu << "] = " << Muon_pt[imu] << " ";
            }
            for (size_t imu=0; imu<nMuon; imu++) {
                std::cout << "Muon_phi[" << imu << "] = " << Muon_phi[imu] << " ";
            }
            for (size_t imu=0; imu<nMuon; imu++) {
                std::cout << "Muon_eta[" << imu << "] = " << Muon_eta[imu] << " ";
            }
            std::cout << std::endl;
        }
    }
    
    return 0;
}
