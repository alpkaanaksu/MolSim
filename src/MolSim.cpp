#include "FileReader.h"
#include "outputWriter/VTKWriter.h"
#include "utils/ArrayUtils.h"
#include "models/ParticleContainer.h"

#include <iostream>
#include <list>
#include <cmath>
#include <boost/program_options.hpp>

#include "Simulation.h"
#include "models/Model.h"

namespace po = boost::program_options;

int main(int argc, char *argv[]) {
    std::string in;
    std::string out = "output";
    double end_time = 1000;
    double delta_t = 0.014;
    int video_duration = 60;
    int fps = 24;

    try {
        po::options_description desc("MolSim program options");

        desc.add_options()
            ("help", "Information about MolSim")
            ("in", po::value<std::string>(), "Provide input file path")
            ("out", po::value<std::string>(), "Provide output file path")
            ("video-duration", po::value<int>(), "Video duration in seconds")
            ("fps", po::value<int>(), "Frames per second")
            ("end-time", po::value<double>(), "End time")
            ("delta-t", po::value<double>(), "Time delta for each iteration")
        ;

        po::positional_options_description p;
        p.add("in", 1);

        po::variables_map vm;
        po::store(po::command_line_parser(argc, argv).options(desc).positional(p).run(), vm);
        po::notify(vm);

        if (vm.count("help")) {
            std::cout << desc << "\n";
            return 1;
        }

        if (vm.count("in")) {
            in = vm["in"].as<std::string>();
        }

        if (vm.count("out")) {
            out = vm["out"].as<std::string>();
        }

        if (vm.count("video-duration")) {
            video_duration = vm["video-duration"].as<int>();
        }

        if (vm.count("fps")) {
            fps = vm["fps"].as<int>();
        }

        if (vm.count("end-time")) {
            end_time = vm["end-time"].as<double>();
        }

        if (vm.count("delta-t")) {
            delta_t = vm["delta-t"].as<double>();
        }
    } catch (std::exception &e) {
        return 1;
    }

    std::cout << "Hello from MolSim for PSE!\n" << std::endl;
    if (argc < 2) {
        std::cout << "Erroneous programme call! " << '\n';
        std::cout << "./MolSim <input filename>\n" << "./MolSim --help" << std::endl;
        return 1;
    }

    Simulation simulation{Model::basicModel(delta_t), end_time, delta_t, video_duration, fps, in, out};

    std::cout << "\nStarting simulation...\n" << simulation << std::endl;

    simulation.run();
}