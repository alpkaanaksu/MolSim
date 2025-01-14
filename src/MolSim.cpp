#include "io/reader/FileReader.h"
#include "utils/ArrayUtils.h"
#include "models/ParticleContainer.h"

#include <iostream>
#include <list>
#include <cmath>
#include <boost/program_options.hpp>
#include <spdlog/spdlog.h>

#include "Simulation.h"
#include "models/Model.h"

namespace po = boost::program_options;

/**
 * @brief Entry point to the program. Takes CLI arguments.
 *
 * @param argc Argument count
 * @param argv Argument vector
 *
 * Entry point to the program. Takes CLI arguments.
 *
 */
int main(int argc, char *argv[]) {
    std::string in;
    std::string out = "output";
    double end_time = 1000;
    double delta_t = 0.014;
    int video_duration = 60;
    int fps = 24;
    outputWriter::OutputType outputType = outputWriter::OutputType::VTK;
    spdlog::level::level_enum logLevel = spdlog::level::info;
    bool legacyFile = false;

    try {
        po::options_description desc("MolSim program options");

        desc.add_options()
            ("help", "Information about MolSim")
            ("legacy-file", "Use old file format instead of JSON")
            ("in", po::value<std::string>(), "Provide input file path")
            ("out", po::value<std::string>(), "Provide output file path")
            ("video-duration", po::value<int>(), "Video duration in seconds")
            ("fps", po::value<int>(), "Frames per second")
            ("end-time", po::value<double>(), "End time")
            ("delta-t", po::value<double>(), "Time delta for each iteration")
            ("output-type", po::value<std::string>(), "Output type (vtk or xyz)")
            ("log", po::value<std::string>(), "Log level (off, info, debug)")
            ("checkpoint", po::value<int>(), "Iteration for checkpoint (a checkpoint will be saved when the iteration is reached")
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

        if (vm.count("legacy-file")) {
            legacyFile = true;
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


        if (vm.count("output-type")) {
            std::string val = vm["output-type"].as<std::string>();

            outputType = outputWriter::stringToOutputType(val);
        }

        if (vm.count("log")) {
            std::string val = vm["log"].as<std::string>();
            if(val == "off")
                logLevel = spdlog::level::off;
            else if(val == "debug")
                logLevel = spdlog::level::debug;
            else if(val == "info")
                logLevel = spdlog::level::info;
            else if(val == "error")
                logLevel = spdlog::level::err;
            else if(val == "critical")
                logLevel = spdlog::level::critical;
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

    spdlog::set_level(logLevel);

    // TODO Fix code duplication
    if (legacyFile) {
        Simulation simulation{Model::gravityModel(delta_t), end_time, delta_t, video_duration, fps, in, out, outputType};

        spdlog::info("Starting simulation...");

        std::cout << simulation << std::endl;

        simulation.run();
    } else {
        Simulation simulation{in};

        spdlog::info("Starting simulation...");

        std::cout << simulation << std::endl;

        simulation.run();
    }
}