#include <filesystem>
#include <iostream>
#include <string>

#include "surge/brady1984/monolayer.hpp"

namespace {

void usage(const char* exe) {
  std::cerr
      << "usage: " << exe << " [run=A] [output_dir=results/brady_A] [steps_override]\n"
      << "supported force-superposition runs: A, D, E, F, G, H\n";
}

}  // namespace

int main(int argc, char** argv) {
  try {
    char run = 'A';
    if (argc > 1) {
      const std::string s = argv[1];
      if (s.empty()) {
        usage(argv[0]);
        return 2;
      }
      run = s[0];
    }

    std::string outDir = "results/brady_A";
    if (argc > 2) outDir = argv[2];

    surge::brady1984::RunConfig config = surge::brady1984::presetRun(run);
    if (argc > 3) {
      config.steps = std::stoi(argv[3]);
      config.warmupSteps = std::min(config.warmupSteps, config.steps / 5);
    }

    std::filesystem::create_directories(outDir);
    const std::string prefix = outDir + "/run_" + config.name;

    surge::brady1984::HistogramOptions hist;
    surge::brady1984::runSimulation(
        config,
        hist,
        prefix + "_trajectory.csv",
        prefix + "_angular.csv",
        prefix + "_radial.csv",
        prefix + "_y.csv");

    std::cout << "Brady 1984 run " << config.name << " complete\n";
    std::cout << "particles=" << config.particles
              << " phi=" << config.areaFraction
              << " gamma*=" << config.gammaStar
              << " tau=" << (config.repulsion ? std::to_string(config.tau) : std::string("no-force"))
              << " dt=" << config.dt
              << " steps=" << config.steps << "\n";
    std::cout << "wrote " << outDir << "\n";
  } catch (const std::exception& e) {
    std::cerr << "error: " << e.what() << "\n";
    usage(argv[0]);
    return 1;
  }
  return 0;
}

