#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "surge/sd3d/core.hpp"

namespace surge::brady1984 {

struct RunConfig {
  std::string name = "A";
  int particles = 25;
  double areaFraction = 0.4;
  double gammaStar = 0.5;
  double tau = 227.0;
  bool repulsion = true;
  double dt = 4e-3;
  int steps = 20000;
  int warmupSteps = 2000;
  std::uint64_t seed = 5141;
  int maxPrintedM = 11;
};

struct HistogramOptions {
  int angularBins = 10;
  double nearMinR = 2.0;
  double nearMaxR = 2.03;
  double radialMaxR = 7.0;
  double radialDr = 0.01;
  double yBinWidth = 0.05;
};

struct Particle2 {
  double x = 0.0;
  double y = 0.0;
};

struct MonolayerState {
  double L = 1.0;
  double H = 1.0;
  std::vector<Particle2> particles;
};

struct BradyStepResult {
  surge::sd3d::Vector qdot;      // [Ux*, Uy*, Omega_z*] per particle
  surge::sd3d::Vector force;     // [Fx, Fy, Tz] nonhydrodynamic force/torque
};

struct BradyStatistics {
  int samples = 0;
  std::vector<double> angularCounts;
  std::vector<double> radialCounts;
  std::vector<double> yCounts;
};

RunConfig presetRun(char runName);

MonolayerState makeInitialState(const RunConfig& config);

BradyStepResult solveProjectedVelocities(
    const MonolayerState& state,
    double time,
    const RunConfig& config);

void rk4Step(MonolayerState& state, double time, const RunConfig& config);

void runSimulation(
    const RunConfig& config,
    const HistogramOptions& histOptions,
    const std::string& trajectoryCsv,
    const std::string& angularCsv,
    const std::string& radialCsv,
    const std::string& yCsv);

}  // namespace surge::brady1984

