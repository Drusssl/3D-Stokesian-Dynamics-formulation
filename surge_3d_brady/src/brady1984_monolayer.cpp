#include "surge/brady1984/monolayer.hpp"

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <random>
#include <stdexcept>

#include "surge/sd3d/near_field.hpp"
#include "surge/sd3d/projection.hpp"

namespace surge::brady1984 {
namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;

struct PairImage {
  double dx = 0.0;
  double dy = 0.0;
  double r = 0.0;
};

double cellLength(int n, double phi) {
  if (n <= 0) throw std::invalid_argument("particle count must be positive");
  if (phi <= 0.0) throw std::invalid_argument("area fraction must be positive");
  return std::sqrt(static_cast<double>(n) * kPi / phi);
}

double wrapPeriodic(double x, double L) {
  x = std::fmod(x, L);
  if (x < 0.0) x += L;
  return x;
}

void wrapShearPeriodic(Particle2& p, const MonolayerState& state, double gammaT) {
  while (p.y >= state.H) {
    p.y -= state.H;
    p.x -= gammaT * state.H;
  }
  while (p.y < 0.0) {
    p.y += state.H;
    p.x += gammaT * state.H;
  }
  p.x = wrapPeriodic(p.x, state.L);
}

PairImage shearMinimumImage(
    const Particle2& from,
    const Particle2& to,
    const MonolayerState& state,
    double gammaT) {
  PairImage best;
  best.r = 1e300;
  for (int ny = -1; ny <= 1; ++ny) {
    const double imageShiftX = gammaT * static_cast<double>(ny) * state.H;
    const double dy = (to.y - from.y) + static_cast<double>(ny) * state.H;
    double dx0 = (to.x - from.x) + imageShiftX;
    dx0 -= state.L * std::round(dx0 / state.L);
    const double r = std::hypot(dx0, dy);
    if (r < best.r) best = {dx0, dy, r};
  }
  return best;
}

void addIsolatedSelf(surge::sd3d::DenseBlocks& full, int particle) {
  const double trans = 6.0 * kPi;
  const double rot = 8.0 * kPi;
  for (int c = 0; c < 3; ++c) {
    full.R(surge::sd3d::dof(particle, c, false), surge::sd3d::dof(particle, c, false)) += trans;
    full.R(surge::sd3d::dof(particle, c, true), surge::sd3d::dof(particle, c, true)) += rot;
  }
}

void subtractLocalSelf(surge::sd3d::DenseBlocks& pairBlocks, int localParticle) {
  const double trans = 6.0 * kPi;
  const double rot = 8.0 * kPi;
  for (int c = 0; c < 3; ++c) {
    const int u = surge::sd3d::dof(localParticle, c, false);
    const int w = surge::sd3d::dof(localParticle, c, true);
    pairBlocks.R(u, u) -= trans;
    pairBlocks.R(w, w) -= rot;
  }
}

surge::sd3d::DenseBlocks assembleMonolayerPairwiseBlocks(
    const MonolayerState& state,
    double time,
    const RunConfig& config) {
  const int n = static_cast<int>(state.particles.size());
  surge::sd3d::DenseBlocks full;
  full.n = n;
  full.R = surge::sd3d::Matrix::Zero(6 * n, 6 * n);
  full.Phi = surge::sd3d::Matrix::Zero(6 * n, 9);

  for (int i = 0; i < n; ++i) addIsolatedSelf(full, i);

  const double gammaT = config.gammaStar * time;
  surge::sd3d::NearFieldOptions near;
  near.maxPrintedM = config.maxPrintedM;

  for (int i = 0; i < n; ++i) {
    for (int j = i + 1; j < n; ++j) {
      const PairImage image = shearMinimumImage(state.particles[static_cast<std::size_t>(i)],
                                                state.particles[static_cast<std::size_t>(j)],
                                                state,
                                                gammaT);
      const double evalR = std::max(image.r, 2.0 + 1e-8);
      const double invR = image.r > 1e-14 ? 1.0 / image.r : 1.0;
      const double evalDx = image.r > 1e-14 ? image.dx * evalR * invR : evalR;
      const double evalDy = image.r > 1e-14 ? image.dy * evalR * invR : 0.0;
      surge::sd3d::System pairSystem;
      pairSystem.viscosity = 1.0;
      pairSystem.particles = {
          {{0.0, 0.0, 0.0}, 1.0},
          {{evalDx, evalDy, 0.0}, 1.0},
      };
      surge::sd3d::DenseBlocks pairBlocks = surge::sd3d::assemblePairwiseResistance3D(pairSystem, near);
      subtractLocalSelf(pairBlocks, 0);
      subtractLocalSelf(pairBlocks, 1);

      const int ids[2] = {i, j};
      for (int la = 0; la < 2; ++la) {
        for (int lc = 0; lc < 6; ++lc) {
          const int gr = 6 * ids[la] + lc;
          const int lr = 6 * la + lc;
          for (int lb = 0; lb < 2; ++lb) {
            for (int ld = 0; ld < 6; ++ld) {
              const int gc = 6 * ids[lb] + ld;
              const int localCol = 6 * lb + ld;
              full.R(gr, gc) += pairBlocks.R(lr, localCol);
            }
          }
          full.Phi.row(gr) += pairBlocks.Phi.row(lr);
        }
      }
    }
  }

  return surge::sd3d::projectBlocksToMonolayer(full);
}

surge::sd3d::Vector repulsiveForces(
    const MonolayerState& state,
    double time,
    const RunConfig& config) {
  const int n = static_cast<int>(state.particles.size());
  surge::sd3d::Vector force = surge::sd3d::Vector::Zero(3 * n);
  if (!config.repulsion) return force;

  const double gammaT = config.gammaStar * time;
  for (int i = 0; i < n; ++i) {
    for (int j = i + 1; j < n; ++j) {
      const PairImage image = shearMinimumImage(state.particles[static_cast<std::size_t>(i)],
                                                state.particles[static_cast<std::size_t>(j)],
                                                state,
                                                gammaT);
      const double r = std::max(image.r, 2.0 + 1e-8);
      const double h = std::max(r - 2.0, 1e-8);
      const double e = std::exp(-config.tau * h);
      const double mag = e / std::max(1.0 - e, 1e-12);
      const double ex = image.dx / r;
      const double ey = image.dy / r;
      force[3 * i + 0] -= mag * ex;
      force[3 * i + 1] -= mag * ey;
      force[3 * j + 0] += mag * ex;
      force[3 * j + 1] += mag * ey;
    }
  }
  return force;
}

surge::sd3d::Vector velocityRhs(
    const MonolayerState& state,
    double time,
    const RunConfig& config) {
  BradyStepResult solved = solveProjectedVelocities(state, time, config);
  surge::sd3d::Vector rhs = surge::sd3d::Vector::Zero(2 * static_cast<int>(state.particles.size()));
  for (int i = 0; i < static_cast<int>(state.particles.size()); ++i) {
    const double y = state.particles[static_cast<std::size_t>(i)].y;
    rhs[2 * i + 0] = solved.qdot[3 * i + 0] + config.gammaStar * y;
    rhs[2 * i + 1] = solved.qdot[3 * i + 1];
  }
  return rhs;
}

void applyIncrement(MonolayerState& state, const surge::sd3d::Vector& rhs, double scale, double newTime, const RunConfig& config) {
  for (int i = 0; i < static_cast<int>(state.particles.size()); ++i) {
    Particle2& p = state.particles[static_cast<std::size_t>(i)];
    p.x += scale * rhs[2 * i + 0];
    p.y += scale * rhs[2 * i + 1];
    wrapShearPeriodic(p, state, config.gammaStar * newTime);
  }
}

void writeTrajectoryHeader(std::ofstream& out, int n) {
  out << "time";
  for (int i = 0; i < n; ++i) out << ",x" << i << ",y" << i;
  out << '\n';
}

void writeTrajectoryFrame(std::ofstream& out, const MonolayerState& state, double time) {
  out << time;
  for (const Particle2& p : state.particles) out << ',' << p.x << ',' << p.y;
  out << '\n';
}

BradyStatistics makeStats(const HistogramOptions& opt, const MonolayerState& state) {
  BradyStatistics stats;
  stats.angularCounts.assign(static_cast<std::size_t>(opt.angularBins), 0.0);
  const int radialBins = static_cast<int>(std::ceil(opt.radialMaxR / opt.radialDr));
  stats.radialCounts.assign(static_cast<std::size_t>(radialBins), 0.0);
  const int yBins = static_cast<int>(std::ceil(state.H / opt.yBinWidth));
  stats.yCounts.assign(static_cast<std::size_t>(yBins), 0.0);
  return stats;
}

void accumulateStats(BradyStatistics& stats, const HistogramOptions& opt, const MonolayerState& state, double time, const RunConfig& config) {
  const int n = static_cast<int>(state.particles.size());
  const double gammaT = config.gammaStar * time;
  for (int i = 0; i < n; ++i) {
    for (int j = i + 1; j < n; ++j) {
      const PairImage image = shearMinimumImage(state.particles[static_cast<std::size_t>(i)],
                                                state.particles[static_cast<std::size_t>(j)],
                                                state,
                                                gammaT);
      double theta = std::atan2(image.dy, image.dx);
      if (theta < 0.0) theta += kPi;
      if (theta >= kPi) theta -= kPi;
      if (image.r >= opt.nearMinR && image.r < opt.nearMaxR) {
        int bin = static_cast<int>(std::floor(theta / kPi * static_cast<double>(opt.angularBins)));
        bin = std::clamp(bin, 0, opt.angularBins - 1);
        stats.angularCounts[static_cast<std::size_t>(bin)] += 2.0;
      }
      if (image.r < opt.radialMaxR) {
        int bin = static_cast<int>(std::floor(image.r / opt.radialDr));
        if (bin >= 0 && bin < static_cast<int>(stats.radialCounts.size())) {
          stats.radialCounts[static_cast<std::size_t>(bin)] += 2.0;
        }
      }
      const double ySep = std::abs(image.dy);
      int yBin = static_cast<int>(std::floor(ySep / opt.yBinWidth));
      if (yBin >= 0 && yBin < static_cast<int>(stats.yCounts.size())) {
        stats.yCounts[static_cast<std::size_t>(yBin)] += 2.0;
      }
    }
  }
  ++stats.samples;
}

void writeAngularCsv(const std::string& path, const BradyStatistics& stats, const HistogramOptions& opt) {
  std::ofstream out(path);
  out << std::setprecision(16);
  out << "theta_deg,count_per_sample\n";
  for (int b = 0; b < opt.angularBins; ++b) {
    const double theta = (static_cast<double>(b) + 0.5) * 180.0 / static_cast<double>(opt.angularBins);
    const double v = stats.samples > 0 ? stats.angularCounts[static_cast<std::size_t>(b)] / stats.samples : 0.0;
    out << theta << ',' << v << '\n';
  }
}

void writeRadialCsv(const std::string& path, const BradyStatistics& stats, const HistogramOptions& opt) {
  std::ofstream out(path);
  out << std::setprecision(16);
  out << "r,count_per_sample\n";
  for (int b = 0; b < static_cast<int>(stats.radialCounts.size()); ++b) {
    const double r = (static_cast<double>(b) + 0.5) * opt.radialDr;
    const double v = stats.samples > 0 ? stats.radialCounts[static_cast<std::size_t>(b)] / stats.samples : 0.0;
    out << r << ',' << v << '\n';
  }
}

void writeYCsv(const std::string& path, const BradyStatistics& stats, const HistogramOptions& opt) {
  std::ofstream out(path);
  out << std::setprecision(16);
  out << "y,count_per_sample\n";
  for (int b = 0; b < static_cast<int>(stats.yCounts.size()); ++b) {
    const double y = (static_cast<double>(b) + 0.5) * opt.yBinWidth;
    const double v = stats.samples > 0 ? stats.yCounts[static_cast<std::size_t>(b)] / stats.samples : 0.0;
    out << y << ',' << v << '\n';
  }
}

}  // namespace

RunConfig presetRun(char runName) {
  RunConfig c;
  c.name = std::string(1, runName);
  switch (runName) {
    case 'A':
      c = {"A", 25, 0.4, 0.5, 227.0, true, 4e-3, 20000, 2000, 5141, 11};
      break;
    case 'D':
      c = {"D", 25, 0.4, 0.5, 0.0, false, 4e-3, 20000, 2000, 5144, 11};
      break;
    case 'E':
      c = {"E", 25, 0.4, 0.5, 10.0, true, 4e-3, 20000, 2000, 5145, 11};
      break;
    case 'F':
      c = {"F", 25, 0.57, 0.5, 10.0, true, 2e-3, 2500, 500, 5146, 11};
      break;
    case 'G':
      c = {"G", 25, 0.57, 0.5, 50.0, true, 4e-3, 16000, 2000, 5147, 11};
      break;
    case 'H':
      c = {"H", 25, 0.1, 1.0, 227.0, true, 4e-3, 20000, 2000, 5148, 11};
      break;
    default:
      throw std::invalid_argument("supported Brady force-superposition presets are A, D, E, F, G, H");
  }
  return c;
}

MonolayerState makeInitialState(const RunConfig& config) {
  MonolayerState state;
  state.L = cellLength(config.particles, config.areaFraction);
  state.H = state.L;
  state.particles.reserve(static_cast<std::size_t>(config.particles));

  const int nx = static_cast<int>(std::ceil(std::sqrt(static_cast<double>(config.particles))));
  const int ny = static_cast<int>(std::ceil(static_cast<double>(config.particles) / nx));
  const double sx = state.L / nx;
  const double sy = state.H / ny;
  std::mt19937_64 rng(config.seed);
  std::uniform_real_distribution<double> jitter(-0.15, 0.15);

  for (int j = 0; j < ny && static_cast<int>(state.particles.size()) < config.particles; ++j) {
    for (int i = 0; i < nx && static_cast<int>(state.particles.size()) < config.particles; ++i) {
      Particle2 p;
      p.x = (i + 0.5 + jitter(rng)) * sx;
      p.y = (j + 0.5 + jitter(rng)) * sy;
      wrapShearPeriodic(p, state, 0.0);
      state.particles.push_back(p);
    }
  }
  return state;
}

BradyStepResult solveProjectedVelocities(
    const MonolayerState& state,
    double time,
    const RunConfig& config) {
  const surge::sd3d::DenseBlocks blocks = assembleMonolayerPairwiseBlocks(state, time, config);
  surge::sd3d::Mat3 E = surge::sd3d::Mat3::Zero();
  E(0, 1) = 0.5;
  E(1, 0) = 0.5;

  BradyStepResult result;
  result.force = repulsiveForces(state, time, config);
  const surge::sd3d::Vector rhs =
      -(config.gammaStar * blocks.Phi * surge::sd3d::flattenStrain(E) + result.force);
  Eigen::LDLT<surge::sd3d::Matrix> ldlt(blocks.R);
  if (ldlt.info() == Eigen::Success) {
    result.qdot = ldlt.solve(rhs);
  } else {
    result.qdot = blocks.R.colPivHouseholderQr().solve(rhs);
  }
  return result;
}

void rk4Step(MonolayerState& state, double time, const RunConfig& config) {
  const double h = config.dt;
  const surge::sd3d::Vector k1 = velocityRhs(state, time, config);

  MonolayerState s2 = state;
  applyIncrement(s2, k1, 0.5 * h, time + 0.5 * h, config);
  const surge::sd3d::Vector k2 = velocityRhs(s2, time + 0.5 * h, config);

  MonolayerState s3 = state;
  applyIncrement(s3, k2, 0.5 * h, time + 0.5 * h, config);
  const surge::sd3d::Vector k3 = velocityRhs(s3, time + 0.5 * h, config);

  MonolayerState s4 = state;
  applyIncrement(s4, k3, h, time + h, config);
  const surge::sd3d::Vector k4 = velocityRhs(s4, time + h, config);

  surge::sd3d::Vector dx = (h / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
  applyIncrement(state, dx, 1.0, time + h, config);
}

void runSimulation(
    const RunConfig& config,
    const HistogramOptions& histOptions,
    const std::string& trajectoryCsv,
    const std::string& angularCsv,
    const std::string& radialCsv,
    const std::string& yCsv) {
  MonolayerState state = makeInitialState(config);
  BradyStatistics stats = makeStats(histOptions, state);

  std::ofstream traj(trajectoryCsv);
  traj << std::setprecision(16);
  writeTrajectoryHeader(traj, static_cast<int>(state.particles.size()));
  writeTrajectoryFrame(traj, state, 0.0);

  double time = 0.0;
  const int stride = std::max(1, config.steps / 1000);
  for (int step = 0; step < config.steps; ++step) {
    rk4Step(state, time, config);
    time += config.dt;
    if (step % stride == 0) writeTrajectoryFrame(traj, state, time);
    if (step >= config.warmupSteps) accumulateStats(stats, histOptions, state, time, config);
  }

  writeAngularCsv(angularCsv, stats, histOptions);
  writeRadialCsv(radialCsv, stats, histOptions);
  writeYCsv(yCsv, stats, histOptions);
}

}  // namespace surge::brady1984
