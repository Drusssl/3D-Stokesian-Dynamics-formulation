#include "surge/sd3d/near_field.hpp"

#include <stdexcept>

#include "surge_pair_resistance.hpp"

namespace surge::sd3d {
namespace {

void addSelf(DenseBlocks& out, int p, double a, double mu) {
  const double trans = 6.0 * surge::kPi * mu * a;
  const double rot = 8.0 * surge::kPi * mu * a * a * a;
  for (int c = 0; c < 3; ++c) {
    out.R(dof(p, c, false), dof(p, c, false)) += trans;
    out.R(dof(p, c, true), dof(p, c, true)) += rot;
  }
}

void subtractLocalSelf(surge::PairBlocks& pair, int local, double a, double mu) {
  const double trans = 6.0 * surge::kPi * mu * a;
  const double rot = 8.0 * surge::kPi * mu * a * a * a;
  for (int c = 0; c < 3; ++c) {
    const int u = local * 6 + c;
    const int w = local * 6 + 3 + c;
    pair.R[static_cast<std::size_t>(u) * 12 + u] -= trans;
    pair.R[static_cast<std::size_t>(w) * 12 + w] -= rot;
  }
}

void scatterPair(DenseBlocks& out, int i, int j, const surge::PairBlocks& pair) {
  const int gp[2] = {i, j};
  for (int la = 0; la < 2; ++la) {
    for (int lc = 0; lc < 6; ++lc) {
      const int gr = gp[la] * 6 + lc;
      const int lr = la * 6 + lc;
      for (int lb = 0; lb < 2; ++lb) {
        for (int ld = 0; ld < 6; ++ld) {
          const int gc = gp[lb] * 6 + ld;
          const int lcol = lb * 6 + ld;
          out.R(gr, gc) += pair.R[static_cast<std::size_t>(lr) * 12 + lcol];
        }
      }
      for (int a = 0; a < 3; ++a) {
        for (int b = 0; b < 3; ++b) {
          out.Phi(gr, strainIndex(a, b)) += pair.Phi[(static_cast<std::size_t>(lr) * 3 + a) * 3 + b];
        }
      }
    }
  }
}

}  // namespace

DenseBlocks assemblePairwiseResistance3D(const System& system, const NearFieldOptions& options) {
  if (system.viscosity <= 0.0) throw std::invalid_argument("viscosity must be positive");
  if (options.highOrderTable != nullptr) {
    // The hook is here deliberately; surge_pair_resistance.hpp still uses the
    // printed coefficients internally. The next scientific step is replacing
    // those scalar evaluators with table/recurrence-backed evaluators.
    (void)options.highOrderTable;
  }

  const int n = static_cast<int>(system.particles.size());
  DenseBlocks out;
  out.n = n;
  out.R = Matrix::Zero(6 * n, 6 * n);
  out.Phi = Matrix::Zero(6 * n, 9);

  for (int p = 0; p < n; ++p) {
    if (system.particles[p].radius <= 0.0) throw std::invalid_argument("particle radius must be positive");
    addSelf(out, p, system.particles[p].radius, system.viscosity);
  }

  surge::EvalOptions pairOptions;
  pairOptions.maxPrintedM = options.maxPrintedM;

  for (int i = 0; i < n; ++i) {
    for (int j = i + 1; j < n; ++j) {
      const Vec3 dx = minimumImage(system.particles[j].x - system.particles[i].x, system.box);
      surge::PairInput pairInput;
      pairInput.x1 = {0.0, 0.0, 0.0};
      pairInput.x2 = {dx.x(), dx.y(), dx.z()};
      pairInput.a1 = system.particles[i].radius;
      pairInput.a2 = system.particles[j].radius;
      pairInput.mu = system.viscosity;

      surge::PairBlocks pair = surge::computePairResistance(pairInput, pairOptions);
      subtractLocalSelf(pair, 0, system.particles[i].radius, system.viscosity);
      subtractLocalSelf(pair, 1, system.particles[j].radius, system.viscosity);
      scatterPair(out, i, j, pair);
    }
  }
  return out;
}

}  // namespace surge::sd3d
