#include "surge/sd3d/ewald.hpp"

#include <cmath>
#include <stdexcept>

#include "surge_pair_resistance.hpp"

namespace surge::sd3d {
namespace {

Mat3 openSpaceRPYBlock(Vec3 rvec, double a, double mu) {
  const double r = rvec.norm();
  if (r <= 0.0) throw std::invalid_argument("zero separation");
  const Vec3 d = rvec / r;
  if (r >= 2.0 * a) {
    const double c = 1.0 / (8.0 * surge::kPi * mu * r);
    const double alpha = 1.0 + 2.0 * a * a / (3.0 * r * r);
    const double beta = 1.0 - 2.0 * a * a / (r * r);
    return c * (alpha * Mat3::Identity() + beta * d * d.transpose());
  }
  const double c = 1.0 / (6.0 * surge::kPi * mu * a);
  const double rho = r / a;
  return c * ((1.0 - 9.0 * rho / 32.0) * Mat3::Identity()
      + (3.0 * rho / 32.0) * d * d.transpose());
}

Mat3 realSpaceStokesletEwald(Vec3 rvec, double xi) {
  const double r = rvec.norm();
  if (r <= 0.0) return Mat3::Zero();
  const Vec3 d = rvec / r;
  const Mat3 rr = d * d.transpose();
  const double xr = xi * r;
  const double erfcTerm = std::erfc(xr);
  const double expTerm = std::exp(-xr * xr);
  const double c = 2.0 * xi / std::sqrt(surge::kPi) * expTerm;

  // Hasimoto-style real-space stokeslet split.
  return erfcTerm * (Mat3::Identity() + rr) / r
      - c * (Mat3::Identity() - rr);
}

}  // namespace

Matrix assembleOpenSpaceRPYMobilityTT(const System& system) {
  const int n = static_cast<int>(system.particles.size());
  Matrix M = Matrix::Zero(3 * n, 3 * n);
  for (int i = 0; i < n; ++i) {
    const double a = system.particles[i].radius;
    M.block<3, 3>(3 * i, 3 * i) =
        (1.0 / (6.0 * surge::kPi * system.viscosity * a)) * Mat3::Identity();
  }
  for (int i = 0; i < n; ++i) {
    for (int j = i + 1; j < n; ++j) {
      const double a = 0.5 * (system.particles[i].radius + system.particles[j].radius);
      const Mat3 block = openSpaceRPYBlock(system.particles[j].x - system.particles[i].x, a, system.viscosity);
      M.block<3, 3>(3 * i, 3 * j) = block;
      M.block<3, 3>(3 * j, 3 * i) = block.transpose();
    }
  }
  return M;
}

Matrix assemblePeriodicStokesletEwaldMobilityTT(const System& system, const EwaldOptions& options) {
  const int n = static_cast<int>(system.particles.size());
  if (!(system.box.periodicX && system.box.periodicY && system.box.periodicZ)) {
    throw std::invalid_argument("periodic Stokeslet Ewald requires a triply periodic box");
  }
  const double V = boxVolume(system.box);
  if (V <= 0.0) throw std::invalid_argument("box volume must be positive");

  Matrix M = Matrix::Zero(3 * n, 3 * n);
  const Vec3 L = system.box.length;
  const double xi = options.xi;

  for (int i = 0; i < n; ++i) {
    if (options.includeSelfMobility) {
      const double a = system.particles[i].radius;
      M.block<3, 3>(3 * i, 3 * i) +=
          (1.0 / (6.0 * surge::kPi * system.viscosity * a)) * Mat3::Identity();
    }
  }

  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      Mat3 J = Mat3::Zero();
      const Vec3 base = system.particles[i].x - system.particles[j].x;
      for (int nx = -options.realCutImages; nx <= options.realCutImages; ++nx) {
        for (int ny = -options.realCutImages; ny <= options.realCutImages; ++ny) {
          for (int nz = -options.realCutImages; nz <= options.realCutImages; ++nz) {
            if (i == j && nx == 0 && ny == 0 && nz == 0) continue;
            const Vec3 image(nx * L.x(), ny * L.y(), nz * L.z());
            J += realSpaceStokesletEwald(base + image, xi);
          }
        }
      }

      for (int hx = -options.kMax; hx <= options.kMax; ++hx) {
        for (int hy = -options.kMax; hy <= options.kMax; ++hy) {
          for (int hz = -options.kMax; hz <= options.kMax; ++hz) {
            if (hx == 0 && hy == 0 && hz == 0) continue;
            const Vec3 kvec(
                2.0 * surge::kPi * hx / L.x(),
                2.0 * surge::kPi * hy / L.y(),
                2.0 * surge::kPi * hz / L.z());
            const double k2 = kvec.squaredNorm();
            const double phase = kvec.dot(base);
            const double factor = (8.0 * surge::kPi / V)
                * (1.0 + k2 / (4.0 * xi * xi))
                * std::exp(-k2 / (4.0 * xi * xi))
                * std::cos(phase) / k2;
            J += factor * (Mat3::Identity() - (kvec * kvec.transpose()) / k2);
          }
        }
      }

      M.block<3, 3>(3 * i, 3 * j) += (1.0 / (8.0 * surge::kPi * system.viscosity)) * J;
    }
  }

  return M;
}

DenseBlocks farFieldResistanceFromMobilityTT(const System& system, const Matrix& Mtt) {
  const int n = static_cast<int>(system.particles.size());
  DenseBlocks out;
  out.n = n;
  out.R = Matrix::Zero(6 * n, 6 * n);
  out.Phi = Matrix::Zero(6 * n, 9);

  Eigen::LDLT<Matrix> ldlt(Mtt);
  Matrix Rtt;
  if (ldlt.info() == Eigen::Success) {
    Rtt = ldlt.solve(Matrix::Identity(3 * n, 3 * n));
  } else {
    Rtt = Mtt.colPivHouseholderQr().solve(Matrix::Identity(3 * n, 3 * n));
  }

  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      out.R.block<3, 3>(6 * i, 6 * j) = Rtt.block<3, 3>(3 * i, 3 * j);
    }
    const double a = system.particles[i].radius;
    const double rot = 8.0 * surge::kPi * system.viscosity * a * a * a;
    out.R.block<3, 3>(6 * i + 3, 6 * i + 3) = rot * Mat3::Identity();
  }
  return out;
}

}  // namespace surge::sd3d
