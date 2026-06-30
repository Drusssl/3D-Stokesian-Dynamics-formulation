#pragma once

#include "core.hpp"

namespace surge::sd3d {

struct EwaldOptions {
  double xi = 0.25;       // splitting parameter
  int realCutImages = 2;  // image shells in each direction
  int kMax = 6;           // reciprocal vector cutoff in each direction
  bool includeSelfMobility = true;
};

Matrix assembleOpenSpaceRPYMobilityTT(const System& system);

// Triply-periodic Stokeslet Ewald translational mobility.
// This is the proper far-field route for periodic 3D SD, but still only the
// force->translation block. Higher multipoles/stresslet blocks are separate work.
Matrix assemblePeriodicStokesletEwaldMobilityTT(const System& system, const EwaldOptions& options);

DenseBlocks farFieldResistanceFromMobilityTT(const System& system, const Matrix& Mtt);

}  // namespace surge::sd3d
