#pragma once

#include <array>
#include <cmath>
#include <stdexcept>

namespace surge {

constexpr double kPi = 3.141592653589793238462643383279502884;

struct Vec3 {
  double x = 0.0;
  double y = 0.0;
  double z = 0.0;
};

struct PairInput {
  Vec3 x1;
  Vec3 x2;
  double a1 = 1.0;
  double a2 = 1.0;
  double mu = 1.0;
};

struct PairBlocks {
  // Order per particle is [Fx,Fy,Fz,Lx,Ly,Lz].
  // Full vector order is particle 1 first, particle 2 second.
  std::array<double, 12 * 12> R{};

  // Phi[row, j, k] contracts with E[j,k] to give force/torque row.
  std::array<double, 12 * 3 * 3> Phi{};

  double r = 0.0;
  double s = 0.0;
  double lambda = 0.0;
  Vec3 d{};
};

struct EvalOptions {
  // The paper prints coefficients through m=11. Higher accuracy needs longer
  // coefficient tables or recurrence-generated coefficients.
  int maxPrintedM = 11;
};

namespace detail {

enum class Family { XA, YA, YB, XC, YC, XG, YG, YH };

struct Scalars {
  double XA[2][2]{};
  double YA[2][2]{};
  double YB[2][2]{};
  double XC[2][2]{};
  double YC[2][2]{};
  double XG[2][2]{};
  double YG[2][2]{};
  double YH[2][2]{};
};

inline double& scalarRef(Scalars& s, Family f, int a, int b) {
  switch (f) {
    case Family::XA: return s.XA[a][b];
    case Family::YA: return s.YA[a][b];
    case Family::YB: return s.YB[a][b];
    case Family::XC: return s.XC[a][b];
    case Family::YC: return s.YC[a][b];
    case Family::XG: return s.XG[a][b];
    case Family::YG: return s.YG[a][b];
    case Family::YH: return s.YH[a][b];
  }
  return s.XA[a][b];
}

inline double sqr(double x) { return x * x; }
inline double cube(double x) { return x * x * x; }
inline double powi(double x, int n) {
  if (n < 0) return 1.0 / powi(x, -n);
  double v = 1.0;
  for (int i = 0; i < n; ++i) v *= x;
  return v;
}

inline Vec3 sub(Vec3 a, Vec3 b) { return {a.x - b.x, a.y - b.y, a.z - b.z}; }
inline double dot(Vec3 a, Vec3 b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
inline double norm(Vec3 a) { return std::sqrt(dot(a, a)); }
inline Vec3 scale(Vec3 a, double c) { return {c * a.x, c * a.y, c * a.z}; }

inline double eps(int i, int j, int k) {
  if (i == j || j == k || i == k) return 0.0;
  if ((i == 0 && j == 1 && k == 2) ||
      (i == 1 && j == 2 && k == 0) ||
      (i == 2 && j == 0 && k == 1)) {
    return 1.0;
  }
  return -1.0;
}

inline double dcomp(Vec3 d, int i) {
  return i == 0 ? d.x : (i == 1 ? d.y : d.z);
}

inline double delta(int i, int j) { return i == j ? 1.0 : 0.0; }

inline int m1(int m) { return m == 2 ? -2 : (m - 2); }

inline double safeLog(double x) {
  if (x <= 0.0) throw std::domain_error("log argument must be positive");
  return std::log(x);
}

inline double fXA(int m, double l) {
  switch (m) {
    case 0: return 1.0;
    case 1: return 3.0 * l;
    case 2: return 9.0 * l;
    case 3: return -4.0 * l + 27.0 * powi(l, 2) - 4.0 * powi(l, 3);
    case 4: return -24.0 * l + 81.0 * powi(l, 2) + 36.0 * powi(l, 3);
    case 5: return 72.0 * powi(l, 2) + 243.0 * powi(l, 3) + 72.0 * powi(l, 4);
    case 6: return 16.0 * l + 108.0 * powi(l, 2) + 281.0 * powi(l, 3) + 648.0 * powi(l, 4) + 144.0 * powi(l, 5);
    case 7: return 288.0 * powi(l, 2) + 1620.0 * powi(l, 3) + 1515.0 * powi(l, 4) + 1620.0 * powi(l, 5) + 288.0 * powi(l, 6);
    case 8: return 576.0 * powi(l, 2) + 4848.0 * powi(l, 3) + 5409.0 * powi(l, 4) + 4524.0 * powi(l, 5) + 3888.0 * powi(l, 6) + 576.0 * powi(l, 7);
    case 9: return 1152.0 * powi(l, 2) + 9072.0 * powi(l, 3) + 14752.0 * powi(l, 4) + 26163.0 * powi(l, 5) + 14752.0 * powi(l, 6) + 9072.0 * powi(l, 7) + 1152.0 * powi(l, 8);
    case 10: return 2304.0 * powi(l, 2) + 20736.0 * powi(l, 3) + 42804.0 * powi(l, 4) + 115849.0 * powi(l, 5) + 76176.0 * powi(l, 6) + 39264.0 * powi(l, 7) + 20736.0 * powi(l, 8) + 2304.0 * powi(l, 9);
    case 11: return 4608.0 * powi(l, 2) + 46656.0 * powi(l, 3) + 108912.0 * powi(l, 4) + 269100.0 * powi(l, 5) + 319899.0 * powi(l, 6) + 269100.0 * powi(l, 7) + 108912.0 * powi(l, 8) + 46656.0 * powi(l, 9) + 4608.0 * powi(l, 10);
    default: return 0.0;
  }
}

inline double fYA(int m, double l) {
  switch (m) {
    case 0: return 1.0;
    case 1: return 1.5 * l;
    case 2: return 2.25 * l;
    case 3: return 2.0 * l + 27.0 / 8.0 * powi(l, 2) + 2.0 * powi(l, 3);
    case 4: return 6.0 * l + 81.0 / 16.0 * powi(l, 2) + 18.0 * powi(l, 3);
    case 5: return 63.0 / 2.0 * powi(l, 2) + 243.0 / 32.0 * powi(l, 3) + 63.0 / 2.0 * powi(l, 4);
    default: return 0.0;
  }
}

inline double fYB(int m, double l) {
  switch (m) {
    case 0:
    case 1: return 0.0;
    case 2: return -6.0 * l;
    case 3: return -9.0 * l;
    case 4: return -27.0 / 2.0 * powi(l, 2);
    case 5: return -12.0 * l - 81.0 / 4.0 * powi(l, 2) - 36.0 * powi(l, 3);
    case 6: return -108.0 * powi(l, 2) - 243.0 / 8.0 * powi(l, 3) - 72.0 * powi(l, 4);
    default: return 0.0;
  }
}

inline double fXC(int m, double l) {
  switch (m) {
    case 0: return 1.0;
    case 1:
    case 2: return 0.0;
    case 3: return 8.0 * powi(l, 3);
    case 4:
    case 5: return 0.0;
    case 6: return 64.0 * powi(l, 3);
    case 7: return 0.0;
    case 8: return 768.0 * powi(l, 5);
    case 9: return 512.0 * powi(l, 6);
    case 10: return 6144.0 * powi(l, 7);
    case 11: return 6144.0 * powi(l, 6) + 6144.0 * powi(l, 8);
    default: return 0.0;
  }
}

inline double fYC(int m, double l) {
  switch (m) {
    case 0: return 1.0;
    case 1:
    case 2: return 0.0;
    case 3: return 4.0 * powi(l, 3);
    case 4: return 12.0 * l;
    case 5: return 18.0 * powi(l, 4);
    case 6: return 27.0 * powi(l, 2) + 256.0 * powi(l, 3);
    case 7: return 72.0 * powi(l, 4) + 81.0 / 2.0 * powi(l, 5) + 72.0 * powi(l, 6);
    default: return 0.0;
  }
}

inline double fXG(int m, double l) {
  switch (m) {
    case 0:
    case 1: return 0.0;
    case 2: return 15.0 * l;
    case 3: return 45.0 * l;
    case 4: return -36.0 * l + 135.0 * powi(l, 2) - 60.0 * powi(l, 3);
    case 5: return -168.0 * l + 405.0 * powi(l, 2) + 360.0 * powi(l, 3);
    case 6: return 216.0 * powi(l, 2) + 1215.0 * powi(l, 3) + 900.0 * powi(l, 4);
    case 7: return 144.0 * l + 108.0 * powi(l, 2) - 1251.0 * powi(l, 3) + 4860.0 * powi(l, 4) + 2160.0 * powi(l, 5);
    case 8: return 864.0 * powi(l, 2) + 6804.0 * powi(l, 3) + 8679.0 * powi(l, 4) + 12960.0 * powi(l, 5) + 5040.0 * powi(l, 6);
    case 9: return 1728.0 * powi(l, 2) + 27072.0 * powi(l, 3) + 25893.0 * powi(l, 4) - 3888.0 * powi(l, 5) + 38880.0 * powi(l, 6) + 11520.0 * powi(l, 7);
    case 10: return 3456.0 * powi(l, 2) + 34992.0 * powi(l, 3) + 66384.0 * powi(l, 4) + 155007.0 * powi(l, 5) + 75900.0 * powi(l, 6) + 97200.0 * powi(l, 7) + 25920.0 * powi(l, 8);
    case 11: return 6912.0 * powi(l, 2) + 77760.0 * powi(l, 3) + 178260.0 * powi(l, 4) + 790845.0 * powi(l, 5) + 500580.0 * powi(l, 6) + 14880.0 * powi(l, 7) + 259200.0 * powi(l, 8) + 57600.0 * powi(l, 9);
    default: return 0.0;
  }
}

inline double fYG(int m, double l) {
  switch (m) {
    case 0:
    case 1:
    case 2:
    case 3: return 0.0;
    case 4: return 12.0 * l + 20.0 * powi(l, 3);
    case 5: return 18.0 * l + 90.0 * powi(l, 3);
    case 6: return 27.0 * powi(l, 2) + 135.0 * powi(l, 4);
    case 7: return 24.0 * l + 81.0 / 2.0 * powi(l, 2) - 336.0 * powi(l, 3) + 405.0 / 2.0 * powi(l, 4) + 600.0 * powi(l, 5);
    case 8: return 216.0 * powi(l, 2) + 243.0 / 4.0 * powi(l, 3) - 1008.0 * powi(l, 4) + 1215.0 / 4.0 * powi(l, 5) + 1080.0 * powi(l, 6);
    case 9: return 378.0 * powi(l, 2) + 21209.0 / 8.0 * powi(l, 3) - 972.0 * powi(l, 4) - 62147.0 / 8.0 * powi(l, 5) + 2970.0 * powi(l, 6) + 3360.0 * powi(l, 7);
    case 10: return 864.0 * powi(l, 2) + 972.0 * powi(l, 3) + 42123.0 / 16.0 * powi(l, 4) + 648.0 * powi(l, 5) - 96073.0 / 16.0 * powi(l, 6) + 4860.0 * powi(l, 7) + 6240.0 * powi(l, 8);
    case 11: return 1728.0 * powi(l, 2) + 3159.0 / 2.0 * powi(l, 3) + 103329.0 / 32.0 * powi(l, 4) + 69387.0 * powi(l, 5) + 763173.0 / 32.0 * powi(l, 6) - 166737.0 / 2.0 * powi(l, 7) + 24840.0 * powi(l, 8) + 17280.0 * powi(l, 9);
    default: return 0.0;
  }
}

inline double fYH(int m, double l) {
  switch (m) {
    case 0:
    case 1:
    case 2: return 0.0;
    case 3: return 10.0 * powi(l, 3);
    case 4:
    case 5: return 0.0;
    case 6: return 24.0 * l - 120.0 * powi(l, 3);
    case 7: return 36.0 * powi(l, 4) + 180.0 * powi(l, 6);
    case 8: return 54.0 * powi(l, 2) + 1280.0 * powi(l, 3) + 270.0 * powi(l, 4) - 1280.0 * powi(l, 5);
    case 9: return 144.0 * powi(l, 4) + 81.0 * powi(l, 5) + 5248.0 * powi(l, 6) + 405.0 * powi(l, 7) + 1200.0 * powi(l, 8);
    case 10: return 432.0 * powi(l, 2) + 243.0 / 2.0 * powi(l, 3) - 1872.0 * powi(l, 4) + 46015.0 / 2.0 * powi(l, 5) + 2880.0 * powi(l, 6) - 9600.0 * powi(l, 7);
    case 11: return 576.0 * powi(l, 4) + 972.0 * powi(l, 5) + 153049.0 / 4.0 * powi(l, 6) - 864.0 * powi(l, 7) + 186173.0 / 4.0 * powi(l, 8) + 5940.0 * powi(l, 9) + 6720.0 * powi(l, 10);
    default: return 0.0;
  }
}

using FFunc = double (*)(int, double);

inline double sumCorrected(int start, int end, int parity, double l, double s, double g1, double g2, double g3, bool includeG1, FFunc f) {
  double sum = 0.0;
  const int maxM = end;
  for (int m = start; m <= maxM; ++m) {
    if ((m & 1) != parity) continue;
    const double fmTilde = std::ldexp(f(m, l), -m);
    double term = fmTilde / powi(1.0 + l, m);
    if (includeG1) term -= g1;
    term -= 2.0 * g2 / static_cast<double>(m);
    term += 4.0 * g3 / (static_cast<double>(m) * static_cast<double>(m1(m)));
    sum += term * powi(2.0 / s, m);
  }
  return sum;
}

inline void eval11_12(Family fam, double l, double s, const EvalOptions& opt, double& v11, double& v12) {
  if (s <= 2.0) throw std::domain_error("s must be greater than 2");
  const double one = 1.0 + l;
  const double q = 1.0 - 4.0 / (s * s);
  const double logPm = safeLog((s + 2.0) / (s - 2.0));
  const double logQ = safeLog(q);
  const double logS2 = safeLog((s * s) / (s * s - 4.0));
  const int maxM = opt.maxPrintedM;

  switch (fam) {
    case Family::XA: {
      const double g1 = 2.0 * sqr(l) / powi(one, 3);
      const double g2 = (1.0 / 5.0) * l * (1.0 + 7.0 * l + sqr(l)) / powi(one, 3);
      const double g3 = (1.0 / 42.0) * (1.0 + 18.0 * l - 29.0 * sqr(l) + 18.0 * cube(l) + powi(l, 4)) / powi(one, 3);
      v11 = g1 / q - g2 * logQ - g3 * q * logQ + fXA(0, l) - g1
          + sumCorrected(2, maxM, 0, l, s, g1, g2, g3, true, fXA);
      const double scaled = 2.0 * g1 / s / q + g2 * logPm + g3 * q * logPm + 4.0 * g3 / s
          + sumCorrected(1, maxM, 1, l, s, g1, g2, g3, true, fXA);
      v12 = -2.0 * scaled / one;
      break;
    }
    case Family::YA: {
      const double g2 = (4.0 / 15.0) * l * (2.0 + l + 2.0 * sqr(l)) / powi(one, 3);
      const double g3 = (2.0 / 375.0) * (16.0 - 45.0 * l + 58.0 * sqr(l) - 45.0 * cube(l) + 16.0 * powi(l, 4)) / powi(one, 3);
      v11 = -g2 * logQ - g3 * q * logQ + fYA(0, l)
          + sumCorrected(2, maxM, 0, l, s, 0.0, g2, g3, false, fYA);
      const double scaled = g2 * logPm + g3 * q * logPm + 4.0 * g3 / s
          + sumCorrected(1, maxM, 1, l, s, 0.0, g2, g3, false, fYA);
      v12 = -2.0 * scaled / one;
      break;
    }
    case Family::YB: {
      const double g2 = -(1.0 / 5.0) * l * (4.0 + l) / sqr(one);
      const double g3 = -(1.0 / 250.0) * (32.0 - 33.0 * l + 83.0 * sqr(l) + 43.0 * cube(l)) / sqr(one);
      v11 = g2 * logPm + g3 * q * logPm + 4.0 * g3 / s
          + sumCorrected(1, maxM, 1, l, s, 0.0, g2, g3, false, fYB);
      const double scaled = -g2 * logQ - g3 * q * logQ
          + sumCorrected(2, maxM, 0, l, s, 0.0, g2, g3, false, fYB);
      v12 = -4.0 * scaled / sqr(one);
      break;
    }
    case Family::XC: {
      // Accelerated form from the documentation. Contact zeta constants are
      // intentionally not used in this all-separation evaluator.
      v11 = -sqr(l) / (2.0 * one) * logQ + sqr(l) / one / s * logPm + 1.0;
      for (int k = 1; 2 * k <= maxM; ++k) {
        const int m = 2 * k;
        const double corr = powi(one, -m) * fXC(m, l)
            - std::pow(2.0, 2 * k + 1) / (static_cast<double>(k) * (2.0 * k - 1.0)) * sqr(l) / (4.0 * one);
        v11 += corr * std::pow(s, -m);
      }
      v12 = 4.0 * sqr(l) / powi(one, 4) * logPm + 8.0 * sqr(l) / powi(one, 4) / s * logQ;
      double series = 0.0;
      for (int k = 1; 2 * k + 1 <= maxM; ++k) {
        const int m = 2 * k + 1;
        const double corr = powi(one, -m) * fXC(m, l)
            - std::pow(2.0, 2 * k + 2) / (static_cast<double>(k) * (2.0 * k + 1.0)) * sqr(l) / one;
        series += corr * std::pow(s, -m);
      }
      v12 += -8.0 / powi(one, 3) * series;
      break;
    }
    case Family::YC: {
      const double g2 = (2.0 / 5.0) * l / one;
      const double g3 = (1.0 / 125.0) * (8.0 + 6.0 * l + 33.0 * sqr(l)) / one;
      const double g4 = (4.0 / 5.0) * sqr(l) / powi(one, 4);
      const double g5 = (4.0 / 125.0) * l * (43.0 - 24.0 * l + 43.0 * sqr(l)) / powi(one, 4);
      v11 = -g2 * logQ - g3 * q * logQ + fYC(0, l)
          + sumCorrected(2, maxM, 0, l, s, 0.0, g2, g3, false, fYC);
      const double scaled = g4 * logPm + g5 * q * logPm + 4.0 * g5 / s
          + sumCorrected(1, maxM, 1, l, s, 0.0, g4, g5, false, fYC);
      v12 = 8.0 * scaled / powi(one, 3);
      break;
    }
    case Family::XG: {
      const double g1 = 3.0 * sqr(l) / powi(one, 3);
      const double g2 = (3.0 / 10.0) * (l + 12.0 * sqr(l) - 4.0 * cube(l)) / powi(one, 3);
      const double g3 = (5.0 + 181.0 * l - 453.0 * sqr(l) + 566.0 * cube(l) - 65.0 * powi(l, 4)) / (140.0 * powi(one, 3));
      v11 = g1 * 2.0 * s / (s * s - 4.0) + (g2 + g3 * (s * s / 4.0 - 1.0)) * logPm - g3 * s
          + sumCorrected(1, maxM, 1, l, s, g1, g2, g3, true, fXG);
      const double scaled = -g1 * 4.0 / (s * s - 4.0) - (g2 + g3 * (s * s / 4.0 - 1.0)) * logS2 + g3
          - sumCorrected(2, maxM, 0, l, s, g1, g2, g3, true, fXG);
      v12 = 4.0 * scaled / sqr(one);
      break;
    }
    case Family::YG: {
      const double g2 = (1.0 / 10.0) * (4.0 * l - sqr(l) + 7.0 * cube(l)) / powi(one, 3);
      const double g3 = (32.0 - 179.0 * l + 532.0 * sqr(l) - 356.0 * cube(l) + 221.0 * powi(l, 4)) / (500.0 * powi(one, 3));
      v11 = (g2 + g3 * (s * s / 4.0 - 1.0)) * logPm - g3 * s
          + sumCorrected(1, maxM, 1, l, s, 0.0, g2, g3, false, fYG);
      const double scaled = -(g2 + g3 * (s * s / 4.0 - 1.0)) * logS2 + g3
          - sumCorrected(2, maxM, 0, l, s, 0.0, g2, g3, false, fYG);
      v12 = 4.0 * scaled / sqr(one);
      break;
    }
    case Family::YH: {
      const double g2 = (1.0 / 10.0) * (2.0 * l - sqr(l)) / sqr(one);
      const double g3 = (16.0 - 61.0 * l + 180.0 * sqr(l) + 2.0 * cube(l)) / (500.0 * sqr(one));
      const double g5 = (1.0 / 20.0) * (sqr(l) + 7.0 * cube(l)) / sqr(one);
      const double g6 = (43.0 * l + 147.0 * sqr(l) - 185.0 * cube(l) + 221.0 * powi(l, 4)) / (1000.0 * sqr(one));
      v11 = (g2 + g3 * (s * s / 4.0 - 1.0)) * logS2 - g3
          + sumCorrected(2, maxM, 0, l, s, 0.0, g2, g3, false, fYH);
      const double scaled = (g5 + g6 * (s * s / 4.0 - 1.0)) * logPm - g6 * s
          + sumCorrected(1, maxM, 1, l, s, 0.0, g5, g6, false, fYH);
      v12 = 8.0 * scaled / powi(one, 3);
      break;
    }
  }
}

inline void fillFamily(Scalars& sc, Family fam, double l, double s, const EvalOptions& opt) {
  double v11 = 0.0;
  double v12 = 0.0;
  eval11_12(fam, l, s, opt, v11, v12);
  scalarRef(sc, fam, 0, 0) = v11;
  scalarRef(sc, fam, 0, 1) = v12;

  const double li = 1.0 / l;
  double inv11 = 0.0;
  double inv12 = 0.0;
  eval11_12(fam, li, s, opt, inv11, inv12);

  switch (fam) {
    case Family::XA:
    case Family::YA:
    case Family::XC:
    case Family::YC:
      scalarRef(sc, fam, 1, 1) = inv11;
      scalarRef(sc, fam, 1, 0) = inv12;
      break;
    case Family::YB:
    case Family::XG:
    case Family::YG:
      scalarRef(sc, fam, 1, 1) = -inv11;
      scalarRef(sc, fam, 1, 0) = -inv12;
      break;
    case Family::YH:
      scalarRef(sc, fam, 1, 1) = inv11;
      scalarRef(sc, fam, 1, 0) = inv12;
      break;
  }
}

inline Scalars evaluateScalars(double lambda, double s, const EvalOptions& opt) {
  Scalars sc;
  fillFamily(sc, Family::XA, lambda, s, opt);
  fillFamily(sc, Family::YA, lambda, s, opt);
  fillFamily(sc, Family::YB, lambda, s, opt);
  fillFamily(sc, Family::XC, lambda, s, opt);
  fillFamily(sc, Family::YC, lambda, s, opt);
  fillFamily(sc, Family::XG, lambda, s, opt);
  fillFamily(sc, Family::YG, lambda, s, opt);
  fillFamily(sc, Family::YH, lambda, s, opt);
  return sc;
}

inline int row(int particle, int component, bool torque) {
  return particle * 6 + (torque ? 3 : 0) + component;
}

inline int col(int particle, int component, bool angular) {
  return particle * 6 + (angular ? 3 : 0) + component;
}

inline double& RAt(PairBlocks& out, int i, int j) {
  return out.R[static_cast<std::size_t>(i) * 12 + j];
}

inline double& PhiAt(PairBlocks& out, int rowIdx, int j, int k) {
  return out.Phi[(static_cast<std::size_t>(rowIdx) * 3 + j) * 3 + k];
}

}  // namespace detail

inline PairBlocks computePairResistance(const PairInput& in, const EvalOptions& opt = {}) {
  if (in.a1 <= 0.0 || in.a2 <= 0.0 || in.mu <= 0.0) {
    throw std::invalid_argument("radii and viscosity must be positive");
  }

  PairBlocks out;
  const Vec3 rvec = detail::sub(in.x2, in.x1);
  out.r = detail::norm(rvec);
  if (out.r == 0.0) throw std::invalid_argument("sphere centers must be distinct");
  out.d = detail::scale(rvec, 1.0 / out.r);
  out.s = 2.0 * out.r / (in.a1 + in.a2);
  if (out.s <= 2.0) throw std::domain_error("spheres overlap or touch; this evaluator requires s > 2");
  out.lambda = in.a2 / in.a1;

  const auto sc = detail::evaluateScalars(out.lambda, out.s, opt);
  const double a[2] = {in.a1, in.a2};
  const double d[3] = {out.d.x, out.d.y, out.d.z};

  for (int alpha = 0; alpha < 2; ++alpha) {
    for (int beta = 0; beta < 2; ++beta) {
      const double apb = a[alpha] + a[beta];
      const double prefA = 3.0 * kPi * in.mu * apb;
      const double prefB = kPi * in.mu * apb * apb;
      const double prefC = kPi * in.mu * apb * apb * apb;
      const double prefG = kPi * in.mu * apb * apb;
      const double prefH = kPi * in.mu * apb * apb * apb;

      for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
          const double ppar = d[i] * d[j];
          const double pperp = detail::delta(i, j) - ppar;
          const double Aij = prefA * (sc.XA[alpha][beta] * ppar + sc.YA[alpha][beta] * pperp);
          const double Cij = prefC * (sc.XC[alpha][beta] * ppar + sc.YC[alpha][beta] * pperp);
          double Bij = 0.0;
          for (int k = 0; k < 3; ++k) Bij += detail::eps(i, j, k) * d[k];
          Bij *= prefB * sc.YB[alpha][beta];

          detail::RAt(out, detail::row(alpha, i, false), detail::col(beta, j, false)) = Aij;
          detail::RAt(out, detail::row(alpha, i, true), detail::col(beta, j, true)) = Cij;
          detail::RAt(out, detail::row(alpha, i, true), detail::col(beta, j, false)) = Bij;

          // Reciprocal force-rotation block: Btilde_{alpha beta} = B^T_{beta alpha}.
          double Btij = 0.0;
          for (int k = 0; k < 3; ++k) Btij += detail::eps(j, i, k) * d[k];
          Btij *= prefB * sc.YB[beta][alpha];
          detail::RAt(out, detail::row(alpha, i, false), detail::col(beta, j, true)) = Btij;
        }
      }

      for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
          for (int k = 0; k < 3; ++k) {
            const double Ghat =
                sc.XG[alpha][beta] * d[i] * (d[j] * d[k] - detail::delta(j, k) / 3.0)
                + sc.YG[alpha][beta] * (d[j] * detail::delta(i, k) + d[k] * detail::delta(i, j) - 2.0 * d[i] * d[j] * d[k]);
            double Hhat = 0.0;
            for (int ell = 0; ell < 3; ++ell) {
              Hhat += sc.YH[alpha][beta] *
                  (detail::eps(i, ell, j) * d[ell] * d[k] + detail::eps(i, ell, k) * d[ell] * d[j]);
            }
            detail::PhiAt(out, detail::row(alpha, i, false), j, k) += prefG * Ghat;
            detail::PhiAt(out, detail::row(alpha, i, true), j, k) += prefH * Hhat;
          }
        }
      }
    }
  }

  return out;
}

inline std::array<double, 12> contractPhi(const PairBlocks& blocks, const std::array<double, 9>& E) {
  std::array<double, 12> out{};
  for (int row = 0; row < 12; ++row) {
    double v = 0.0;
    for (int j = 0; j < 3; ++j) {
      for (int k = 0; k < 3; ++k) {
        v += blocks.Phi[(static_cast<std::size_t>(row) * 3 + j) * 3 + k] * E[j * 3 + k];
      }
    }
    out[row] = v;
  }
  return out;
}

inline std::array<double, 12> hydrodynamicForceTorque(
    const PairBlocks& blocks,
    const std::array<double, 12>& Ustar,
    const std::array<double, 9>& E) {
  std::array<double, 12> out = contractPhi(blocks, E);
  for (int i = 0; i < 12; ++i) {
    for (int j = 0; j < 12; ++j) {
      out[i] += blocks.R[static_cast<std::size_t>(i) * 12 + j] * Ustar[j];
    }
  }
  return out;
}

}  // namespace surge
