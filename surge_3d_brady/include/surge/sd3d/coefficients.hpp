#pragma once

#include <string>
#include <unordered_map>
#include <vector>

namespace surge::sd3d {

enum class ScalarFamily { XA, YA, YB, XC, YC, XG, YG, YH };

struct CoefficientTable {
  // Optional high-order coefficient store. Keys are family names like "XA".
  // For now this is an extension point: production agreement needs recurrence-
  // generated or tabulated coefficients beyond the printed f_0...f_11.
  std::unordered_map<std::string, std::vector<double>> values;

  bool has(const std::string& family) const;
  double get(const std::string& family, int m) const;
};

CoefficientTable loadCoefficientCsv(const std::string& path);

}  // namespace surge::sd3d
