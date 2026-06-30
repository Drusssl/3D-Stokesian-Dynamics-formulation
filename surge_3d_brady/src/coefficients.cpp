#include "surge/sd3d/coefficients.hpp"

#include <fstream>
#include <sstream>
#include <stdexcept>

namespace surge::sd3d {

bool CoefficientTable::has(const std::string& family) const {
  return values.find(family) != values.end();
}

double CoefficientTable::get(const std::string& family, int m) const {
  const auto it = values.find(family);
  if (it == values.end() || m < 0 || static_cast<std::size_t>(m) >= it->second.size()) {
    throw std::out_of_range("coefficient not available");
  }
  return it->second[static_cast<std::size_t>(m)];
}

CoefficientTable loadCoefficientCsv(const std::string& path) {
  std::ifstream in(path);
  if (!in) throw std::runtime_error("could not open coefficient CSV");

  CoefficientTable table;
  std::string line;
  while (std::getline(in, line)) {
    if (line.empty() || line[0] == '#') continue;
    std::stringstream ss(line);
    std::string family;
    std::string mText;
    std::string valueText;
    if (!std::getline(ss, family, ',')) continue;
    if (!std::getline(ss, mText, ',')) continue;
    if (!std::getline(ss, valueText, ',')) continue;
    const int m = std::stoi(mText);
    const double value = std::stod(valueText);
    auto& coeffs = table.values[family];
    if (m >= static_cast<int>(coeffs.size())) coeffs.resize(static_cast<std::size_t>(m + 1), 0.0);
    coeffs[static_cast<std::size_t>(m)] = value;
  }
  return table;
}

}  // namespace surge::sd3d
