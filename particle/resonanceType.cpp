#include "resonanceType.hpp"

void ResonanceType::print() const {
  ParticleType::print();
  std::cout << "Particle width: " << width_ << std::endl;
}

ResonanceType::ResonanceType(const char* name, const double mass,
                             const int charge, const double width)
    : ParticleType(name, mass, charge), width_(width) {
  if (width < 0) {
    std::cout << "Warning! Negative particle width provided.\n";
  };
}
