#include "particleType.hpp"

void ParticleType::print() const {
  std::cout << "Particle name: " << name_ << std::endl;
  std::cout << "Particle mass: " << mass_ << std::endl;
  std::cout << "Particle charge: " << charge_ << std::endl;
};

ParticleType::ParticleType(const char* name, const double mass,
                           const int charge)
    : name_(name), mass_(mass), charge_(charge) {
  if (mass < 0.) std::cout << "Warning! Negative particle mass provided.\n";
}
