#ifndef PARTICLETYPE_HPP
#define PARTICLETYPE_HPP

#include <iostream>

class ParticleType {  // Particle type representation class
 public:
  // getters
  virtual const char* getName() const { return name_; };
  virtual double getMass() const { return mass_; };
  virtual int getCharge() const { return charge_; };
  virtual double getWidth() const { return 0.; }

  virtual void print() const;

  ParticleType(const char* name, const double mass, const int charge);

 private:
  const char* name_;
  const double mass_;
  const int charge_;
};

#endif
