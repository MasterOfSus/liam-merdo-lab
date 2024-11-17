#ifndef RESONANCETYPE_HPP
#define RESONANCETYPE_HPP

#include "particleType.hpp"

class ResonanceType : public ParticleType {
 public:
  double getWidth() const { return width_; };

  void print() const;

  ResonanceType(const char* name, const double mass, const int charge,
                const double width);

 private:
  const double width_;
};

#endif
