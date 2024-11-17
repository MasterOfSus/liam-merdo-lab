#ifndef PARTICLE_HPP
#define PARTICLE_HPP

#include <cmath>
#include <cstdlib> 

#include "resonanceType.hpp"

struct Momentum {

double px;
double py;
double pz;

double norm2() const { return px*px + py*py + pz*pz; }
double norm() const { return sqrt(norm2()); }

Momentum operator+(const Momentum p) {
	return {px + p.px, py + p.py, pz+p.pz};
}

}; // struct momentum

class Particle {
public:
	Particle(const int index = -1, const Momentum momentum = {0., 0., 0.}) : index_(index), momentum_(momentum) {};
	Particle(const char* name, const Momentum momentum = {0., 0., 0.});

	static void addParticleType(const char* name, double mass, int charge, double width = 0.);
	
	int getIndex() const { return index_; }
	Momentum getMomentum() const { return momentum_; }
	int getCharge() const { return particleTypes_[getIndex()]->getCharge(); };
	double getPx() const { return momentum_.px; }
	double getPy() const { return momentum_.py; }
	double getPz() const { return momentum_.pz; }
	double getMass() const { return particleTypes_[index_]->getMass(); }
	double getEnergy() const { return sqrt(getMass()*getMass() + getMomentum().norm2()); }

	void setIndex(const int index);
	void setIndex(const char* name);
	void setMomentum(const Momentum p) { momentum_ = p; }

	static void printTypes();
	void print() const;

	int decayToBody(Particle &dau1,Particle &dau2) const;

private:
	void Boost(double bx, double by, double bz);
	static int findParticle(const char* name);
	static const int maxNumParticleTypes_{10};
	static ParticleType* particleTypes_[maxNumParticleTypes_];
	static int nParticleTypes_; 
	int index_;
	Momentum momentum_;
}; // class Particle

double invMass(const Particle& p1, const Particle& p2);

int tailIndex(const Particle* array, const int arraySize, const int start);

#endif
