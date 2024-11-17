#include "particle/particle.hpp"
#include "particle/resonanceType.hpp"

void test() {

	// ParticleType and ResonanceType minimal tests
	const int particleN {4};

	const char* names[particleN] {
		"freaky",
		"sussy",
		"boring"
	};

	double masses[particleN] {
		-1.,
		3.14,
		1.,
	};

	int charges[particleN] {
		-1,
		2,
		5
	};

	double width {1.618};

	ParticleType* susPart = new ParticleType("bruh", 1., -1);

	ParticleType* particles[particleN + 1] {susPart, susPart, susPart};

	for (int i{}; i < particleN - 1; ++i) {
		ParticleType* particle = new ParticleType(names[i], masses[i], charges[i]);
		particles[i] = particle;
		Particle::addParticleType(particle->getName(), particle->getMass(), particle->getCharge());
	}
	{
		ResonanceType* resonance = new ResonanceType("sussonance", 2., 1., width);
		particles[3] = resonance;
		Particle::addParticleType(resonance->getName(), resonance->getMass(), resonance->getCharge(), resonance->getWidth());
	}

	for (int i{}; i < particleN; ++i) {
		particles[i]->print();
	}

	Particle::printTypes();

	// Particle minimal tests
}
