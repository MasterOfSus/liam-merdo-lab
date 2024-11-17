#include "particle.hpp"

// setting static variable
int Particle::nParticleTypes_ = 0;
ParticleType* Particle::particleTypes_[]; 

int Particle::findParticle(const char* name) {
	for (int i{}; i < nParticleTypes_; ++i) {
		if (particleTypes_[i]->getName() == name) return i;
	}
	return -1;
}

Particle::Particle(const char* name, Momentum momentum) : index_(findParticle(name)), momentum_(momentum) {
	if (findParticle(name) == -1) std::cout << "Invalid name. Initialized particle as default particle with index -1.\n";
}

void Particle::addParticleType(const char* name, double mass, int charge, double width) {
	int index = findParticle(name);
	if (nParticleTypes_ < maxNumParticleTypes_) {
		if (index != -1) std::cout << "Particle with name \"" << name << "\" has already been added at index " << index << ".\n";
		else {
			if (width) {
				particleTypes_[nParticleTypes_] = new ResonanceType(name, mass, charge, width);
				nParticleTypes_++;
				std::cout << "Added resonance type with name " << name << " at index " << nParticleTypes_ - 1 << ".\n";
			} else {
				particleTypes_[nParticleTypes_] = new ParticleType(name, mass, charge);
				nParticleTypes_++;
				std::cout << "Added particle type with name " << name << " at index " << nParticleTypes_ - 1 << ".\n";
			}
		}
	} else std::cout << "Maximum number of particles reached.\n";
}

void Particle::setIndex(const int index) {
	if (0 <= index && index < nParticleTypes_) {
		index_ = index;
	} else {
		std::cout << "Invalid index provided: " << index << std::endl;
	}
}

void Particle::setIndex(const char* name) {
	setIndex(findParticle(name));
}

void Particle::printTypes() {
	std::cout << "Available particle types:\n";
	int i {};
	for (ParticleType* particleType: particleTypes_) {
		std::cout << "Particle at index " << i << " attributes:\n";
		particleType->print();
		++i;
		if (i >= nParticleTypes_) break;
	}
}

double invMass(const Particle& p1, const Particle& p2) {
	double totEnergy { p1.getEnergy() + p2.getEnergy() };
	Momentum totP { p1.getMomentum() + p2.getMomentum() };
	return sqrt(totEnergy*totEnergy - totP.norm2());
}

int Particle::decayToBody(Particle &dau1,Particle &dau2) const {
  if(getMass() == 0.0){
    printf("Decayment cannot be preformed if mass is zero\n");
    return 1;
  }
  
  double massMot = getMass();
  double massDau1 = dau1.getMass();
  double massDau2 = dau2.getMass();

  if(index_ > -1){ // add width effect

    // gaussian random numbers

    float x1, x2, w, y1;
    
    double invnum = 1./RAND_MAX;
    do {
      x1 = 2.0 * rand()*invnum - 1.0;
      x2 = 2.0 * rand()*invnum - 1.0;
      w = x1 * x1 + x2 * x2;
    } while ( w >= 1.0 );
    
    w = sqrt( (-2.0 * log( w ) ) / w );
    y1 = x1 * w;

    massMot += particleTypes_[index_]->getWidth() * y1;

  }

  if(massMot < massDau1 + massDau2){
    printf("Decayment cannot be preformed because mass is too low in this channel\n");
    return 2;
  }
  
  double pout = sqrt((massMot*massMot - (massDau1+massDau2)*(massDau1+massDau2))*(massMot*massMot - (massDau1-massDau2)*(massDau1-massDau2)))/massMot*0.5;

  double norm = 2*M_PI/RAND_MAX;

  double phi = rand()*norm;
  double theta = rand()*norm*0.5 - M_PI/2.;
  dau1.setMomentum({pout*sin(theta)*cos(phi),pout*sin(theta)*sin(phi),pout*cos(theta)});
  dau2.setMomentum({-pout*sin(theta)*cos(phi),-pout*sin(theta)*sin(phi),-pout*cos(theta)});

  double energy = sqrt(momentum_.norm2() + massMot*massMot);

  double bx = getPx()/energy;
  double by = getPx()/energy;
  double bz = getPz()/energy;

  dau1.Boost(bx,by,bz);
  dau2.Boost(bx,by,bz);

  return 0;
}

void Particle::Boost(double bx, double by, double bz)
{
  double energy = getEnergy();

  //Boost this Lorentz vector
  double b2 = bx*bx + by*by + bz*bz;
  double gamma = 1.0 / sqrt(1.0 - b2);
  double bp = bx*getPx() + by*getPy() + bz*getPz();
  double gamma2 = b2 > 0 ? (gamma - 1.0)/b2 : 0.0;

  momentum_.px += gamma2*bp*bx + gamma*bx*energy;
  momentum_.py += gamma2*bp*by + gamma*by*energy;
  momentum_.pz += gamma2*bp*bz + gamma*bz*energy;
}

int tailIndex(const Particle* array, const int arraySize,int start) {
	for (int i {start - 1}; i < arraySize; ++i) {
		if (array[i].getIndex() == -1) {
			return i;
		}
	}
	return -1;
}
