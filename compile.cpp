void compile(const char option) {

	if (option == 's') {
		gROOT->LoadMacro("particle/particleType.cpp+");
		gROOT->LoadMacro("particle/resonanceType.cpp+");
		gROOT->LoadMacro("particle/particle.cpp+");
		gROOT->LoadMacro("simulate.cpp+");
	}
	else if (option == 't') {
		gROOT->LoadMacro("particle/particleType.cpp+");
		gROOT->LoadMacro("particle/resonanceType.cpp+");
		gROOT->LoadMacro("particle/particle.cpp+");
		gROOT->LoadMacro("test.cpp+");
	}
	else if (option == 'a') gROOT->LoadMacro("analyze.cpp+");
	else std::cout << "Invalid option.";
}

