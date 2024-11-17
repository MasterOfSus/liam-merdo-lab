void compile(const char option) {
	gROOT->LoadMacro("particle/particleType.cpp+");
	gROOT->LoadMacro("particle/resonanceType.cpp+");
	gROOT->LoadMacro("particle/particle.cpp+");
	if (option == 's') gROOT->LoadMacro("simulate.cpp+");
	else if (option == 't') gROOT->LoadMacro("test.cpp+");
	else if (option == 'a') gROOT->LoadMacro("analyze.cpp+");
	else std::cout << "Invalid option.";
}

void recompile(const char& option) {
	gROOT->LoadMacro("particle/particleType.cpp++");
	gROOT->LoadMacro("particle/resonanceType.cpp++");
	gROOT->LoadMacro("particle/particle.cpp++");
	if (option == 's') gROOT->LoadMacro("simulate.cpp++");
	else if (option == 't') gROOT->LoadMacro("test.cpp++");
	else if (option == 'a') gROOT->LoadMacro("analyze.cpp++");
	else std::cout << "Invalid option.";
}
