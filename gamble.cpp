void gamble() {
  std::string amogus{""};
  while (amogus != "no") {
    simulate();
    analyze();
    std::cout << "Again? ";
    std::cin >> amogus;
    std::cout << std::endl;
  }
}
