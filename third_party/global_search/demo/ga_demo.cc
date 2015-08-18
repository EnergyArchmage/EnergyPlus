#include <iostream>
#include "genetic_algorithm.hh"

int main(int argc, char *argv[])
{
  auto ga = global_search::GeneticAlgorithm(50, 100);
  std::cout << "Hello World!" << std::endl;
  return 0;
}
