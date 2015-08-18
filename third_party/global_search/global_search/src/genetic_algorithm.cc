#include "genetic_algorithm.hh"
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <ctime>
#include <cstring>

using namespace std;

global_search::GeneticAlgorithm::GeneticAlgorithm(int popsize, int maxgens)
{
  this->population_size = popsize;
  this->maximum_generations = maxgens;
}

double global_search::GeneticAlgorithm::run()
{
  return 0.0;
}
