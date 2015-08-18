#ifndef GLOBAL_SEARCH__GENETIC_ALGORITHM_H_
#define GLOBAL_SEARCH__GENETIC_ALGORITHM_H_

namespace global_search
{
class GeneticAlgorithm
{
    int population_size, maximum_generations;
  public:
    GeneticAlgorithm(
      int popsize,
      int maxgens);
    double run();
};
}
#endif // GLOBAL_SEARCH__GENETIC_ALGORITHM_H_
