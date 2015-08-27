#include <iostream>
#include <iomanip>
#include <vector>
#include "global_search/genetic_algorithm.hh"

int main(int argc, char *argv[])
{
	auto ga = globalSearch::GeneticAlgorithm();
	auto out = ga.run();
	std::cout << "\n\nFittest: " << std::setprecision(6) << out.fitness << std::endl;
	std::cout << "[";
	int idx(0);
	int lastIdx = out.genes.size() - 1;
	for (const auto &gene : out.genes) {
		std::cout << std::setw(10) << gene;
		if (idx == lastIdx) {
			std::cout << "]\n";
		} else {
			std::cout << ", ";
		}
		++idx;
	}
	return 0;
}
