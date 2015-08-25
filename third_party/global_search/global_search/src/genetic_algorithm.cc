#include "genetic_algorithm.hh"
#include <algorithm>
#include <assert.h>
#include <numeric>
#include <cstdlib>
#include <chrono>
#include <random>
#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <ctime>
#include <cstring>
#include <array>

double
globalSearch::ExampleObjectiveFunction::call(const std::vector<double>& xs)
{
	double x0 = xs[0];
	double x1 = xs[1];
	double x2 = xs[2];
	return (x0 * x0) - (x0 * x1) + x2;
}

globalSearch::GeneticAlgorithm::GeneticAlgorithm(
	globalSearch::ObjectiveFunction* objFn,
	std::vector<double> defaultGenes,
	std::vector<bool> continuousFlags,
	std::vector<bool> variableFlags,
	std::vector<double> upperBounds,
	std::vector<double> lowerBounds,
	std::size_t populationSize,
	std::size_t maximumGenerations,
	double probabilityOfCrossover,
	double probabilityOfMutation
) :
	objFn_(objFn),
	defaultGenes_(defaultGenes),
	continuousFlags_(continuousFlags),
	variableFlags_(variableFlags),
	upperBounds_(upperBounds),
	lowerBounds_(lowerBounds),
	populationSize_(populationSize),
	maximumGenerations_(maximumGenerations),
	probabilityOfCrossover_(probabilityOfCrossover),
	probabilityOfMutation_(probabilityOfMutation),
	rng_()
{
	numberOfGenes_ = defaultGenes_.size();
	assert(numberOfGenes_ == continuousFlags_.size());
	assert(numberOfGenes_ == variableFlags_.size());
	assert(numberOfGenes_ == upperBounds_.size());
	assert(numberOfGenes_ == lowerBounds_.size());
}

globalSearch::GeneticAlgorithm::GeneticAlgorithm() :
	globalSearch::GeneticAlgorithm::GeneticAlgorithm{
		new globalSearch::ExampleObjectiveFunction(),
		{2.5, 2.5, 0.0},
		{true, true, true},
		{true, true, true},
		{5.0, 5.0, 2.0},
		{0.0, 0.0, -2.0},
		50,
		1000,
		0.8,
		0.15}
{}

globalSearch::Individual
globalSearch::GeneticAlgorithm::run()
{
	startReport();
	std::vector<globalSearch::Individual> population =
		generateInitialPopulation();
	std::vector<globalSearch::Individual> nextPopulation;
	globalSearch::Individual fittest;
	globalSearch::Individual globalFittest;
	evaluatePopulation(population);
	fittest = fittestMemberOf(population);
	globalFittest = fittest;
	std::size_t generation;
	for (generation=0;
			generation < maximumGenerations_;
		 	generation++) {
		reportProgress(generation, population, fittest);
		nextPopulation = selectSurvivors(population);
		assert(populationSize_ == nextPopulation.size());
		nextPopulation = crossover(nextPopulation);
		assert(populationSize_ == nextPopulation.size());
		mutatePopulation(nextPopulation);
		assert(populationSize_ == nextPopulation.size());
		population = copyPopulation(nextPopulation);
		assert(populationSize_ == population.size());
		evaluatePopulation(population);
		assert(populationSize_ == population.size());
		insertElite(population, fittest);
		fittest = fittestMemberOf(population);
		if (fittest.fitness > globalFittest.fitness) {
			globalFittest = fittest;
		}
		assert(populationSize_ == population.size());
	}
	return globalFittest;
}

double
globalSearch::GeneticAlgorithm::randomRange(
	double lower,
	double upper
)
{
	std::uniform_real_distribution<double> distribution(lower, upper);
	return distribution(rng_);
}

double
globalSearch::GeneticAlgorithm::nthRandomGene(int n) {
	double offset = 0.49999;
	bool continuous = continuousFlags_[n];
 	bool variable = variableFlags_[n];
	double newValue, lowerBound, upperBound;
	if (variable) {
		lowerBound = lowerBounds_[n];
		upperBound = upperBounds_[n];
		if (continuous) {
			newValue = randomRange(lowerBound, upperBound);
		} else {
			newValue = std::round(
					randomRange(lowerBound - offset, upperBound + offset));
		}
	} else {
		newValue = defaultGenes_[n];
	}
	return newValue;
}

std::vector<globalSearch::Individual>
globalSearch::GeneticAlgorithm::generateInitialPopulation()
{
	std::vector<globalSearch::Individual> inds;
	std::size_t i;
	for (i=0; i < populationSize_; i++) {
		inds.push_back(initIndividual());
	}
	return inds;
}

globalSearch::Individual
globalSearch::GeneticAlgorithm::initIndividual()
{
	std::vector<double> vars;
	std::size_t i;
	for (i = 0; i < numberOfGenes_; ++i) {
		vars.push_back(nthRandomGene(i));
	}
	globalSearch::Individual indy;
	indy.genes = vars;
	indy.fitness = 0.0;
	indy.relativeFitness = 0.0;
	indy.cumulativeFitness = 0.0;
	return indy;
}

void
globalSearch::GeneticAlgorithm::evaluatePopulation(std::vector<globalSearch::Individual>& pop)
{
	for (auto &member : pop) {
		member.fitness = objFn_->call(member.genes);
	}
	return;
}

globalSearch::Individual
globalSearch::GeneticAlgorithm::fittestMemberOf(const std::vector<globalSearch::Individual>& pop)
{
	double bestFitnessThusFar;
	globalSearch::Individual fittest;
	bool first = true;
	for (const auto &member : pop) {
		if (first) {
			fittest = member;
			bestFitnessThusFar = fittest.fitness;
			first = false;
			continue;
		}
		if (member.fitness > bestFitnessThusFar) {
			fittest = member;
			bestFitnessThusFar = fittest.fitness;
		}
	}
	return fittest;
}

std::vector<globalSearch::Individual>
globalSearch::GeneticAlgorithm::selectSurvivors(std::vector<globalSearch::Individual>& pop)
{
	std::vector<globalSearch::Individual> newPop;
	std::size_t i, j;
	double p, cf0, cfj, cfjj;
	std::size_t lastIdx = populationSize_ - 1;
	double cumulativeFitness = 0.0;
	double totalFitnessOfPopulation = 0.0;
	for (const auto &member : pop) {
		totalFitnessOfPopulation += member.fitness;
	}
	for (auto &member : pop) {
		member.relativeFitness = member.fitness / totalFitnessOfPopulation;
	}
	for (auto &member : pop) {
		cumulativeFitness += member.relativeFitness;
		member.cumulativeFitness = cumulativeFitness;
	}
	// Roulette Wheel Selection for the New Population
	cf0 = pop[0].cumulativeFitness;
	for (i = 0; i < populationSize_; ++i) {
		p = randomRange(0.0, 1.0);
		if (p < cf0) {
			newPop.push_back(pop[0]);
		} else {
			for (j = 0; j < populationSize_; ++j) {
				if (j == lastIdx) {
					newPop.push_back(pop[j]);
					break;
				} else {
					cfj = pop[j].cumulativeFitness;
					cfjj = pop[j+1].cumulativeFitness;
					if ((p >= cfj) && (p < cfjj)) {
						newPop.push_back(pop[j+1]);
						break;
					}
				}
			}
		}
	}
	return newPop;
}

std::vector<globalSearch::Individual>
globalSearch::GeneticAlgorithm::crossover(const std::vector<globalSearch::Individual>& pop)
{
	globalSearch::Individual mother;
	globalSearch::Individual father;
	std::array<globalSearch::Individual,2> children;
	bool has_mother = false;
	double p;
	std::size_t idx = 0;
	std::size_t lastIdx = populationSize_ - 1;
	std::vector<globalSearch::Individual> nextPop;
	for (const auto &member : pop) {
		p = randomRange(0.0, 1.0);
		if (((p < probabilityOfCrossover_) && (idx != lastIdx)) ||
				(has_mother && idx == lastIdx)) {
			if (has_mother) {
				father = member;
				children = makeChildren(mother, father);
				nextPop.push_back(children[0]);
				nextPop.push_back(children[1]);
				has_mother = false;
			} else {
				mother = member;
				has_mother = true;
			}
		} else {
			nextPop.push_back(member);
		}
		++idx;
	}
	return nextPop;
}

// Note that it is irrelevant if the genes are variable or not
// as genes that are not variable (i.e., constant) will have
// the same value for mother and father. Similarly, continuous
// or discrete is also irrelevant for making children
std::array<globalSearch::Individual,2>
globalSearch::GeneticAlgorithm::makeChildren(
		const globalSearch::Individual& mom,
		const globalSearch::Individual& dad
)
{
	std::array<globalSearch::Individual,2> children;
	double v1, v2, p;
	auto first = globalSearch::Individual();
	auto second = globalSearch::Individual();
	std::size_t gIdx;
	for (gIdx = 0; gIdx < numberOfGenes_; gIdx++) {
		v1 = mom.genes[gIdx];
		v2 = dad.genes[gIdx];
		p = randomRange(0.0, 1.0);
		if (p < 0.5) {
			first.genes.push_back(v1);
			second.genes.push_back(v2);
		} else {
			first.genes.push_back(v2);
			second.genes.push_back(v1);
		}
	}
	children[0] = first;
	children[1] = second;
	return children;
}

void
globalSearch::GeneticAlgorithm::mutatePopulation(std::vector<globalSearch::Individual>& pop)
{
	double p;
	std::size_t geneIdx;
	for (auto &member : pop) {
		for (geneIdx = 0; geneIdx < numberOfGenes_; ++geneIdx) {
			p = randomRange(0.0, 1.0);
			if (p < probabilityOfMutation_) {
				member.genes[geneIdx] = nthRandomGene(geneIdx);
			}
		}
	}
}

void
globalSearch::GeneticAlgorithm::startReport()
{
	std::size_t i;
	std::cout << "#globalSearch/GeneticAlgorithm\n";
	std::cout << "{:population-size " << populationSize_ << "\n";
	std::cout << " :maximum-generations " << maximumGenerations_ << "\n";
	std::cout << " :probability-of-crossover " << std::setprecision(3) << probabilityOfCrossover_ << "\n";
	std::cout << " :probability-of-mutation " << std::setprecision(3) << probabilityOfMutation_ << "\n";
	std::cout << " :number-of-genes " << numberOfGenes_ << "\n";
	std::cout << " :genes\n";
	for (i = 0; i < numberOfGenes_; ++i) {
		if (i == 0) {
			std::cout << " [{idx " << i << ", ";
		} else {
			std::cout << "  {idx " << i << ", ";
		}
		std::cout << "default " << defaultGenes_[i] << ", ";
		std::cout << "continuous? " << continuousFlags_[i] << ", ";
		std::cout << "variable? " << variableFlags_[i] << ", ";
		std::cout << "bounds [" << lowerBounds_[i] << ", " << upperBounds_[i] << "]}";
		if (i == numberOfGenes_ - 1) {
			std::cout << "]}\n";
		} else {
			std::cout << "\n";
		}
	}
	std::cout << "\n";
	// line 1
	std::cout << std::setw(14) << "Generation No.";
	std::cout << "  ";
	std::cout << std::setw(14) << "Best Fitness";
	std::cout << "  ";
	std::cout << std::setw(14) << "Pop Fitness";
	std::cout << "  ";
	std::cout << std::setw(14) << "Avg Fitness";
	std::cout << "  ";
	std::cout << std::setw(14) << "Fittest Genes\n";
	// line 2
	std::cout << std::setfill('-') << std::setw(14) << "";
	std::cout << "  ";
	std::cout << std::setfill('-') << std::setw(14) << "";
	std::cout << "  ";
	std::cout << std::setfill('-') << std::setw(14) << "";
	std::cout << "  ";
	std::cout << std::setfill('-') << std::setw(14) << "";
	std::cout << "  ";
	std::cout << std::setfill('-') << std::setw(14) << "" << "\n";
	std::cout << std::setfill(' ');
}

void
globalSearch::GeneticAlgorithm::reportProgress(
		int generationNumber,
		const std::vector<globalSearch::Individual>& pop,
		const globalSearch::Individual& fittest
)
{
	double totalFitness = std::accumulate(
			pop.begin(),
			pop.end(),
			0.0,
			[](double sum, globalSearch::Individual ind) { return sum + ind.fitness; });
	double avgFitness = totalFitness / ((double)populationSize_);
	double bestFitness = fittest.fitness;
	double g;
	std::size_t i;
	std::cout << std::setw(14) << std::setprecision(6) << generationNumber << "  ";
	std::cout << std::setw(14) << std::setprecision(6) << bestFitness << "  ";
	std::cout << std::setw(14) << std::setprecision(6) << totalFitness << "  ";
	std::cout << std::setw(14) << std::setprecision(6) << avgFitness << "  ";
	std::cout << std::setprecision(4) << "[";
	for (i=0; i < fittest.genes.size(); ++i) {
		g = fittest.genes[i];
		if (i == fittest.genes.size() - 1) {
			std::cout << std::setw(4) << g << "]\n";
		} else {
			std::cout << std::setw(4) << g << ", ";
		}
	}
}

std::vector<globalSearch::Individual>
globalSearch::GeneticAlgorithm::copyPopulation(const std::vector<globalSearch::Individual>& pop)
{
	std::vector<globalSearch::Individual> out;
	for (const auto &member : pop) {
		out.push_back(member);
	}
	return out;
}

void
globalSearch::GeneticAlgorithm::insertElite(
		std::vector<globalSearch::Individual>& population,
		const globalSearch::Individual& fittest)
{
	int bestMemberIdx = 0;
	int worstMemberIdx = 0;
  double bestFitness = population[bestMemberIdx].fitness;
	double worstFitness = population[worstMemberIdx].fitness;
	globalSearch::Individual member, newFittest;
	std::size_t i;
  for (i = 0; i < population.size(); ++i) {
		member = population[i];
    if (member.fitness > bestFitness ) {
			bestFitness = member.fitness;
			bestMemberIdx = i;
		}
		if (member.fitness < worstFitness) {
			worstFitness = member.fitness;
			worstMemberIdx = i;
		}
  }
  if (bestFitness < fittest.fitness) {
    population[worstMemberIdx] = fittest;
  }
  return;
}
