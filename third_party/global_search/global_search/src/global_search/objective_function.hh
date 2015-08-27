#ifndef GLOBAL_SEARCH__OBJECTIVE_FUNCTION_H_
#define GLOBAL_SEARCH__OBJECTIVE_FUNCTION_H_
#include <vector>

namespace globalSearch {

class ObjectiveFunction {
	public: // Methods
		virtual double call(const std::vector<double>& values) = 0;
};

} // global_search
#endif // GLOBAL_SEARCH__GENETIC_ALGORITHM_H_
