set(SRC
  src/global_search/genetic_algorithm.hh
  src/global_search/genetic_algorithm.cc
  src/global_search/objective_function.hh
)

#if(CMAKE_COMPILER_IS_GNUCXX)
#  # Suppress warnings about enum comparisons (and everything else)
#  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -w")
#endif()

add_library(global_search STATIC ${SRC})
target_include_directories(global_search PUBLIC "src")

set_property(TARGET global_search PROPERTY CXX_STANDARD 11)
set_property(TARGET global_search PROPERTY CXX_STANDARD_REQUIRED ON)

