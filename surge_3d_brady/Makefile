CXX ?= c++
CXXFLAGS ?= -std=c++17 -O2 -Wall -Wextra -pedantic
EIGEN3_INCLUDE ?= /opt/homebrew/include/eigen3

INCLUDES = -Iinclude -I. -I.. -I$(EIGEN3_INCLUDE)
SRC = src/coefficients.cpp src/ewald.cpp src/near_field.cpp src/projection.cpp src/solver.cpp src/time_integrator.cpp

.PHONY: all clean

all: two_sphere_trajectory smoke_checks

two_sphere_trajectory: $(SRC) examples/two_sphere_trajectory.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(SRC) examples/two_sphere_trajectory.cpp -o two_sphere_trajectory

smoke_checks: $(SRC) examples/smoke_checks.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(SRC) examples/smoke_checks.cpp -o smoke_checks

clean:
	rm -f two_sphere_trajectory smoke_checks
