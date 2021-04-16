// nsgaii.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <boost/mpi.hpp>
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/filesystem.hpp>
//#include <random>
//#include <chrono>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/timer/timer.hpp>
#include <boost/chrono.hpp>

#include "nsgaii_options.h"
#include "nsgaii_core.h"
#include "nsgaii_helper.h"
#include "nsgaii_testFunctions.h"
#include "c2vsim_io.h"
#include "c2vsim_ga.h"

NSGAII::SingletonRealGenerator* NSGAII::SingletonRealGenerator::_instance = nullptr;

int main(int argc, char* argv[])
{

	boost::mpi::environment env(argc, argv);
	boost::mpi::communicator world;
	
	NSGAII::SingletonRealGenerator *RG = RG->getInstance();
	RG->printSeed();


	NSGAII::options opt;
	bool bCheckInputs = false;
	C2VSIM::OPTIONS::options cvopt;
	if (world.rank() == 0) {
		std::cout << "Current path:" << std::endl;
		std::cout << boost::filesystem::current_path() << std::endl;
		bCheckInputs = NSGAII::readInputParameters(argc, argv, opt);
	}
	boost::mpi::broadcast(world, bCheckInputs, 0);
	if (!bCheckInputs)
		return 0;
	boost::mpi::broadcast(world, opt.Nobjectives, 0);
	boost::mpi::broadcast(world, opt.MaxGenerations, 0);
	boost::mpi::broadcast(world, opt.bUseModel, 0);

	if (opt.bUseModel) {
		bCheckInputs = C2VSIM::OPTIONS::readConfigFile(argc, argv, cvopt);
	}

	
	// Start timing
	//auto start = std::chrono::high_resolution_clock::now();
	boost::timer::cpu_timer timer;

	// Generate and distribute initial population
	NSGAII::Population pop(opt);
	C2VSIM::c2vsimData CVD(cvopt);
	if (opt.bUseModel)
		CVD.readInputFiles();
	world.barrier();
	//std::cout << "Processor " << world.rank() << " has " << cvopt.Nsteps << std::endl;
	//CVD.debugMsg();
	
	
	if (world.rank() == 0) {
		pop.initializePopulation();
		if (opt.bWriteHistory)
			pop.openHistoryFile();
		//pop.printPopulation();
	}
	world.barrier();

	int currentGeneration = 0;
	while (currentGeneration < opt.MaxGenerations) {
		if (world.rank() == 0){
            std::cout << "Generation: " << currentGeneration << ", # Pareto Solutions: " << pop.ParetoSize() << " (Tabu size: " << pop.TabuSize() << ")" << std::endl;
        }

		pop.broadcast(world);
		//if (world.rank() == 3) {
		//	std::cout << "Processor " << world.rank() << " has " << pop.population.size() << std::endl;
		//}

		// Evaluate solutions
		std::map<int, std::vector<double> > solutions;
		std::map<int, NSGAII::Individual>::iterator itind;
		if (pop.population.size() % world.size() != 0){
		    std::cout << "The population size must be divided by the number of processors" << std::endl;
		    return 0;
		}

        int Nloops = pop.population.size() / world.size();
        std::cout << "Nloops " << Nloops << std::endl;

		for (int i = 0; i < Nloops; ++i){
            world.barrier();
		    int solutionID = world.rank() + i*world.size();
		    std::cout << "Rank " <<  world.rank() << " evaluates " << solutionID << std::endl;

            std::vector<double> ObjectiveFunctionValues;
            itind = pop.population.find(solutionID);
            if (itind != pop.population.end()){
                //NSGAII::Kursawe(itind->second.decisionVariables, ObjectiveFunctionValues);
                C2VSIM::OF::maxGWSTminCost(itind->second.decisionVariables, ObjectiveFunctionValues, CVD, world.rank());
                ObjectiveFunctionValues[0] += world.rank();
                ObjectiveFunctionValues[1] += world.rank();
                std::cout <<  "Rank " <<  world.rank() << " OF: [" << ObjectiveFunctionValues[0] << "," << ObjectiveFunctionValues[1] << "]" << std::endl;
                //C2VSIM::OF::maxWTminArea(itind->second.decisionVariables, ObjectiveFunctionValues, CVD);
                //solutions.insert(std::pair<int, std::vector<double > >(solutionID, ObjectiveFunctionValues));
            }
            else{
                std::cerr << "Rank " << world.rank() << " didnt find solution with id " << solutionID << std::endl;
            }
            world.barrier();


            // processor 0 will gather the solutions from all processors
            if (world.rank() == 0){
                std::vector<int> ids;
                // Receive the solution ids from the other processors
                boost::mpi::gather(world, solutionID, ids, 0);

                // Receive the objective function values
                std::vector<std::vector<double> > allF;
                for (int j = 0; j < opt.Nobjectives; ++j){
                    std::vector<double> f;
                    boost::mpi::gather(world, ObjectiveFunctionValues[j], f, 0);
                    allF.push_back(f);
                }
                for (int j = 0; j < world.size(); ++j){
                    itind = pop.population.find(ids[j]);
                    if (itind != pop.population.end()){
                        std::vector<double> objf;
                        for (int k = 0; k < opt.Nobjectives; ++k) {
                            itind->second.objectiveFunctions.push_back(allF[k][j]);
                        }
                    }
                    else{
                        std::cerr << "I can't find the solution id " << ids[j] << " in population" << std::endl;
                    }
                }

            }
            else{
                // Send the solution id to processor 0
                boost::mpi::gather(world, solutionID, 0);

                // Send th objective function values
                for (int j = 0; j < opt.Nobjectives; ++j){
                    boost::mpi::gather(world, ObjectiveFunctionValues[j], 0);
                }
            }
		}
        world.barrier();
		if (world.rank() == 0){
            for (itind  = pop.population.begin(); itind != pop.population.end(); ++itind){
                std::cout << itind->second.objectiveFunctions[0] << ", " << itind->second.objectiveFunctions[0] << std::endl;
            }
		}
        world.barrier();
        return 0;

        /*
		for (unsigned int i = world.rank(); i < pop.population.size(); i = i + world.size()) {
			itind = pop.population.find(i);
			if (itind != pop.population.end()) {
				std::vector<double> ObjectiveFunctionValues;
				//NSGAII::Kursawe(itind->second.decisionVariables, ObjectiveFunctionValues);
				C2VSIM::OF::maxGWSTminCost(itind->second.decisionVariables, ObjectiveFunctionValues, CVD, world.rank());
				//C2VSIM::OF::maxWTminArea(itind->second.decisionVariables, ObjectiveFunctionValues, CVD);
				solutions.insert(std::pair<int, std::vector<double > >(itind->first, ObjectiveFunctionValues));
			}
		}

		world.barrier();
		

		// Find out the maximum number of solutions that a processor has to send
		int MaxSolutions = 0;
		if (world.rank() == 0) {
			std::vector<int> SolutionsPerProcessor;
			boost::mpi::gather(world, static_cast<int>(solutions.size()), SolutionsPerProcessor, 0);

			for (unsigned int i = 0; i < SolutionsPerProcessor.size(); ++i) {
				if (SolutionsPerProcessor[i] > MaxSolutions) {
					MaxSolutions = SolutionsPerProcessor[i];
				}
			}
			//std::cout << "There are maximum " << MaxSolutions << std::endl;
		}
		else {
			boost::mpi::gather(world, static_cast<int>(solutions.size()), 0);
		}
		// Broadcast the MaxSolutions to all processors
		boost::mpi::broadcast(world, MaxSolutions, 0);
		//std::cout << "There are maximum " << MaxSolutions << std::endl;

		// gather the solutions from all processors
		std::map<int, std::vector<double> >::iterator itsol = solutions.begin();
		for (int i = 0; i < MaxSolutions; ++i) {
			if (world.rank() == 0) {
				// get the solution id
				std::vector<int> ids;
				boost::mpi::gather(world, itsol->first, ids, 0);

				//for (unsigned int k = 0; k < ids.size(); ++k)
				//	std::cout << ids[k] << " ";
				//std::cout << std::endl;


				std::vector<std::vector<double> > allF;
				for (int j = 0; j < opt.Nobjectives; ++j) {
					std::vector<double> f;
					boost::mpi::gather(world, itsol->second[j], f, 0);
					//for (unsigned int k = 0; k < f.size(); ++k)
					//	std::cout << f[k] << " ";
					//std::cout << std::endl;

					allF.push_back(f);
				}

				for (int j = 0; j < world.size(); ++j) {
					if (ids[j] < 0)
						continue;

					//std::cout << ids[j] << " ";
					itind = pop.population.find(ids[j]);
					if (itind != pop.population.end()) {
						std::vector<double> objf;
						for (int k = 0; k < opt.Nobjectives; ++k) {
							itind->second.objectiveFunctions.push_back(allF[k][j]);
						}
					}
				}
			}
			else {
				int id = -9;
				double f = 0;
				if (itsol != solutions.end()) {
					id = itsol->first;
				}
				boost::mpi::gather(world, id, 0);
				for (int j = 0; j < opt.Nobjectives; ++j) {
					if (itsol != solutions.end()) {
						f = itsol->second[j];
					}
					boost::mpi::gather(world, f, 0);
				}
			}
			if (itsol != solutions.end())
				++itsol;
		}
		*/

		world.barrier();
		

		if (world.rank() == 0) {
			//pop.printPopulation();
			//pop.printObjectives();
			pop.nonDominatingSorting();
			pop.selectFitest();
			pop.selectParents();
			pop.crossOver();
			pop.mutation();

			if (opt.bWriteHistory)
				pop.printCurrentPareto(currentGeneration);

			if (opt.bPrintRestartFile)
			    pop.printResumeFile();
		}
		else {
			pop.population.clear();
		}
		world.barrier();
		currentGeneration++;
	}

	boost::chrono::duration<double> seconds = boost::chrono::nanoseconds(timer.elapsed().user);
	//std::cout << seconds.count() << std::endl;
	//auto finish = std::chrono::high_resolution_clock::now();
	//std::chrono::duration<double> elapsed = finish - start;
	if (world.rank() == 0) {
		std::cout << "Optimization executed in " << seconds.count() << std::endl;
		if (opt.bWriteHistory)
			pop.closeHistoryFile();
	}
	pop.printPareto();
	return 0;
}
