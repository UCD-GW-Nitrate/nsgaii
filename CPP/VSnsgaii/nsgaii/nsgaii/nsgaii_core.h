#pragma once

#include <iostream>
//#include <random>
//#include <chrono>
#include <boost/mpi.hpp>

#include "nsgaii_options.h"
#include "nsgaii_helper.h"


namespace NSGAII {

	struct Rank {
		int rank = -1;
		double accum = -1;
	};

	class Individual {
	public:
		Individual() {};
		int compareObjectives(Individual other);
		int compareRanking(Individual other);
		std::vector<double> decisionVariables;
		std::vector<double> objectiveFunctions;
		Rank rank;
	};

	int Individual::compareObjectives(Individual other) {
		unsigned int win_this = 0;
		unsigned int win_other = 0;
		for (unsigned int i = 0; i < objectiveFunctions.size(); ++i) {
			if (objectiveFunctions[i] <= other.objectiveFunctions[i])
				win_this++;
			else
				win_other++;
		}
		if (win_this == objectiveFunctions.size())
			return 1;
		else if (win_other == objectiveFunctions.size())
			return -1;
		else
			return 0;
	}

	int Individual::compareRanking(Individual other) {
		if (rank.rank < other.rank.rank)
			return 1;
		else if (rank.rank > other.rank.rank)
			return 0;
		else if (rank.rank == other.rank.rank) {
			if (rank.accum < 0)
				return 1;
			if (other.rank.rank < 0)
				return 0;

			if (rank.accum > other.rank.accum)
				return 1;
			else if (rank.accum < other.rank.accum)
				return 0;
			else if (rank.accum == other.rank.accum)
				return 1;
		}
	}


	class Population {
	public:
		Population(NSGAII::options& opt);
		void initializePopulation();
		void broadcast(boost::mpi::communicator& world);

		void nonDominatingSorting();
		void selectFitest();
		void selectParents();
		void crossOver();
		void mutation();
		void printPareto();

		// Debug methods
		void printObjectives();
		void printPopulation();
		int ParetoSize();
		void openHistoryFile();
		void closeHistoryFile();
		void printCurrentPareto(int gen);

		NSGAII::options options;
		std::map<int, Individual> population;
		
	private:
		std::vector<std::vector<int> > Fronts;
		std::map<int, Individual> ParetoSolutions;
		std::vector<int> selectedIds;
		std::vector<std::pair<int, int> > selectedParentIds;
		void appendParetoSolutions();
		int tournament();
		std::list<std::vector<double> > tabuList;
		bool isTabu(std::vector<double>& v);
		std::ofstream historyfile;
	};

	Population::Population(NSGAII::options& opt)
		:
		options(opt)
	{}

	void Population::initializePopulation() {
		
		NSGAII::SingletonRealGenerator* RG = RG->getInstance();
		//RG->printSeed();
		//std::cout << "Initialize Population" << std::endl;
		//std::map<int, std::vector<double> >::iterator it;
		std::map<int, Individual>::iterator it;

		for (int i = 0; i < options.PopulationSize; ++i) {
			Individual individual;
			while (true) {
				individual.decisionVariables.clear();
				for (int ivar = 0; ivar < options.ProblemSize; ++ivar) {
					if (options.populationType.compare("REAL") == 0) {
						individual.decisionVariables.push_back(RG->randomNumber(options.LowerBound[ivar], options.UpperBound[ivar]));
					}
					else if (options.populationType.compare("BINARY") == 0)
					{	
						double test_r = static_cast<double>(i)/ options.PopulationSize;
						if (i == 0)
							test_r = -1;
						else if (i == options.PopulationSize - 1)
							test_r = 2;
						double r = RG->randomNumber();
						if (r > test_r)
							individual.decisionVariables.push_back(1);
						else
							individual.decisionVariables.push_back(0);
					}
				}
				if (options.useTabu) {
					if (!isTabu(individual.decisionVariables)) 
						break;
					else {
						// The following applies only for the custom initial population implementation.
						// If any of the previously generated individuals have all zero the code will
						// loop forever as the last individual is forced to have all zeros
						if (i == options.PopulationSize - 1) {
							break;
						}
					}
				}
				else
					break;
			}
			population.insert(std::pair<int, Individual>(i, individual));
			tabuList.push_front(individual.decisionVariables);
		}
	}

	void Population::broadcast(boost::mpi::communicator& world) {
		int PopulationSize;
		if (world.rank() == 0)
			PopulationSize = population.size();
		boost::mpi::broadcast(world, PopulationSize, 0);

		for (int i = 0; i < PopulationSize; ++i) {
			std::vector<double> temp;
			if (world.rank() == 0) {
				temp = population[i].decisionVariables;
			}
			broad_cast_vector<double>(temp, 0, world);
			if (world.rank() != 0) {
				Individual individual;
				individual.decisionVariables = temp;
				population.insert(std::pair<int, Individual>(i, individual));
			}
			world.barrier();
		}
	}

	void Population::nonDominatingSorting() {
		appendParetoSolutions();
		std::map<int, Individual>::iterator itp, itq;
		std::vector<std::vector<int> > Sp;
		std::vector<int> np;
		
		std::vector<int> tempF;
		std::vector<double> fmin(options.Nobjectives, 99999999999);
		std::vector<double> fmax(options.Nobjectives, -99999999999);
		for (itp = population.begin(); itp != population.end(); ++itp) {
			for (int iobj = 0; iobj < options.Nobjectives; ++iobj) {
				if (itp->second.objectiveFunctions[iobj] < fmin[iobj])
					fmin[iobj] = itp->second.objectiveFunctions[iobj];
				if (itp->second.objectiveFunctions[iobj] > fmax[iobj])
					fmax[iobj] = itp->second.objectiveFunctions[iobj];
			}

			std::vector<int> tempSp;
			int tempnp = 0;
			for (itq = population.begin(); itq != population.end(); ++itq) {
				if (itp->first == itq->first)
					continue;
				int result = itp->second.compareObjectives(itq->second);
				if (result == 1)
					tempSp.push_back(itq->first);
				else if (result == -1)
					tempnp++;
			}
			if (tempnp == 0) 
				tempF.push_back(itp->first);
			

			Sp.push_back(tempSp);
			np.push_back(tempnp);
		}
		Fronts.push_back(tempF);
		int cntF = 0;
		while (true) {
			std::vector<int> Q;
			for (unsigned int i = 0; i < Fronts[cntF].size(); ++i) {
				//std::cout << F[cntF][i] << std::endl;
				for (unsigned int j = 0; j < Sp[Fronts[cntF][i]].size(); ++j) {
					//std::cout << Sp[F[cntF][i]][j] << std::endl;
					np[Sp[Fronts[cntF][i]][j]] = np[Sp[Fronts[cntF][i]][j]] - 1;
					if (np[Sp[Fronts[cntF][i]][j]] == 0)
						Q.push_back(Sp[Fronts[cntF][i]][j]);
				}
			}
			if (Q.size() == 0)
				break;
			Fronts.push_back(Q);
			cntF++;
		}
		// Accumulator measure
		std::vector<std::vector<double>> A;
		for (unsigned int ifront = 0; ifront < Fronts.size(); ++ifront) {
			std::vector<double> a(Fronts[ifront].size(), 0);
			for (int iobj = 0; iobj < options.Nobjectives; ++iobj) {
				std::vector<std::pair<double,int> > f;
				double Df = fmax[iobj] - fmin[iobj];
				for (unsigned int i = 0; i < Fronts[ifront].size(); ++i) {
					itp = population.find(Fronts[ifront][i]);
					if (itp != population.end()) 
						f.push_back(std::pair<double, int> (itp->second.objectiveFunctions[iobj],i));
				}
				
				std::sort(f.begin(), f.end());
				for (unsigned int i = 0; i < f.size(); ++i) {
					if (i == 0 || i == f.size() - 1) {
						a[f[i].second] = -1.0;
					}
					else {
						a[f[i].second] += (f[i + 1].first - f[i - 1].first)/Df;
					}
				}
			}
			A.push_back(a);
		}
		// Set the Ranking and the accumulator value to population
		for (unsigned int ifront = 0; ifront < Fronts.size(); ++ifront) {
			for (unsigned int i = 0; i < Fronts[ifront].size(); ++i) {
				itp = population.find(Fronts[ifront][i]);
				if (itp != population.end()) {
					itp->second.rank.rank = ifront;
					itp->second.rank.accum = A[ifront][i];
				}
			}
		}

		// Set the solutions with 0 ranking as new pareto front
		for (unsigned int i = 0; i < Fronts[0].size(); ++i) {
			itp = population.find(Fronts[0][i]);
			if (itp != population.end()) {
				ParetoSolutions.insert(std::pair<int, Individual>(i, itp->second));
			}
		}
	}

	void Population::printObjectives() {
		std::map<int, Individual >::iterator it;
		for (it = population.begin(); it != population.end(); ++it) {
			for (unsigned int i = 0; i < it->second.objectiveFunctions.size(); ++i)
				std::cout << it->second.objectiveFunctions[i] << " ";
			std::cout << std::endl;
		}
	}

	void Population::printPopulation() {
		std::cout << "Population :--------------------------------------" << std::endl;
		std::map<int, Individual>::iterator it;
		for (it = population.begin(); it != population.end(); ++it) {
			for (unsigned int i = 0; i < it->second.decisionVariables.size(); ++i)
				std::cout << it->second.decisionVariables[i] << " ";
			std::cout << std::endl;
		}
		std::cout << "----------------------------------------------------" << std::endl;
	}

	void Population::selectFitest() {
		selectedIds.clear();
		std::map<int, Individual>::iterator it;
		for (unsigned int ifront = 0; ifront < Fronts.size(); ++ifront) {
			if (selectedIds.size() + Fronts[ifront].size() <= options.PopulationSize) {
				for (unsigned int i = 0; i < Fronts[ifront].size(); ++i) {
					selectedIds.push_back(Fronts[ifront][i]);
				}
			}
			else {
				std::vector<std::pair<double, int> >accum;
				for (unsigned int i = 0; i < Fronts[ifront].size(); ++i) {
					it = population.find(Fronts[ifront][i]);
					if (it != population.end()) {
						double x = it->second.rank.accum;
						if (x < 0)
							x = 100000000000;
						accum.push_back(std::pair<double, int>(x, i));
					}
				}
				std::sort(accum.begin(), accum.end());
				std::reverse(accum.begin(), accum.end());

				for (unsigned int i = 0; i < accum.size(); ++i) {
					selectedIds.push_back(Fronts[ifront][accum[i].second]);
					if (static_cast<int>(selectedIds.size()) >= options.PopulationSize)
						break;
				}
			}
			if (static_cast<int>(selectedIds.size()) >= options.PopulationSize)
				break;
		}
	}

	int Population::tournament() {
		std::vector<std::pair<int, int> > players;
		NSGAII::SingletonRealGenerator* RG = RG->getInstance();
		std::map<int, Individual>::iterator itb, it;
		std::list<int> tournament_list;

		for (int i = 0; i < options.TournamentSize; ++i) {
			while (true) {
				int id = RG->randomNumber(0, static_cast<int>(selectedIds.size() - 1));
				it = population.find(id);
				if (it != population.end()) {
					bool idexists = false;
					for (std::list<int>::iterator itl = tournament_list.begin(); itl != tournament_list.end(); ++itl) {
						if (*itl == it->first) {
							idexists = true;
							break;
						}
					}
					if (idexists)
						continue;

					tournament_list.push_back(it->first);
					if (i == 0) {
						itb = it;
						break;
					}
					else {
						int out = it->second.compareRanking(itb->second);
						if (out == 1) 
							itb = it;
						break;
					}
				}
			}
		}
		return itb->first;
	}

	void Population::selectParents() {
		//unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
		//std::default_random_engine generator(seed);
		//std::uniform_int_distribution<int> distribution(0, static_cast<int>(selectedIds.size()-1));
		std::map<int, Individual>::iterator ita, itb;

		for (int i = 0; i < options.PopulationSize; ++i) {
			int parentA = tournament();
			int parentB;
			while (true) {
				parentB = tournament();
				if (parentA != parentB)
					break;
			}
			selectedParentIds.push_back(std::pair<int, int>(parentA, parentB));
		}
	}

	void Population::crossOver() {
		NSGAII::SingletonRealGenerator* RG = RG->getInstance();
		//RG->printSeed();

		std::map<int, Individual>::iterator ita, itb;
		std::vector<std::vector<double> > newPopulation;

		for (unsigned ichild = 0; ichild < selectedParentIds.size(); ++ichild) {
			std::vector<double> child;
			ita = population.find(selectedParentIds[ichild].first);
			itb = population.find(selectedParentIds[ichild].second);
			if (ita == population.end() || itb == population.end())
				continue;

			if (options.XoverType.compare("ONECUT") == 0) {
				int cut = RG->randomNumber(1, options.ProblemSize-1);
				for (int ivar = 0; ivar < options.ProblemSize; ++ivar) {
					if (ivar <= cut)
						child.push_back(ita->second.decisionVariables[ivar]);
					else
						child.push_back(itb->second.decisionVariables[ivar]);
				}
			}
			else if (options.XoverType.compare("TWOCUT") == 0) {
				double one_third = static_cast<double>(options.ProblemSize) / 3.0;
				int cut1 = RG->randomNumber(1, static_cast<int>(2*one_third));
				int cut2 = RG->randomNumber(cut1+1, options.ProblemSize - 1);
				std::vector<double> pa, pb;
				if (RG->randomNumber() > 0.5) {
					pa = ita->second.decisionVariables;
					pb = itb->second.decisionVariables;
				}
				else {
					pa = itb->second.decisionVariables;
					pb = ita->second.decisionVariables;
				}

				for (int ivar = 0; ivar < options.ProblemSize; ++ivar) {
					if (ivar <= cut1)
						child.push_back(pa[ivar]);
					else if ((ivar > cut1) & (ivar <= cut2)){
						child.push_back(pb[ivar]);
					}
					else {
						child.push_back(pa[ivar]);
					}
				}
			}
			else if (options.XoverType.compare("UNIBIN") == 0) {
				for (int ivar = 0; ivar < options.ProblemSize; ++ivar) {
					if (RG->randomNumber() > 0.5) 
						child.push_back(ita->second.decisionVariables[ivar]);
					else
						child.push_back(itb->second.decisionVariables[ivar]);
				}
			}
			else if (options.XoverType.compare("HEURISTIC") == 0) {
				int out = ita->second.compareRanking(itb->second);
				double best, worst;
				// find which parent is the fittests

				for (int ivar = 0; ivar < options.ProblemSize; ++ivar) {
					double r = RG->normalRandom(0, options.SDheuristic);
					double x;
					if (out == 1) {
						best = ita->second.decisionVariables[ivar];
						worst = itb->second.decisionVariables[ivar];
					
					}
					else {
						best = itb->second.decisionVariables[ivar];
						worst = ita->second.decisionVariables[ivar];
					}

					x = best + r * (best - worst);

					if (x < options.LowerBound[ivar] || x > options.UpperBound[ivar]) {
						double reduseS = 0.95;
						while (true) {
							r = RG->normalRandom(0, options.SDheuristic) * reduseS;
							x = best + r * (best - worst);
							if (x < options.LowerBound[ivar] || x > options.UpperBound[ivar]) {
								reduseS = reduseS * 0.95;
							}
							else
								break;
						}
					}
					child.push_back(x);
				}
			}
			
			newPopulation.push_back(child);
		}

		// Replace old population
		population.clear();
		Fronts.clear();
		selectedIds.clear();
		selectedParentIds.clear();
		for (unsigned int i = 0; i < newPopulation.size(); ++i) {
			Individual indivudual;
			indivudual.decisionVariables = newPopulation[i];
			population.insert(std::pair<int, Individual>(static_cast<int>(i), indivudual));
		}
	}

	void Population::mutation() {

		NSGAII::SingletonRealGenerator* RG = RG->getInstance();

		bool binaryMutation = options.populationType.compare("BINARY") == 0;
		
		std::map<int, Individual>::iterator it;
		for (it = population.begin(); it != population.end(); ++it) {
			// make a copy of the individual and work on that
			std::vector<double> individual_copy;
			double mut_mult = 1;
			while (true) {
				individual_copy = it->second.decisionVariables;
				for (unsigned int ivar = 0; ivar < individual_copy.size(); ++ivar) {
					double r = RG->randomNumber();
					if (r < options.MutationProbability* mut_mult) {
						if (binaryMutation) {
							if (individual_copy[ivar] > 0.5)
								individual_copy[ivar] = 0;
							else
								individual_copy[ivar] = 1;
						}
						else
							individual_copy[ivar] = RG->randomNumber(options.LowerBound[ivar], options.UpperBound[ivar]);
					}

					if (!binaryMutation) {
						if (individual_copy[ivar] < options.LowerBound[ivar])
							individual_copy[ivar] = options.LowerBound[ivar];
						if (individual_copy[ivar] > options.UpperBound[ivar])
							individual_copy[ivar] = options.UpperBound[ivar];
					}
				}
				if (options.useTabu) {
					if (isTabu(individual_copy))
						mut_mult = mut_mult * 1.5;
					else
						break;

				}	
				else
					break;
			}
			it->second.decisionVariables = individual_copy;
			tabuList.push_front(individual_copy);
		}	
	}

	void Population::appendParetoSolutions() {
		std::map<int, Individual>::iterator it, it1;
		int index = population.size();
		for (it = ParetoSolutions.begin(); it != ParetoSolutions.end(); ++it) {
			while (true) {
				it1 = population.find(index);
				if (it1 == population.end())
					break;
				else
					index++;
			}
			population.insert(std::pair<int, Individual>(index, it->second));
			index++;
		}
		ParetoSolutions.clear();
	}

	int Population::ParetoSize() {
		return ParetoSolutions.size();
	}

	void Population::printPareto() {
		std::map<int, Individual>::iterator it;
		std::ofstream outfile;
		outfile.open(options.OutputFile);
		outfile << ParetoSize() << " " << options.ProblemSize << " " << options.Nobjectives << std::endl;
		for (it = ParetoSolutions.begin(); it != ParetoSolutions.end(); ++it) {
			for (unsigned int i = 0; i < it->second.decisionVariables.size(); ++i) {
				outfile << it->second.decisionVariables[i] << " ";
			}
			outfile << std::endl;
		}
		for (it = ParetoSolutions.begin(); it != ParetoSolutions.end(); ++it) {
			for (unsigned int i = 0; i < it->second.objectiveFunctions.size(); ++i) {
				outfile << it->second.objectiveFunctions[i] << " ";
			}
			outfile << std::endl;
		}


		outfile.close();
	}

	void Population::openHistoryFile() {
		historyfile.open(options.HistoryFile);
	}

	void Population::closeHistoryFile() {
		historyfile.close();
	}

	void Population::printCurrentPareto(int gen) {
		std::map<int, Individual>::iterator it;
		historyfile << gen << " " << options.Nobjectives << " " << ParetoSize() << std::endl;
		for (it = ParetoSolutions.begin(); it != ParetoSolutions.end(); ++it) {
			for (unsigned int i = 0; i < it->second.objectiveFunctions.size(); ++i) {
				historyfile << it->second.objectiveFunctions[i] << " ";
			}
			historyfile << std::endl;
		}
	}

	bool Population::isTabu(std::vector<double>& v) {
		bool exists = false;
		std::list<std::vector<double> >::iterator it;
		for (it = tabuList.begin(); it != tabuList.end(); ++it){
			double d = 0;
			for (unsigned int j = 0; j < it->size(); ++j) {
				d += std::abs(it->at(j) - v[j]);
				if (d > options.tabuThreshold)
					break;
			}
			if (d < options.tabuThreshold) {
				std::cout << "Hey I found one here" << std::endl;
				exists = true;
				break;
			}
		}
		return exists;
	}
}
