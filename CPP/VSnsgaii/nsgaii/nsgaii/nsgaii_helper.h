#pragma once

#include <iostream>
//#include <random>
//#include <chrono>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>

#include <boost/mpi.hpp>

namespace NSGAII {
	//void test_random_Numbers() {
	//	unsigned int seed = std::chrono::system_clock::now().time_since_epoch().count();
	//	std::cout << "Seed: " << seed << std::endl;
	//	std::default_random_engine generator(seed);
	//	std::uniform_real_distribution<double> distribution(0.0, 10.0);
	//	for (int i = 0; i < 20; ++i) {
	//		double x = distribution(generator);
	//		std::cout << std::fixed << x << std::endl;
	//	}
	//}


	template<typename T>
	void broad_cast_vector(std::vector<T>& v, int proc, boost::mpi::communicator& world) {
		T value;
		int VectorSize;

		if (world.rank() == proc) {
			VectorSize = static_cast<int>(v.size());
		}

		boost::mpi::broadcast(world, VectorSize, proc);
		if (world.rank() != proc) {
			v.clear();
		}

		for (int i = 0; i < VectorSize; ++i) {
			if (world.rank() == proc) {
				value = v[i];
			}
			boost::mpi::broadcast(world, value, proc);
			if (world.rank() != proc) {
				v.push_back(value);
			}
		}
	}

	class SingletonRealGenerator {
		static SingletonRealGenerator* _instance;
		unsigned int seed;

		// private constructor
		SingletonRealGenerator() {
			//seed = std::chrono::system_clock::now().time_since_epoch().count();
			generator.seed(std::time(0));
		}

		boost::random::uniform_real_distribution<double> uniformDistribution;
		//std::uniform_real_distribution<double> uniformDistribution;
		boost::random::normal_distribution<double> normalDistribution;
		//std::normal_distribution<double> normalDistribution;

		//std::default_random_engine generator;
		//std::mt19937 generator;
		boost::mt19937 generator;

	public:

		void printSeed() {
			std::cout << "Seed: " << seed << std::endl;
			//std::cout << "a,b" << uniformDistribution.a() << "," << uniformDistribution.b() << std::endl;
			//std::cout << "m,s" << normalDistribution.mean() << "," << normalDistribution.sigma() << std::endl;
		}
		static SingletonRealGenerator* getInstance() {
			if (!_instance)
				_instance = new SingletonRealGenerator;
			return _instance;
		}
		double randomNumber() {
			return uniformDistribution(generator);
		}
		double randomNumber(double min, double max) {
			return min + (max - min) * randomNumber();
		}
		double normalRandom() {
			return normalDistribution(generator);
		}
		double normalRandom(double m, double s) {
			return m + normalRandom() * s;
		}
		int randomNumber(int min, int max) {
			double r = randomNumber(static_cast<double>(min), static_cast<double>(max));
			return static_cast<int>(r + 0.5);
		}
	};

	long long convertBinaryToDecimal(std::vector<int>& v) {
		long long n = 0;
		for (unsigned int i = 0; i < v.size(); ++i) {
			if (v[i]>0)
				n += std::pow(2, i);
		}
		return n;
	}
}
