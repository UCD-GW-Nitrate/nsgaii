#pragma once

#include <iostream>
#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>


namespace po = boost::program_options;

namespace NSGAII {
	void splitString(std::string& inp, std::vector<double>& out) {
		std::vector<std::string> temp;
		boost::split(temp, inp, [](char c) {return c == ' '; });
		double x;
		for (unsigned int i = 0; i < temp.size(); ++i) {
			x = boost::lexical_cast<double>(temp[i]);
			out.push_back(x);
		}

	}

	struct options {
		int ProblemSize;
		int PopulationSize;
		int Nobjectives;
		int TournamentSize;
		int MaxGenerations;
		double SDheuristic;
		double MutationProbability;
		std::string OutputFile;
		std::vector<double> LowerBound;
		std::vector<double> UpperBound;

	};

	bool readInputParameters(int argc, char* argv[], options& opt) {
		// Command line options
		po::options_description commandLineOptions("Command line options");
		commandLineOptions.add_options()
			("version,v", "print version information")
			("help,h", "Get a list of options in the configuration file")
			("config,c", po::value<std::string >(), "Set configuration file")
			;

		po::variables_map vm_cmd;

		po::store(po::parse_command_line(argc, argv, commandLineOptions), vm_cmd);

		if (vm_cmd.size() == 0) {
			std::cout << " To run nsgaii specify the configuration file as" << std::endl;
			std::cout << "-c configfilename" << std::endl << std::endl;;
			std::cout << "Other command line options are:" << std::endl;
			std::cout << commandLineOptions << std::endl;
			return false;
		}

		if (vm_cmd.count("version")) {
			std::cout << "|------------------|" << std::endl;
			std::cout << "|      NSGAII      |" << std::endl;
			std::cout << "| Version : 1.0.00 |" << std::endl;
			std::cout << "|    by  giorgk    |" << std::endl;
			std::cout << "|------------------|" << std::endl;
			return false;
		}

		// Configuration file options
		po::options_description config_options("Configuration file options");
		config_options.add_options()
			("OutputFile", po::value<std::string>(), "The file where the Output will be written")
			("ProblemSize", po::value<int>()->default_value(10), "Number of desicion variables")
			("Nobjectives", po::value<int>()->default_value(2), "Number of conflicting objectives")
			("PopulationSize", po::value<int>()->default_value(100), "Population size of GA")
			("LowerBound", po::value<std::vector<std::string> >()->multitoken(), "Lower bound of decision variables")
			("UpperBound", po::value<std::vector<std::string> >()->multitoken(), "Lower bound of decision variables")
			("TournamentSize", po::value<int>()->default_value(4), "Selection tournament size")
			("SDheuristic", po::value<double>()->default_value(0.25), "Standard deviation for heuristic crossover")
			("MutationProbability", po::value<double>(), "Mutation probability (default: 1/ProblemSize)")
			("MaxGenerations", po::value<int>()->default_value(200), "Maximum number of generations")
			;

		if (vm_cmd.count("help")) {
			std::cout << "Copy paste the following list of options" << std::endl;
			std::cout << "into a configuration file." << std::endl;
			std::cout << "The options without default values are mandatory" << std::endl;
			std::cout << "All options are case sensitive" << std::endl;
			std::cout << "------------------------------" << std::endl;
			std::cout << config_options << std::endl;
			return false;
		}

		po::variables_map vm_cfg;

		if (vm_cmd.count("config")) {
			std::cout << vm_cmd["config"].as<std::string>().c_str() << std::endl;
			po::store(po::parse_config_file<char>(vm_cmd["config"].as<std::string>().c_str(), config_options), vm_cfg);
			opt.Nobjectives = vm_cfg["Nobjectives"].as<int>();
			opt.ProblemSize = vm_cfg["ProblemSize"].as<int>();
			opt.PopulationSize = vm_cfg["PopulationSize"].as<int>();
			opt.TournamentSize = vm_cfg["TournamentSize"].as<int>();
			opt.SDheuristic = vm_cfg["SDheuristic"].as<double>();
			opt.MaxGenerations = vm_cfg["MaxGenerations"].as<int>();

			if (vm_cfg.count("MutationProbability"))
				opt.MutationProbability = vm_cfg["MutationProbability"].as<double>();
			else
				opt.MutationProbability = 1.0 / static_cast<double>(opt.ProblemSize);

			if (vm_cfg.count("OutputFile"))
				opt.OutputFile = vm_cfg["OutputFile"].as<std::string>();
			else
				opt.OutputFile = "nsgaiiOutput.dat";
			
			std::vector<std::string> temp = vm_cfg["LowerBound"].as<std::vector<std::string> >();
			splitString(temp[0], opt.LowerBound);
			if (opt.ProblemSize != opt.LowerBound.size()) {
				std::cout << "The Problem size ~= LowerBound size" << std::endl;
				return false;
			}

			temp.clear();
			temp = vm_cfg["UpperBound"].as<std::vector<std::string> >();
			splitString(temp[0], opt.UpperBound);
			if (opt.ProblemSize != opt.UpperBound.size()) {
				std::cout << "The Problem size ~= UpperBound" << std::endl;
				return false;
			}
		}
		return true;
	}
}
