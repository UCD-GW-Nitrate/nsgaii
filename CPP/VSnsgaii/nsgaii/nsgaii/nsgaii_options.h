#ifndef NSGAII_OPTIONS_H
#define NSGAII_OPTIONS_H

#include <iostream>
#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include "c2vsim_io.h"


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

	//enum crossoverType {ONECUT, TWOCUT, UNIFORM,HEURISTIC};

	struct options {
		std::string populationType;
		int ProblemSize;
		int PopulationSize;
		int Nobjectives;
		int TournamentSize;
		int MaxGenerations;
		double SDheuristic;
		double MutationProbability;
		std::string OutputFile;
		std::string HistoryFile;
		std::vector<double> LowerBound;
		std::vector<double> UpperBound;
		bool bUseModel = false;
		bool bWriteHistory = false;
		std::string XoverType;
		bool useTabu = false;
		double tabuThreshold;
		bool bResumeOptimization = false;
		bool bPrintRestartFile = false;
		std::string RestartOutFile;
		std::string RestartInFile;
	};

	bool readInputParameters(int argc, char* argv[], options& opt) {
		// Command line options
		po::options_description commandLineOptions("Command line options");
		commandLineOptions.add_options()
			("version,v", "print version information")
			("help,h", "Get a list of options in the configuration file")
			("config,c", po::value<std::string >(), "Set configuration file")
			("model,m", po::value<std::string >(), "Set configuration file for the C2Vsim")
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
		    // GAoptions
			("GAoptions.PopulationType", po::value<std::string>(), "Type of Population: REAL or BINARY")
			("GAoptions.ProblemSize", po::value<int>()->default_value(10), "Number of desicion variables")
			("GAoptions.Nobjectives", po::value<int>()->default_value(2), "Number of conflicting objectives")
			("GAoptions.PopulationSize", po::value<int>()->default_value(100), "Population size of GA")
			("GAoptions.LowerBound", po::value<std::vector<std::string> >()->multitoken(), "Lower bound of decision variables (Used only for REAL)")
			("GAoptions.UpperBound", po::value<std::vector<std::string> >()->multitoken(), "Lower bound of decision variables (Used only for REAL)")
			("GAoptions.TournamentSize", po::value<int>()->default_value(4), "Selection tournament size")
			("GAoptions.SDheuristic", po::value<double>()->default_value(0.25), "Standard deviation for heuristic crossover")
			("GAoptions.MutationProbability", po::value<double>(), "Mutation probability (default: 1/ProblemSize)")
			("GAoptions.MaxGenerations", po::value<int>()->default_value(200), "Maximum number of generations")
			("GAoptions.XoverType", po::value<std::string>(), "Crossover: ONECUT,TWOCUT,UNIREAL,UNIBIN,HEURISTIC")
			("GAoptions.tabuThres", po::value<double>()->default_value(0.0), "Set tabu threshold to enable it")

            // Outputfiles
            ("Output.SolutionFile", po::value<std::string>(), "The file where the Output will be written")
            ("Output.HistoryFile", po::value<std::string>(), "The file where the fronts for each generation are written")
            ("Output.RestartInFile", po::value<std::string>(), "Provide the file if you need to resume the optimization")
            ("Output.RestartOutFile", po::value<std::string>(), "The file where the last state of the optimization will be printed")
			;

		if (vm_cmd.count("help")) {
			std::cout << " To run nsgaii specify the configuration file as" << std::endl;
			std::cout << "-c configfilename" << std::endl << std::endl;;
			std::cout << "Other command line options are:" << std::endl;
			std::cout << commandLineOptions << std::endl;

			std::cout << "NSGAII configuration file options:" << std::endl;
			std::cout << "The options without default values are mandatory" << std::endl;
			std::cout << "(All options are case sensitive)" << std::endl;
			std::cout << "------------------------------" << std::endl;
			std::cout << config_options << std::endl;

			C2VSIM::OPTIONS::options opt;
			C2VSIM::OPTIONS::readConfigFile(argc, argv, opt);

			return false;
		}

		if (vm_cmd.count("model")) {
			opt.bUseModel = true;
		}

		po::variables_map vm_cfg;

		if (vm_cmd.count("config")) {
			std::cout << "Configuration file: " << vm_cmd["config"].as<std::string>().c_str() << std::endl;
			po::store(po::parse_config_file<char>(vm_cmd["config"].as<std::string>().c_str(), config_options), vm_cfg);

            { // GA options
                opt.populationType = vm_cfg["GAoptions.PopulationType"].as<std::string>();
                opt.ProblemSize = vm_cfg["GAoptions.ProblemSize"].as<int>();
                opt.Nobjectives = vm_cfg["GAoptions.Nobjectives"].as<int>();
                opt.PopulationSize = vm_cfg["GAoptions.PopulationSize"].as<int>();
                opt.TournamentSize = vm_cfg["GAoptions.TournamentSize"].as<int>();
                opt.SDheuristic = vm_cfg["GAoptions.SDheuristic"].as<double>();
                opt.MaxGenerations = vm_cfg["GAoptions.MaxGenerations"].as<int>();
                opt.XoverType = vm_cfg["GAoptions.XoverType"].as<std::string>();

                if (vm_cfg.count("GAoptions.MutationProbability"))
                    opt.MutationProbability = vm_cfg["GAoptions.MutationProbability"].as<double>();
                else
                    opt.MutationProbability = 1.0 / static_cast<double>(opt.ProblemSize);

                opt.tabuThreshold = vm_cfg["GAoptions.tabuThres"].as<double>();
                if (opt.tabuThreshold > 0)
                    opt.useTabu = true;

                if (opt.populationType.compare("REAL") == 0) {
                    std::vector<std::string> temp = vm_cfg["GAoptions.LowerBound"].as<std::vector<std::string> >();
                    splitString(temp[0], opt.LowerBound);
                    if (opt.ProblemSize != opt.LowerBound.size()) {
                        if (opt.LowerBound.size() == 1) {
                            opt.LowerBound.resize(opt.ProblemSize, opt.LowerBound[0]);
                        }
                        else {
                            std::cout << "The Problem size ~= LowerBound size" << std::endl;
                            return false;
                        }
                    }

                    temp.clear();
                    temp = vm_cfg["GAoptions.UpperBound"].as<std::vector<std::string> >();
                    splitString(temp[0], opt.UpperBound);
                    if (opt.ProblemSize != opt.UpperBound.size()) {
                        if (opt.UpperBound.size() == 1) {
                            opt.UpperBound.resize(opt.ProblemSize, opt.UpperBound[0]);
                        }
                        else {
                            std::cout << "The Problem size ~= UpperBound" << std::endl;
                            return false;
                        }
                    }
                }
            }

            { // Output files options
                if (vm_cfg.count("Output.SolutionFile"))
                    opt.OutputFile = vm_cfg["Output.SolutionFile"].as<std::string>();
                else
                    opt.OutputFile = "nsgaiiOutput.dat";

                if (vm_cfg.count("Output.HistoryFile")) {
                    opt.HistoryFile = vm_cfg["Output.HistoryFile"].as<std::string>();
                    opt.bWriteHistory = true;
                }

                if (vm_cfg.count("Output.RestartInFile")){
                    opt.RestartInFile  = vm_cfg["Output.RestartInFile"].as<std::string>();
                    if (!opt.RestartInFile.empty())
                        opt.bResumeOptimization = true;
                }

                if (vm_cfg.count("Output.RestartOutFile")){
                    opt.RestartOutFile = vm_cfg["Output.RestartOutFile"].as<std::string>();
                    if (!opt.RestartOutFile.empty())
                        opt.bPrintRestartFile = true;
                }
            }
		}
		return true;
	}
}

#endif // NSGAII_OPTIONS_H
