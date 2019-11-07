#pragma once
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <boost/algorithm/string/trim.hpp>
#include <boost/program_options.hpp>
namespace po = boost::program_options;

namespace C2VSIM {

	struct ElemInfo {
		double area;
	};

	struct Diversion {
		int IRDV;
		int ICDVMAX;
		double FDVMAX;
		int ICOLRL;
		double FRACRL;
		int ICOLNL;
		double FRACNL;
		int NDLDV;
		int IRGDL;
		int ICOLDL;
		double FRACDL;
		int ICFSIRIG;
		int ICADJ;
		std::vector<int> IERELS;
		std::vector<double> FERELS;
		void readDivlist(std::stringstream& is) {
			is >> IRDV;
			is >> ICDVMAX;
			is >> FDVMAX;
			is >> ICOLRL;
			is >> FRACRL;
			is >> ICOLNL;
			is >> FRACNL;
			is >> NDLDV;
			is >> IRGDL;
			is >> ICOLDL;
			is >> FRACDL;
			is >> ICFSIRIG;
			is >> ICADJ;
		}
		void readdivElem(std::stringstream& ss, std::ifstream& is) {
			int Nel, ie;
			double fe;
			ss >> Nel;
			for (unsigned int i = 0; i < Nel; ++i) {
				if (i == 0) {
					ss >> ie;
					ss >> fe;
				}
				else {
					is >> ie;
					is >> fe;
				}
				IERELS.push_back(ie);
				FERELS.push_back(fe);
			}
		}
		void writeDivList(std::ofstream& os) {
			os << IRDV << " ";
			os << ICDVMAX << " ";
			os << FDVMAX << " ";
			os << ICOLRL << " ";
			os << FRACRL << " ";
			os << ICOLNL << " ";
			os << FRACNL << " ";
			os << NDLDV << " ";
			os << IRGDL << " ";
			os << ICOLDL << " ";
			os << FRACDL << " ";
			os << ICFSIRIG << " ";
			os << ICADJ << std::endl;;
		}
		void writeDivElem(std::ofstream& os) {
			os << IERELS.size() << " ";
			if (IERELS.size() == 0) {
				os << "0 0" << std::endl;;
			}
			else {
				for (unsigned int i = 0; i < IERELS.size(); ++i)
					os << IERELS[i] << " " << FERELS[i] << std::endl;
			}
		}
	};
	struct ByPass {
		int IA;
		int IDIVT;
		int IDIVC;
		double DIVRL;
		double DIVNL;
		std::vector<double> DIVX;
		std::vector<double> DIVY;
		std::vector<int> IERELS;
		std::vector<double> FERELS;
		void readBPlist(std::stringstream& ss, std::ifstream& is) {
			ss >> IA;
			ss >> IDIVT;
			ss >> IDIVC;
			ss >> DIVRL;
			ss >> DIVNL;
			if (IDIVC < 0) {
				double x, y;
				for (int i = 0; i < std::abs(IDIVC); ++i) {
					is >> x;
					is >> y;
					DIVX.push_back(x);
					DIVY.push_back(y);
				}
			}
		}

		void writeBPList(std::ofstream& os) {
			os << IA << " ";
			os << IDIVT << " ";
			os << IDIVC << " ";
			os << DIVRL << " ";
			os << DIVNL << std::endl;;
			if (IDIVC < 0) {
				for (unsigned int i = 0; i < DIVX.size(); ++i) {
					os << std::fixed << DIVX[i] << " " << DIVY[i] << std::endl;
				}
			}
		}

		void readSeepageLocations(std::stringstream& ss, std::ifstream& is) {
			int Nel, ie;
			double fe;
			ss >> Nel;
			if (Nel == 0) {
				ss >> ie;
				ss >> fe;
				IERELS.push_back(ie);
				FERELS.push_back(fe);
			}
			else {
				for (int i = 0; i < Nel; ++i) {
					if (i == 0) {
						ss >> ie;
						ss >> fe;
					}
					else {
						is >> ie;
						is >> fe;
					}
					IERELS.push_back(ie);
					FERELS.push_back(fe);
				}
			}
		}

		void writeSL(std::ofstream& os) {
			if (IERELS.size() == 0) {
				os << "0 0 0" << std::endl;
			}
			else {
				if (IERELS[0] == 0)
					os << "0 0 0" << std::endl;
				else {
					os << IERELS.size() << " ";
					for (unsigned int i = 0; i < IERELS.size(); ++i)
						os << IERELS[i] << " " << FERELS[i] << std::endl;
				}
			}			
		}
	};

	struct DiversionData {
		int NRDV; // Number of diversions
		std::map<int, Diversion> DIVS;
		int NDIVS; // number of bypasses
		int FACTX;
		int TUNITX;
		int FACTY;
		int TUNITY;
		std::map<int, ByPass> BPS;
		std::vector<std::string> time;
		std::vector<std::vector<double> > Data;
		int NCOLDV;
		double FACTDV;
		int NSPDV;
		int NFQDV;
		std::string DSSFL = "";

		void appendDiversion(Diversion& div, std::vector<double> TS) {
			NRDV++;
			NCOLDV++;
			for (unsigned int i = 0; i < Data.size(); ++i)
				Data[i].push_back(TS[i]);
			int icol = Data[0].size();
			div.ICOLRL = icol;
			div.ICOLNL = icol;
			div.ICOLDL = icol;
			DIVS.insert(std::pair<int, Diversion>(NRDV, div));
		}
	};

	struct GWBUD {
		double DP; //Deep Percolation
		double S; //Ending Storage
		double NDP; // Net Deep Percolation
		double GS; // Gain from stream
		double R; //Recharge
		double L; // Lake
		double BI; // Boundary Inflow
		double SS; // Subsidence
		double SI; //Subsurface Irrigation
		double TDO; // Tile Drain Outflow
		double P; // Pumping
		double NSI; //Net subsurface Inflow
		double D; // Discrepancy
		double CS; // Cumulative Subsidence
		void read(std::ifstream& ifs) {
			ifs >> DP;
			ifs >> S; ifs >> S;
			ifs >> NDP;
			ifs >> GS;
			ifs >> R;
			ifs >> L;
			ifs >> BI;
			ifs >> SS;
			ifs >> SI;
			ifs >> TDO;
			ifs >> P;
			ifs >> NSI;
			ifs >> D;
			ifs >> CS;
		}

		GWBUD(){
			DP = 0.0; S = 0.0; NDP = 0.0; GS = 0.0; R = 0.0; L = 0.0; BI = 0.0;
			SS = 0.0; SI = 0.0; TDO = 0.0; P = 0.0; NSI = 0.0; D = 0.0; CS = 0.0;
		}

		GWBUD operator+(const GWBUD& a) const {
			GWBUD b;
			b.DP = DP + a.DP;
			b.S = S + a.S;
			b.NDP = NDP + a.NDP;
			b.GS = GS + a.GS;
			b.R = R + a.R;
			b.L = L + a.L;
			b.BI = BI + a.BI;
			b.SS = SS + a.SS;
			b.SI = SI + a.SI;
			b.TDO = TDO + a.TDO;
			b.P = P + a.P;
			b.NSI = NSI + a.NSI;
			b.D = D + a.D;
			b.CS = CS + a.CS;
			return b;
		}
	};

	class GWbudTimeSeries {
	public:
		void add(GWbudTimeSeries other);
		std::vector<GWBUD> GwBudTS;
	};

	void GWbudTimeSeries::add(GWbudTimeSeries other) {
		if (GwBudTS.size() < other.GwBudTS.size())
			GwBudTS.resize(other.GwBudTS.size(), GWBUD());
		for (unsigned int i = 0; i < GwBudTS.size(); ++i) {
			GwBudTS[i] = GwBudTS[i] + other.GwBudTS[i];
		}
	}

	namespace READERS {
		bool readGWHydOut(std::string filename, std::map<int, std::vector<double> >& GWH) {
			std::ifstream GWfile;
			GWfile.open(filename);
			if (!GWfile.is_open()) {
				std::cout << "Cant open file: " << filename << std::endl;
				return false;
			}
			std::string temp;
			for (unsigned int i = 0; i < 7; ++i) {
				std::getline(GWfile, temp);
			}
			std::string s;
			double x;
			std::vector<std::vector<double> > values(1393, std::vector<double>(1056, 0.0));

			for (unsigned int i = 0; i < 1056; ++i) {
				GWfile >> s;
				for (unsigned int j = 0; j < 1393; ++j) {
					GWfile >> x;
					values[j][i] = x;
				}

			}
			GWfile.close();

			for (unsigned int j = 0; j < 1393; ++j) {
				GWH[j] = values[j];
			}
			return true;
		}

		bool readGWBud(std::string filename, GWbudTimeSeries& GWTS) {
			std::ifstream gwBudfile;
			gwBudfile.open(filename);
			if (!gwBudfile.is_open()) {
				std::cout << "Cant open file: " << filename << std::endl;
				return false;
			}

			std::string temp;
			for (unsigned int isub = 0; isub < 21; ++isub) {
				GWbudTimeSeries gwts;

				int cnt_dash = 0;
				while (true) {
					std::getline(gwBudfile, temp);
					std::string::iterator its = temp.begin();
					if (its != temp.end()) {
						if (std::string(1,*its).compare("-") == 0)
							cnt_dash++;
					}
					if (cnt_dash >= 2)
						break;
				}

				for (unsigned int i = 0; i < 1056; ++i) {
					gwBudfile >> temp;
					GWBUD gw;
					gw.read(gwBudfile);
					gwts.GwBudTS.push_back(gw);
				}
				GWTS.add(gwts);
			}
			gwBudfile.close();
			return true;
		}

		bool readDivSpec(std::string filename, DiversionData& divSpec) {

			std::ifstream divfile;
			divfile.open(filename);
			if (!divfile.is_open()) {
				std::cout << "Cant open file: " << filename << std::endl;
				return false;
			}

			int section = 1;
			std::string temp;
			std::string::iterator its;
			std::map<int, Diversion>::iterator it;
			std::map<int, ByPass>::iterator itbp;
			std::stringstream ss;
			int id;
			int cnt_div = 0;
			while (true) {
				std::getline(divfile, temp);
				boost::algorithm::trim(temp);
				
				its = temp.begin();
				if (its == temp.end())
					continue;
				if (temp.compare(0, 1, "C") == 0)
					continue;

				ss.clear();
				ss.str(temp);
				if (section == 1) {
					ss >> divSpec.NRDV;
					section++;

				}
				else if (section == 2) {
					ss >> id;
					Diversion div;
					div.readDivlist(ss);
					divSpec.DIVS.insert(std::pair<int, Diversion>(id, div));
					cnt_div++;
					if (cnt_div >= divSpec.NRDV) {
						section++;
						cnt_div = 0;
					}
					
				}
				else if (section == 3) {
					ss >> id;
					it = divSpec.DIVS.find(id);
					if (it != divSpec.DIVS.end()) {
						it->second.readdivElem(ss, divfile);
						cnt_div++;
					}
					if (cnt_div >= divSpec.NRDV) {
						section++;
						cnt_div = 0;
					}
				}
				else if (section == 4) {
					ss >> divSpec.NDIVS;
					std::getline(divfile, temp);
					ss.clear();
					ss.str(temp);
					ss >> divSpec.FACTX;
					std::getline(divfile, temp);
					ss.clear();
					ss.str(temp);
					ss >> divSpec.TUNITX;
					std::getline(divfile, temp);
					ss.clear();
					ss.str(temp);
					ss >> divSpec.FACTY;
					std::getline(divfile, temp);
					ss.clear();
					ss.str(temp);
					ss >> divSpec.TUNITY;
					section++;
				}
				else if (section == 5) {
					ss >> id;
					ByPass bp;
					bp.readBPlist(ss, divfile);
					divSpec.BPS.insert(std::pair<int, ByPass>(id, bp));
					cnt_div++;
					if (cnt_div >= divSpec.NDIVS) {
						section++;
						cnt_div = 0;
					}
				}
				else if (section == 6) {
					ss >> id;
					itbp = divSpec.BPS.find(id);
					if (itbp != divSpec.BPS.end()) {
							itbp->second.readSeepageLocations(ss, divfile);
						cnt_div++;
					}
					if (cnt_div >= divSpec.NDIVS) {
						break;
					}
				}
			}
			divfile.close();
			return true;
		}

		bool readDivData(std::string filename, DiversionData& divdata, int NtimeSteps) {
			std::ifstream divfile;
			divfile.open(filename);
			if (!divfile.is_open()) {
				std::cout << "Cant open file: " << filename << std::endl;
				return false;
			}

			std::string temp, tm;
			std::string::iterator its;
			std::stringstream ss;
			double x;
			int section = 1;
			int itm = 0;
			while (true) {
				std::getline(divfile, temp);
				boost::algorithm::trim(temp);
				its = temp.begin();
				if (its == temp.end())
					continue;
				if (temp.compare(0, 1, "C") == 0)
					continue;
				if (section == 1){
					ss.clear();
					ss.str(temp);
					ss >> divdata.NCOLDV;
					divdata.Data.resize(NtimeSteps, std::vector<double>(divdata.NCOLDV, 0.0));
					std::getline(divfile, temp);
					ss.clear();
					ss.str(temp);
					ss >> divdata.FACTDV;
					std::getline(divfile, temp);
					ss.clear();
					ss.str(temp);
					ss >> divdata.NSPDV;
					std::getline(divfile, temp);
					ss.clear();
					ss.str(temp);
					ss >> divdata.NFQDV;
					std::getline(divfile, temp);
					section++;
				}
				else if (section == 2) {
					ss.clear();
					ss.str(temp);
					
					ss >> tm;
					divdata.time.push_back(tm);
					for (int i = 0; i < divdata.NCOLDV; ++i) {
						ss >> x;
						divdata.Data[itm][i] = x;
					}
					itm++;
					if (itm >= NtimeSteps)
						break;
				}
			}
			divfile.close();
			return true;
		}

		bool readDivTimeSeries(std::string filename, std::map<int, std::vector<double> >& DTS) {
			std::ifstream dtsfile;
			dtsfile.open(filename);
			if (!dtsfile.is_open()) {
				std::cout << "Cant open file: " << filename << std::endl;
				return false;
			}
			std::string temp;
			std::stringstream ss;
			int Ndivs, Nsteps, id;
			std::getline(dtsfile, temp);
			ss.clear();
			ss.str(temp);
			ss >> Ndivs;
			ss >> Nsteps;
			double x;
			for (int i = 0; i < Ndivs; ++i) {
				std::getline(dtsfile, temp);
				ss.clear();
				ss.str(temp);
				std::vector<double> Xv;
				ss >> id;
				for (int j = 0; j < Nsteps; ++j) {
					ss >> x;
					Xv.push_back(x);
				}
				DTS.insert(std::pair<int, std::vector<double > >(id, Xv));
			}
			dtsfile.close();
			return true;
		}

		bool readDivElems(std::string filename, std::map<int, std::vector<int> >& divElem) {
			std::ifstream dvfile;
			dvfile.open(filename);
			if (!dvfile.is_open()) {
				std::cout << "Cant open file: " << filename << std::endl;
				return false;
			}
			std::string temp;
			std::stringstream ss;
			int Ndivs, Nelem, id, k;
			std::getline(dvfile, temp);
			ss.clear();
			ss.str(temp);
			ss >> Ndivs;
			for (int i = 0; i < Ndivs; ++i) {
				std::getline(dvfile, temp);
				ss.clear();
				ss.str(temp);
				ss >> id;
				ss >> Nelem;
				std::vector<int> elids;
				for (int j = 0; j < Nelem; ++j) {
					ss >> k;
					elids.push_back(k);
				}
				divElem.insert(std::pair<int, std::vector<int> >(id, elids));
			}
			dvfile.close();
			return true;
		}

		bool readElemInfo(std::string filename, std::map<int, ElemInfo>& elemInfoMap) {
			std::ifstream elfile;
			elfile.open(filename);
			if (!elfile.is_open()) {
				std::cout << "Cant open file: " << filename << std::endl;
				return false;
			}
			int Nelem, id;
			C2VSIM::ElemInfo v;
			elfile >> Nelem;
			for (int i = 0; i < Nelem; ++i) {
				elfile >> id;
				elfile >> v.area;
				elemInfoMap[id] = v;
			}
			elfile.close();
			return true;
		}
	}

	namespace WRITERS {
		void writeDivSpec(std::string filename, DiversionData& divSpec) {
			std::ofstream outfile;
			outfile.open(filename);
			outfile << divSpec.NRDV << std::endl;
			std::map<int, Diversion>::iterator itd;
			for (itd = divSpec.DIVS.begin(); itd != divSpec.DIVS.end(); ++itd) {
				outfile << itd->first << " ";
				itd->second.writeDivList(outfile);
			}

			for (itd = divSpec.DIVS.begin(); itd != divSpec.DIVS.end(); ++itd) {
				outfile << itd->first << " ";
				itd->second.writeDivElem(outfile);
			}

			outfile << divSpec.NDIVS << std::endl;
			outfile << divSpec.FACTX << std::endl;
			outfile << divSpec.TUNITX << "min" << std::endl;
			outfile << divSpec.FACTY << std::endl;
			outfile << divSpec.TUNITY << "min" << std::endl;
			outfile << "C" << std::endl;
			std::map<int, ByPass>::iterator itb;
			for (itb = divSpec.BPS.begin(); itb != divSpec.BPS.end(); ++itb) {
				outfile << itb->first << " ";
				itb->second.writeBPList(outfile);
			}

			for (itb = divSpec.BPS.begin(); itb != divSpec.BPS.end(); ++itb) {
				outfile << itb->first << " ";
				itb->second.writeSL(outfile);
			}
			outfile.close();
		}

		void writeDivData(std::string filename, DiversionData& divdata) {
			std::ofstream outfile;
			outfile.open(filename);
			outfile << divdata.NCOLDV << std::endl;
			outfile << std::fixed << divdata.FACTDV << std::endl;
			outfile << divdata.NSPDV << std::endl;
			outfile << divdata.NFQDV << std::endl;
			outfile << divdata.DSSFL << std::endl;

			for (unsigned int i = 0; i < divdata.time.size(); ++i) {
				outfile << divdata.time[i] << " ";
				for (unsigned int j = 0; j < divdata.Data[i].size(); ++j) {
					outfile << divdata.Data[i][j] << " ";
				}
				outfile << std::endl;
			}
			outfile.close();
		}
	}

	namespace OPTIONS {
		struct options {
			int Nsteps = 1056;
			std::string divSpecFile; 
			std::string divDataFile;
			std::string divElemFile;
			std::string divTimeSeriesFile;
			std::string BaseWTfile;
			std::string BaseGWbudFile;
			std::string ElementInfofile;
			std::string SimulationExe;
			std::string BudgetExe;
		};

		bool readConfigFile(int argc, char* argv[], C2VSIM::OPTIONS::options& opt) {
			po::options_description commandLineOptions("Command line options");
			commandLineOptions.add_options()
				("help,h", "Get a list of options in the configuration file")
				("config,c", po::value<std::string >(), "Set configuration file")
				("model,m", po::value<std::string >(), "Set configuration file for the C2Vsim")
				;

			po::variables_map vm_cmd;
			po::store(po::parse_command_line(argc, argv, commandLineOptions), vm_cmd);

			po::options_description config_options("Configuration file options");
			config_options.add_options()
				("Nsteps", po::value<int>()->default_value(1056), "Number of monthly steps")
				("divSpecFile", po::value<std::string>(), "Diversion Specification file from C2Vsim")
				("divDataFile", po::value<std::string>(), "Diversion Data file from C2Vsim")
				("divElemFile", po::value<std::string>(), "File with the elements to recieve diversions")
				("divTimeSeriesFile", po::value<std::string>(), "Timeseries of the diversions")
				("BaseWTfile", po::value<std::string>(), "Groundwater table for the base scenario")
				("ElementInfofile", po::value<std::string>(), "Element info (id area)")
				("SimulationExe", po::value<std::string>(), "Path including the simulation executable")
				("BudgetExe", po::value<std::string>(), "Path including the budget executable")
				("BaseGWbudFile", po::value<std::string>(), "Groundwater Storage for the base scenario")
				;

			if (vm_cmd.count("help")) {
				std::cout << "If you are coupling optimization" << std::endl;
				std::cout << "with C2Vsim then create an additional file" << std::endl;
				std::cout << "with the following options:" << std::endl;
				std::cout << "------------------------------" << std::endl;
				std::cout << config_options << std::endl;
				std::cout << "To run nsgaii with c2vsim pass both config and model arguments:" << std::endl;
				std::cout << "nsgaii.exe -c configfile -m modelfile" << std::endl;
				return false;
			}

			po::variables_map vm_cfg;

			if (vm_cmd.count("model")) {
				std::cout << "Configuration model file: " << vm_cmd["model"].as<std::string>().c_str() << std::endl;
				po::store(po::parse_config_file<char>(vm_cmd["model"].as<std::string>().c_str(), config_options), vm_cfg);

				opt.Nsteps = vm_cfg["Nsteps"].as<int>();
				opt.divSpecFile = vm_cfg["divSpecFile"].as<std::string>();
				opt.divDataFile = vm_cfg["divDataFile"].as<std::string>();
				opt.divElemFile = vm_cfg["divElemFile"].as<std::string>();
				opt.divTimeSeriesFile = vm_cfg["divTimeSeriesFile"].as<std::string>();
				opt.BaseWTfile = vm_cfg["BaseWTfile"].as<std::string>();
				opt.ElementInfofile = vm_cfg["ElementInfofile"].as<std::string>();
				opt.SimulationExe = vm_cfg["SimulationExe"].as<std::string>();
				opt.BudgetExe = vm_cfg["BudgetExe"].as<std::string>();
				opt.BaseGWbudFile = vm_cfg["BaseGWbudFile"].as<std::string>();
			}
			return true;
		}
	}
}
