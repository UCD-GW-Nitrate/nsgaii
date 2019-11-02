#pragma once
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>

namespace C2VSIM {

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

		bool readGWBud(std::string filename) {
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
			}
			gwBudfile.close();
			return true;
		}
	}
}
