#pragma once
#include <map>
#include <vector>
#include <string>

#include <boost/filesystem.hpp>

#include "c2vsim_io.h"
#include "nsgaii_options.h"
#include "nsgaii_helper.h"
#include "nsgaii_core.h"

namespace C2VSIM {
	class c2vsimData {
	public:
		c2vsimData(C2VSIM::OPTIONS::options& opt);
		void readInputFiles();
		void debugMsg();
		// For a given location in the encoding return the diversion node and the element 
		bool getNodeElemId(unsigned int i, int& node, int& elem);
		C2VSIM::DiversionData getDiversionData();
		std::vector<double> getDTS(int node);
		double getArea(int elemId);
		std::string simulationExe();
		std::string budgetExe();
		std::vector<double> GWHbaseValue(int node);

	private:
		C2VSIM::OPTIONS::options options;
		std::map<int, std::vector<int> > divElem;
		std::map<int, int> ElemDiv;
		std::map<int, int> idElem;
		std::map<int, std::vector<double> > DTS;
		C2VSIM::DiversionData divData;
		C2VSIM::GWbudTimeSeries GWBUD;
		std::map<int, std::vector<double> > GWH;
		std::map<int, ElemInfo> elemInfoMap;

		void makeElemDiv();
	};

	c2vsimData::c2vsimData(C2VSIM::OPTIONS::options& opt)
		:
		options(opt)
	{}

	void c2vsimData::readInputFiles() {
		C2VSIM::READERS::readDivElems(options.divElemFile, divElem);
		makeElemDiv();
		C2VSIM::READERS::readDivTimeSeries(options.divTimeSeriesFile, DTS);
		C2VSIM::READERS::readDivSpec(options.divSpecFile, divData);
		C2VSIM::READERS::readDivData(options.divDataFile, divData, options.Nsteps);
		C2VSIM::READERS::readGWHydOut(options.BaseWTfile, GWH);
		C2VSIM::READERS::readElemInfo(options.ElementInfofile, elemInfoMap);
		C2VSIM::READERS::readGWBud(options.BaseGWbudFile, GWBUD);
	}

	void c2vsimData::makeElemDiv() {
		std::map<int, std::vector<int> >::iterator it;
		int cnt = 0;
		for (it = divElem.begin(); it != divElem.end(); ++it) {
			for (unsigned int i = 0; i < it->second.size(); ++i) {
				ElemDiv.insert(std::pair<int, int>(it->second[i], it->first));
				idElem.insert(std::pair<int, int>(cnt++, it->second[i]));
			}
		}
	}

	bool c2vsimData::getNodeElemId(unsigned int i, int& node, int& elem) {
		std::map<int, int>::iterator it;
		it = idElem.find(i);
		if (it != idElem.end()) {
			elem = it->second;
			it = ElemDiv.find(elem);
			if (it != ElemDiv.end()) {
				node = it->second;
				return true;
			}
			else {
				node = -9;
				return false;
			}
		}
		else {
			elem = -9;
			node = -9;
			return false;
		}
	}

	C2VSIM::DiversionData c2vsimData::getDiversionData() {
		return divData;
	}

	std::vector<double> c2vsimData::getDTS(int node) {
		std::map<int, std::vector<double> >::iterator it = DTS.find(node);
		if (it != DTS.end())
			return it->second;
		else
			return std::vector<double>();
	}

	double c2vsimData::getArea(int elemId) {
		std::map<int, ElemInfo>::iterator it;
		it = elemInfoMap.find(elemId);
		if (it != elemInfoMap.end()) {
			return it->second.area;
		}
		else {
			return -1;
		}
	}

	std::string c2vsimData::simulationExe() {
		return options.SimulationExe;
	}
	std::string c2vsimData::budgetExe() {
		return options.BudgetExe;
	}

	std::vector<double> c2vsimData::GWHbaseValue(int node) {
		std::map<int, std::vector<double> >::iterator it = GWH.find(node);
		if (it != GWH.end())
			return it->second;
		else
			std::vector<double>();
	}

	void c2vsimData::debugMsg() {
		std::cout << GWH.size() << std::endl;
	}

	namespace OF {
		void maxWTminArea(std::vector<double>& var, std::vector<double>& fun, C2VSIM::c2vsimData& cvd) {
			NSGAII::SingletonRealGenerator* RG = RG->getInstance();
			
			// find the diversion nodes and the elements to apply the water
			std::map<int, std::vector<int> > nodeElemMap;
			std::map<int, std::vector<int> >::iterator it;
			int node, elem;
			for (unsigned int i = 0; i < var.size(); ++i) {
				if (var[i] > 0.5) {
					if (cvd.getNodeElemId(i, node, elem)) {
						it = nodeElemMap.find(node);
						if (it == nodeElemMap.end()) {
							std::vector<int> temp;
							temp.push_back(elem);
							nodeElemMap.insert(std::pair<int, std::vector<int> >(node, temp));
						}
						else {
							it->second.push_back(elem);
						}
					}
				}
			}


			// get a copy of the existing diversion data
			C2VSIM::DiversionData divData = cvd.getDiversionData();
			double area = 0;
			for (it = nodeElemMap.begin(); it != nodeElemMap.end(); ++it) {
				std::vector<double> TS = cvd.getDTS(it->first);
				C2VSIM::Diversion div;
				div.IRDV = it->first;
				div.ICDVMAX = 265;
				div.FDVMAX = 1;
				div.FRACRL = 0.95; //Fraction for recoverable loss = recharge
				div.FRACNL = 0.05; //Fraction for non recoverable loss 
				div.NDLDV = 1;
				div.IRGDL = 1;
				div.FRACDL = 0; // Fraction for irrigation
				div.ICFSIRIG = 1;
				div.ICADJ = 1;
				div.IERELS = it->second;
				for (unsigned int i = 0; i < it->second.size(); ++i) {
					double a = cvd.getArea(it->second[i]);
					area += a;
					div.FERELS.push_back(a);
				}
				divData.appendDiversion(div, TS);
			}

			fun.clear();
			fun.push_back(RG->randomNumber(1000000, 2000000));
			fun.push_back(area);
			return;


			C2VSIM::WRITERS::writeDivSpec("d:/giorgk/Documents/GitHub/C2VsimCG/RunC2Vsim/tempSpec.dat", divData);
			C2VSIM::WRITERS::writeDivData("d:/giorgk/Documents/GitHub/C2VsimCG/RunC2Vsim/tempData.dat", divData);
			std::string main_dir = boost::filesystem::current_path().string();
			// Enter into the simulation path
			boost::filesystem::current_path("d:/giorgk/Documents/GitHub/C2VsimCG/RunC2Vsim");
			std::cout << boost::filesystem::current_path() << std::endl;
			
			std::string sim_command = cvd.simulationExe();
			std::string bud_command = cvd.budgetExe();
			sim_command.append(" CVsim.in");
			bud_command.append(" CVBudget.in");


			system(sim_command.c_str());
			system(bud_command.c_str());
			//Read the output while in the path
			std::map<int, std::vector<double> > simGWhyd;
			C2VSIM::GWbudTimeSeries simGBbud;
			C2VSIM::READERS::readGWHydOut("d:/giorgk/Documents/GitHub/C2VsimCG/RunC2Vsim/Results/CVGWhyd.out", simGWhyd);
			C2VSIM::READERS::readGWBud("d:/giorgk/Documents/GitHub/C2VsimCG/RunC2Vsim/Results/CVground.BUD", simGBbud);

			std::cout << "Model finished" << std::endl;
			boost::filesystem::current_path(main_dir);
			std::cout << boost::filesystem::current_path() << std::endl;
			bool tf = true;
			std::cout << "all is" << tf << std::endl;

			std::map<int, std::vector<double> >::iterator itwt;
			double f = 0;
			double cntf = 0;
			double fmax = -999999;
			double x;
			for (itwt = simGWhyd.begin(); itwt != simGWhyd.end(); ++itwt) {
				std::vector<double> baseHead = cvd.GWHbaseValue(itwt->first);
				for (unsigned int i = 0; i < itwt->second.size(); ++i) {
					x = itwt->second[i] - baseHead[i];
					f += x;
					if (x > fmax)
						fmax = x;
					cntf = cntf + 1.0;
				}
			}
			std::cout << "Mean: " << f / cntf << ", Max: " << fmax << ", Sum: " << f << std::endl;
			fun.clear();
			fun.push_back(f);
			fun.push_back(area);
		}
	}
}
