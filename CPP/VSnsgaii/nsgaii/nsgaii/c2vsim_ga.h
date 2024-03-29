#ifndef C2VSIM_GA_H
#define C2VSIM_GA_H

#define AQUA false

#include <map>
#include <vector>
#include <string>

#include <boost/filesystem.hpp>

#include "c2vsim_io.h"
//#include "nsgaii_options.h"
//#include "nsgaii_helper.h"
//#include "nsgaii_core.h"

namespace C2VSIM {

	struct COST {
		COST() {
			Land = 0;
			Capital = 0;
			Water = 0;
			Lift = 0;
			Conveyance = 0;
		}
		double Land;
		double Capital;
		double Water;
		double Lift;
		double Conveyance;
		double Total() {
			return Land + Capital + Water + Lift + Conveyance;
		}
		void setCost(double pland, double maxQ, double totQ, double x_lift, double x_distance) {
			Land = pland * (maxQ*1000.0 / 75);
			Capital = 5000 * (maxQ*1000.0 / 75);
			Water = 0.0;
			Lift = std::abs(0.17 * 1.45 * x_lift * totQ*1000);
			Conveyance = 0.02 * x_distance * totQ*1000;
		}

		void setCost(double pland, double maxQ, std::vector<double>& TS, double x_lift, double x_distance, std::vector<double>& YearlyDiscount, int startIdx/*, int nElem*/) {
			Land = pland * (maxQ * 1000.0 / 75);
			Capital = 5000 * (maxQ * 1000.0 / 75);
			Water = 0.0;
            Lift = 0.0;
            Conveyance = 0.0;
            int idx = 0;
            int count_months = 0;
			for (int j = startIdx; j < TS.size(); ++j) {
			    if (idx >= YearlyDiscount.size()){
			        idx = YearlyDiscount.size() - 1;
			    }
			    Lift += std::abs(0.17 * 1.45 * x_lift * TS[j]*YearlyDiscount[idx] * 1000);
			    Conveyance += 0.02 * x_distance * TS[j]*YearlyDiscount[idx]*1000;
                count_months++;
                if (count_months >= 12){
                    count_months = 0;
                    idx++;
                }
			}
		}

		COST operator+(const COST& a) const {
			COST b;
			b.Land = Land + a.Land;
			b.Capital = Capital + a.Capital;
			b.Water = Water + a.Water;
			b.Lift = Lift + a.Lift;
			b.Conveyance = Conveyance + a.Conveyance;

			return b;
		}
	};



	class c2vsimData {
	public:
		c2vsimData(C2VSIM::OPTIONS::options& opt);
		void readInputFiles();
		void debugMsg();
		// For a given location in the encoding return the diversion node and the element 
		bool getNodeElemId(int i, int& node, int& elem);
		C2VSIM::DiversionData getDiversionData();
		void getDTS(int node, std::vector<double>& TS, double& maxq, double& totq);
        DVAR_info getValue(int elemId);
		std::string simulationExe();
		std::string budgetExe();
		std::string simulationPath() {return options.SimulationPath; }
		std::string BudgetOutputFile() { return options.BudgetOutputFile; }
		std::string DivSpecOpt() { return options.DivSpecOpt; }
		std::string DivDataOpt() { return options.DivDataOpt; }
		std::vector<double> GWHbaseValue(int node);
		int nsteps();
		double calcStorageChange(C2VSIM::GWbudTimeSeries& newGWBUB);
		void getDiscountFactors(std::vector<double>& DF){DF = DiscountFactors;}
        int getStartDivStep(){return options.StartDivStep;}

	private:
		C2VSIM::OPTIONS::options options;
		std::map<int, std::vector<int> > divElem;
		std::map<int, int> ElemDiv;
		std::map<int, int> idElem;
		std::map<int, DTSdata > DTS;
		C2VSIM::DiversionData divData;
		C2VSIM::GWbudTimeSeries GWBUD;
		std::map<int, std::vector<double> > GWH;
		std::map<int, ElemInfo> elemInfoMap;
		std::vector<double> DiscountFactors;
		std::map<int, DVAR_info> decisionVariablesMap;

		void makeElemDiv();
		void printElementMapping();
	};

	c2vsimData::c2vsimData(C2VSIM::OPTIONS::options& opt)
		:
		options(opt)
	{}

	void c2vsimData::readInputFiles() {
	    C2VSIM::READERS::readDiscountFactors(options.DiscountFile, DiscountFactors);
        C2VSIM::READERS::readDivTimeSeries(options.divTimeSeriesFile, DTS, options.StartDivStep);
        C2VSIM::READERS::readDecisionVariableMaps(options.DVMfile, decisionVariablesMap);
        //C2VSIM::READERS::readDivElems(options.divElemFile, divElem);
        //makeElemDiv();
        //printElementMapping();
		C2VSIM::READERS::readDivSpec(options.divSpecFile, divData);
		C2VSIM::READERS::readDivData(options.divDataFile, divData, 1056);
		C2VSIM::READERS::readGWHydOut(options.BaseWTfile, GWH, options.Nsteps);
		//C2VSIM::READERS::readElemInfo(options.ElementInfofile, elemInfoMap);
		C2VSIM::READERS::readGWBud(options.BaseGWbudFile, GWBUD, options.Nsteps);
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
	void c2vsimData::printElementMapping() {
		std::map<int, int >::iterator it, it1;
		std::ofstream outfile;
		outfile.open(options.ElemMapFile);
		for (it = idElem.begin(); it != idElem.end(); ++it) {
			it1 = ElemDiv.find(it->second);
			outfile << it->first << " " << it->second << " " << it1->second << std::endl;
		}
		outfile.close();
	}

	bool c2vsimData::getNodeElemId(int i, int& node, int& elem) {
        std::map<int, DVAR_info>::iterator itdvm;

        itdvm = decisionVariablesMap.find(i);
        if (itdvm != decisionVariablesMap.end()){
            node = itdvm->second.divID;
            elem = itdvm->second.elemId;
            return true;
        }
        else{
            elem = -9;
            node = -9;
            return false;
        }

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

	double c2vsimData::calcStorageChange(C2VSIM::GWbudTimeSeries& newGWBUB) {
		double out = 0.0;
		int n = newGWBUB.GwBudTS.size() - 1;
		// The budget files are in acft so we return TAF
		out = (newGWBUB.GwBudTS[n].S - GWBUD.GwBudTS[n].S)/1000.0;
		// Cumulative change
		//for (unsigned int i = 0; i < newGWBUB.GwBudTS.size(); ++i) {
		//	out += newGWBUB.GwBudTS[i].S - GWBUD.GwBudTS[i].S;
		//}
		return out;
	}

	C2VSIM::DiversionData c2vsimData::getDiversionData() {
		return divData;
	}

	void c2vsimData::getDTS(int node, std::vector<double>& TS, double& maxq, double& totq) {
		std::map<int, DTSdata >::iterator it = DTS.find(node);
		if (it != DTS.end()) {
			TS = it->second.TS;
			maxq = it->second.max;
			totq = it->second.Total;
		}
	}

    DVAR_info c2vsimData::getValue(int elemId) {
        std::map<int, DVAR_info>::iterator it;
		//std::map<int, ElemInfo>::iterator it;
		it = decisionVariablesMap.find(elemId);
		if (it != decisionVariablesMap.end()) {
			return it->second;
		}
		else {
			return DVAR_info();
		}
	}

	std::string c2vsimData::simulationExe() {
		return options.SimulationExe;
	}
	std::string c2vsimData::budgetExe() {
		return options.BudgetExe;
	}
	int c2vsimData::nsteps() {
		return options.Nsteps;
	}

	std::vector<double> c2vsimData::GWHbaseValue(int node) {
		// This code is not currently in use but if we used it then we have to make
		// sure that it will work if returns an empty vector
		std::map<int, std::vector<double> >::iterator it = GWH.find(node);
		if (it != GWH.end())
			return it->second;
		else
			return std::vector<double>();
	}

	void c2vsimData::debugMsg() {
		std::cout << GWH.size() << std::endl;
	}

	namespace OF {
	    struct ElemDV{
	        int elem; // This is the element id
	        int id; // This is the decision variable id
	    };

		double setupInputFiles(std::vector<double>& var, C2VSIM::c2vsimData& cvd) {
			// find the diversion nodes and the elements to apply the water
			std::vector<double> DF;
			cvd.getDiscountFactors(DF);
			std::map<int, std::vector<ElemDV> > nodeElemMap;
			std::map<int, std::vector<ElemDV> >::iterator it;
			int node, elem;
			for (unsigned int i = 0; i < var.size(); ++i) {
				if (var[i] > 0.5) {
					if (cvd.getNodeElemId(i+1, node, elem)) {
						it = nodeElemMap.find(node);
						if (it == nodeElemMap.end()) {
                            ElemDV ed;
                            ed.elem = elem;
                            ed.id = i+1;
							std::vector<ElemDV> temp;
							temp.push_back(ed);
							nodeElemMap.insert(std::pair<int, std::vector<ElemDV> >(node, temp));
						}
						else {
                            ElemDV ed;
                            ed.elem = elem;
                            ed.id = i+1;
							it->second.push_back(ed);
						}
					}
				}
			}

			// get a copy of the existing diversion data
			C2VSIM::DiversionData divData = cvd.getDiversionData();
			COST Totalcost;
			for (it = nodeElemMap.begin(); it != nodeElemMap.end(); ++it) {
				std::vector<double> TS;
				double maxq = 0;
				double totq = 0;
				cvd.getDTS(it->first, TS, maxq, totq);
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
				//div.IERELS = it->second;
				double NreceivElem = static_cast<double>(it->second.size());
				for (unsigned int i = 0; i < it->second.size(); ++i) {
					COST cost;
                    DVAR_info v = cvd.getValue(it->second[i].id);
					//cost.setCost(v.p_land, maxq / NreceivElem, totq / NreceivElem, v.x_lift, v.x_distance);
					cost.setCost(v.p_land, maxq / NreceivElem, TS, v.x_lift, v.x_dist,DF, cvd.getStartDivStep()/*,  NreceivElem*/);
                    div.IERELS.push_back(it->second[i].elem);
					div.FERELS.push_back(1);
					Totalcost = Totalcost + cost;
				}
				divData.appendDiversion(div, TS);
			}

			C2VSIM::WRITERS::writeDivSpec(cvd.DivSpecOpt(), divData);
			C2VSIM::WRITERS::writeDivData(cvd.DivDataOpt(), divData);
			return Totalcost.Total()/1000000;
		}

		void maxGWSTminCost(std::vector<double>& var, std::vector<double>& fun, C2VSIM::c2vsimData& cvd, int rank) {
            fun.clear();
            fun.push_back(sqrt((double)rand() / RAND_MAX));
            fun.push_back(sqrt((double)rand() / RAND_MAX));
            return;

#if AQUA
            // make sure these files do not exists
            boost::filesystem::remove("Results/CVground.bin");
            boost::filesystem::remove("Results/CVground.BUD");
            // keep the main path
            std::string main_dir = "/nfs/aqua-z0/giorgk/nsgaii";
            std::string run_dir = main_dir + "/RUNS/RUN";
            run_dir.append(std::to_string(rank));
            // Enter into the simulation path
            boost::filesystem::current_path(run_dir);
            //std::cout << "Rank: " << rank << boost::filesystem::current_path() << std::endl;
            std::string sim_command = "/opt/wine-stable/bin/wine ";
            sim_command.append(cvd.simulationExe()).append(" CVsimAqua.in >/dev/null");
#else
 
            std::string main_dir = boost::filesystem::current_path().string();
            // Enter into the simulation path
			boost::filesystem::current_path(cvd.simulationPath());
			std::cout << boost::filesystem::current_path() << std::endl;
			std::string sim_command = cvd.simulationExe();
			sim_command.append(" CVsim.in");
#endif
            double costOF = setupInputFiles(var, cvd);

            if (!system(NULL))
                std::cout << "Rank: " << rank << " System not ready" << std::endl;

			int sys = system(sim_command.c_str());
			if (!boost::filesystem::exists("Results/CVground.bin")) {
                std::cout << "Rank: " << rank << " Can't find the CVground.bin file" << std::endl;
				fun.clear();
				fun.push_back(10000000);
				fun.push_back(10000000);
                boost::filesystem::current_path(main_dir);
				return;
			}

#if AQUA
            std::string bud_command = "/opt/wine-stable/bin/wine ";
            bud_command.append(cvd.budgetExe()).append(" CVBudgetAqua.in >/dev/null"); //
#else
            std::string bud_command = cvd.budgetExe();
			bud_command.append(" CVBudget.in");
#endif

			sys = system(bud_command.c_str());
			if (!boost::filesystem::exists("Results/CVground.BUD")) {
                std::cout << "Rank: " << rank << " Can't find the CVground.BUD file" << std::endl;
				fun.clear();
				fun.push_back(10000000);
				fun.push_back(10000000);
                boost::filesystem::current_path(main_dir);
				return;
			}

			C2VSIM::GWbudTimeSeries simGBbud;
#if AQUA
            C2VSIM::READERS::readGWBud("Results/CVground.BUD", simGBbud, cvd.nsteps());
#else
			C2VSIM::READERS::readGWBud(cvd.BudgetOutputFile(), simGBbud, cvd.nsteps());
#endif
			double envOF = cvd.calcStorageChange(simGBbud);
			boost::filesystem::remove("Results/CVground.bin");
			boost::filesystem::remove("Results/CVground.BUD");
			fun.clear();
			fun.push_back(-envOF);
			fun.push_back(costOF);
			boost::filesystem::current_path(main_dir);
		}

		/*void maxWTminArea(std::vector<double>& var, std::vector<double>& fun, C2VSIM::c2vsimData& cvd) {
			NSGAII::SingletonRealGenerator* RG = RG->getInstance();
			*//*
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
					//double a = cvd.getValue(it->second[i]);
					//area += a;
					//div.FERELS.push_back(a);
				}
				divData.appendDiversion(div, TS);
			}
			//fun.clear();
			//fun.push_back(-RG->randomNumber(100000,600000));
			//fun.push_back(area);
			//return;
			C2VSIM::WRITERS::writeDivSpec("d:/giorgk/Documents/GitHub/C2VsimCG/RunC2Vsim/tempSpec.dat", divData);
			C2VSIM::WRITERS::writeDivData("d:/giorgk/Documents/GitHub/C2VsimCG/RunC2Vsim/tempData.dat", divData);
			*//*
			double costOF = setupInputFiles(var, cvd);

			std::string main_dir = boost::filesystem::current_path().string();
			// Enter into the simulation path
			boost::filesystem::current_path("d:/giorgk/Documents/GitHub/C2VsimCG/RunC2Vsim");
			std::cout << boost::filesystem::current_path() << std::endl;
			
			std::string sim_command = cvd.simulationExe();
			//std::string bud_command = cvd.budgetExe();
			sim_command.append(" CVsim.in");
			//bud_command.append(" CVBudget.in");


			system(sim_command.c_str());
			//system(bud_command.c_str());
			//Read the output while in the path
			std::map<int, std::vector<double> > simGWhyd;
			C2VSIM::GWbudTimeSeries simGBbud;
			C2VSIM::READERS::readGWHydOut("d:/giorgk/Documents/GitHub/C2VsimCG/RunC2Vsim/Results/CVGWhyd.out", simGWhyd, cvd.nsteps());
			//C2VSIM::READERS::readGWBud("d:/giorgk/Documents/GitHub/C2VsimCG/RunC2Vsim/Results/CVground.BUD", simGBbud, cvd.nsteps());

			std::cout << "Model finished" << std::endl;
			boost::filesystem::current_path(main_dir);
			std::cout << boost::filesystem::current_path() << std::endl;
			bool tf = true;
			std::cout << "all is" << tf << std::endl;

			std::map<int, std::vector<double> >::iterator itwt;
			//double f = 0;
			//double cntf = 0;
			//double fmax = -999999;
			//double x;
			//for (itwt = simGWhyd.begin(); itwt != simGWhyd.end(); ++itwt) {
			//	std::vector<double> baseHead = cvd.GWHbaseValue(itwt->first);
			//	for (unsigned int i = 0; i < itwt->second.size(); ++i) {
			//		x = itwt->second[i] - baseHead[i];
			//		//if (x < -50)
			//		//	continue;
			//		f += x;
			//		if (x > fmax)
			//			fmax = x;
			//		cntf = cntf + 1.0;
			//	}
			//}

			double fp = 0;
			double rate = 0.9;
			double dx = 3.28084;
			double x, xp, w, xs, xe;
			for (itwt = simGWhyd.begin(); itwt != simGWhyd.end(); ++itwt) {
				std::vector<double> baseHead = cvd.GWHbaseValue(itwt->first);
				for (unsigned int i = 0; i < itwt->second.size(); ++i) {
					x = itwt->second[i] - baseHead[i];
					if (x < 0)
						fp += x;
					else if (x > 0) {
						if (x > 10)
							bool tf = true;
						xp = 0;
						xs = 0;
						xe = dx;
						w = 1;
						while (xp < x) {
							if (xe > x)
								xe = x;
							fp += (xe - xs) * w;
							xs += dx;
							xe += dx;
							w = w * rate;
							xp += dx;
						}
					}
					//if (x < -60)
					//	continue;
					//if (x > fmax)
					//	fmax = x;
					//cntf = cntf + 1.0;
				}
			}

			//std::cout << "Mean: " << f / cntf << ", Max: " << fmax << ", Sum: " << f << std::endl;
			fun.clear();
			fun.push_back(-fp);
			fun.push_back(costOF);
		}*/
	}
}

#endif // C2VSIM_GA_H