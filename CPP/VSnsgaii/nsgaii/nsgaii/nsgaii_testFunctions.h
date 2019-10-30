#pragma once

#include <vector>
#include <cmath>

namespace NSGAII {
	const double invsq3 = 1 / std::sqrt(3.0);
	void FonsecaFleming(std::vector<double>& var, std::vector<double>& fun) {
		double a = 0.0;
		double b = 0.0;
		for (unsigned int i = 0; i < 3; ++i) {
			a += (var[i] - invsq3) * (var[i] - invsq3);
			b += (var[i] + invsq3) * (var[i] + invsq3);
		}
		fun.push_back(1 - std::exp(-a));
		fun.push_back(1 - std::exp(-b));
	}

	void Kursawe(std::vector<double>& var, std::vector<double>& fun) {
		//std::cout << "x=[" << var[0] << ";" << var[1] << ";" << var[2] << "]" << std::endl;
		double a = 0;
		double b = 0;
		for (int i = 0; i < 3; ++i) {
			if (i < 2) {
				a += -10.0 * std::exp(-0.2 * std::sqrt(var[i] * var[i] + var[i + 1] * var[i + 1]));
			}
			b += std::pow(std::abs(var[i]), 0.8) + 5 * std::sin(var[i]*var[i]*var[i]);
		}
		fun.push_back(a);
		fun.push_back(b);
	}
}
