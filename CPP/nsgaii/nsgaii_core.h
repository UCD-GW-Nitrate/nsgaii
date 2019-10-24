#ifndef NSGAII_CORE_H
#define NSGAII_CORE_H

#include <vector>
#include <map>
#include <boost/random/linear_congruential.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/generator_iterator.hpp>

#include "nsgaii_options.h"

typedef boost::minstd_rand base_generator_type;
namespace NSGAII {

    class Population{
    public:
        Population(NSGAII::options &opt);
        void initializePopulation();

    private:
        NSGAII::options options;
        std::map<int, std::vector<double> > members;




    };



    Population::Population(NSGAII::options &opt)
        :
          options(opt)
    {}

    void Population::initializePopulation(){
        std::cout << "Initialize Population" << std::endl;
        std::map<int, std::vector<double> >::iterator it;

        for (int ivar = 0; ivar < options.ProblemSize; ++ivar){
            boost::uniform_real<>uniform(options.LowerBound[ivar], options.UpperBound[ivar]);
            base_generator_type generator(42);
            boost::variate_generator<base_generator_type&, boost::uniform_real<> > uni(generator, uniform);
            for (int i = 0 ; i < options.PopulationSize; ++i){
                it = members.find(i);
                if (it == members.end()){
                    std::vector<double> x;
                    x.push_back(uni());
                    members.insert(std::pair<int,std::vector<double> >(i, x));
                }
                else {
                    it->second.push_back(uni());
                }
            }
        }

        std::cout << "Done" << std::endl;
    }


}
#endif // NSGAII_CORE_H
