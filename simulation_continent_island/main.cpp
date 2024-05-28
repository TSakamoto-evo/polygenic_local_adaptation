#include "parameter.hpp"
#include "population.hpp"
#include <unordered_map>
#include <chrono>

int main(int argc, char* argv[]){
  int pop_size_island, pop_size_continent;
  double mut_rate, mig_rate, theta, gamma, sd, s;

  if(argc == 9){
    sscanf(argv[1], "%d", &pop_size_island);
    sscanf(argv[2], "%d", &pop_size_continent);
    sscanf(argv[3], "%lf", &mut_rate);
    sscanf(argv[4], "%lf", &mig_rate);
    sscanf(argv[5], "%lf", &theta);
    sscanf(argv[6], "%lf", &gamma);
    sscanf(argv[7], "%lf", &sd);
    sscanf(argv[8], "%lf", &s);
  }else{
    std::cerr << "Error!" << std::endl;
    std::exit(1);
  }

  Parameter para;
  para.pop_size_island = pop_size_island;
  para.pop_size_continent = pop_size_continent;
  para.mut_rate_island = mut_rate;
  para.mut_rate_continent = mut_rate;
  para.mig_rate = mig_rate;

  para.num_loci = 1000;
  para.reco_rate = 10.0;

  para.theta = theta;
  para.sd_limit = 3.0;
  para.sd = sd;
  para.s = s;
  para.gamma = gamma;

  para.bin = 1;

  std::ofstream ofs1("mean_var.txt", std::ios::app);
  std::ofstream ofs2("bin_contribution.txt", std::ios::app);

  int generation = 0;
  Population pop(para, generation);

  std::chrono::system_clock::time_point last_output;
  last_output = std::chrono::system_clock::now();

  std::vector<std::vector<double>> regi;
  std::vector<int> regi_time;
  std::vector<std::vector<double>> bin_contribution;
  std::vector<int> bin_time;

  for(int i = generation; i <= 200000; i++){
    pop.one_generation();

    if(i % 100 == 0){
      // register mean and variance of phenotypic value and fitness every 100 generations
      std::vector<double> ret;
      pop.calculate_mean_var(ret);
      regi.push_back(ret);
      regi_time.push_back(i);
    }

    if(i % 1000 == 0){
      // register difference of mean genotypic value at each locus between populations every 1000 generations
      std::vector<double> ret;
      pop.calculate_bin_contribution(ret);
      bin_contribution.push_back(ret);
      bin_time.push_back(i);
    }

    // back-up the population state every 30 minutes
    std::chrono::system_clock::time_point now = std::chrono::system_clock::now();
    auto time = now - last_output;
    double minute = std::chrono::duration_cast<std::chrono::minutes>(time).count();

    if(minute > 30){
      pop.regi_state();

      // output mean and variance of phenotype and fitness
      for(int j = 0; j < static_cast<int>(regi_time.size()); j++){
        ofs1 << regi_time.at(j);
        for(const auto& k: regi.at(j)){
          ofs1 << "\t" << k;
        }
        ofs1 << std::endl;
      }

      regi.clear();
      regi_time.clear();

      // output locus-based genotypic value
      for(int j = 0; j < static_cast<int>(bin_time.size()); j++){
        ofs2 << bin_time.at(j);
        for(const auto& k: bin_contribution.at(j)){
          ofs2 << "\t" << std::setprecision(4) << k;
        }
        ofs2 << std::endl;
      }

      bin_contribution.clear();
      bin_time.clear();

      std::system("mv regi_state2.txt regi_state.txt");

      // renew the time point of last output
      last_output = now;
    }
  }

  for(int j = 0; j < static_cast<int>(regi_time.size()); j++){
    ofs1 << regi_time.at(j);
    for(const auto& k: regi.at(j)){
      ofs1 << "\t" << k;
    }
    ofs1 << std::endl;
  }

  for(int j = 0; j < static_cast<int>(bin_time.size()); j++){
    ofs2 << bin_time.at(j);
    for(const auto& k: bin_contribution.at(j)){
      ofs2 << "\t" << std::setprecision(4) << k;
    }
    ofs2 << std::endl;
  }

  return(0);
}
