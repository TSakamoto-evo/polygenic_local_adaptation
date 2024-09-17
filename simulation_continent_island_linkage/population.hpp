#ifndef POPULATION
#define POPULATION

#include "parameter.hpp"
#include "individual.hpp"
#include <vector>
#include <algorithm>
#include <random>
#include <iostream>
#include <fstream>
#include <unordered_map>
#include <chrono>
#include <filesystem>
#include <string>
#include <sstream>

class Population{
private:
  Parameter para;
  std::vector<Individual> pop_island;
  std::vector<Individual> pop_continent;

  std::vector<Individual> parent_island;
  std::vector<Individual> parent_continent;

  std::vector<int> order;
  std::vector<double> fitness_island;
  std::vector<double> fitness_continent;

  std::mt19937 mt;
  std::bernoulli_distribution p;
  
  std::exponential_distribution<> det_reco;
  std::poisson_distribution<> det_mut_island;
  std::poisson_distribution<> det_mut_continent;
  std::poisson_distribution<> det_mig;
  std::uniform_int_distribution<> choose_ind_island;
  std::uniform_int_distribution<> choose_ind_continent;
  std::uniform_int_distribution<> choose_locus;
  std::normal_distribution<> choose_gamma_i;

  int gen;

public:
  Population(const Parameter input_para, int& generation);
  void ret_reco_genotype_island(const int parent,
    std::vector<double>& ret_small, double& ret_major);
  void ret_reco_genotype_continent(const int parent,
    std::vector<double>& ret_small, double& ret_major);
  void selection();
  void mutation_island();
  void mutation_continent();
  void one_generation();

  void calculate_mean_var(std::vector<double>& ret) const;
  void calculate_bin_contribution(std::vector<double>& ret) const;
  void regi_state();
  void initialize_from_regi();
};


#endif
