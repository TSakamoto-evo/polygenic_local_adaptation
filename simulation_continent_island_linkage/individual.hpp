#ifndef INDIVIDUAL
#define INDIVIDUAL

#include <vector>
#include <iostream>

class Individual{
private:
  // background region
  std::vector<double> genotype1;
  std::vector<double> genotype2;
  // focal locus
  double effect1;
  double effect2;

public:
  Individual(const std::vector<double>& input_genotype1,
    const std::vector<double>& input_genotype2, 
    const double input_effect1, const double input_effect2);
  void input_genotype(const std::vector<double>& input_genotype, 
    const double input_effect, const bool which);
  double ret_phenotypic_value() const;
  void ret_haplotype(std::vector<double>& ret_small, double& ret_major, 
    const std::vector<double>& switch_pos, const bool which) const;
  void mutation(const bool which, const int locus, const double add);
  double ret_major_value() const{ return(effect1 + effect2); };
  void return_genotype(std::vector<double>& geno1, 
    std::vector<double>& geno2, double& eff1, double& eff2) const;
};

#endif
