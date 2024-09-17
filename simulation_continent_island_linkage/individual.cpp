#include "individual.hpp"

Individual::Individual(const std::vector<double>& input_genotype1,
    const std::vector<double>& input_genotype2,
    const double input_effect1, const double input_effect2){

  genotype1 = input_genotype1;
  genotype2 = input_genotype2;

  effect1 = input_effect1;
  effect2 = input_effect2;
}

void Individual::input_genotype(const std::vector<double>& input_genotype, 
    const double input_effect, const bool which){
  
  if(which){
    genotype1 = input_genotype;
    effect1 = input_effect;
  }else{
    genotype2 = input_genotype;
    effect2 = input_effect;
  }
}

double Individual::ret_phenotypic_value() const{
  double ret = effect1 + effect2;

  for(int i = 0; i < static_cast<int>(genotype1.size()); i++){
    ret += genotype1[i];
    ret += genotype2[i];
  }

  return(ret);
}

void Individual::ret_haplotype(std::vector<double>& ret_small, double& ret_major, 
    const std::vector<double>& switch_pos, const bool which) const{

  // background region
  // determine which allele is inherited at the first locus
  bool now_which = which;
  int locus = 0;

  int total_locus = static_cast<int>(genotype1.size());
  double focal_pos = (total_locus - 1.0) / 2.0;

  for(int i = 0; i < static_cast<int>(switch_pos.size()); i++){
    if(now_which){
      while(locus < static_cast<int>(genotype1.size()) && locus < switch_pos[i]){
        // first genome is inherited
        ret_small[locus] = genotype1[locus];
        locus++;
      }
    }else{
      while(locus < static_cast<int>(genotype1.size()) && locus < switch_pos[i]){
        // second genome is inherited
        ret_small[locus] = genotype2[locus];
        locus++;
      }
    }

    if((i == 0 && focal_pos < switch_pos[i]) || 
      (switch_pos[i - 1] <= focal_pos && focal_pos < switch_pos[i])){
      if(now_which){
        ret_major = effect1;
      }else{
        ret_major = effect2;
      }
    }

    // crossing over
    now_which = !now_which;
  }
}

void Individual::mutation(const bool which, const int locus, const double add){
  if(which){
    genotype1[locus] += add;
  }else{
    genotype2[locus] += add;
  }
}

void Individual::return_genotype(std::vector<double>& geno1, 
    std::vector<double>& geno2, double& eff1, double& eff2) const{
  
  geno1 = genotype1;
  geno2 = genotype2;
  eff1 = effect1;
  eff2 = effect2;
}