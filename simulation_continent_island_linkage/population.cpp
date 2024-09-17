#include "population.hpp"

Population::Population(const Parameter input_para, int& generation){
  para = input_para;

  std::random_device seed;
  mt.seed(seed());

  pop_island.clear();

  std::vector<double> tmp(para.num_loci);
  for(int i = 0; i < para.pop_size_island; i++){
    // start from the state where the focal locus is diverged (mutant allele fixed)
    pop_island.emplace_back(tmp, tmp, para.gamma, para.gamma);
  }

  pop_continent.clear();

  for(int i = 0; i < para.pop_size_continent; i++){
    // start from the state where the focal locus is diverged (wildtype allele fixed)
    pop_continent.emplace_back(tmp, tmp, 0.0, 0.0);
  }

  gen = 0;

  // if there is an unfinished simulation data, use it for initialization and restart
  if(std::filesystem::is_regular_file("regi_state.txt")){
    initialize_from_regi();
  }

  parent_island = pop_island;
  parent_continent = pop_continent;

  // determine distance between adjacent crossing over points
  // (in the unit of distance between adjacent locus in the background region)
  std::exponential_distribution<> tmp_det_reco(1.0 * para.reco_rate / para.num_loci);
  det_reco = tmp_det_reco;

  // determine the number of new mutations in the island population
  std::poisson_distribution<> tmp_det_mut_island(2.0 * para.pop_size_island * para.mut_rate_island);
  det_mut_island = tmp_det_mut_island;

  // determine the number of new mutations in the continent population
  std::poisson_distribution<> tmp_det_mut_continent(2.0 * para.pop_size_continent * para.mut_rate_continent);
  det_mut_continent = tmp_det_mut_continent;

  // determine the number of migrant parents
  std::poisson_distribution<> tmp_det_mig(2.0 * para.pop_size_island * para.mig_rate);
  det_mig = tmp_det_mig;

  std::uniform_int_distribution<> tmp_choose_ind_island(0, para.pop_size_island - 1);
  choose_ind_island = tmp_choose_ind_island;

  std::uniform_int_distribution<> tmp_choose_ind_continent(0, para.pop_size_continent - 1);
  choose_ind_continent = tmp_choose_ind_continent;

  std::uniform_int_distribution<> tmp_choose_locus(0, para.num_loci - 1);
  choose_locus = tmp_choose_locus;

  std::normal_distribution<> tmp_choose_gamma_i(0.0, para.sd);
  choose_gamma_i = tmp_choose_gamma_i;

  std::bernoulli_distribution tmp_p(0.5);
  p = tmp_p;

  std::vector<int> tmp_order(2 * para.pop_size_island);
  order = tmp_order;
  std::iota(order.begin(), order.end(), 0);

  std::vector<double> tmp_fitness_island(para.pop_size_island);
  fitness_island = tmp_fitness_island;

  std::vector<double> tmp_fitness_continent(para.pop_size_continent);
  fitness_continent = tmp_fitness_continent;

  generation = gen;
}

void Population::ret_reco_genotype_island(const int parent,
  std::vector<double>& ret_small, double& ret_major){

  // determine which allele is inherited at the first locus in the background region
  bool which = p(mt);

  // determine recombination breakpoints in the background region
  std::vector<double> reco_pos;
  reco_pos.reserve(5 * para.reco_rate);

  double start = 0.0;
  while(start < para.num_loci){
    double dist = det_reco(mt);
    reco_pos.push_back(start + dist);
    start += dist;
  }

  parent_island[parent].ret_haplotype(ret_small, ret_major, reco_pos, which);
}

void Population::ret_reco_genotype_continent(const int parent,
  std::vector<double>& ret_small, double& ret_major){

  // determine which allele is inherited at the first locus in the background region
  bool which = p(mt);

  // determine recombination breakpoints in the background region
  std::vector<double> reco_pos;
  reco_pos.reserve(5 * para.reco_rate);

  double start = 0.0;
  while(start < para.num_loci){
    double dist = det_reco(mt);
    reco_pos.push_back(start + dist);
    start += dist;
  }

   parent_continent[parent].ret_haplotype(ret_small, ret_major, reco_pos, which);
}

void Population::selection(){
  // calculate fitness of each individual in the island population
  for(int i = 0; i < para.pop_size_island; i++){
    double phenotype = pop_island[i].ret_phenotypic_value();
    fitness_island[i] = std::exp(-para.s * (phenotype - para.theta) * (phenotype - para.theta));
  }

  // calculate fitness of each individual in the continent population
  for(int i = 0; i < para.pop_size_continent; i++){
    double phenotype = pop_continent[i].ret_phenotypic_value();
    fitness_continent[i] = std::exp(-para.s * phenotype * phenotype);
  }

  parent_island = pop_island;
  parent_continent = pop_continent;

  std::discrete_distribution<> choose_parent_island(fitness_island.begin(), fitness_island.end());
  std::discrete_distribution<> choose_parent_continent(fitness_continent.begin(), fitness_continent.end());
  std::shuffle(order.begin(), order.end(), mt);

  std::vector<double> geno_small(para.num_loci);
  double geno_major = 0.0;

  int mig_num = det_mig(mt);

  // parent is a migrant from the continent population 
  for(int i = 0; i < mig_num; i++){
    int ind = order[i] / 2;
    bool which = order[i] % 2;

    int parent = choose_parent_continent(mt);
    ret_reco_genotype_continent(parent, geno_small, geno_major);

    pop_island[ind].input_genotype(geno_small, geno_major, which);
  }

  // parent is a nonmigrant
  for(int i = mig_num; i < 2 * para.pop_size_island; i++){
    int ind = order[i] / 2;
    bool which = order[i] % 2;

    int parent = choose_parent_island(mt);
    ret_reco_genotype_island(parent, geno_small, geno_major);

    pop_island[ind].input_genotype(geno_small, geno_major, which);
  }

  // reproduction in the continent
  for(int i = 0; i < 2 * para.pop_size_continent; i++){
    int ind = i / 2;
    bool which = i % 2;

    int parent = choose_parent_continent(mt);
    ret_reco_genotype_continent(parent, geno_small, geno_major);

    pop_continent[ind].input_genotype(geno_small, geno_major, which);
  }
}

void Population::mutation_island(){
  int num_mut_island = det_mut_island(mt);

  for(int i = 0; i < num_mut_island; i++){
    int ind = choose_ind_island(mt);
    int locus = choose_locus(mt);
    bool which = p(mt);
    double add = choose_gamma_i(mt);
    
    // truncation in effect size
    while(std::abs(add) > para.sd * para.sd_limit){
      add = choose_gamma_i(mt);
    }

    pop_island[ind].mutation(which, locus, add);
  }
}

void Population::mutation_continent(){
  int num_mut_continent = det_mut_continent(mt);

  for(int i = 0; i < num_mut_continent; i++){
    int ind = choose_ind_continent(mt);
    int locus = choose_locus(mt);
    bool which = p(mt);
    double add = choose_gamma_i(mt);
    
    // truncation in effect size
    while(std::abs(add) > para.sd * para.sd_limit){
      add = choose_gamma_i(mt);
    }

    pop_continent[ind].mutation(which, locus, add);
  }
}

void Population::one_generation(){
  selection();
  mutation_island();
  mutation_continent();

  gen++;
}

void Population::calculate_mean_var(std::vector<double>& ret) const{
  // calculate quantities in the island population
  {
    double mean_fit_island = 0.0;
    double var_fit_island = 0.0;

    double mean_pheno_island = 0.0;
    double var_pheno_island = 0.0;

    double freq_island = 0.0;

    for(int i = 0; i < para.pop_size_island; i++){
      double phenotype = parent_island[i].ret_phenotypic_value();
      double fit = std::exp(-para.s * (phenotype - para.theta) * (phenotype - para.theta));
      double major = parent_island[i].ret_major_value();

      mean_fit_island += fit;
      var_fit_island += fit * fit;
      mean_pheno_island += phenotype;
      var_pheno_island += phenotype * phenotype;
      freq_island += major;
    }

    mean_fit_island /= para.pop_size_island;
    var_fit_island /= para.pop_size_island;
    var_fit_island -= mean_fit_island * mean_fit_island;

    mean_pheno_island /= para.pop_size_island;
    var_pheno_island /= para.pop_size_island;
    var_pheno_island -= mean_pheno_island * mean_pheno_island;

    freq_island /= (2.0 * para.pop_size_island * para.gamma);

    ret.clear();
    ret.push_back(mean_fit_island);
    ret.push_back(var_fit_island);
    ret.push_back(mean_pheno_island);
    ret.push_back(var_pheno_island);
    ret.push_back(freq_island);
  }

  // calculate quantities in the continent population
  {
    double mean_fit_continent = 0.0;
    double var_fit_continent = 0.0;

    double mean_pheno_continent = 0.0;
    double var_pheno_continent = 0.0;

    double freq_continent = 0.0;

    for(int i = 0; i < para.pop_size_continent; i++){
      double phenotype = parent_continent[i].ret_phenotypic_value();
      double fit = std::exp(-para.s * phenotype * phenotype);
      double major = parent_continent[i].ret_major_value();

      mean_fit_continent += fit;
      var_fit_continent += fit * fit;
      mean_pheno_continent += phenotype;
      var_pheno_continent += phenotype * phenotype;
      freq_continent += major;
    }

    mean_fit_continent /= para.pop_size_continent;
    var_fit_continent /= para.pop_size_continent;
    var_fit_continent -= mean_fit_continent * mean_fit_continent;

    mean_pheno_continent /= para.pop_size_continent;
    var_pheno_continent /= para.pop_size_continent;
    var_pheno_continent -= mean_pheno_continent * mean_pheno_continent;

    freq_continent /= (2.0 * para.pop_size_continent * para.gamma);

    ret.push_back(mean_fit_continent);
    ret.push_back(var_fit_continent);
    ret.push_back(mean_pheno_continent);
    ret.push_back(var_pheno_continent);
    ret.push_back(freq_continent);
  }
}

void Population::calculate_bin_contribution(std::vector<double>& ret) const{
  std::vector<double> geno1, geno2;
  double eff1, eff2;

  std::vector<double> tmp_island(para.num_loci / para.bin + 1);

  for(int i = 0; i < para.pop_size_island; i++){
    parent_island[i].return_genotype(geno1, geno2, eff1, eff2);

    for(int j = 0; j < para.num_loci; j++){
      tmp_island[j / para.bin] += geno1[j] + geno2[j];
    }
    tmp_island.back() += eff1 + eff2;
  }

  for(auto& i: tmp_island){
    i /= 2.0 * para.pop_size_island;
  }

  std::vector<double> tmp_continent(para.num_loci / para.bin + 1);

  for(int i = 0; i < para.pop_size_continent; i++){
    parent_continent[i].return_genotype(geno1, geno2, eff1, eff2);

    for(int j = 0; j < para.num_loci; j++){
      tmp_continent[j / para.bin] += geno1[j] + geno2[j];
    }
    tmp_continent.back() += eff1 + eff2;
  }

  for(auto& i: tmp_continent){
    i /= 2.0 * para.pop_size_continent;
  }

  ret = tmp_island;

  for(int i = 0; i < static_cast<int>(ret.size()); i++){
    ret[i] -= tmp_continent[i];
  }
}

void Population::regi_state(){
  std::ofstream ofs("regi_state2.txt");
  ofs << gen << std::endl;

  std::vector<double> geno1;
  std::vector<double> geno2;
  double eff1, eff2;

  for(int i = 0; i < para.pop_size_island; i++){
    pop_island[i].return_genotype(geno1, geno2, eff1, eff2);

    ofs << "i\t" << 2 * i << "\t" << eff1;
    for(int j = 0; j < para.num_loci; j++){
      ofs << "\t" << geno1[j];
    }
    ofs << std::endl;

    ofs << "i\t" << 2 * i + 1 << "\t" << eff2;
    for(int j = 0; j < para.num_loci; j++){
      ofs << "\t" << geno2[j];
    }
    ofs << std::endl;
  }

  for(int i = 0; i < para.pop_size_continent; i++){
    pop_continent[i].return_genotype(geno1, geno2, eff1, eff2);

    ofs << "c\t" << 2 * i << "\t" << eff1;
    for(int j = 0; j < para.num_loci; j++){
      ofs << "\t" << geno1[j];
    }
    ofs << std::endl;

    ofs << "c\t" << 2 * i + 1 << "\t" << eff2;
    for(int j = 0; j < para.num_loci; j++){
      ofs << "\t" << geno2[j];
    }
    ofs << std::endl;
  }
}

void Population::initialize_from_regi(){
  std::ifstream ifs("regi_state.txt");
  if(!ifs){
    std::cerr << "Fail to open the file!" << std::endl;
    std::exit(1);
  }

  std::vector<std::string> lines;
  std::string line;
  int num_line = 0;

  while (getline(ifs, line)){
    std::istringstream iss(line);
    std::string tmp_list;
    std::vector<std::string> list;

    while(getline(iss, tmp_list, '\t')){
      list.push_back(tmp_list);
    }

    if(num_line == 0){
      gen = std::stoi(list[0]);
    }else{
      int index = std::stoi(list[1]);
      double effect = std::stod(list[2]);

      std::vector<double> geno;
      for(int i = 3; i < static_cast<int>(list.size()); i++){
        geno.push_back(std::stod(list[i]));
      }

      int ind = index / 2;
      int which = (index + 1) % 2;

      if(list[0] == "i"){
        pop_island[ind].input_genotype(geno, effect, which);
      }else{
        pop_continent[ind].input_genotype(geno, effect, which);
      }
    }
    num_line++;
  }
}
