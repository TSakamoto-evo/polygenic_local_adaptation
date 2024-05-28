#ifndef PARAMETER
#define PARAMETER

class Parameter{
public:
  int pop_size_island;
  int pop_size_continent;
  double mut_rate_island;
  double mut_rate_continent;
  double mig_rate;

  int num_loci;
  double reco_rate;

  double theta;
  double s;
  double sd_limit;
  double sd;
  double gamma;

  int bin;
};

#endif
