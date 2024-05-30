#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

void calculate_beta_vec(const double theta2, const double sigma, 
  const double sd_cut, const double s, const double m,
  const double n1, const double n2, 
  std::vector<double>& vec_z2, std::vector<double>& vec_beta);
double integrate_x_t1(const double s, const double n1,
  const double gamma);
double f_g1(const double a, const double x);
double int_g1(const double a, const double x);
double x_t1(const double a, const double x);
double integrate_gamma_x_t2(const double s, const double gamma, const double z2,
  const double theta2, const double m, const double n2);
double integrate_gamma_two_x_1minusx_t1(const double s, const double gamma, const double n1);
double integrate_gamma_two_x_1minusx_t2(const double s, const double gamma, const double z2,
  const double theta2, const double m, const double n2);
double interpolate_beta(const double z2, const std::vector<double>& vec_z2, 
  const std::vector<double>& vec_beta);
double calculate_variance1(const double sigma, 
  const double sd_cut, const double s, const double n1);
double calculate_variance2(const double theta2, const double sigma, 
  const double sd_cut, const double s, const double m,
  const double n1, const double n2, double z2);

int main(){
  double n1 = 10000.0;
  double n2 = 10000.0;
  double s = 0.025;
  double m = 0.001;
  double gamma_f = 0.2;
  double sigma = 0.02;
  double theta2 = 1.0;

  double sd_cut = 3.0;

  std::vector<double> vec_z2, vec_beta;

  // first calculate y2 as a function of z2
  calculate_beta_vec(theta2, sigma, sd_cut, s, m, n1, n2, vec_z2, vec_beta);

  std::ofstream ofs("equilibrium_value.txt");
  ofs << "mu\tz2_nodiv\tvar_pop1\tvar_pop2_nodiv\tz2_withdiv\tp2_withdiv\ty2_withdiv\tvar_pop2_withdiv" << std::endl;

  for(int i = 0; i <= 250; i++){
    double mu = std::pow(10.0, -i / 50.0);
    ofs << mu;

    // equilibrium in the absense of the divergence at the focal locus
    {
      double z2_min = vec_z2.at(0);
      double z2_max = vec_z2.back();

      for(int j = 0; j <= 50; j++){
        double z2_mid = (z2_min + z2_max) / 2.0;
        double y2 = mu * interpolate_beta(z2_mid, vec_z2, vec_beta);

        if(2.0 * y2 > z2_mid){
          z2_min = z2_mid;
        }else{
          z2_max = z2_mid;
        }

        if(j == 50){
          double variance1 = mu * calculate_variance1(sigma, sd_cut, s, n1);
          double variance2 = mu * calculate_variance2(theta2, sigma, sd_cut, s, m, n1, n2, z2_mid);
          ofs << "\t" << z2_mid << "\t" << variance1 << "\t" << variance2;
        }
      }
    }

    // equilibrium starting from a divergent focal locus

    // check the stability of the cluster divergence
    bool stability;
    {
      double z2_limit, p2_limit;
      if(gamma_f * gamma_f >= m / 2.0 / s){
        z2_limit = theta2 + gamma_f / 2.0 - std::sqrt(2.0 * m / s);
        p2_limit = 1.0 - std::sqrt(m / 2.0 / s / gamma_f / gamma_f);
      }else{
        z2_limit = theta2 - gamma_f / 2.0 - m / 2.0 / s / gamma_f;
        p2_limit = 0.0;
      }

      if(z2_limit < 0.0){
        stability = 0;
      }else if(z2_limit > theta2){
        stability = 1;
      }else{
        double y2 = mu * interpolate_beta(z2_limit, vec_z2, vec_beta);

        if(2.0 * y2 + 2.0 * gamma_f * p2_limit > z2_limit){
          stability = 0;
        }else{
          stability = 1;
        }
      }
    }

    if(stability == 0){
      ofs << "\t" << "nan" << "\t" << "nan" << "\t" << "nan" << "\t" << "nan" << std::endl;
    }else{
      double z2_min = vec_z2.at(0);
      double z2_max;

      if(gamma_f * gamma_f >= m / 2.0 / s){
        z2_max = theta2 + gamma_f / 2.0 - std::sqrt(2.0 * m / s);
      }else{
        z2_max = theta2 - gamma_f / 2.0 - m / 2.0 / s / gamma_f;
      }

      if(z2_max > vec_z2.back()){
        z2_max = vec_z2.back();
      }

      for(int j = 0; j <= 50; j++){
        double z2_mid = (z2_min + z2_max) / 2.0;
        double y2 = mu * interpolate_beta(z2_mid, vec_z2, vec_beta);
        double p2 = 3.0 / 4.0 - (theta2 - z2_mid) / 2.0 / gamma_f +
          std::sqrt((gamma_f + 2.0 * (theta2 - z2_mid)) * 
          (gamma_f + 2.0 * (theta2 - z2_mid)) - 8.0 * m / s) / 4.0 / gamma_f;

        if(2.0 * y2 + 2.0 * gamma_f * p2 > z2_mid){
          z2_min = z2_mid;
        }else{
          z2_max = z2_mid;
        }

        if(j == 50){
          double variance2 = 2.0 * gamma_f * gamma_f * p2 * (1.0 - p2) + 
            mu * calculate_variance2(theta2, sigma, sd_cut, s, m, n1, n2, z2_mid);

          ofs << "\t" << z2_mid << "\t" << p2 << "\t" << y2 << "\t" << variance2 << std::endl;
        }
      }
    }
  }

  return(0);
}

void calculate_beta_vec(const double theta2, const double sigma, 
  const double sd_cut, const double s, const double m,
  const double n1, const double n2, 
  std::vector<double>& vec_z2, std::vector<double>& vec_beta){

  std::ofstream ofs("regi_beta.txt");

  vec_z2.clear();
  vec_beta.clear();

  int z_sep = 1000;
  int int_sep = 500;
  double trunk_correction = -1.0 + (1.0 + std::erf(sd_cut / std::sqrt(2.0)));

  double h_int = sd_cut * sigma / int_sep;

  for(int i = 0; i <= z_sep; i++){
    double z2 = 1.0 * theta2 * i / z_sep;
    double sum_mean = 0.0;

    for(int j = -int_sep; j < int_sep; j++){
      double gamma0 = sd_cut * sigma * j / int_sep;
      double gamma1 = sd_cut * sigma * (j + 0.5) / int_sep;
      double gamma2 = sd_cut * sigma * (j + 1.0) / int_sep;

      double mean0 = integrate_gamma_x_t2(s, gamma0, z2, theta2, m, n2);
      double mean1 = integrate_gamma_x_t2(s, gamma1, z2, theta2, m, n2);
      double mean2 = integrate_gamma_x_t2(s, gamma2, z2, theta2, m, n2);

      double f0 = 1.0 / std::sqrt(2.0 * std::acos(-1.0) * sigma * sigma) *
        std::exp(-gamma0 * gamma0 / 2.0 / sigma / sigma) / trunk_correction;
      double f1 = 1.0 / std::sqrt(2.0 * std::acos(-1.0) * sigma * sigma) *
        std::exp(-gamma1 * gamma1 / 2.0 / sigma / sigma) / trunk_correction;
      double f2 = 1.0 / std::sqrt(2.0 * std::acos(-1.0) * sigma * sigma) *
        std::exp(-gamma2 * gamma2 / 2.0 / sigma / sigma) / trunk_correction;

      double mval0 = 2.0 * n1 * 2.0 * n2 * m * integrate_x_t1(s, n1, gamma0) + 2.0 * n2;
      double mval1 = 2.0 * n1 * 2.0 * n2 * m * integrate_x_t1(s, n1, gamma1) + 2.0 * n2;
      double mval2 = 2.0 * n1 * 2.0 * n2 * m * integrate_x_t1(s, n1, gamma2) + 2.0 * n2;

      sum_mean += h_int / 6.0 * (mean0 * f0 * mval0 +
        4.0 * mean1 * f1 * mval1 + mean2 * f2 * mval2);
    }

    vec_z2.push_back(z2);
    vec_beta.push_back(sum_mean);

    ofs << z2 << "\t" << sum_mean << std::endl;
  }
}

double integrate_gamma_x_t2(const double s, const double gamma, const double z2,
  const double theta2, const double m, const double n2){
  // calculate int_0^1 \gamma x T_2 dx
  // use tanh-sinh quadrature for integration

  double ret = 0.0;

  double h = 0.01;
  double pi = std::acos(-1.0);

  for(int k = -1000; k <= 1000; k++){
    double y = std::tanh(pi * std::sinh(k * h) / 2.0);
    double x = (y + 1.0) / 2.0;
    // g_x = -int_0^x 2M2/V2 when m = 0
    double g_x = -8.0 * n2 * s * gamma * (theta2 - z2) * x +
      4.0 * n2 * s * gamma * gamma * x * (1.0 - x);

    double log_sum = (1.0 - 4.0 * n2 * m) * std::log(2.0);

    log_sum += (1.0 - 4.0 * n2 * m) * pi / 2.0 * std::sinh(k * h);
    log_sum += std::log(pi / 2.0 * std::cosh(k * h));
    log_sum -= g_x;

    double log_value2 = (1.0 + 4.0 * n2 * m) * std::log(std::cosh(pi * std::sinh(k * h) / 2.0));

    if(std::sinh(k * h) < -100.0){
      log_value2 = (1.0 + 4.0 * n2 * m) * (-pi * std::sinh(k * h) / 2.0 - std::log(2.0));
    }
    if(std::sinh(k * h) > 100.0){
      log_value2 = (1.0 + 4.0 * n2 * m) * (pi * std::sinh(k * h) / 2.0 - std::log(2.0));
    }

    log_sum -= log_value2;
    ret += h * std::exp(log_sum);
  }

  return(ret * gamma);
}

double integrate_gamma_two_x_1minusx_t1(const double s, const double gamma, const double n1){

  int sep = 1000;

  // calculate int_0^1 2\gamma^2 x(1-x) T_1 dx
  double ret = 0.0;
  double a = 4.0 * n1 * s * gamma * gamma;

  for(int i = 0; i < sep; i++){
    double x0 = 1.0 * i / sep;
    double x1 = 1.0 * (i + 0.5) / sep;
    double x2 = 1.0 * (i + 1.0) / sep;

    double vx0 = (int_g1(a, 1.0) - int_g1(a, x0)) / (int_g1(a, 1.0) - int_g1(a, 0.0));
    double vx1 = (int_g1(a, 1.0) - int_g1(a, x1)) / (int_g1(a, 1.0) - int_g1(a, 0.0));
    double vx2 = (int_g1(a, 1.0) - int_g1(a, x2)) / (int_g1(a, 1.0) - int_g1(a, 0.0));

    double y0 = 4.0 * gamma * gamma * vx0 / std::exp(4.0 * n1 * s * gamma * gamma * x0 * (1.0 - x0));
    double y1 = 4.0 * gamma * gamma * vx1 / std::exp(4.0 * n1 * s * gamma * gamma * x1 * (1.0 - x1));
    double y2 = 4.0 * gamma * gamma * vx2 / std::exp(4.0 * n1 * s * gamma * gamma * x2 * (1.0 - x2));

    ret += 1.0 / sep / 6.0 * (y0 + 4.0 * y1 + y2);
  }

  return(ret);
}

double integrate_gamma_two_x_1minusx_t2(const double s, const double gamma, const double z2,
  const double theta2, const double m, const double n2){

  int sep = 1000;

  // calculate int_0^1 2\gamma^2 x(1-x) T_2 dx
  double ret = 0.0;

  for(int i = 0; i < sep; i++){
    double x0 = 1.0 * i / sep;
    double x1 = 1.0 * (i + 0.5) / sep;
    double x2 = 1.0 * (i + 1.0) / sep;

    double y0 = 4.0 * gamma * gamma / std::exp(-8.0 * n2 * s * gamma * (theta2 - z2) * x0 +
      4.0 * n2 * s * gamma * gamma * x0 * (1.0 - x0)) * std::pow(1.0 - x0, 4.0 * n2 * m);
    double y1 = 4.0 * gamma * gamma / std::exp(-8.0 * n2 * s * gamma * (theta2 - z2) * x1 +
      4.0 * n2 * s * gamma * gamma * x1 * (1.0 - x1)) * std::pow(1.0 - x1, 4.0 * n2 * m);
    double y2 = 4.0 * gamma * gamma / std::exp(-8.0 * n2 * s * gamma * (theta2 - z2) * x2 +
      4.0 * n2 * s * gamma * gamma * x2 * (1.0 - x2)) * std::pow(1.0 - x2, 4.0 * n2 * m);

    ret += 1.0 / sep / 6.0 * (y0 + 4.0 * y1 + y2);
  }

  return(ret);
}

double f_g1(const double a, const double x){
  double ret = std::exp(a * x * (1.0 - x));
  return(ret);
}

double int_g1(const double a, const double x){
  if(a < 1e-5){
    double ret = 0.5 * (2.0 * x - 1.0) -
      (2.0 * x - 1.0) * (2.0 * x * x - 2.0 * x - 1.0) / 12.0 * a;
    return(ret);
  }else{
    double ret = std::exp(a / 4.0) * std::sqrt(std::acos(-1.0)) *
      std::erf(0.5 * std::sqrt(a) * (2.0 * x - 1.0)) / 2.0 / std::sqrt(a);
    return(ret);
  }
}

double x_t1(const double a, const double x){
  if(1.0 - x < 1e-6){
    double ret = 2.0 * (1.0 - x) / (int_g1(a, 1.0) - int_g1(a, 0.0)) /
      (int_g1(a, 1.0) - int_g1(a, 0.0)) / f_g1(a, x);
    return(ret);
  }else{
    double vx = (int_g1(a, 1.0) - int_g1(a, x)) / (int_g1(a, 1.0) - int_g1(a, 0.0));
    double ret = 2.0 * vx * vx / (1.0 - x) / f_g1(a, x);
    return(ret);
  }
}

double integrate_x_t1(const double s, const double n1,
  const double gamma){

  double sum = 0.0;
  int sep = 1000;
  double h = 1.0 / sep;

  for(int i = 0; i < sep; i++){
    double x0 = 1.0 * i / sep;
    double x1 = 1.0 * (i + 0.5) / sep;
    double x2 = 1.0 * (i + 1.0) / sep;

    double y0 = x_t1(4.0 * n1 * s * gamma * gamma, x0);
    double y1 = x_t1(4.0 * n1 * s * gamma * gamma, x1);
    double y2 = x_t1(4.0 * n1 * s * gamma * gamma, x2);

    sum += h / 6.0 * (y0 + 4.0 * y1 + y2);
  }

  return(sum);
}

double interpolate_beta(const double z2, const std::vector<double>& vec_z2, 
  const std::vector<double>& vec_beta){
  
  int index = 0;
  while(z2 > vec_z2.at(index)){
    index++;
  }

  double z2_upper = vec_z2.at(index);
  double z2_lower = vec_z2.at(index - 1);

  double beta_upper = vec_beta.at(index);
  double beta_lower = vec_beta.at(index - 1);

  // linear interpolation
  double beta = beta_lower + (beta_upper - beta_lower) * (z2 - z2_lower) / 
    (z2_upper - z2_lower);

  return(beta);
}

double calculate_variance1(const double sigma, 
  const double sd_cut, const double s, const double n1){

  int int_sep = 1000;
  double trunk_correction = -1.0 + (1.0 + std::erf(sd_cut / std::sqrt(2.0)));

  double h_int = sd_cut * sigma / int_sep;
  double ret = 0.0;

  for(int i = -int_sep; i < int_sep; i++){
    double gamma0 = sd_cut * sigma * i / int_sep;
    double gamma1 = sd_cut * sigma * (i + 0.5) / int_sep;
    double gamma2 = sd_cut * sigma * (i + 1.0) / int_sep;

    double var0 = integrate_gamma_two_x_1minusx_t1(s, gamma0, n1);
    double var1 = integrate_gamma_two_x_1minusx_t1(s, gamma1, n1);
    double var2 = integrate_gamma_two_x_1minusx_t1(s, gamma2, n1);

    double f0 = 1.0 / std::sqrt(2.0 * std::acos(-1.0) * sigma * sigma) *
      std::exp(-gamma0 * gamma0 / 2.0 / sigma / sigma) / trunk_correction;
    double f1 = 1.0 / std::sqrt(2.0 * std::acos(-1.0) * sigma * sigma) *
      std::exp(-gamma1 * gamma1 / 2.0 / sigma / sigma) / trunk_correction;
    double f2 = 1.0 / std::sqrt(2.0 * std::acos(-1.0) * sigma * sigma) *
      std::exp(-gamma2 * gamma2 / 2.0 / sigma / sigma) / trunk_correction;

    ret += h_int / 6.0 * (var0 * f0 + 4.0 * var1 * f1 + var2 * f2);
  }

  return(2.0 * n1 * ret);
}

double calculate_variance2(const double theta2, const double sigma, 
  const double sd_cut, const double s, const double m,
  const double n1, const double n2, double z2){

  int int_sep = 1000;
  double trunk_correction = -1.0 + (1.0 + std::erf(sd_cut / std::sqrt(2.0)));

  double h_int = sd_cut * sigma / int_sep;
  double ret = 0.0;

  for(int i = -int_sep; i < int_sep; i++){
    double gamma0 = sd_cut * sigma * i / int_sep;
    double gamma1 = sd_cut * sigma * (i + 0.5) / int_sep;
    double gamma2 = sd_cut * sigma * (i + 1.0) / int_sep;

    double var0 = integrate_gamma_two_x_1minusx_t2(s, gamma0, z2, theta2, m, n2);
    double var1 = integrate_gamma_two_x_1minusx_t2(s, gamma1, z2, theta2, m, n2);
    double var2 = integrate_gamma_two_x_1minusx_t2(s, gamma2, z2, theta2, m, n2);

    double f0 = 1.0 / std::sqrt(2.0 * std::acos(-1.0) * sigma * sigma) *
      std::exp(-gamma0 * gamma0 / 2.0 / sigma / sigma) / trunk_correction;
    double f1 = 1.0 / std::sqrt(2.0 * std::acos(-1.0) * sigma * sigma) *
      std::exp(-gamma1 * gamma1 / 2.0 / sigma / sigma) / trunk_correction;
    double f2 = 1.0 / std::sqrt(2.0 * std::acos(-1.0) * sigma * sigma) *
      std::exp(-gamma2 * gamma2 / 2.0 / sigma / sigma) / trunk_correction;

    double mval0 = 2.0 * n1 * 2.0 * n2 * m * integrate_x_t1(s, n1, gamma0) + 2.0 * n2;
    double mval1 = 2.0 * n1 * 2.0 * n2 * m * integrate_x_t1(s, n1, gamma1) + 2.0 * n2;
    double mval2 = 2.0 * n1 * 2.0 * n2 * m * integrate_x_t1(s, n1, gamma2) + 2.0 * n2;

    ret += h_int / 6.0 * (var0 * f0 * mval0 +
      4.0 * var1 * f1 * mval1 + var2 * f2 * mval2);
  }

  return(ret);
}
