#include <cmath>
#include <vector>
#include <fstream>
#include <iostream>

double integrate_gamma_x_t2(const double s, const double gamma, const double z,
  const double theta2, const double m, const double n2);
double calculate_beta(const double z2, const double theta2, const double sigma, 
  const double sd_cut, const double s, const double m, 
  const double n1, const double n2);
double integrate_x_t1(const double s, const double n1,
  const double gamma);
double f_g1(const double a, const double x);
double int_g1(const double a, const double x);
double x_t1(const double a, const double x);

int main(){
  double n1 = 10000.0;
  double n2 = 10000.0;
  double s = 0.025;
  double gamma_f = 0.2;
  double theta2 = 1.0;
  double sigma = 0.02;

  double sd_cut = 3.0;

  if(theta2 < 0.0 || gamma_f < 0.0 || 2.0 * gamma_f > theta2){
    std::cerr << "Parameter error!" << std::endl;
    std::exit(0);
  }

  std::ofstream ofs("boundary.txt");
  ofs << "m\tu" << std::endl;

  int m_sep = 100;
  double m_min = s * gamma_f * gamma_f / 8.0;
  double m_max = 2.0 * s * gamma_f * (theta2 - gamma_f / 2.0);

  std::cout << "min_m: " << m_min << "\t" << "max_m: " << m_max << std::endl;
  ofs << m_min << "\t" << "Inf" << std::endl;

  double log_m_min = std::log(m_min);
  double log_m_max = std::log(m_max);

  for(int i = 1; i < m_sep; i++){
    double log_m = log_m_min + 1.0 * i / m_sep * (log_m_max - log_m_min);
    double m = std::exp(log_m);

    double z2, p2;
    if(gamma_f * gamma_f > m / 2.0 / s){
      z2 = theta2 + gamma_f / 2.0 - std::sqrt(2.0 * m / s);
      p2 = 1.0 - std::sqrt(m / 2.0 / s / gamma_f / gamma_f);
    }else{
      z2 = theta2 - gamma_f / 2.0 - m / 2.0 / s / gamma_f;
      p2 = 0.0;
    }

    double y2 = (z2 - 2.0 * gamma_f * p2) / 2.0;
    double beta = calculate_beta(z2, theta2, sigma, sd_cut, s, m, n1, n2);
    double thre_u = y2 / beta;

    ofs << m << "\t" << thre_u << std::endl;
  }

  ofs << m_max << "\t0.0" << std::endl;

  return(0);
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

double calculate_beta(const double z2, const double theta2, const double sigma, const double sd_cut,
  const double s, const double m, const double n1, const double n2){

  int int_sep = 500;
  double trunk_correction = -1.0 + (1.0 + std::erf(sd_cut / std::sqrt(2.0)));

  double h_int = sd_cut * sigma / int_sep;
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
  return(sum_mean);
}
