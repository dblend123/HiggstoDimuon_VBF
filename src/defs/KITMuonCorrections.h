#pragma once

#include <boost/math/special_functions/erf.hpp>
#include <correction.h>

#include <TRandom3.h>

#include <algorithm>
#include <cstdint>
#include <cmath>
#include <initializer_list>
#include <iostream>
#include <string>
#include <vector>
#include "Muon.h"

namespace kit {

// ------------------------------
// SeedSequence (KIT reference)
// ------------------------------
class SeedSequence {
public:
  explicit SeedSequence(std::initializer_list<uint32_t> seeds)
      : m_seeds(seeds) {}

  template <typename Iter>
  void generate(Iter begin, Iter end) const {
    const size_t n = std::distance(begin, end);
    if (n == 0) return;

    const uint32_t mult = 0x9e3779b9;
    const uint32_t mix_const = 0x85ebca6b;

    std::vector<uint32_t> buffer(n, 0x8b8b8b8b);

    size_t s = m_seeds.size();

    size_t i = 0;
    for (; i < std::min(n, s); ++i) {
      buffer[i] = buffer[i] ^ (m_seeds[i] + mult * (uint32_t)i);
    }
    for (; i < n; ++i) {
      buffer[i] = buffer[i] ^ (mult * (uint32_t)i);
    }

    for (size_t k = 0; k < n; ++k) {
      uint32_t z = buffer[(k + n - 1) % n] ^ (buffer[k] >> 27);
      buffer[k] = (z * mix_const) ^ (buffer[k] << 13);
    }

    std::copy(buffer.begin(), buffer.end(), begin);
  }

private:
  std::vector<uint32_t> m_seeds;
};

// ------------------------------
// CrystalBall (KIT reference)
// ------------------------------
struct CrystalBall {
  double pi = 3.14159;
  double sqrtPiOver2 = std::sqrt(pi / 2.0);
  double sqrt2 = std::sqrt(2.0);
  double m;
  double s;
  double a;
  double n;
  double B;
  double C;
  double D;
  double N;
  double NA;
  double Ns;
  double NC;
  double F;
  double G;
  double k;
  double cdfMa;
  double cdfPa;

  CrystalBall() : m(0), s(1), a(10), n(10) { init(); }
  CrystalBall(double mean, double sigma, double alpha, double nn)
      : m(mean), s(sigma), a(alpha), n(nn) { init(); }

  void init() {
    double fa = std::fabs(a);
    double ex = std::exp(-fa * fa / 2);
    double A  = std::pow(n / fa, n) * ex;
    double C1 = n / fa / (n - 1) * ex;
    double D1 = 2 * sqrtPiOver2 * std::erf(fa / sqrt2);
    B = n / fa - fa;
    C = (D1 + 2 * C1) / C1;
    D = (D1 + 2 * C1) / 2;
    N = 1.0 / s / (D1 + 2 * C1);
    k = 1.0 / (n - 1);
    NA = N * A;
    Ns = N * s;
    NC = Ns * C1;
    F = 1 - fa * fa / n;
    G = s * n / fa;
    cdfMa = cdf(m - a * s);
    cdfPa = cdf(m + a * s);
  }

  double cdf(double x) const {
    double d = (x - m) / s;
    if (d < -a) return NC / std::pow(F - s * d / G, n - 1);
    if (d >  a) return NC * (C - std::pow(F + s * d / G, 1 - n));
    return Ns * (D - sqrtPiOver2 * std::erf(-d / sqrt2));
  }

  double invcdf(double u) const {
    if (u < cdfMa) return m + G * (F - std::pow(NC / u, k));
    if (u > cdfPa) return m - G * (F - std::pow(C - u / NC, -k));
    return m - sqrt2 * s *
           boost::math::erf_inv((D - u / Ns) / sqrtPiOver2);
  }
};

// ------------------------------
// Helper evaluators (avoid global cset)
// ------------------------------
inline double get_rndm(const correction::CorrectionSet& cset,
                       double eta, double phi, float nL,
                       uint32_t evtNumber, uint32_t lumiNumber)
{
  const double mean  = cset.at("cb_params")->evaluate({std::abs(eta), (double)nL, 0});
  const double sigma = cset.at("cb_params")->evaluate({std::abs(eta), (double)nL, 1});
  const double nn    = cset.at("cb_params")->evaluate({std::abs(eta), (double)nL, 2});
  const double alpha = cset.at("cb_params")->evaluate({std::abs(eta), (double)nL, 3});

  CrystalBall cb(mean, sigma, alpha, nn);

  const int64_t phi_seed =
      (int64_t)((phi / M_PI) * ((1LL << 31) - 1)) & 0xFFF;

  SeedSequence seq{evtNumber, lumiNumber, (uint32_t)phi_seed};
  uint32_t seed = 0;
  seq.generate(&seed, &seed + 1);

  TRandom3 rnd(seed);
  const double u = rnd.Rndm();
  return cb.invcdf(u);
}

inline double get_std(const correction::CorrectionSet& cset,
                      double pt, double eta, float nL)
{
  const double p0 = cset.at("poly_params")->evaluate({std::abs(eta), (double)nL, 0});
  const double p1 = cset.at("poly_params")->evaluate({std::abs(eta), (double)nL, 1});
  const double p2 = cset.at("poly_params")->evaluate({std::abs(eta), (double)nL, 2});

  double sigma = p0 + p1 * pt + p2 * pt * pt;
  if (sigma < 0.0) sigma = 0.0;
  return sigma;
}

inline double get_k(const correction::CorrectionSet& cset,
                    double eta, const std::string& var)
{
  const double k_data = cset.at("k_data")->evaluate({std::abs(eta), var});
  const double k_mc   = cset.at("k_mc")->evaluate({std::abs(eta), var});

  double k = 0.0;
  if (k_mc < k_data) k = std::sqrt(k_data * k_data - k_mc * k_mc);
  return k;
}

// ------------------------------
// Public API you call from analyzer
// ------------------------------
inline double pt_resol_num(const correction::CorrectionSet& cset,
                           double pt, double eta, double phi, float nL,
                           uint32_t evtNumber, uint32_t lumiNumber,
                           double low_pt_threshold = 26.0)
{
  const double rndm = get_rndm(cset, eta, phi, nL, evtNumber, lumiNumber);
  const double std  = get_std (cset, pt, eta, nL);
  const double k    = get_k   (cset, eta, "nom");

  double ptc = pt * (1.0 + k * std * rndm);

  if (std::isnan(ptc)) ptc = pt;
  if (pt < low_pt_threshold || pt > 200.0 ||
      ptc / pt > 2.0 || ptc / pt < 0.1 || ptc < 0.0)
    ptc = pt;

  return ptc;
}

inline double basePt(const analysis::core::Muon& mu) {
  if (std::isfinite(mu._bsConstrainedChi2) &&
      mu._bsConstrainedChi2 < 30.0 &&
      std::isfinite(mu._bsConstrainedPt) &&
      mu._bsConstrainedPt > 0.0) {
    return mu._bsConstrainedPt;
  }
  return std::max(1e-3, (double)mu._pt);
}

inline float nLayers(const analysis::core::Muon& mu) {
  return (mu._nTrackerLayers > 0 ? (float)mu._nTrackerLayers : 12.0f);
}

inline double pt_resol(const correction::CorrectionSet& cset,
                       const analysis::core::Muon& mu,
                       uint32_t evtNumber, uint32_t lumiNumber,
                       double low_pt_threshold = 26.0)
{
  const double pt  = basePt(mu);
  const double eta = mu._eta;
  const double phi = mu._phi;
  const float  nL  = nLayers(mu);

  // call the existing numeric overload (the one you pasted earlier)
  return pt_resol_num(cset, pt, eta, phi, nL, evtNumber, lumiNumber, low_pt_threshold);
}

inline double pt_scale(const correction::CorrectionSet& cset,
                       bool is_data, double pt, double eta, double phi,
                       int charge, double low_pt_threshold = 26.0)
{
  const std::string dtmc = is_data ? "data" : "mc";

 const double abseta = std::abs(eta);

 const double a = cset.at("a_" + dtmc)->evaluate({abseta, phi, "nom"});
 const double m = cset.at("m_" + dtmc)->evaluate({abseta, phi, "nom"});

  if (pt < low_pt_threshold) return pt;
  return 1.0 / (m / pt + charge * a);
}

inline double pt_scale_var(const correction::CorrectionSet& cset,
                           double pt, double eta, double phi,
                           int charge, const std::string& updn)
{
  const double abseta = std::abs(eta);

  const double stat_a   = cset.at("a_mc")->evaluate({abseta, phi, std::string("stat")});
  const double stat_m   = cset.at("m_mc")->evaluate({abseta, phi, std::string("stat")});
  const double stat_rho = cset.at("m_mc")->evaluate({abseta, phi, std::string("rho_stat")});

  const double unc =
      pt * pt * std::sqrt(stat_m * stat_m / (pt * pt) +
                          stat_a * stat_a +
                          2.0 * charge * stat_rho * stat_m / pt * stat_a);

  double pt_var = pt;
  if (updn == "up") pt_var = pt + unc;
  if (updn == "dn") pt_var = pt - unc;
  return pt_var;
}

} // namespace kit

