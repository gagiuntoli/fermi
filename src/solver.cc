
#include "solver.h"

#include <cmath>

double solver_keff(const ellpack_t &A, const ellpack_t &B, std::vector<size_t> dirichlet) {
  const size_t MAX_ITERS = 100;
  size_t n = A.nrows;
  std::vector<double> phi(n, 1.0);
  std::vector<double> source(n);
  std::vector<double> source_new(n);

  double keff = 1.0;
  size_t iters = 0;
  while (true) {
    for (const size_t &node : dirichlet) {
      phi[node] = 0.0;
    }

    ellpack_mvp(source, B, phi);
    double norm_source = 0;
    for (size_t i = 0; i < n; i++) {
      norm_source += source[i];
      source[i] /= keff;
    }

    ellpack_solve_cg(phi, A, source);

    ellpack_mvp(source_new, B, phi);

    double norm_source_new = 0;
    for (size_t i = 0; i < n; i++) {
      norm_source_new += source_new[i];
    }
    double keff_new = keff * norm_source_new / norm_source;

    if (iters > MAX_ITERS || (std::abs(keff_new - keff) / keff) < 0.01) return keff_new;
    keff = keff_new;

    double power = norm(phi, n);
    for (size_t i = 0; i < n; i++) {
      phi[i] /= power;
    }

    iters++;
  }
}
