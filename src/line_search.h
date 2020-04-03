#ifndef SNPLIB_LINE_SEARCH_H
#define SNPLIB_LINE_SEARCH_H

#include <cmath>

template <class T>
class LineSearch {
 private:
  T &worker_;
  double f_new;
  double f_old, g_old;
  double a_h, a_l;
  double f_h, f_l;
  double g_l;

  double Zoom() {
    double a_j;
    for (size_t i = 0; i < 10; ++i) {
      auto t = a_h - a_l;
      auto a = (f_h - f_l) - t * g_l;
      a /= t * t;
      a_j = a_l - g_l / a / 2.0;
      worker_.UpdateVars(a_j);
      double f_j = worker_.CalcLikelihood();
      if (a_j == a_l || i == 9) {
        worker_.UpdateGradients();
        f_new = f_j;
        return a_j;
      }
      if ((f_j - f_old) > 1e-4 * a_j * g_old || f_j >= f_l) {
        a_h = a_j;
        f_h = f_j;
      } else {
        worker_.UpdateGradients();
        double g_j = worker_.CalcLineGradient();
        if (std::fabs(g_j) <= -0.9 * g_old) {
          f_new = f_j;
          return a_j;
        }
        if (g_j * (a_h - a_l) >= 0) {
          a_h = a_l;
          f_h = f_l;
        }
        a_l = a_j;
        f_l = f_j;
        g_l = g_j;
      }
    }
    return a_j;
  }

 public:
  LineSearch(T &worker) : worker_(worker) {}
  ~LineSearch() {}
  double Search(double f) {
    g_old = worker_.CalcLineGradient();
    f_old = f;
    double a_1 = 1.0;
    double a_0 = 0.0;
    double f_0 = f_old;
    double g_0 = g_old;
    size_t num_iter = 0;
    while (true) {
      worker_.UpdateVars(a_1);
      double f_1 = worker_.CalcLikelihood();
      if ((f_1 - f_old) > 1e-4 * a_1 * g_old || (num_iter > 0 && f_1 >= f_0)) {
        a_h = a_1;
        f_h = f_1;
        a_l = a_0;
        f_l = f_0;
        g_l = g_0;
        return Zoom();
      }
      worker_.UpdateGradients();
      double g_1 = worker_.CalcLineGradient();
      if (std::fabs(g_1) <= -0.9 * g_old) {
        f_new = f_1;
        return a_1;
      }
      if (g_1 >= 0) {
        a_h = a_0;
        a_l = a_1;
        f_h = f_0;
        f_l = f_1;
        g_l = g_1;
        return Zoom();
      }
      a_0 = a_1;
      f_0 = f_1;
      g_0 = g_1;
      a_1 *= 1.618;
      num_iter++;
    }
  }
  double GetFNew() const { return f_new; }
};

#endif  // SNPLIB_LINE_SEARCH_H