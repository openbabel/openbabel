// CppNumericalSolver
// based on:
// L-BFGS-B: A LIMITED MEMORY ALGORITHM FOR BOUND CONSTRAINED OPTIMIZATION
// Richard H. Byrd, Peihuang Lu, Jorge Nocedal and Ciyou Zhu
#include <iostream>
#include <list>
#include <Eigen/LU>
#include "isolver.h"
#include "../boundedproblem.h"
#include "../linesearch/morethuente.h"
#ifndef LBFGSBSOLVER_H
#define LBFGSBSOLVER_H
namespace cppoptlib {
template<typename TProblem>
class LbfgsbSolver : public ISolver<TProblem, 1> {
  public:
    using Superclass = ISolver<TProblem, 1>;
    using typename Superclass::Scalar;
    using typename Superclass::TVector;
    using MatrixType = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;
    using VariableTVector = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
  protected:
  // workspace matrices
  MatrixType W, M;
  Scalar theta;
  int DIM;
  int m_historySize = 5;

  /**
   * @brief sort pairs (k,v) according v ascending
   * @details [long description]
   *
   * @param v [description]
   * @return [description]
   */
  std::vector<int> sort_indexes(const std::vector< std::pair<int, Scalar> > &v) {
    std::vector<int> idx(v.size());
    for (size_t i = 0; i != idx.size(); ++i)
      idx[i] = v[i].first;
    sort(idx.begin(), idx.end(), [&v](size_t i1, size_t i2) {
      return v[i1].second < v[i2].second;
    });
    return idx;
  }
  /**
   * @brief Algorithm CP: Computation of the generalized Cauchy point
   * @details PAGE 8
   *
   * @param c [description]
   */
  void getGeneralizedCauchyPoint(const TProblem &problem, const TVector &x, const TVector &g, TVector &x_cauchy, VariableTVector &c) {
    const int DIM = x.rows();
    // Given x,l,u,g, and B = \theta I-WMW
    // {all t_i} = { (idx,value), ... }
    // TODO: use "std::set" ?
    std::vector<std::pair<int, Scalar> > SetOfT;
    // the feasible set is implicitly given by "SetOfT - {t_i==0}"
    TVector d = -g;
    // n operations
    for (int j = 0; j < DIM; j++) {
      if (g(j) == 0) {
        SetOfT.push_back(std::make_pair(j, std::numeric_limits<Scalar>::max()));
      } else {
        Scalar tmp = 0;
        if (g(j) < 0) {
          tmp = (x(j) - problem.upperBound()(j)) / g(j);
        } else {
          tmp = (x(j) - problem.lowerBound()(j)) / g(j);
        }
        SetOfT.push_back(std::make_pair(j, tmp));
      if (tmp == 0) d(j) = 0;
      }
    }
    // sortedindices [1,0,2] means the minimal element is on the 1-st entry
    std::vector<int> sortedIndices = sort_indexes(SetOfT);
    x_cauchy = x;
    // Initialize
    // p :=     W^Scalar*p
    VariableTVector p = (W.transpose() * d);                     // (2mn operations)
    // c :=     0
    c = VariableTVector::Zero(W.cols());
    // f' :=    g^Scalar*d = -d^Td
    Scalar f_prime = -d.dot(d);                         // (n operations)
    // f'' :=   \theta*d^Scalar*d-d^Scalar*W*M*W^Scalar*d = -\theta*f' - p^Scalar*M*p
    Scalar f_doubleprime = (Scalar)(-1.0 * theta) * f_prime - p.dot(M * p); // (O(m^2) operations)
    f_doubleprime = std::max<Scalar>(std::numeric_limits<Scalar>::epsilon(), f_doubleprime);
    Scalar f_dp_orig = f_doubleprime;
    // \delta t_min :=  -f'/f''
    Scalar dt_min = -f_prime / f_doubleprime;
    // t_old :=     0
    Scalar t_old = 0;
    // b :=     argmin {t_i , t_i >0}
    int i = 0;
    for (int j = 0; j < DIM; j++) {
      i = j;
      if (SetOfT[sortedIndices[j]].second > 0)
        break;
    }
    int b = sortedIndices[i];
    // see below
    // t                    :=  min{t_i : i in F}
    Scalar t = SetOfT[b].second;
    // \delta Scalar             :=  t - 0
    Scalar dt = t ;
    // examination of subsequent segments
    while ((dt_min >= dt) && (i < DIM)) {
      if (d(b) > 0)
        x_cauchy(b) = problem.upperBound()(b);
      else if (d(b) < 0)
        x_cauchy(b) = problem.lowerBound()(b);
      // z_b = x_p^{cp} - x_b
      Scalar zb = x_cauchy(b) - x(b);
      // c   :=  c +\delta t*p
      c += dt * p;
      // cache
      VariableTVector wbt = W.row(b);
      f_prime += dt * f_doubleprime + (Scalar) g(b) * g(b) + (Scalar) theta * g(b) * zb - (Scalar) g(b) *
      wbt.transpose() * (M * c);
      f_doubleprime += (Scalar) - 1.0 * theta * g(b) * g(b)
                       - (Scalar) 2.0 * (g(b) * (wbt.dot(M * p)))
                       - (Scalar) g(b) * g(b) * wbt.transpose() * (M * wbt);
      f_doubleprime = std::max<Scalar>(std::numeric_limits<Scalar>::epsilon() * f_dp_orig, f_doubleprime);
      p += g(b) * wbt.transpose();
      d(b) = 0;
      dt_min = -f_prime / f_doubleprime;
      t_old = t;
      ++i;
      if (i < DIM) {
        b = sortedIndices[i];
        t = SetOfT[b].second;
        dt = t - t_old;
      }
    }
    dt_min = std::max<Scalar>(dt_min, (Scalar)0.0);
    t_old += dt_min;
    #pragma omp parallel for
    for (int ii = i; ii < x_cauchy.rows(); ii++) {
      x_cauchy(sortedIndices[ii]) = x(sortedIndices[ii]) + t_old * d(sortedIndices[ii]);
    }
    c += dt_min * p;
  }
  /**
   * @brief find alpha* = max {a : a <= 1 and  l_i-xc_i <= a*d_i <= u_i-xc_i}
   * @details [long description]
   *
   * @param FreeVariables [description]
   * @return [description]
   */
  Scalar findAlpha(const TProblem &problem, TVector &x_cp, VariableTVector &du, std::vector<int> &FreeVariables) {
    Scalar alphastar = 1;
    const unsigned int n = FreeVariables.size();
    assert(du.rows() == n);
    for (unsigned int i = 0; i < n; i++) {
      if (du(i) > 0) {
        alphastar = std::min<Scalar>(alphastar, (problem.upperBound()(FreeVariables[i]) - x_cp(FreeVariables[i])) / du(i));
      } else {
        alphastar = std::min<Scalar>(alphastar, (problem.lowerBound()(FreeVariables[i]) - x_cp(FreeVariables[i])) / du(i));
      }
    }
    return alphastar;
  }
  /**
   * @brief solving unbounded probelm
   * @details [long description]
   *
   * @param SubspaceMin [description]
   */
  void SubspaceMinimization(const TProblem &problem, TVector &x_cauchy, TVector &x, VariableTVector &c, TVector &g,
  TVector &SubspaceMin) {
    Scalar theta_inverse = 1 / theta;
    std::vector<int> FreeVariablesIndex;
    for (int i = 0; i < x_cauchy.rows(); i++) {
      if ((x_cauchy(i) != problem.upperBound()(i)) && (x_cauchy(i) != problem.lowerBound()(i))) {
        FreeVariablesIndex.push_back(i);
      }
    }
    const int FreeVarCount = FreeVariablesIndex.size();
    MatrixType WZ = MatrixType::Zero(W.cols(), FreeVarCount);
    for (int i = 0; i < FreeVarCount; i++)
      WZ.col(i) = W.row(FreeVariablesIndex[i]);
    TVector rr = (g + theta * (x_cauchy - x) - W * (M * c));
    // r=r(FreeVariables);
    MatrixType r = MatrixType::Zero(FreeVarCount, 1);
    for (int i = 0; i < FreeVarCount; i++)
      r.row(i) = rr.row(FreeVariablesIndex[i]);
    // STEP 2: "v = w^T*Z*r" and STEP 3: "v = M*v"
    VariableTVector v = M * (WZ * r);
    // STEP 4: N = 1/theta*W^T*Z*(W^T*Z)^T
    MatrixType N = theta_inverse * WZ * WZ.transpose();
    // N = I - MN
    N = MatrixType::Identity(N.rows(), N.rows()) - M * N;
    // STEP: 5
    // v = N^{-1}*v
    if (v.size() > 0)
      v = N.lu().solve(v);
    // STEP: 6
    // HERE IS A MISTAKE IN THE ORIGINAL PAPER!
    VariableTVector du = -theta_inverse * r - theta_inverse * theta_inverse * WZ.transpose() * v;
    // STEP: 7
    Scalar alpha_star = findAlpha(problem, x_cauchy, du, FreeVariablesIndex);
    // STEP: 8
    VariableTVector dStar = alpha_star * du;
    SubspaceMin = x_cauchy;
    for (int i = 0; i < FreeVarCount; i++) {
      SubspaceMin(FreeVariablesIndex[i]) = SubspaceMin(FreeVariablesIndex[i]) + dStar(i);
    }
  }
 public:
  void setHistorySize(const int hs) { m_historySize = hs; }

  void minimize(TProblem &problem, TVector &x0) {
    if(!problem.isValid(x0))
      std::cerr << "start with invalid x0" << std::endl;
    DIM = x0.rows();
    theta = 1.0;
    W = MatrixType::Zero(DIM, 0);
    M = MatrixType::Zero(0, 0);
    MatrixType yHistory = MatrixType::Zero(DIM, 0);
    MatrixType sHistory = MatrixType::Zero(DIM, 0);
    TVector x = x0, g = x0;
    Scalar f = problem.value(x);
    problem.gradient(x, g);
    // conv. crit.
    auto noConvergence =
    [&](TVector &x, TVector &g)->bool {
      return (((x - g).cwiseMax(problem.lowerBound()).cwiseMin(problem.upperBound()) - x).template lpNorm<Eigen::Infinity>() >= 1e-4);
    };
    this->m_current.reset();
    this->m_status = Status::Continue;
    while (problem.callback(this->m_current, x) && noConvergence(x, g) && (this->m_status == Status::Continue)) {
      Scalar f_old = f;
      TVector x_old = x;
    TVector g_old = g;
      // STEP 2: compute the cauchy point
      TVector CauchyPoint = TVector::Zero(DIM);
      VariableTVector c = VariableTVector::Zero(W.cols());
    getGeneralizedCauchyPoint(problem, x, g, CauchyPoint, c);
      // STEP 3: compute a search direction d_k by the primal method for the sub-problem
      TVector SubspaceMin;
    SubspaceMinimization(problem, CauchyPoint, x, c, g, SubspaceMin);
      // STEP 4: perform linesearch and STEP 5: compute gradient
      Scalar alpha_init = 1.0;
      const Scalar rate = MoreThuente<TProblem, 1>::linesearch(x,  SubspaceMin-x ,  problem, alpha_init);
      // update current guess and function information
      x = x - rate*(x-SubspaceMin);
      f = problem.value(x);
      problem.gradient(x, g);
      // prepare for next iteration
      TVector newY = g - g_old;
      TVector newS = x - x_old;
      // STEP 6:
      Scalar test = newS.dot(newY);
      test = (test < 0) ? -1.0 * test : test;
      if (test > 1e-7 * newY.squaredNorm()) {
        if (yHistory.cols() < m_historySize) {
          yHistory.conservativeResize(DIM, yHistory.cols() + 1);
          sHistory.conservativeResize(DIM, sHistory.cols() + 1);
        } else {
          yHistory.leftCols(m_historySize - 1) = yHistory.rightCols(m_historySize - 1).eval();
          sHistory.leftCols(m_historySize - 1) = sHistory.rightCols(m_historySize - 1).eval();
        }
        yHistory.rightCols(1) = newY;
        sHistory.rightCols(1) = newS;
    // STEP 7:
        theta = (Scalar)(newY.transpose() * newY) / (newY.transpose() * newS);
        W = MatrixType::Zero(yHistory.rows(), yHistory.cols() + sHistory.cols());
        W << yHistory, (theta * sHistory);
        MatrixType A = sHistory.transpose() * yHistory;
        MatrixType L = A.template triangularView<Eigen::StrictlyLower>();
        MatrixType MM(A.rows() + L.rows(), A.rows() + L.cols());
        MatrixType D = -1 * A.diagonal().asDiagonal();
        MM << D, L.transpose(), L, ((sHistory.transpose() * sHistory) * theta);
        M = MM.inverse();
      }
      if (fabs(f_old - f) < 1e-8) {
        // successive function values too similar
        break;
      }
      ++this->m_current.iterations;
      this->m_current.gradNorm = g.norm();
      this->m_status = checkConvergence(this->m_stop, this->m_current);
    }
    x0 = x;
    if (this->m_debug > DebugLevel::None) {
        std::cout << "Stop status was: " << this->m_status << std::endl;
        std::cout << "Stop criteria were: " << std::endl << this->m_stop << std::endl;
        std::cout << "Current values are: " << std::endl << this->m_current << std::endl;
    }
  }
};
}
/* namespace cppoptlib */
#endif /* LBFGSBSOLVER_H_ */
