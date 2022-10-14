#include <RcppArmadillo.h>
#include <RcppEigen.h>
//#include <Rcpp.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo, RcppEigen)]]

//typedef arma::mat MATTYPE;
//typedef arma::vec VECTYPE;
//typedef arma::fmat MATTYPE;
//typedef arma::fvec VECTYPE;


// [[Rcpp::export]]
List norm_S_batches(const arma::vec& x, const arma::vec& p, const arma::vec& i, int ncol, const arma::vec& groups, int N) {
    arma::mat exp_sum = arma::zeros<arma::mat>(ncol, N);
    std::vector<int> jj(i.size());
    std::vector<double> xx(i.size());
    int ii = 0;
    
    for (int c = 0; c < ncol; c++) {
        for (int j = p[c]; j < p[c + 1]; j++) {
            exp_sum(c, groups[i[j]]) += std::exp(x[j]);
            jj[ii] = c + 1;
            ii++;
        }
    }
    for (int c = 0; c < ncol; c++) {
        for (int j = p[c]; j < p[c + 1]; j++) {
            xx[j] = std::exp(x[j])/exp_sum(c, groups[i[j]]);
        }
    }
    return Rcpp::List::create(Named("i") = i + 1, Named("j") = jj, Named("x") = xx);//Rcpp::IntegerVector(i(), i.end()),
}

// [[Rcpp::export]]
List upper_tri(const arma::vec& x, const arma::vec& p, const arma::vec& i, int ncol) {
    int l = i.size() / 2;
    std::vector<int> ii(l);
    std::vector<int> jj(l);
    std::vector<double> xx(l);
    int tt = 0;

    for (int c = 0; c < ncol; c++) {
        for (int j = p[c]; j < p[c + 1]; j++) {
            if (i[j] < c) {
                ii[tt] = i[j] + 1;
                jj[tt] = c + 1;
                xx[tt] = x[j];
                tt++;
            }
        }
    }

    return Rcpp::List::create(Named("i") = ii, Named("j") = jj, Named("x") = xx);
}

// [[Rcpp::export]]
List filter_T(const arma::vec& x, const arma::vec& p, const arma::vec& i, int ncol, double T_th, const arma::mat& A, const arma::mat& B) {
    int l = i.size();
    std::vector<int> ii(l);
    std::vector<int> jj(l);
    std::vector<double> xx(l);
    int tt = 0;

    for (int c = 0; c < ncol; c++) {
        for (int j = p[c]; j < p[c + 1]; j++) {
            arma::vec a_j = A.col(i[j]);
            arma::vec b_c = B.col(c);
            arma::mat ip = a_j.t() * b_c;
            if (ip(0) < T_th) {
                xx[tt] = 0;
            } else {
                xx[tt] = x[j];
            }
            ii[tt] = i[j] + 1;
            jj[tt] = c + 1;
            tt++;
        }
    }

    return Rcpp::List::create(Named("i") = ii, Named("j") = jj, Named("x") = xx);
}

// [[Rcpp::export]]
SEXP eigenMapMatMult(const Eigen::Map<Eigen::MatrixXd> A, Eigen::Map<Eigen::MatrixXd> B){
    Eigen::MatrixXd C = A * B;

    return Rcpp::wrap(C);
}

// [[Rcpp::export]]
SEXP eigenMapMatcrossprod(const Eigen::Map<Eigen::MatrixXd> A, Eigen::Map<Eigen::MatrixXd> B){
  Eigen::MatrixXd C = A.adjoint() * B;
  
  return Rcpp::wrap(C);
}
// [[Rcpp::export]]
SEXP eigenMapMattcrossprod(const Eigen::Map<Eigen::MatrixXd> A, Eigen::Map<Eigen::MatrixXd> B){
  Eigen::MatrixXd C = A * B.adjoint();
  
  return Rcpp::wrap(C);
}
