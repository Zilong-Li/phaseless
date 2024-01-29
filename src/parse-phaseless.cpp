#include "phaseless.hpp"
#include <Rcpp.h>
#include <RcppEigen.h>
#include <alpaca/alpaca.h>

using namespace Rcpp;
using namespace std;

//' parse parameters in joint model
//' @param filename path to binary file from joint command
//' @export
// [[Rcpp::export]]
List parse_joint_par(std::string filename)
{
    auto filesize = std::filesystem::file_size(filename);
    std::error_code ec;
    std::ifstream ifs(filename, std::ios::in | std::ios::binary);
    constexpr auto OPTIONS = alpaca::options::fixed_length_encoding;
    std::unique_ptr<Pars> par = std::make_unique<Pars>(alpaca::deserialize<OPTIONS, Pars>(ifs, filesize, ec));
    ifs.close();
    assert((bool)ec == false);
    return List::create(Named("C") = par->C, Named("K") = par->K, Named("M") = par->M, Named("N") = par->N,
                        Named("er") = par->er, Named("P") = par->P, Named("F") = par->F, Named("Q") = par->Q);
}

//' parse posterior probabilty in joint model
//' @param filename path to binary file from joint command
//' @param chunk which chunk to extract
//' @export
// [[Rcpp::export]]
List parse_joint_post(std::string filename, int chunk = 0)
{
    auto filesize = std::filesystem::file_size(filename);
    std::error_code ec;
    std::ifstream ifs(filename, std::ios::in | std::ios::binary);
    constexpr auto OPTIONS = alpaca::options::fixed_length_encoding;
    std::unique_ptr<Pars> par = std::make_unique<Pars>(alpaca::deserialize<OPTIONS, Pars>(ifs, filesize, ec));
    ifs.close();
    assert((bool)ec == false);
    int nchunks = par->gls.size();
    const int C = par->C;
    const int K = par->K;
    const int CC = par->C * par->C;
    Eigen::Map<const MyArr2D> Q(par->Q.data(), par->K, par->N);
    Eigen::Map<const MyArr2D> P(par->P.data(), par->M, par->C);
    Eigen::Map<const MyArr1D> er(par->er.data(), par->M);
    auto R = er2R(er);
    MyFloat2D ret_gamma(par->N), ret_post_y(par->N), ret_post_zy(par->N);
    int m{0}, s{0}, z1{0}, z2{0}, y1{0}, zz{0}, S{0};
    for(int ic = 0, pos_chunk = 0; ic < nchunks; ic++)
    {
        S = par->pos[ic].size();
        pos_chunk += S;
        if(ic != chunk) continue;
        MyArr2D ind_post_zg1(C, S), ind_post_zg2(C, S);
        MyArr2D ind_post_zy(C * K, S);
        MyArr1D gamma_div_emit(CC);
        MyArr2D ind_post_y = MyArr2D::Zero(K, S);
        for(int ind = 0; ind < par->N; ind++)
        {
            Eigen::Map<const MyArr2D> gli(par->gls[ic].data() + ind * S * 3, S, 3);
            MyArr2D emit = get_emission_by_gl(gli, P.middleRows(pos_chunk - S, S)).transpose(); // CC x S
            // first get H ie old PI in fastphase
            MyArr2D H = MyArr2D::Zero(C, S);
            int z1, y1, s; // m * C + z1
            for(s = 0; s < S; s++)
                for(z1 = 0; z1 < C; z1++)
                    for(y1 = 0; y1 < par->K; y1++) H(z1, s) += Q(y1, ind) * par->F[y1][(pos_chunk - S + s) * C + z1];
            // cs is 1 / colsum(alpha)
            const auto [alpha, beta, cs] = forward_backwards_diploid(emit, R.middleCols(pos_chunk - S, S), H);
            MyArr2D gamma = alpha * beta;
            ret_gamma[ind] = MyFloat1D(gamma.data(), gamma.data() + gamma.size());
            ind_post_zg1.setZero(), ind_post_zg2.setZero(), ind_post_y.setZero(), ind_post_zy.setZero();
            for(s = 0; s < S; s++)
            {
                m = pos_chunk - S + s;
                gamma_div_emit = (alpha.col(s) * beta.col(s)) / emit.col(s); // what if emit is 0
                for(z1 = 0; z1 < C; z1++)
                {
                    ind_post_zg1(z1, s) = (gamma_div_emit(Eigen::seqN(z1, C, C)) * (1 - P(m, z1))
                                           * (gli(s, 0) * (1 - P.row(m)) + gli(s, 1) * P.row(m)).transpose())
                                              .sum();
                    ind_post_zg2(z1, s) = (gamma_div_emit(Eigen::seqN(z1, C, C)) * (P(m, z1))
                                           * (gli(s, 1) * (1 - P.row(m)) + gli(s, 2) * P.row(m)).transpose())
                                              .sum();
                    if(s == 0)
                    {
                        auto tmp = (alpha.col(0) * beta.col(0)).segment(z1 * C, C).sum();
                        for(y1 = 0; y1 < K; y1++)
                        {
                            ind_post_zy(y1 * C + z1, 0) = tmp * Q(y1, ind);
                            ind_post_y(y1, 0) += ind_post_zy(y1 * C + z1, 0);
                        }
                    }
                }
                if(s == 0) continue;
                MyArr1D alphaprev(C); // previous alpha colsums
                for(z1 = 0; z1 < C; z1++) alphaprev(z1) = alpha(Eigen::seqN(z1, C, C), s - 1).sum();
                for(z1 = 0; z1 < C; z1++)
                {
                    double tmp{0};
                    for(z2 = 0; z2 < C; z2++)
                    {
                        zz = z1 * C + z2;
                        double eb = emit(zz, s) * beta(zz, s);
                        tmp += eb * (R(1, m) * alphaprev(z2) + R(2, m) * H(z2, s));
                        for(y1 = 0; y1 < K; y1++)
                        {
                            ind_post_y(y1, s) +=
                                eb * cs(s) * Q(y1, ind)
                                * (R(0, m) * alpha(zz, s - 1)
                                   + R(1, m) * (alphaprev(z2) * par->F[y1][m * C + z1] + alphaprev(z1) * H(z2, s))
                                   + R(2, m) * par->F[y1][m * C + z1] * H(z2, s));
                        }
                    }
                    for(y1 = 0; y1 < K; y1++)
                        ind_post_zy(y1 * C + z1, s) = tmp * Q(y1, ind) * par->F[y1][m * C + z1] * cs(s);
                }
            }
            ret_post_y[ind] = MyFloat1D(ind_post_y.data(), ind_post_y.data() + ind_post_y.size());
            ret_post_zy[ind] = MyFloat1D(ind_post_zy.data(), ind_post_zy.data() + ind_post_zy.size());
        }
    }
    return List::create(Named("C") = par->C, Named("K") = par->K, Named("S") = S, Named("N") = par->N,
                        Named("gamma") = ret_gamma, Named("ancestry") = ret_post_y,
                        Named("clusterancestry") = ret_post_zy);
}

//' parse options in the fastphase model
//' @param filename path to binary file from impute command
//' @export
// [[Rcpp::export]]
List parse_impute_opt(std::string filename) {
    auto filesize = std::filesystem::file_size(filename);
    std::error_code ec;
    std::ifstream ifs(filename, std::ios::in | std::ios::binary);
    constexpr auto OPTIONS = alpaca::options::fixed_length_encoding;
    std::unique_ptr<BigAss> genome = std::make_unique<BigAss>(alpaca::deserialize<OPTIONS, BigAss>(ifs, filesize, ec));
    ifs.close();
    assert((bool)ec == false);
    return List::create(Named("C") = genome->C,
                        Named("B") = genome->B,
                        Named("G") = genome->G,
                        Named("chunksize") = genome->chunksize,
                        Named("nsamples") = genome->nsamples,
                        Named("nsnps") = genome->nsnps,
                        Named("nchunks") = genome->nchunks);
    
}


//' parse parameters in the fastphase model
//' @param filename path to binary file from impute command
//' @export
// [[Rcpp::export]]
List parse_impute_par(std::string filename, int ic = -1)
{
    auto filesize = std::filesystem::file_size(filename);
    std::error_code ec;
    std::ifstream ifs(filename, std::ios::in | std::ios::binary);
    constexpr auto OPTIONS = alpaca::options::fixed_length_encoding;
    std::unique_ptr<BigAss> genome = std::make_unique<BigAss>(alpaca::deserialize<OPTIONS, BigAss>(ifs, filesize, ec));
    ifs.close();
    assert((bool)ec == false);
    MyArr2D alpha, beta, ae;
    Int1D ids;
    for(int ind = 0; ind < genome->nsamples; ind++) ids.push_back(ind);
    int nchunks = ic < 0 ? genome->nchunks : 1;
    int N = ids.size();
    List ret(N);
    for(auto ind : ids)
    {
        List gammaI(nchunks), aeI(nchunks);
        for(int c = 0; c < nchunks; c++) {
            ic = nchunks > 1 ? c : std::max(ic, c);
            const int iM = genome->pos[ic].size();
            const int nGrids = genome->B > 1 ? (iM + genome->B - 1) / genome->B : iM;
            alpha.setZero(genome->C * genome->C, nGrids);
            beta.setZero(genome->C * genome->C, nGrids);
            get_cluster_probability(ind, iM, alpha, beta, genome->gls[ic], genome->R[ic], genome->PI[ic], genome->F[ic]);
            if(!((1 - (alpha * beta).colwise().sum()).abs() < 1e-6).all()) cao.error("gamma sum is not 1.0!\n");
            ae.setZero(genome->C * genome->C, nGrids);
            get_cluster_frequency(ae, genome->R[ic], genome->PI[ic]);
            gammaI[c] = alpha * beta;
            aeI[c] = ae;
        }
        ret[ind] =  List::create(Named("gamma") = gammaI,
                                 Named("ae") = aeI);
    }
    return ret;
}
