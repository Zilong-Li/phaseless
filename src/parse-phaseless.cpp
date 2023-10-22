#include "phaseless.hpp"
#include <Rcpp.h>
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

//' parse posterior probabilty of Z in joint model
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
    const int CC = par->C * par->C;
    Eigen::Map<const MyArr2D> Q(par->Q.data(), par->K, par->N);
    Eigen::Map<const MyArr2D> P(par->P.data(), par->M, par->C);
    Eigen::Map<const MyArr1D> er(par->er.data(), par->M);
    auto R = er2R(er);
    MyFloat2D res(par->N);
    int S{0};
    for(int ind = 0; ind < par->N; ind++)
    {
        for(int ic = 0, pos_chunk = 0; ic < nchunks; ic++)
        {
            S = par->pos[ic].size();
            pos_chunk += S;
            if(ic == chunk)
            {
                Eigen::Map<const MyArr2D> gli(par->gls[ic].data() + ind * S * 3, S, 3);
                MyArr2D emit = get_emission_by_gl(gli, P.middleRows(pos_chunk - S, S)).transpose(); // CC x S
                MyArr2D alpha(CC, S), beta(CC, S);
                // first get H ie old PI in fastphase
                MyArr2D H = MyArr2D::Zero(C, S);
                int z1, y1, s; // m * C + z1
                for(s = 0; s < S; s++)
                    for(z1 = 0; z1 < C; z1++)
                        for(y1 = 0; y1 < par->K; y1++)
                            H(z1, s) += Q(y1, ind) * par->F[y1][(pos_chunk - S + s) * C + z1];
                forward_backwards_diploid(alpha, beta, emit, R.middleCols(pos_chunk - S, S),
                                          H); // cs is 1 / colsum(alpha)
                MyArr2D gamma = alpha * beta;
                res[ind] = MyFloat1D(gamma.data(), gamma.data() + gamma.size());
                break;
            }
        }
    }
    return List::create(Named("C") = par->C, Named("K") = par->K, Named("S") = S, Named("N") = par->N,
                        Named("gamma") = res);
}
