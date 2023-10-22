#include <Rcpp.h>
#include "common.hpp"
#include "phaseless.hpp"
#include <alpaca/alpaca.h>

using namespace Rcpp;
using namespace std;

//' parse parameters in joint model
//' @param filename path to binary file from joint command
//' @export
// [[Rcpp::export]]
List parse_joint(std::string filename)
{
    auto filesize = std::filesystem::file_size(filename);
    std::error_code ec;
    std::ifstream ifs(filename, std::ios::in | std::ios::binary);
    constexpr auto OPTIONS = alpaca::options::fixed_length_encoding;
    std::unique_ptr<Pars> par = std::make_unique<Pars>(alpaca::deserialize<OPTIONS, Pars>(ifs, filesize, ec));
    ifs.close();
    assert((bool)ec == false);
    return List::create(Named("F") = par->F);
}
