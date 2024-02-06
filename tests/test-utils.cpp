#define _DECLARE_TOOLBOX_HERE

#include "../src/common.hpp"
#include "catch.hh"

using namespace std;
using namespace Eigen;

TEST_CASE("calc_grid_distance", "[test-utils]")
{
    // Int1D pos{1,2,3,4,5,6,7,8,9,10,11};
    Int1D pos{1, 2, 3, 4, 5, 6, 7, 8, 9};
    const auto b = find_grid_to_collapse(pos, 4);
    
    Char1D cb(b.size());
    for(int i = 0; i < b.size();i++) cb[i] = (b(i)==true);
    for(auto i : cb) cao.cerr(i==1);
    auto d = calc_grid_distance(pos, b);
    // cao.cerr("size:", d.size());
}
