#include "../external/argparse/argparse.hpp"
#include "../src/io.hpp"
#include "catch.hh"

using namespace std;

TEST_CASE("test parse negative integer", "[test-cli]")
{
    argparse::ArgumentParser program;
    program.add_argument("--verbose", "-v")
        .help("enable verbose logging")
        .default_value(false)
        .implicit_value(true);

    program.add_argument("-n","--number").help("Input number").scan<'i', int>();

    program.parse_args({"./main", "-n", "-1"});
    REQUIRE(program.get<int>("--number") == -1);
}
