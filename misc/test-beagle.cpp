// -*- compile-command: "g++ beagle.cpp -o beagle -std=c++17 -g -O3 -Wall -lz && ./beagle ../data/bgl.gz > t"; -*-
#include <iomanip>
#include <iostream>
#include <vector>
#include <zlib.h>

using namespace std;

using FloatVec1D = std::vector<float>;
using FloatVec2D = std::vector<FloatVec1D>;
using DoubleVec1D = std::vector<double>;
using DoubleVec2D = std::vector<DoubleVec1D>;

/*
** @param gz   file hander returned by gzopen
** @param buf  buffer used for storing data
** @param size buffer size for realloc buffer
** @return length and buffer of current line
*/
int zgets(gzFile gz, char** buf, uint64_t* size)
{
    int rlen = 0;
    char* tok = gzgets(gz, *buf + rlen, *size - rlen); // return buf or NULL
    if (!tok)
        return rlen;
    int tmp = tok ? strlen(tok) : 0;
    if (tok[tmp - 1] != '\n')
    {
        // expand buf size if no end-of-line found
        rlen += tmp;
        *size *= 2;
        *buf = (char*)realloc(*buf, *size);
    }
    rlen += tmp;
    return rlen;
}

FloatVec2D read_genotype_likelihoods(const std::string& beagle)
{
    // VARIBLES
    gzFile fp = nullptr;
    char *original, *buffer, *tok;
    uint64_t bufsize = (uint64_t)128 * 1024 * 1024;
    int nsamples, nsnps, i, j;
    const char* delims = "\t \n";
    FloatVec2D gls;
    FloatVec1D gli;

    // PARSE BEAGLE FILE
    fp = gzopen(beagle.c_str(), "r");
    original = buffer = (char*)calloc(bufsize, sizeof(char));
    zgets(fp, &buffer, &bufsize);
    if (buffer != original)
        original = buffer;
    strtok_r(buffer, delims, &buffer);
    int nCol = 1;
    while ((tok = strtok_r(NULL, delims, &buffer)))
    {
        nCol++;
        if (nCol > 3)
            gli.push_back(strtod(tok, NULL));
    }
    if (nCol % 3)
        throw std::runtime_error("Number of columns should be a multiple of 3.\n");
    nsamples = nCol / 3 - 1;
    nsnps = 1;
    assert(gli.size() == nsamples * 3);
    gls.push_back(gli);
    buffer = original;
    while (zgets(fp, &buffer, &bufsize))
    {
        if (buffer != original)
            original = buffer;
        tok = strtok_r(buffer, delims, &buffer); // id
        tok = strtok_r(NULL, delims, &buffer);   // ref
        tok = strtok_r(NULL, delims, &buffer);   // alt
        for (i = 0; i < nsamples; i++)
        {
            tok = strtok_r(NULL, delims, &buffer);
            gli[3 * i + 0] = strtod(tok, NULL);
            tok = strtok_r(NULL, delims, &buffer);
            gli[3 * i + 1] = strtod(tok, NULL);
            tok = strtok_r(NULL, delims, &buffer);
            gli[3 * i + 2] = strtod(tok, NULL);
        }
        gls.push_back(gli);
        nsnps++;
    }
    gzclose(fp);

    cerr << nsamples << "," << nsnps << endl;
    cerr << gli.size() << "," << gls.size() << endl;
    for (i = 0; i < nsnps; i++)
    {
        for (j = 0; j < nsamples * 3; j++)
        {
            cout << std::fixed << std::setprecision(6) << gls[i][j] << "\t";
        }
        cout << endl;
    }

    return gls;
}

int main(int argc, char* argv[])
{
    read_genotype_likelihoods(argv[1]);
    return 0;
}
