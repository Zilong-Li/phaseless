#include "utils.hpp"

#include "io.hpp"
#include "phaseless.hpp"
#include "threadpool.hpp"
#include <alpaca/alpaca.h>

using namespace std;

std::string make_beagle_header(std::string fam)
{
    std::ifstream ifs(fam);
    std::string line;
    std::string hdr{"marker\tallele1\tallele2"};
    while(getline(ifs, line))
    {
        std::stringstream ss(line);
        std::vector<std::string> token(std::istream_iterator<std::string>{ss}, std::istream_iterator<std::string>{});
        hdr += "\t" + token[1] + "\t" + token[1] + "\t" + token[1];
    }
    hdr += "\n";
    return hdr;
}

int run_parse_main(Options & opts)
{
    cao.cao.open(opts.out + ".log");
    cao.is_screen = !opts.noscreen;
    cao.print(opts.opts_in_effect);
    if(!opts.in_joint.empty())
    {
        opts.in_bin = opts.in_joint;
        auto filesize = std::filesystem::file_size(opts.in_bin);
        std::error_code ec;
        std::ifstream ifs(opts.in_bin, std::ios::in | std::ios::binary);
        constexpr auto OPTIONS = alpaca::options::fixed_length_encoding;
        std::unique_ptr<Pars> par = std::make_unique<Pars>(alpaca::deserialize<OPTIONS, Pars>(ifs, filesize, ec));
        ifs.close();
        assert((bool)ec == false);
        Phaseless faith(par->K, par->C, par->N, par->M, opts.seed);
        faith.setFlags(opts.ptol, opts.ftol, opts.qtol, opts.debug, opts.nQ, opts.nP, opts.nF, opts.nR);
        faith.initRecombination(par->pos);
        faith.setStartPoint(par);
        faith.setStartPoint(opts.in_qfile, opts.in_pfile);
        double loglike, diff, prevlike{std::numeric_limits<double>::lowest()};
        Eigen::IOFormat fmt(6, Eigen::DontAlignCols, " ", "\n");
        vector<future<double>> res;
        std::ofstream oanc(opts.out + ".Q");
        ThreadPool pool(opts.nthreads);
        if(opts.noaccel)
        {
            for(int it = 0; SIG_COND && it <= opts.nimpute; it++)
            {
                tim.clock();
                faith.initIteration();
                for(int i = 0; i < faith.N; i++)
                {
                    if(it == opts.nimpute)
                        res.emplace_back(pool.enqueue(&Phaseless::runBigass, &faith, i, std::ref(par->gls), true));
                    else
                        res.emplace_back(pool.enqueue(&Phaseless::runBigass, &faith, i, std::ref(par->gls), false));
                }
                loglike = 0;
                for(auto && ll : res) loglike += ll.get();
                res.clear(); // clear future and renew
                faith.updateIteration();
                diff = it ? loglike - prevlike : NAN;
                prevlike = loglike;
                cao.print(tim.date(), "run whole genome, iteration", it, ", likelihoods =", loglike, ", diff =", diff,
                          ", time", tim.reltime(), " sec");
                if(diff < opts.ltol)
                {
                    cao.print(tim.date(), "hit stopping criteria, diff =", std::scientific, diff, " <", opts.ltol);
                    break;
                }
            }
        }
        else
        {
            MyArr2D Q0, Q1, Q2, Qt;
            MyArr2D F0, F1, F2, Ft;
            const int istep{4};
            double alpha{0}, stepMax{4}, alphaMax{1280}, logcheck{0};
            for(int it = 0; SIG_COND && (it < opts.nimpute / 4); it++)
            {
                // first normal iter
                faith.initIteration();
                Q0 = faith.Q;
                F0 = cat_stdvec_of_eigen(faith.F);
                for(int i = 0; i < faith.N; i++)
                    res.emplace_back(pool.enqueue(&Phaseless::runBigass, &faith, i, std::ref(par->gls), false));
                loglike = 0;
                for(auto && ll : res) loglike += ll.get();
                res.clear(); // clear future and renew
                faith.updateIteration();
                // second normal iter
                tim.clock();
                faith.initIteration();
                Q1 = faith.Q;
                F1 = cat_stdvec_of_eigen(faith.F);
                for(int i = 0; i < faith.N; i++)
                    res.emplace_back(pool.enqueue(&Phaseless::runBigass, &faith, i, std::ref(par->gls), false));
                loglike = 0;
                for(auto && ll : res) loglike += ll.get();
                res.clear(); // clear future and renew
                faith.updateIteration();
                diff = it ? loglike - prevlike : NAN;
                prevlike = loglike;
                cao.print(tim.date(), "SqS3 iteration", it * 4 + 1, ", alpha=", alpha, ", likelihoods =", std::fixed,
                          loglike, ", diff =", diff, ", time", tim.reltime(), " sec");
                if(diff < opts.ltol)
                {
                    cao.print(tim.date(), "hit stopping criteria, diff =", std::scientific, diff, " <", opts.ltol);
                    break;
                }
                // save for later comparison
                Qt = faith.Q;
                Ft = cat_stdvec_of_eigen(faith.F);
                // calculate alpha based on first two pars
                // alpha = ((Q1 - Q0).square().sum()) / ((faith.Q - 2 * Q1 + Q0).square().sum());
                alpha = ((F1 - F0).square().sum() + (Q1 - Q0).square().sum())
                        / ((cat_stdvec_of_eigen(faith.F) - 2 * F1 + F0).square().sum()
                           + (faith.Q - 2 * Q1 + Q0).square().sum());
                alpha = max(1.0, sqrt(alpha));
                if(alpha >= stepMax)
                {
                    alpha = min(stepMax, alphaMax);
                    stepMax = min(stepMax * istep, alphaMax);
                }
                // third accel iter
                // update Q and F using the second em iter
                faith.Q = Q0 + 2 * alpha * (Q1 - Q0) + alpha * alpha * (faith.Q - 2 * Q1 + Q0);
                for(int k = 0; k < faith.K; k++)
                    faith.F[k] =
                        F0.middleRows(k * faith.C, faith.C)
                        + 2 * alpha * (F1.middleRows(k * faith.C, faith.C) - F0.middleRows(k * faith.C, faith.C))
                        + alpha * alpha
                              * (faith.F[k] - 2 * F1.middleRows(k * faith.C, faith.C)
                                 + F0.middleRows(k * faith.C, faith.C));
                faith.protectPars();
                faith.initIteration();
                for(int i = 0; i < faith.N; i++)
                    res.emplace_back(pool.enqueue(&Phaseless::runBigass, &faith, i, std::ref(par->gls), false));
                loglike = 0;
                for(auto && ll : res) loglike += ll.get();
                res.clear(); // clear future and renew
                faith.updateIteration();
                // save current pars
                Q2 = faith.Q;
                F2 = cat_stdvec_of_eigen(faith.F);
                // check if normal third iter is better
                faith.Q = Qt;
                for(int k = 0; k < faith.K; k++) faith.F[k] = Ft.middleRows(k * faith.C, faith.C);
                faith.initIteration();
                for(int i = 0; i < faith.N; i++)
                    res.emplace_back(pool.enqueue(&Phaseless::runBigass, &faith, i, std::ref(par->gls), false));
                logcheck = 0;
                for(auto && ll : res) logcheck += ll.get();
                res.clear(); // clear future and renew
                faith.updateIteration();
                if(logcheck - loglike > 0.1)
                {
                    if(faith.debug)
                        cao.warn(tim.date(), "normal EM yields better likelihoods than the accelerated EM.", logcheck,
                                 " -", loglike, "> 0.1");
                }
                else
                {
                    faith.Q = Q2;
                    for(int k = 0; k < faith.K; k++) faith.F[k] = F2.middleRows(k * faith.C, faith.C);
                }
            }
        }
        faith.Q = (faith.Q * 1e6).round() / 1e6;
        oanc << std::fixed << faith.Q.transpose().format(fmt) << "\n";
        par->init(faith.K, faith.C, faith.M, faith.N, faith.er, faith.P, faith.Q, faith.F);
        std::ofstream opar(opts.out + ".pars.bin", std::ios::out | std::ios::binary);
        auto bytes_written = alpaca::serialize<OPTIONS, Pars>(*par, opar);
        opar.close();
        assert(std::filesystem::file_size(opts.out + ".pars.bin") == bytes_written);
        cao.done(tim.date(), "parse done and outputting.", bytes_written, " bytes written to file");
        return 0;
    }
    if(!opts.in_impute.empty())
    {
        opts.in_bin = opts.in_impute;
        auto filesize = std::filesystem::file_size(opts.in_bin);
        std::error_code ec;
        std::ifstream ifs(opts.in_bin, std::ios::in | std::ios::binary);
        constexpr auto OPTIONS = alpaca::options::fixed_length_encoding;
        std::unique_ptr<BigAss> genome =
            std::make_unique<BigAss>(alpaca::deserialize<OPTIONS, BigAss>(ifs, filesize, ec));
        ifs.close();
        assert((bool)ec == false);
        int ic = opts.ichunk;
        if(ic >= genome->nchunks)
        {
            cao.error("the chunk ", ic, " (0-based) to be extracted is not less than total chunks", genome->nchunks);
            return 1;
        }
        else if(ic < 0)
        {
            cao.warn("the chunk ", ic, " to be extracted is less than 0. all chunks will be extracted!");
        }
        Int1D ids;
        if(opts.samples.empty())
        {
            for(int ind = 0; ind < genome->nsamples; ind++) ids.push_back(ind);
        }
        else
        {
            UMapStringInt ids_m;
            int i = 0;
            for(auto & id : genome->sampleids) ids_m[id] = i++;
            std::ifstream ifs(opts.samples);
            if(!ifs.is_open()) cao.error("can not open the file: ", opts.samples);
            std::string line;
            while(getline(ifs, line)) ids.push_back(ids_m[line]);
        }
        // haplike is p(Z_is |X_is , theta) = alpha * beta / p(X|theta) = gamma
        std::ofstream ofs_alpha(opts.out + ".alpha.bin", std::ios::binary);
        std::ofstream ofs_beta(opts.out + ".beta.bin", std::ios::binary);
        int nsamples = ids.size();
        ofs_alpha.write((char *)&genome->C, 4);
        ofs_alpha.write((char *)&nsamples, 4);
        ofs_beta.write((char *)&genome->C, 4);
        ofs_beta.write((char *)&nsamples, 4);
        MyArr2D alpha, beta;
        if(ic < 0)
        {
            ofs_alpha.write((char *)&genome->G, 4);
            ofs_beta.write((char *)&genome->G, 4);
            for(auto ind : ids)
            {
                for(ic = 0; ic < genome->nchunks; ic++)
                {
                    const int iM = genome->pos[ic].size();
                    const int nGrids = genome->B > 1 ? (iM + genome->B - 1) / genome->B : iM;
                    alpha.setZero(genome->C * genome->C, nGrids);
                    beta.setZero(genome->C * genome->C, nGrids);
                    get_cluster_likelihood(ind, iM, alpha, beta, genome->gls[ic], genome->R[ic], genome->PI[ic],
                                           genome->F[ic]);
                    if(!((1 - (alpha * beta).colwise().sum()).abs() < 1e-6).all()) cao.error("gamma sum is not 1.0!\n");
                    ofs_alpha.write((char *)alpha.data(), genome->C * genome->C * nGrids * sizeof(MyFloat));
                    ofs_beta.write((char *)beta.data(), genome->C * genome->C * nGrids * sizeof(MyFloat));
                }
            }
        }
        else
        {
            const int iM = genome->pos[ic].size();
            const int nGrids = genome->B > 1 ? (iM + genome->B - 1) / genome->B : iM;
            ofs_alpha.write((char *)&nGrids, 4);
            ofs_beta.write((char *)&nGrids, 4);
            std::ofstream ofs_ae(opts.out + ".cluster.freq");
            Eigen::IOFormat fmt(6, Eigen::DontAlignCols, "\t", "\n");
            MyArr2D ae = MyArr2D::Zero(genome->C, nGrids);
            for(auto ind : ids)
            {
                alpha.setZero(genome->C * genome->C, nGrids);
                beta.setZero(genome->C * genome->C, nGrids);
                get_cluster_likelihood(ind, iM, alpha, beta, genome->gls[ic], genome->R[ic], genome->PI[ic],
                                       genome->F[ic]);
                if(!((1 - (alpha * beta).colwise().sum()).abs() < 1e-6).all()) cao.error("gamma sum is not 1.0!\n");
                ofs_alpha.write((char *)alpha.data(), genome->C * genome->C * nGrids * sizeof(MyFloat));
                ofs_beta.write((char *)beta.data(), genome->C * genome->C * nGrids * sizeof(MyFloat));
                for(int i = 0; i < nGrids; i++)
                    ae.col(i) += (alpha.col(i) * beta.col(i)).reshaped(genome->C, genome->C).colwise().sum();
            }
            ae.rowwise() /= ae.colwise().sum();
            ofs_ae << ae.transpose().format(fmt) << "\n";
        }
        return 0;
    }
    return 0;
}

int run_convert_main(Options & opts)
{
    cao.cao.open(opts.out + ".log");
    cao.is_screen = !opts.noscreen;
    cao.print(opts.opts_in_effect);
    std::setlocale(LC_ALL, "C"); // use minial locale
    std::ios_base::sync_with_stdio(false); // don't sync
    uint64_t nsamples = count_lines(opts.in_plink + ".fam");
    uint64_t nsnps = count_lines(opts.in_plink + ".bim");
    cao.print("nsamples :", nsamples, ", nsnps :", nsnps);
    std::string hdr = make_beagle_header(opts.in_plink + ".fam");
    int nchunks = (nsnps + opts.chunksize - 1) / opts.chunksize;
    std::ifstream ifs_bed(opts.in_plink + ".bed", std::ios::in | std::ios::binary);
    uint8_t header[3];
    ifs_bed.read(reinterpret_cast<char *>(&header[0]), 3);
    if((header[0] != 0x6c) || (header[1] != 0x1b) || (header[2] != 0x01))
        cao.error("Incorrect magic number in plink bed file.\n");
    ThreadPool pool(opts.nthreads);
    vector<future<string>> res;
    std::ifstream ifs_bim(opts.in_plink + ".bim", std::ios::in);
    for(int ic = 0; ic < nchunks; ic++)
    {
        const auto [bed, marker] = read_plink_bed(ifs_bed, ifs_bim, nsamples, opts.chunksize);
        res.emplace_back(pool.enqueue(convert_geno2like, bed, marker, nsamples));
    }
    string out{opts.out + ".gz"};
    gzFile gzfp = gzopen(out.c_str(), "wb");
    gzwrite(gzfp, hdr.c_str(), hdr.size());
    for(auto && ll : res)
    {
        auto gl = ll.get();
        gzwrite(gzfp, gl.c_str(), gl.size());
    }
    gzclose(gzfp);
    return 0;
}
