/*******************************************************************************
 * @file        https://github.com/Zilong-Li/phaseless/src/main.cpp
 * @author      Zilong Li
 * Copyright (C) 2023. The use of this code is governed by the LICENSE file.
 ******************************************************************************/
#define _DECLARE_TOOLBOX_HERE

#include "admixture.hpp"
#include "fastphase.hpp"
#include "phaseless.hpp"
#include "utils.hpp"
#include <argparse/argparse.hpp>
#include <signal.h>

using namespace argparse;

int main(int argc, char * argv[])
{
    // ========= helper message and parameters parsing ===========================

    const std::string VERSION{"0.6.0"};

    // below for catching ctrl+c, and dumping files
    struct sigaction sa;
    sigemptyset(&sa.sa_mask);
    sa.sa_flags = 0;
    sa.sa_handler = handler;
    sigaction(SIGPIPE, &sa, 0);
    sigaction(SIGINT, &sa, 0);

    // clang-format off
    ArgumentParser program("phaseless", VERSION, default_arguments::version);
    program.add_epilog("This project is still under development!\n"
                       "Contact: zilong.dk@gmail.com");
    program.add_argument("-D","--debug")
        .help("enable debug mode")
        .flag();
    program.add_argument("-S", "--no-stdout")
        .help("disable print log to screen")
        .flag();
    program.add_argument("-a", "--no-accel")
        .help("disable accelerated EM")
        .flag();
    program.add_argument("-q","--NQ")
        .help("disable updating Q")
        .flag();
    program.add_argument("-p", "--NP")
        .help("disable updating P")
        .flag();
    program.add_argument("-r", "--NR")
        .help("disable updating R")
        .flag();
    program.add_argument("-f","--NF")
        .help("disable updating F")
        .flag();
    program.add_argument("-F","--write-F")
        .help("output F")
        .flag();
    program.add_argument("--ltol")
        .help("convergence tolerance for difference in log likelihoods")
        .default_value(1e-2)
        .scan<'g', double>();
    program.add_argument("--ptol")
        .help("lower boundary for P")
        .default_value(1e-6)
        .scan<'g', double>();
    program.add_argument("--ftol")
        .help("lower boundary for F")
        .default_value(1e-9)
        .scan<'g', double>();
    program.add_argument("--qtol")
        .help("lower boundary for Q")
        .default_value(1e-9)
        .scan<'g', double>();
    program.add_argument("--qfile")
        .help("read Q file as the start point")
        .default_value(std::string{""});
    program.add_argument("--pfile")
        .help("read P file as the start point")
        .default_value(std::string{""});
    program.add_argument("--rfile")
        .help("read R file as the start point")
        .default_value(std::string{""});

    argparse::ArgumentParser cmd_joint("joint", VERSION, default_arguments::help);
    cmd_joint.add_description("run phasing and admixture inference in one goal");
    cmd_joint.add_argument("-c", "--cluster")
        .help("number of haplotype clusters")
        .default_value(8)
        .scan<'i', int>();
    cmd_joint.add_argument("-k", "--ancestry")
        .help("number of ancestries (required)")
        .required()
        .scan<'i', int>();
    cmd_joint.add_argument("-g", "--beagle")
        .help("gziped beagle format as input")
        .default_value(std::string{""});
    cmd_joint.add_argument("-i", "--iterations")
        .help("number of EM iterations")
        .default_value(1000)
        .scan<'i', int>();
    cmd_joint.add_argument("-n", "--threads")
        .help("number of threads; -1 selects all available CPUs")
        .default_value(-1)
        .scan<'i', int>();
    cmd_joint.add_argument("--gpu")
        .help("run the joint-model E step on an NVIDIA CUDA GPU")
        .flag();
    cmd_joint.add_argument("-o", "--out")
        .help("output prefix")
        .default_value(std::string{"joint"});
    cmd_joint.add_argument("-s", "--chunksize")
        .help("size of each chunk in sites unit ")
        .default_value(32000)
        .scan<'i', int>();
    cmd_joint.add_argument("-S", "--single-chunk")
        .help("treat input as big single chunk")
        .flag();
    cmd_joint.add_argument("-V", "--vcf")
        .help("output the VCF file")
        .flag();
    cmd_joint.add_argument("-Q", "--aQ")
        .help("aphla is accelarated with Q only")
        .flag();
    cmd_joint.add_argument("-d","--seed")
        .help("seed for reproducibility")
        .default_value(996)
        .scan<'i', int>();
    cmd_joint.add_argument("--conv-gap-tol")
        .help("observed log likelihood improvement per observation")
        .default_value(1e-5)
        .scan<'g', double>();
    cmd_joint.add_argument("--conv-relative-tol")
        .help("relative log likelihood convergence tolerance")
        .default_value(2e-6)
        .scan<'g', double>();
    cmd_joint.add_argument("--conv-parameter-tol")
        .help("joint parameter stability tolerance")
        .default_value(1e-3)
        .scan<'g', double>();
    cmd_joint.add_argument("--conv-stable-iterations")
        .help("consecutive stable accepted iterations required")
        .default_value(3)
        .scan<'i', int>();
    cmd_joint.add_argument("--stitch-heuristics")
        .help("enable STITCH-inspired alignment and revival inside shared-haplotype initialization")
        .flag();
    cmd_joint.add_argument("--random-init")
        .help("skip posterior-driven joint initialization and retain the legacy random start")
        .flag();
    cmd_joint.add_argument("--init-haplotype-iterations")
        .help("maximum ordinary EM scans used to learn the shared haplotype start")
        .default_value(40)
        .scan<'i', int>();
    cmd_joint.add_argument("--init-haplotype-min-iterations")
        .help("minimum shared-haplotype scans before adaptive stopping")
        .default_value(20)
        .scan<'i', int>();
    cmd_joint.add_argument("--init-haplotype-relative-tol")
        .help("relative likelihood tolerance for adaptive shared-haplotype stopping")
        .default_value(5e-4)
        .scan<'g', double>();
    cmd_joint.add_argument("--init-haplotype-profile-tol")
        .help("posterior cluster-profile RMS tolerance for adaptive shared-haplotype stopping")
        .default_value(5e-3)
        .scan<'g', double>();
    cmd_joint.add_argument("--init-haplotype-stable-iterations")
        .help("consecutive stable shared-haplotype scans required")
        .default_value(3)
        .scan<'i', int>();
    cmd_joint.add_argument("--init-profile-pruning")
        .help("use low-cost informative-SNP pruning in the posterior profile used to initialize Q and F")
        .flag();
    cmd_joint.add_argument("--init-profile-block-size")
        .help("SNP block size for low-cost informative-site selection in the initialization profile")
        .default_value(100)
        .scan<'i', int>();
    cmd_joint.add_argument("--init-profile-information-fraction")
        .help("fraction of cluster-allele mutual information retained within each profile block")
        .default_value(0.95)
        .scan<'g', double>();
    cmd_joint.add_argument("--init-profile-min-snps")
        .help("minimum informative SNPs retained per initialization-profile block")
        .default_value(5)
        .scan<'i', int>();
    cmd_joint.add_argument("--init-ancestry-iterations")
        .help("ordinary EM scans used to refine each posterior-driven ancestry start")
        .default_value(15)
        .scan<'i', int>();
    cmd_joint.add_argument("--init-noise")
        .help("relative jitter and soft-assignment temperature for posterior-driven ancestry starts")
        .default_value(0.05)
        .scan<'g', double>();
    cmd_joint.add_argument("--init-restarts")
        .help("posterior-driven ancestry starts to try, retaining the highest likelihood")
        .default_value(3)
        .scan<'i', int>();
    cmd_joint.add_argument("--q-pseudocount")
        .help("symmetric pseudocount added to each ancestry when updating Q; 0 disables it")
        .default_value(0.5)
        .scan<'g', double>();
    cmd_joint.add_argument("--p-shrinkage")
        .help("site-frequency-centered prior weight for P updates; 0 disables it")
        .default_value(0.5)
        .scan<'g', double>();
    cmd_joint.add_argument("--heuristic-block-size")
        .help("SNP block size used by STITCH-inspired heuristics")
        .default_value(100)
        .scan<'i', int>();
    cmd_joint.add_argument("--heuristic-reset-radius")
        .help("number of SNPs reset on either side of a relabelled boundary")
        .default_value(20)
        .scan<'i', int>();
    cmd_joint.add_argument("--heuristic-min-usage")
        .help("posterior chromosome-copy usage below which a cluster interval is revived")
        .default_value(0.01)
        .scan<'g', double>();
    cmd_joint.add_argument("--heuristic-donor-weight")
        .help("weight of the sampled donor emissions when reviving a cluster")
        .default_value(0.8)
        .scan<'g', double>();
    // cmd_joint.add_parents(program);

    argparse::ArgumentParser cmd_impute("impute", VERSION, default_arguments::help);
    cmd_impute.add_description("run imputation for low coverage sequencing data");
    cmd_impute.add_argument("-c", "--cluster")
        .help("number of ancestral haplotype clusters")
        .default_value(10)
        .scan<'i', int>();
    cmd_impute.add_argument("-C", "--collapse")
        .help("collapse SNPs in a reasonable window")
        .flag();
    cmd_impute.add_argument("-B", "--grid-size")
        .help("number of SNPs (>=3) in each grid. 1 disables collapsing")
        .default_value(1)
        .scan<'i', int>();
    cmd_impute.add_argument("-f", "--vcf")
        .help("vcf/bcf format with GL/PL tag as input")
        .default_value(std::string{""});
    cmd_impute.add_argument("-g", "--beagle")
        .help("gziped beagle format as input")
        .required()
        .default_value(std::string{""});
    cmd_impute.add_argument("-i", "--iterations")
        .help("number of EM iterations")
        .default_value(100)
        .scan<'i', int>();
    cmd_impute.add_argument("-n", "--threads")
        .help("number of threads")
        .default_value(10)
        .scan<'i', int>();
    cmd_impute.add_argument("-o", "--out")
        .help("output prefix")
        .default_value(std::string{"impute"});
    cmd_impute.add_argument("-r", "--region")
        .help("region in vcf/bcf to subset")
        .default_value(std::string{""});
    cmd_impute.add_argument("-s", "--chunksize")
        .help("size of each chunk in sites unit ")
        .default_value(10000)
        .scan<'i', int>();
    cmd_impute.add_argument("-S", "--single-chunk")
        .help("treat input as big single chunk")
        .flag();
    cmd_impute.add_argument("-d","--seed")
        .help("seed for reproducibility")
        .default_value(996)
        .scan<'i', int>();
    cmd_impute.add_argument("--write-hapsum")
        .help("write Hapsum instead of AE into parse.bin")
        .flag();
    cmd_impute.add_argument("--refill-haps")
        .help("refill infrequently used haplotype clusters.\n"
              "1: reset P to min allele emission probability for that haplotype cluster\n"
              "2: re-sample P by copying from haplotype with the highest probability\n"
              "3: re-sample P by copying from others with respect to their probability\n"
              "0: disable this")
        .default_value(2)
        .scan<'i', int>();
    cmd_impute.add_argument("--minRecombRate")
        .help("min recombination rate to determine if a SNP should be collapsed")
        .default_value(1e-4)
        .scan<'g', double>();
    // cmd_impute.add_parents(program);

    argparse::ArgumentParser cmd_admix("admix", VERSION, default_arguments::help);
    cmd_admix.add_description("run admixture with cluster likelihoods as input");
    cmd_admix.add_argument("-b", "--bin")
        .help("binary format from impute command as input")
        .default_value(std::string{""});
    cmd_admix.add_argument("-k", "--ancestry")
        .help("number of ancestry in admixture assumption")
        .default_value(2)
        .scan<'i', int>();
    cmd_admix.add_argument("-i", "--iterations")
        .help("number of maximun EM iterations")
        .default_value(2000)
        .scan<'i', int>();
    cmd_admix.add_argument("-n", "--threads")
        .help("number of threads")
        .default_value(10)
        .scan<'i', int>();
    cmd_admix.add_argument("-o", "--out")
        .help("output prefix")
        .default_value(std::string{"admix"});
    cmd_admix.add_argument("-d","--seed")
        .help("seed for reproducibility")
        .default_value(996)
        .scan<'i', int>();
    cmd_admix.add_argument("-f", "--force-accept")
        .help("always accept the acceleration solution")
        .flag();
    cmd_admix.add_argument("-F", "--constrain-F")
        .help("apply constraint on F so that it is not smaller than cluster frequency in fastphase model")
        .flag();
    cmd_admix.add_argument("-P", "--min-P")
        .help("set cluster likelihood to zeros if P (in fastphase) < min-P")
        .default_value(0.0)
        .scan<'g', double>();

    argparse::ArgumentParser cmd_convert("convert", VERSION, default_arguments::help);
    cmd_convert.add_description("different file format converter");
    cmd_convert.add_argument("-i","--input")
        .help("input file to be converted")
        .default_value(std::string{""});
    cmd_convert.add_argument("-o","--output")
        .help("output file prefix")
        .default_value(std::string{"convert"});
    cmd_convert.add_argument("-p", "--plink2beagle")
        .help("use plink1 file as input without .bed")
        .flag();
    cmd_convert.add_argument("-n", "--threads")
        .help("number of threads")
        .default_value(4)
        .scan<'i', int>();
    cmd_convert.add_argument("-s", "--chunksize")
        .help("size of each chunk in sites unit ")
        .default_value(10000)
        .scan<'i', int>();
    // cmd_convert.add_parents(program);

    program.add_subparser(cmd_impute);
    program.add_subparser(cmd_admix);
    program.add_subparser(cmd_convert);
    program.add_subparser(cmd_joint);
    // clang-format on

    try
    {
        Options opts;
        for(int i = 0; i < argc; i++) opts.opts_in_effect += " " + std::string{argv[i]};
        opts.opts_in_effect += "\nVersion: " + VERSION + "\n" + get_machine();
        program.parse_args(argc, argv);
        opts.debug = program.get<bool>("--debug");
        opts.noscreen = program.get<bool>("--no-stdout");
        opts.in_rfile.assign(program.get("--rfile"));
        opts.in_qfile.assign(program.get("--qfile"));
        opts.in_pfile.assign(program.get("--pfile"));
        opts.ptol = program.get<double>("--ptol");
        opts.ftol = program.get<double>("--ftol");
        opts.qtol = program.get<double>("--qtol");
        if(opts.ptol <= 0 || opts.ptol >= 0.5 || opts.ftol <= 0 || opts.ftol >= 0.5
           || opts.qtol <= 0 || opts.qtol >= 0.5)
            throw std::invalid_argument("P, F, and Q lower boundaries must be between zero and 0.5");
        opts.nQ = program.get<bool>("--NQ");
        opts.nP = program.get<bool>("--NP");
        opts.nR = program.get<bool>("--NR");
        opts.nF = program.get<bool>("--NF");
        opts.oF = program.get<bool>("--write-F");
        opts.ltol = program.get<double>("--ltol");
        opts.noaccel = program.get<bool>("--no-accel");

        if(program.is_subcommand_used(cmd_joint))
        {
            opts.in_beagle.assign(cmd_joint.get("--beagle"));
            opts.out.assign(cmd_joint.get("--out"));
            opts.C = cmd_joint.get<int>("--cluster");
            opts.K = cmd_joint.get<int>("--ancestry");
            opts.nthreads = cmd_joint.get<int>("--threads");
            opts.gpu = cmd_joint.get<bool>("--gpu");
            opts.nimpute = cmd_joint.get<int>("--iterations");
            opts.seed = cmd_joint.get<int>("--seed");
            opts.conv_gap_tol = cmd_joint.get<double>("--conv-gap-tol");
            opts.conv_relative_tol = cmd_joint.get<double>("--conv-relative-tol");
            opts.conv_parameter_tol = cmd_joint.get<double>("--conv-parameter-tol");
            opts.conv_stable_iterations = cmd_joint.get<int>("--conv-stable-iterations");
            opts.stitch_heuristics = cmd_joint.get<bool>("--stitch-heuristics");
            opts.random_init = cmd_joint.get<bool>("--random-init");
            opts.init_haplotype_iterations = cmd_joint.get<int>("--init-haplotype-iterations");
            opts.init_haplotype_min_iterations = cmd_joint.get<int>("--init-haplotype-min-iterations");
            opts.init_haplotype_relative_tol = cmd_joint.get<double>("--init-haplotype-relative-tol");
            opts.init_haplotype_profile_tol = cmd_joint.get<double>("--init-haplotype-profile-tol");
            opts.init_haplotype_stable_iterations = cmd_joint.get<int>("--init-haplotype-stable-iterations");
            opts.init_profile_pruning = cmd_joint.get<bool>("--init-profile-pruning");
            opts.init_profile_block_size = cmd_joint.get<int>("--init-profile-block-size");
            opts.init_profile_information_fraction =
                cmd_joint.get<double>("--init-profile-information-fraction");
            opts.init_profile_min_snps = cmd_joint.get<int>("--init-profile-min-snps");
            opts.init_ancestry_iterations = cmd_joint.get<int>("--init-ancestry-iterations");
            opts.init_noise = cmd_joint.get<double>("--init-noise");
            opts.init_restarts = cmd_joint.get<int>("--init-restarts");
            opts.q_pseudocount = cmd_joint.get<double>("--q-pseudocount");
            opts.p_shrinkage = cmd_joint.get<double>("--p-shrinkage");
            opts.heuristic_block_size = cmd_joint.get<int>("--heuristic-block-size");
            opts.heuristic_reset_radius = cmd_joint.get<int>("--heuristic-reset-radius");
            opts.heuristic_min_usage = cmd_joint.get<double>("--heuristic-min-usage");
            opts.heuristic_donor_weight = cmd_joint.get<double>("--heuristic-donor-weight");
            if(opts.conv_gap_tol <= 0 || opts.conv_relative_tol <= 0 || opts.conv_parameter_tol <= 0
               || opts.conv_stable_iterations < 1)
                throw std::invalid_argument("joint convergence tolerances and stable iterations must be positive");
            if(opts.heuristic_block_size < 1 || opts.heuristic_reset_radius < 0
               || opts.heuristic_min_usage < 0 || opts.heuristic_min_usage >= 1
               || opts.heuristic_donor_weight < 0 || opts.heuristic_donor_weight > 1)
                throw std::invalid_argument("invalid STITCH heuristic configuration");
            if(opts.init_haplotype_iterations < 1 || opts.init_haplotype_min_iterations < 1
               || opts.init_haplotype_min_iterations > opts.init_haplotype_iterations
               || (opts.stitch_heuristics && !opts.random_init && opts.in_qfile.empty()
                   && opts.init_haplotype_iterations < opts.init_haplotype_min_iterations + 8)
               || opts.init_haplotype_relative_tol <= 0 || opts.init_haplotype_profile_tol <= 0
               || opts.init_haplotype_stable_iterations < 1 || opts.init_ancestry_iterations < 1
               || opts.init_restarts < 1 || opts.init_noise < 0 || opts.init_noise > 1
               || opts.init_profile_block_size < 1 || opts.init_profile_min_snps < 1
               || opts.init_profile_min_snps > opts.init_profile_block_size
               || !std::isfinite(opts.init_profile_information_fraction)
               || opts.init_profile_information_fraction <= 0
               || opts.init_profile_information_fraction > 1
               || !std::isfinite(opts.q_pseudocount) || opts.q_pseudocount < 0
               || !std::isfinite(opts.p_shrinkage) || opts.p_shrinkage < 0)
                throw std::invalid_argument("invalid joint initialization configuration");
            opts.chunksize = cmd_joint.get<int>("--chunksize");
            opts.single_chunk = cmd_joint.get<bool>("--single-chunk");
            opts.oVCF = cmd_joint.get<bool>("--vcf");
            opts.aQ = cmd_joint.get<bool>("--aQ");
            if(opts.single_chunk) opts.chunksize = INT_MAX;
            if((opts.in_beagle.empty() && opts.in_vcf.empty()) || cmd_joint.get<bool>("--help"))
                throw std::runtime_error(cmd_joint.help().str());
            run_phaseless_main(opts);
        }
        else if(program.is_subcommand_used(cmd_impute))
        {
            opts.in_beagle.assign(cmd_impute.get("--beagle"));
            opts.in_vcf.assign(cmd_impute.get("--vcf"));
            opts.out.assign(cmd_impute.get("--out"));
            opts.C = cmd_impute.get<int>("--cluster");
            opts.gridsize = cmd_impute.get<int>("--grid-size");
            opts.nthreads = cmd_impute.get<int>("--threads");
            opts.nimpute = cmd_impute.get<int>("--iterations");
            opts.seed = cmd_impute.get<int>("--seed");
            opts.chunksize = cmd_impute.get<int>("--chunksize");
            opts.single_chunk = cmd_impute.get<bool>("--single-chunk");
            opts.eHap = cmd_impute.get<bool>("--write-hapsum");
            opts.collapse = cmd_impute.get<bool>("--collapse");
            opts.refillHaps = cmd_impute.get<int>("--refill-haps");
            opts.tol_r = cmd_impute.get<double>("--minRecombRate");
            if(opts.single_chunk) opts.chunksize = INT_MAX;
            if((opts.in_beagle.empty() && opts.in_vcf.empty()) || cmd_impute.get<bool>("--help"))
                throw std::runtime_error(cmd_impute.help().str());
            run_impute_main(opts);
        }
        else if(program.is_subcommand_used(cmd_admix))
        {
            opts.ptol = cmd_admix.get<double>("--min-P");
            opts.in_bin.assign(cmd_admix.get("--bin"));
            opts.out.assign(cmd_admix.get("--out"));
            opts.seed = cmd_admix.get<int>("--seed");
            opts.K = cmd_admix.get<int>("-k");
            opts.nthreads = cmd_admix.get<int>("--threads");
            opts.nadmix = cmd_admix.get<int>("--iterations");
            opts.cF = cmd_admix.get<bool>("--constrain-F");
            opts.force = cmd_admix.get<bool>("--force-accept");
            if(opts.in_bin.empty() || cmd_admix.get<bool>("--help")) throw std::runtime_error(cmd_admix.help().str());
            run_admix_main(opts);
        }
        else if(program.is_subcommand_used(cmd_convert))
        {
            opts.nthreads = cmd_convert.get<int>("--threads");
            opts.in_plink = cmd_convert.get("--input");
            opts.out = cmd_convert.get("--output");
            opts.chunksize = cmd_convert.get<int>("--chunksize");
            if(opts.in_plink.empty() || cmd_convert.get<bool>("--help"))
                throw std::runtime_error(cmd_convert.help().str());
            if(cmd_convert.get<bool>("--plink2beagle")) run_convert_main(opts);
        }
        else
        {
            cao.cerr(program.help().str());
            std::exit(1);
        }
    }
    catch(const std::exception & err)
    {
        std::cerr << err.what() << std::endl;
        std::exit(1);
    }

    return 0;
}
