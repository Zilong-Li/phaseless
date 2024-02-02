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

    const std::string VERSION{"0.4.1"};

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
    program.add_argument("-l", "--ltol")
        .help("convergence tolerance of difference in log likelihoods")
        .default_value(1e-1)
        .scan<'g', double>();
    program.add_argument("-P", "--ptol")
        .help("lower boundary for P")
        .default_value(1e-4)
        .scan<'g', double>();
    program.add_argument("-F", "--ftol")
        .help("lower boundary for F")
        .default_value(1e-4)
        .scan<'g', double>();
    program.add_argument("-Q", "--qtol")
        .help("lower boundary for Q")
        .default_value(1e-6)
        .scan<'g', double>();
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
        .default_value(10)
        .scan<'i', int>();
    cmd_joint.add_argument("-k", "--ancestry")
        .help("number of ancestry")
        .default_value(3)
        .scan<'i', int>();
    cmd_joint.add_argument("-g", "--beagle")
        .help("gziped beagle format as input")
        .default_value(std::string{""});
    cmd_joint.add_argument("-i", "--iterations")
        .help("number of EM iterations")
        .default_value(1000)
        .scan<'i', int>();
    cmd_joint.add_argument("-n", "--threads")
        .help("number of threads")
        .default_value(1)
        .scan<'i', int>();
    cmd_joint.add_argument("-o", "--out")
        .help("output prefix")
        .default_value(std::string{"joint"});
    cmd_joint.add_argument("-s", "--chunksize")
        .help("size of each chunk in sites unit ")
        .default_value(10000)
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
        .default_value(999)
        .scan<'i', int>();
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
        .help("number of SNPs (>1) in each grid. 1 disables collapsing")
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
        .default_value(999)
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
        .default_value(0)
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
        .default_value(1000)
        .scan<'i', int>();
    cmd_admix.add_argument("-n", "--threads")
        .help("number of threads")
        .default_value(4)
        .scan<'i', int>();
    cmd_admix.add_argument("-o", "--out")
        .help("output prefix")
        .default_value(std::string{"admix"});
    cmd_admix.add_argument("-d","--seed")
        .help("seed for reproducibility")
        .default_value(999)
        .scan<'i', int>();
    // cmd_admix.add_parents(program);

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
        opts.nQ = program.get<bool>("--NQ");
        opts.nP = program.get<bool>("--NP");
        opts.nR = program.get<bool>("--NR");
        opts.nF = program.get<bool>("--NF");
        opts.ltol = program.get<double>("--ltol");
        opts.noaccel = program.get<bool>("--no-accel");

        if(program.is_subcommand_used(cmd_joint))
        {
            opts.in_beagle.assign(cmd_joint.get("--beagle"));
            opts.out.assign(cmd_joint.get("--out"));
            opts.C = cmd_joint.get<int>("--cluster");
            opts.K = cmd_joint.get<int>("--ancestry");
            opts.nthreads = cmd_joint.get<int>("--threads");
            opts.nimpute = cmd_joint.get<int>("--iterations");
            opts.seed = cmd_joint.get<int>("--seed");
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
            opts.in_bin.assign(cmd_admix.get("--bin"));
            opts.out.assign(cmd_admix.get("--out"));
            opts.seed = cmd_admix.get<int>("--seed");
            opts.K = cmd_admix.get<int>("-k");
            opts.nthreads = cmd_admix.get<int>("--threads");
            opts.nadmix = cmd_admix.get<int>("--iterations");
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
