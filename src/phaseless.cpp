/*******************************************************************************
 * @file        https://github.com/Zilong-Li/phaseless/src/phaseless.cpp
 * @author      Zilong Li
 * Copyright (C) 2023. The use of this code is governed by the LICENSE file.
 ******************************************************************************/
#define _DECLARE_TOOLBOX_HERE

#include "admix.cpp"
#include "convert.cpp"
#include "impute.cpp"
#include "parse.cpp"
#include <argparse/argparse.hpp>

int main(int argc, char * argv[])
{
    // ========= helper message and parameters parsing ===========================

    const std::string VERSION{"0.2.1"};

    // clang-format off
    argparse::ArgumentParser program("phaseless", VERSION);
    program.add_argument("--debug")
        .help("enable debug mode")
        .default_value(false)
        .implicit_value(true);


    argparse::ArgumentParser cmd_impute("impute", VERSION, argparse::default_arguments::help);
    cmd_impute.add_parents(program);
    cmd_impute.add_description("run imputation for low coverage sequencing data");
    cmd_impute.add_argument("-c", "--cluster")
        .help("number of ancestral haplotype clusters")
        .default_value(10)
        .scan<'i', int>();
    cmd_impute.add_argument("-f", "--vcf")
        .help("vcf/bcf format with GL/PL tag as input")
        .default_value(std::string{""});
    cmd_impute.add_argument("-g", "--beagle")
        .help("gziped beagle format as input")
        .default_value(std::string{""});
    cmd_impute.add_argument("-i", "--iterations")
        .help("number of EM iterations")
        .default_value(40)
        .scan<'i', int>();
    cmd_impute.add_argument("-n", "--threads")
        .help("number of threads")
        .default_value(1)
        .scan<'i', int>();
    cmd_impute.add_argument("-o", "--out")
        .help("output prefix")
        .default_value(std::string{"impute"});
    cmd_impute.add_argument("-p", "--no-print")
        .help("disable print log to screen")
        .default_value(false)
        .implicit_value(true);
    cmd_impute.add_argument("-r", "--region")
        .help("region in vcf/bcf to subset")
        .default_value(std::string{""});
    cmd_impute.add_argument("-s", "--chunksize")
        .help("size of each chunk in sites unit ")
        .default_value(10000)
        .scan<'i', int>();
    cmd_impute.add_argument("-S", "--single-chunk")
        .help("treat input as big single chunk")
        .default_value(false)
        .implicit_value(true);
    cmd_impute.add_argument("--seed")
        .help("seed for reproducing results")
        .default_value(999)
        .scan<'i', int>();

    argparse::ArgumentParser cmd_admix("admix", VERSION, argparse::default_arguments::help);
    cmd_admix.add_parents(program);
    cmd_admix.add_description("run admixture with cluster likelihoods as input");
    cmd_admix.add_argument("-a", "--no-accel")
        .help("disable accelerated EM")
        .default_value(false)
        .implicit_value(true);
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
    cmd_admix.add_argument("-p", "--no-print")
        .help("disable print log to screen")
        .default_value(false)
        .implicit_value(true);
    cmd_admix.add_argument("-q", "--qtol")
        .help("tolerance of stopping criteria for diff(Q)")
        .default_value(1e-6)
        .scan<'g', double>();
    cmd_admix.add_argument("-l", "--ltol")
        .help("tolerance of stopping criteria for diff(loglikelihood)")
        .default_value(1e-1)
        .scan<'g', double>();
    cmd_admix.add_argument("-s","--seed")
        .help("seed for reproducing results")
        .default_value(999)
        .scan<'i', int>();

    argparse::ArgumentParser cmd_convert("convert", VERSION, argparse::default_arguments::help);
    cmd_convert.add_parents(program);
    cmd_convert.add_description("different file format converter");
    cmd_convert.add_argument("input")
        .help("input file to be converted");
    cmd_convert.add_argument("output")
        .help("output file prefix")
        .default_value(std::string{"convert"});
    cmd_convert.add_argument("-p", "--plink2beagle")
        .help("use plink1 file as input without .bed")
        .default_value(false)
        .implicit_value(true);
    cmd_convert.add_argument("-n", "--threads")
        .help("number of threads")
        .default_value(4)
        .scan<'i', int>();
    cmd_convert.add_argument("-s", "--chunksize")
        .help("size of each chunk in sites unit ")
        .default_value(10000)
        .scan<'i', int>();

    argparse::ArgumentParser cmd_parse("parse", VERSION, argparse::default_arguments::help);
    cmd_parse.add_parents(program);
    cmd_parse.add_description("manipulate pars.bin file outputted by impute command");
    cmd_parse.add_argument("-b", "--bin")
        .help("binary format from impute command as input")
        .default_value(std::string{""});
    cmd_parse.add_argument("-c", "--chunk")
        .help("which chunk to extract (0-based), negative means all chunks")
        .default_value(0)
        .scan<'i', int>();
    cmd_parse.add_argument("-s", "--samples")
        .help("extract samples in the file, one sample id per line")
        .default_value(std::string{""});
    cmd_parse.add_argument("-o", "--out")
        .help("output prefix")
        .default_value(std::string{"parse"});

    program.add_subparser(cmd_impute);
    program.add_subparser(cmd_admix);
    program.add_subparser(cmd_parse);
    program.add_subparser(cmd_convert);
    // clang-format on

    try
    {
        Options opts;
        for(int i = 0; i < argc; i++) opts.opts_in_effect += " " + std::string{argv[i]};
        program.parse_args(argc, argv);

        if(program.is_subcommand_used(cmd_impute))
        {
            opts.in_beagle.assign(cmd_impute.get("--beagle"));
            opts.in_vcf.assign(cmd_impute.get("--vcf"));
            opts.out.assign(cmd_impute.get("--out"));
            opts.C = cmd_impute.get<int>("--cluster");
            opts.nthreads = cmd_impute.get<int>("--threads");
            opts.nimpute = cmd_impute.get<int>("--iterations");
            opts.seed = cmd_impute.get<int>("--seed");
            opts.chunksize = cmd_impute.get<int>("--chunksize");
            opts.single_chunk = cmd_impute.get<bool>("--single-chunk");
            opts.noscreen = cmd_impute.get<bool>("--no-print");
            opts.debug = cmd_impute.get<bool>("--debug");
            if(opts.single_chunk) opts.chunksize = INT_MAX;
            if(opts.in_beagle.empty() && opts.in_vcf.empty()) cao.error(cmd_impute.help().str());
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
            opts.qtol = cmd_admix.get<double>("--qtol");
            opts.ltol = cmd_admix.get<double>("--ltol");
            opts.noaccel = cmd_admix.get<bool>("--no-accel");
            opts.noscreen = cmd_admix.get<bool>("--no-print");
            opts.debug = cmd_admix.get<bool>("--debug");
            if(opts.in_bin.empty()) cao.error(cmd_admix.help().str());
            run_admix_main(opts);
        }
        else if(program.is_subcommand_used(cmd_parse))
        {
            opts.in_bin.assign(cmd_parse.get("--bin"));
            opts.out.assign(cmd_parse.get("--out"));
            opts.samples = cmd_parse.get("--samples");
            opts.ichunk = cmd_parse.get<int>("--chunk");
            if(opts.in_bin.empty()) cao.error(cmd_parse.help().str());
            run_parse_main(opts);
        }
        else if(program.is_subcommand_used(cmd_convert))
        {
            opts.nthreads = cmd_convert.get<int>("--threads");
            opts.in_plink = cmd_convert.get("input");
            opts.out = cmd_convert.get("output");
            opts.chunksize = cmd_convert.get<int>("--chunksize");
            if(cmd_convert.get<bool>("--plink2beagle"))
                run_convert_main(opts);
            else
                cao.error(cmd_convert.help().str());
        }
        else
        {
            std::cerr << program.help().str();
            std::exit(1);
        }
    }
    catch(const std::runtime_error & err)
    {
        std::cerr << err.what() << std::endl;
        std::cerr << program;
        std::exit(1);
    }

    return 0;
}
