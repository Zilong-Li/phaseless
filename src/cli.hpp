#ifndef CLI_H_
#define CLI_H_

#include <argparse/argparse.hpp>
#include <filesystem>

struct Options
{
    int chunksize{10000}, K{2}, C{10}, nadmix{1000}, nphase{40}, nthreads{1}, seed{999};
    double qtol{1e-6}, info{0};
    bool noaccel{0}, noscreen{0}, run_impute{0}, run_admix{0}, run_pars{0},single_chunk{0};
    std::filesystem::path out, in_beagle, in_vcf, in_bin;
    std::string samples{"-"}, region{""};
    std::string opts_in_effect{"Options in effect:\n   "};
};

inline auto parsecli(int argc, char * argv[])
{
    const std::string VERSION{"0.1.9"};

    // clang-format off
    argparse::ArgumentParser program("phaseless", VERSION, argparse::default_arguments::help);

    argparse::ArgumentParser impute_command("impute", VERSION, argparse::default_arguments::help);
    impute_command.add_description("run imputation for low coverage sequencing data");
    impute_command.add_argument("-c", "--cluster")
        .help("number of ancestral haplotype clusters")
        .default_value(10)
        .scan<'i', int>();
    impute_command.add_argument("-f", "--vcf")
        .help("vcf/bcf format with GL/PL tag as input")
        .default_value(std::string{""});
    impute_command.add_argument("-g", "--beagle")
        .help("gziped beagle format as input")
        .default_value(std::string{""});
    impute_command.add_argument("-n", "--threads")
        .help("number of threads")
        .default_value(1)
        .scan<'i', int>();
    impute_command.add_argument("-o", "--out")
        .help("output prefix")
        .default_value(std::string{"impute."});
    impute_command.add_argument("-p", "--no-print")
        .help("disable print log to screen")
        .default_value(false)
        .implicit_value(true);
    impute_command.add_argument("-r", "--region")
        .help("region in vcf/bcf to subset")
        .default_value(std::string{""});
    impute_command.add_argument("-s", "--chunksize")
        .help("size of each chunk in sites unit ")
        .default_value(10000)
        .scan<'i', int>();
    impute_command.add_argument("-S", "--single-chunk")
        .help("treat input as big single chunk")
        .default_value(false)
        .implicit_value(true);
    impute_command.add_argument("--info")
        .help("filter and re-impute sites with low info")
        .default_value(0.0)
        .scan<'g', double>();
    impute_command.add_argument("--seed")
        .help("seed for reproducing results and different inits")
        .default_value(999)
        .scan<'i', int>();

    argparse::ArgumentParser admix_command("admix", VERSION, argparse::default_arguments::help);
    admix_command.add_description("run admixture with cluster likelihoods as input");
    admix_command.add_argument("-a", "--no-accel")
        .help("disable accelerated EM")
        .default_value(false)
        .implicit_value(true);
    admix_command.add_argument("-b", "--bin")
        .help("binary format from impute command as input")
        .default_value(std::string{""});
    admix_command.add_argument("-k", "--ancestry")
        .help("number of ancestry in admixture assumption")
        .default_value(2)
        .scan<'i', int>();
    admix_command.add_argument("-n", "--threads")
        .help("number of threads")
        .default_value(4)
        .scan<'i', int>();
    admix_command.add_argument("-o", "--out")
        .help("output prefix")
        .default_value(std::string{"admix."});
    admix_command.add_argument("-p", "--no-print")
        .help("disable print log to screen")
        .default_value(false)
        .implicit_value(true);
    admix_command.add_argument("-q", "--qtol")
        .help("tolerance of stopping criteria for diff(Q)")
        .default_value(1e-6)
        .scan<'g', double>();
    admix_command.add_argument("-s","--seed")
        .help("seed for reproducing results and different inits")
        .default_value(999)
        .scan<'i', int>();

    argparse::ArgumentParser pars_command("parse", VERSION, argparse::default_arguments::help);
    pars_command.add_description("manipulate pars.bin file outputted from impute command");
    pars_command.add_argument("-b", "--bin")
        .help("binary format from impute command as input")
        .default_value(std::string{""});
    pars_command.add_argument("-o", "--out")
        .help("output prefix")
        .default_value(std::string{"parse."});

    // clang-format on

    program.add_subparser(impute_command);
    program.add_subparser(admix_command);
    program.add_subparser(pars_command);

    try
    {
        Options opts;
        for(int i = 0; i < argc; i++) opts.opts_in_effect += " " + std::string{argv[i]};
        program.parse_args(argc, argv);
        if(program.is_subcommand_used(impute_command))
        {
            opts.run_impute = true;
            opts.in_beagle.assign(impute_command.get("--beagle"));
            opts.in_vcf.assign(impute_command.get("--vcf"));
            opts.out.assign(impute_command.get("--out"));
            opts.C = impute_command.get<int>("--cluster");
            opts.nthreads = impute_command.get<int>("--threads");
            opts.seed = impute_command.get<int>("--seed");
            opts.chunksize = impute_command.get<int>("--chunksize");
            opts.info = impute_command.get<double>("--info");
            opts.single_chunk = impute_command.get<bool>("--single-chunk");
            opts.noscreen = impute_command.get<bool>("--no-print");
            if(opts.in_beagle.empty() && opts.in_vcf.empty())
            {
                std::cerr << impute_command.help().str();
                std::exit(1);
            }
        }
        else if(program.is_subcommand_used(admix_command))
        {
            opts.run_admix = true;
            opts.in_bin.assign(admix_command.get("--bin"));
            opts.out.assign(admix_command.get("--out"));
            opts.seed = admix_command.get<int>("--seed");
            opts.K = admix_command.get<int>("-k");
            opts.nthreads = admix_command.get<int>("--threads");
            opts.qtol = admix_command.get<double>("--qtol");
            opts.noaccel = admix_command.get<bool>("--no-accel");
            opts.noscreen = admix_command.get<bool>("--no-print");
            if(opts.in_bin.empty())
            {
                std::cerr << admix_command.help().str();
                std::exit(1);
            }
        }
        else if(program.is_subcommand_used(pars_command))
        {
            opts.run_pars = true;
            opts.in_bin.assign(pars_command.get("--bin"));
            opts.out.assign(pars_command.get("--out"));
            if(opts.in_bin.empty())
            {
                std::cerr << pars_command.help().str();
                std::exit(1);
            }
        }
        else
        {
            std::cerr << program.help().str();
            std::exit(1);
        }
        return opts;
    }
    catch(const std::runtime_error & err)
    {
        std::cerr << err.what() << std::endl;
        std::cerr << program;
        std::exit(1);
    }
}

#endif // CLI_H_
