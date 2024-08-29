/*
 * An application to evaluate the Deep Learning energy estimator,
 * that is using Particle Flow information.
 */

#include "config.h"
#if (HAVE_VLNEVAL_INC == 1) && (HAVE_VLNEVAL_LIB == 1)

#include <cmath>
#include <chrono>
#include <iostream>
#include <string>
#include <stdexcept>
#include <unordered_map>

#include "TFile.h"
#include "TTree.h"

#include <boost/format.hpp>
#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include <vlneval/struct/VarDict.h>
#include <vlneval/zoo/VLNEnergyModel.h>

/*
 * TODO: 30k is a value extracted from
 *      https://github.com/HaiwangYu/wcp-mc-eval/blob/main/muon_energy.C
 *      need to double check that this value is correct
 */
constexpr size_t MAX_TRACKS = 30000;
constexpr double MASK_VALUE = 0.;       // Value to replace NaNs with

struct Config
{
    std::string input;
    std::string branch;
    std::string model;
    bool isNue;
    bool isNumu;
    bool fullContain;
};

struct Trees
{
    TTree *bdt;
    TTree *eval;
    TTree *pf;
    TTree *ke;
};

struct RootVars
{
    int run;
    int subrun;
    int event;

    float numu_score;
    float nue_score;
    float numu_cc_flag;

    int   reco_Ntrack;
    int   reco_pdg[MAX_TRACKS];
    float reco_startXYZT[MAX_TRACKS][4];
    float reco_endXYZT[MAX_TRACKS][4];
    float reco_startMomentum[MAX_TRACKS][4];

    void setup(Trees &trees);
};

class ProgressBar
{
private:
    size_t total;
    size_t current;
    size_t redraw;

    float  momentum;
    float  speed;

    std::chrono::time_point<std::chrono::system_clock> start_time;

    void update_speed(int time_taken)
    {
        if ((current == 0) || (time_taken == 0)) {
            return;
        }

        const float curr_speed = static_cast<float>(current) / time_taken;
        speed = (1 - momentum) * curr_speed + momentum * speed;
    }

public:
    ProgressBar(size_t total, size_t redraw = 100, float momentum = 0.999)
      : total(total), current(0), redraw(redraw), momentum(momentum),
        speed(0), start_time(std::chrono::system_clock::now())
    {
        update(0);
    }

    ~ProgressBar()
    {
        std::cerr << std::endl;
    }

    static void pprintDuration(int seconds)
    {
        if (seconds < 0) {
            return;
        }

        static boost::format fmt_short("%1$02d:%2$02d");
        static boost::format fmt_long ("%1$02d:%2$02d:%3$02d");

        int minutes = seconds / 60;
        seconds     = seconds % 60;

        int hours = minutes / 60;
        minutes   = minutes % 60;

        if (hours == 0) {
            std::cerr << fmt_short % minutes % seconds;
        }
        else {
            std::cerr << fmt_long % hours % minutes % seconds;
        }
    }

    int get_elapsed_time_in_seconds()
    {
        const auto curr_time = std::chrono::system_clock::now();

        return std::chrono::duration_cast<std::chrono::seconds>(
            curr_time - start_time
        ).count();
    }

    int estimate_remaining_time_in_seconds(int time_taken)
    {
        if (speed <= 0) {
            return -1;
        }

        return static_cast<float>(total - current) / speed;
    }

    void update(size_t step = 1)
    {
        current += step;

        const int time_taken = get_elapsed_time_in_seconds();
        update_speed(time_taken);

        if ((current % redraw != 0) && (current != total)) {
            return;
        }

        const int time_remain = estimate_remaining_time_in_seconds(time_taken);

        std::cerr << boost::format("Progress: %1% / %2% [%3$.1f %%]")
            % current % total % (100. * current / total);;

        std::cerr << ". Time: ";
        ProgressBar::pprintDuration(time_taken);

        std::cerr << ". ETA: ";
        ProgressBar::pprintDuration(time_remain);

        std::cerr << "    " << '\r';
    }
};

void RootVars::setup(Trees &trees)
{
    trees.eval->SetBranchAddress("run",    &run);
    trees.eval->SetBranchAddress("subrun", &subrun);
    trees.eval->SetBranchAddress("event",  &event);

    trees.bdt->SetBranchAddress("numu_score",   &numu_score);
    trees.bdt->SetBranchAddress("nue_score",    &nue_score);
    trees.bdt->SetBranchAddress("numu_cc_flag", &numu_cc_flag);

    trees.pf->SetBranchAddress("reco_Ntrack",        &reco_Ntrack);
    trees.pf->SetBranchAddress("reco_pdg",           &reco_pdg);
    trees.pf->SetBranchAddress("reco_startXYZT",     &reco_startXYZT);
    trees.pf->SetBranchAddress("reco_endXYZT",       &reco_endXYZT);
    trees.pf->SetBranchAddress("reco_startMomentum", &reco_startMomentum);
}

Trees getTrees(TFile &file)
{
    Trees result;

    result.pf   = static_cast<TTree*>(file.Get("wcpselection/T_PFeval"));
    result.eval = static_cast<TTree*>(file.Get("wcpselection/T_eval"));
    result.bdt  = static_cast<TTree*>(file.Get("wcpselection/T_BDTvars"));
    result.ke   = static_cast<TTree*>(file.Get("wcpselection/T_KINEvars"));

    result.pf->AddFriend(result.eval);
    result.pf->AddFriend(result.bdt);
    result.pf->AddFriend(result.ke);

    return result;
}

void extractAddressVars(const RootVars &rootVars, VarDict &vars)
{
    vars.scalar["addr.run"]    = rootVars.run;
    vars.scalar["addr.subrun"] = rootVars.subrun;
    vars.scalar["addr.event"]  = rootVars.event;
}

void extractBDTVars(const RootVars &rootVars, VarDict &vars)
{
    vars.scalar["event.numu_score"]   = rootVars.numu_score;
    vars.scalar["event.nue_score"]    = rootVars.nue_score;
    vars.scalar["event.numu_cc_flag"] = rootVars.numu_cc_flag;
}

void initVectorVar(const std::string &var, size_t size, VarDict &vars)
{
    auto it = vars.vector.find(var);

    if (it != vars.vector.end()) {
        it->second.resize(size);
        std::fill(it->second.begin(), it->second.end(), 0);
    }
    else {
        vars.vector[var] = std::vector<double>(size, 0);
    }
}

void extractPFPDGVars(const RootVars &rootVars, VarDict &vars)
{
    const std::string PDG_PREFIX = "particle.pdg.";
    const std::unordered_map<int, std::string> PDG_MAP({
        { 11,   PDG_PREFIX + "electron" },
        { 13,   PDG_PREFIX + "muon"     },
        { 22,   PDG_PREFIX + "gamma"    },
        { 2212, PDG_PREFIX + "proton"   },
        { 2112, PDG_PREFIX + "neutron"  },
        { 211,  PDG_PREFIX + "pion"     },
        { 111,  PDG_PREFIX + "pizero"   },
        { 0,    PDG_PREFIX + "other"    }
    });

    const size_t nTracks = rootVars.reco_Ntrack;

    /* NOTE: this is very inefficient. Maybe optimize */
    for (const auto &item : PDG_MAP) {
        initVectorVar(item.second, nTracks, vars);
    }

    for (size_t idx = 0; idx < nTracks; idx++) {
        int pdg = std::abs(rootVars.reco_pdg[idx]);
        auto it = PDG_MAP.find(pdg);

        if (it != PDG_MAP.end()) {
            vars.vector[it->second][idx] = 1;
        }
        else {
            vars.vector[PDG_MAP.at(0)][idx] = 1;
        }
    }
}

void unpackXYZTArray(
    const float array[MAX_TRACKS][4],
    size_t n,
    const std::string &prefix,
    VarDict &vars
)
{
    const std::vector<std::string> COORDS({
        prefix + "x", prefix + "y", prefix + "z", prefix + "t"
    });

    for (const auto &var : COORDS) {
        initVectorVar(var, n, vars);
    }

    for (size_t idx = 0; idx < n; idx++) {
        for (size_t coord_idx = 0; coord_idx < COORDS.size(); coord_idx++) {
            vars.vector[COORDS[coord_idx]][idx] = array[idx][coord_idx];
        }
    }
}

void extractPFVars(const RootVars &rootVars, VarDict &vars)
{
    size_t nTracks = rootVars.reco_Ntrack;

    extractPFPDGVars(rootVars, vars);
    unpackXYZTArray(
        rootVars.reco_startXYZT, nTracks, "particle.start.", vars
    );
    unpackXYZTArray(
        rootVars.reco_endXYZT, nTracks, "particle.end.", vars
    );
    unpackXYZTArray(
        rootVars.reco_startMomentum, nTracks, "particle.startMomentum.", vars
    );
}

void extractVars(const RootVars &rootVars, VarDict &vars)
{
    extractAddressVars(rootVars, vars);
    extractBDTVars(rootVars, vars);
    extractPFVars(rootVars, vars);
}

void maskNans(VarDict &vars, double mask)
{
    for (auto &item : vars.scalar) {
        if (! std::isfinite(item.second)) {
            item.second = mask;
        }
    }

    for (auto &item : vars.vector) {
        for (auto &value : item.second) {
            if (! std::isfinite(value)) {
                value = mask;
            }
        }
    }

}

void usage(const char *name, const po::options_description &options)
{
    std::cout << "USAGE: " << name << " [OPTION..] INPUT" << std::endl;
    std::cout << options << std::endl;
}

void parseVMFlavor(Config &config, const po::variables_map &vm)
{
    const std::string flavor = vm["flavor"].as<std::string>();

    config.isNue  = false;
    config.isNumu = false;

    if (flavor == "numu") {
        config.isNue  = false;
        config.isNumu = true;
    }
    else if (flavor == "nue") {
        config.isNue  = true;
        config.isNumu = false;
    }
    else {
        throw std::runtime_error("Unknown flavor: '" + flavor + "'");
    }
}

void parseVMContainment(Config &config, const po::variables_map &vm)
{
    const std::string contain = vm["contain"].as<std::string>();

    config.fullContain = true;

    if (contain == "full") {
        config.fullContain = true;
    }
    else if (contain == "partial") {
        config.fullContain = false;
    }
    else {
        throw std::runtime_error("Unknown containment: '" + contain + "'");
    }
}

Config parseVM(po::variables_map vm)
{
    Config result;

    result.input  = vm["input"].as<std::string>();
    result.branch = vm["branch"].as<std::string>();
    result.model  = vm["model"].as<std::string>();

    parseVMFlavor(result, vm);
    parseVMContainment(result, vm);

    return result;
}

Config parseArgs(int argc, char** argv)
{
    po::options_description options("Allowed options");
    options.add_options()
        ("help,h", "print this help message")
        (
            "branch,b",
            po::value<std::string>()->default_value("vlne"),
            "branch to save energy"
        )
        (
            "contain",
            po::value<std::string>()->default_value("full"),
            "event containment [full|partial]"
        )
        (
            "flavor",
            po::value<std::string>()->default_value("numu"),
            "flavor of inputs [numu|nue]"
        )
        ("input,i",  po::value<std::string>(), "input file")
        ("model,m",  po::value<std::string>()->required(), "model directory")
        ;

    po::positional_options_description positional_args;
    positional_args.add("input",  1);

    po::variables_map vm;
    po::store(
        po::command_line_parser(argc, argv)
            .options(options)
            .positional(positional_args)
            .run(),
        vm
    );

    if (vm.count("help")) {
        usage(argv[0], options);
        exit(0);
    }

    po::notify(vm);

    if (vm.count("input") == 0) {
        usage(argv[0], options);
        std::cerr << "ERROR: Must provide INPUT file" << std::endl;
        exit(1);
    }

    return parseVM(vm);
}

VarDict setupVarDict(const Config &config)
{
    VarDict result;

    result.scalar["config.flavor.numu"] = config.isNumu;
    result.scalar["config.flavor.nue"]  = config.isNue;
    result.scalar["config.contain.full"] = config.fullContain;

    return result;
}

TBranch* createEnergyBranch(
    Trees &trees, const Config &config, const std::string &label, float *value
)
{
    const std::string name = config.branch + "_" + label;
    const std::string leaf = name + "/F";

    TBranch *branch = trees.ke->GetBranch(name.c_str());
    if (branch != nullptr) {
        throw std::runtime_error(
            "DL EE branch '" + name + "' already exists."
            " Refusing to override."
        );
    }

    return trees.ke->Branch(name.c_str(), value, leaf.c_str());
}

void processFile(TFile &file, const Config &config, VLN::VLNEnergyModel &model)
{
    Trees trees = getTrees(file);

    RootVars rootVars;
    rootVars.setup(trees);

    VLN::VLNEnergy energy;

    TBranch *primaryEBranch = createEnergyBranch(
        trees, config, "primaryE", &energy.primaryE
    );
    TBranch *totalEBranch = createEnergyBranch(
        trees, config, "totalE", &energy.totalE
    );

    VarDict vars = setupVarDict(config);
    ProgressBar pbar(trees.pf->GetEntries());

    for (int idx = 0; idx < trees.pf->GetEntries(); idx++) {
        trees.pf->GetEntry(idx);

        extractVars(rootVars, vars);
        maskNans(vars, MASK_VALUE);

        energy = model.predict(vars);

        primaryEBranch->Fill();
        totalEBranch->Fill();

        pbar.update();
    }

    trees.ke->Write("wcpselection/T_KINEvars", TObject::kOverwrite);
}

int main(int argc, char** argv)
{
    const Config config = parseArgs(argc, argv);
    TFile file(config.input.c_str(), "update");

    VLN::VLNEnergyModel model(config.model);

    std::cout
        << "Starting energy evaluation"             << std::endl
        << "  - File: "         << config.input     << std::endl
        << "  - Model: "        << config.model     << std::endl
        << "  - Save Branch: "  << config.branch    << std::endl;

    processFile(file, config, model);
    file.Write();

    std::cout << "Finished energy evaluation" << std::endl;
}

#else  /* HAVE_VLNEVAL != 1 */

#include <iostream>
int main(int argc, char** argv)
{ std::cout << "Stub Application" << std::endl; }

#endif /* HAVE_VLNEVAL */

