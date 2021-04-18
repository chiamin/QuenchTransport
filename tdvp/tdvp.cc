#include <csignal>
#include "itensor/all.h"
#include "Timer.h"
Timers timer;
#include "ReadInput.h"
#include "IUtility.h"
#include "MyObserver.h"
#include "tdvp.h"
#include "MPSUtility.h"
#include "Entanglement.h"
#include "GeneralUtility.h"
#include "MixedBasis.h"
#include "SubCorr.h"
#include "MeaCurrent.h"
using namespace itensor;
using namespace std;

vector<Index> get_site_inds (const MPS& psi)
{
    int N = length(psi);
    vector<Index> sites (N);
    for(int i = 1; i <= N; i++)
    {
        auto ii = findIndex (psi(i), "Site");
        sites.at(i-1) = ii;
    }
    return sites;
}

void writeAll (const string& filename,
               const MPS& psi, const MPO& H,
               int step)
{
    ofstream ofs (filename);
    write (ofs, psi);
    write (ofs, H);
    write (ofs, step);
}

void readAll (const string& filename,
              MPS& psi, MPO& H,
              int& step)
{
    ifstream ifs = open_file (filename);
    read (ifs, psi);
    read (ifs, H);
    read (ifs, step);
}

int main(int argc, char* argv[])
{
timer.start();
    string infile = argv[1];
    InputGroup input (infile,"basic");
    auto mu_leadL   = input.getReal("mu_leadL");
    auto mu_leadR   = input.getReal("mu_leadR");
    auto mu_device  = input.getReal("mu_device");

    auto SubCorrN   = input.getInt("SubCorrN",-1);
    auto corr_cutoff = input.getReal("corr_cutoff");

    auto input_dir       = input.getString("input_dir");
    auto psi_file       = input.getString("psi_file");
    auto H_file         = input.getString("H_file");

    auto dt            = input.getReal("dt");
    auto time_steps    = input.getInt("time_steps");
    auto NumCenter     = input.getInt("NumCenter");
    auto ConserveQNs   = input.getYesNo("ConserveQNs",false);
    auto sweeps        = Read_sweeps (infile);

    auto UseSVD        = input.getYesNo("UseSVD",true);
    auto SVDmethod     = input.getString("SVDMethod","gesdd");  // can be also "ITensor"
    auto WriteDim      = input.getInt("WriteDim");

    auto write         = input.getYesNo("write",false);
    auto write_dir     = input.getString("write_dir",".");
    auto write_file    = input.getString("write_file","");
    auto read          = input.getYesNo("read",false);
    auto read_dir      = input.getString("read_dir",".");
    auto read_file     = input.getString("read_file","");

    auto out_dir     = input.getString("outdir",".");
    if (write_dir == "." && out_dir != ".")
        write_dir = out_dir;

    auto system = make_shared <WireSystem> ();
    system->read (input_dir+"/"+H_file);

    // Declare variables
    MPS psi;
    MPO H;
    Fermion sites;
    int step = 1;

    if (!read)
    {
        // Read MPS in the window
        readFromFile (input_dir+"/"+psi_file, psi);
        int N = length (psi);
        // Site indices
        auto site_inds = get_site_inds (psi);
        sites = Fermion (site_inds);
        // Hamiltonian
        system->add_mu (SegL, mu_leadL);
        system->add_mu (SegR, mu_leadR);
        system->add_mu (SegS, mu_device);
        auto ampo = system->ampo (sites);
        H = toMPO (ampo);
    }
    else
    {
        readAll (read_dir+"/"+read_file, psi, H, step);
        auto site_inds = get_site_inds (psi);
        sites = Fermion (site_inds);
    }
    cout << "H MPO dim = " << maxLinkDim(H) << endl;
    cout << "device site = " << system->idevL() << " " << system->idevR() << endl;

    // Args parameters
    if (SubCorrN == -1) SubCorrN = length(sites);
    Args args_obs   = {"ConserveQNs",ConserveQNs};
    Args args_tdvp  = {"Quiet",true,"NumCenter",NumCenter,"DoNormalize",true,
                       "UseSVD",UseSVD,"SVDmethod",SVDmethod,"WriteDim",WriteDim};

    // Observer
    auto obs = MyObserver (sites, psi, args_obs);
    // Current
    int N = length (psi);
    // MPO
    vector<int> spec_links = {system->idevL()-1, system->idevR()};
    vector<MPO> JMPOs (N);
    for(int i : spec_links)
    {
        auto mpo = get_current_mpo2 (sites, *system, i);
        JMPOs.at(i) = mpo;
    }
    // Correlation
    int l = (N - SubCorrN) / 2;
    int ibeg = l+1,
        iend = l+SubCorrN;
    auto sub_corr = SubCorr (system, ibeg, iend);


    // Time evolution
    cout << "Start time evolution" << endl;
    cout << sweeps << endl;
    psi.position(1);
    Real en, err;
    for(int i = 0; i < time_steps; i++)
    {
        cout << "step = " << step << endl;

        timer["tdvp"].start();
        tdvp (psi, H, 1_i*dt, sweeps, obs, args_tdvp);
        timer["tdvp"].stop();

        // Measure currents by correlations

        timer["current corr"].start();
        sub_corr.measure (psi, {"Cutoff",corr_cutoff});
        timer["current corr"].stop();
        for(int j = 1; j < N; j++)
        {
            timer["current corr"].start();
            auto J = sub_corr.get_current (j);
            timer["current corr"].stop();
            cout << "\t*current " << j << " " << j+1 << " = " << J << endl;
        }
        for(int j : spec_links)
        {
            timer["current mps"].start();
            auto J = get_current (JMPOs.at(j), psi);
            timer["current mps"].stop();
            cout << "\t*current spec " << j << " " << j+1 << " = " << J << endl;
        }

        if (write)
        {            
            writeAll (write_dir+"/"+write_file, psi, H, step);
        }
        step++;
    }
        timer.print();

    return 0;
}
