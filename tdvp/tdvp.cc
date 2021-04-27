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
#include "ReadWriteFile.h"
#include "ContainerUtility.h"
#include "basisextension.h"
using namespace vectool;
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
    itensor::write (ofs, psi);
    itensor::write (ofs, H);
    itensor::write (ofs, step);
}

void readAll (const string& filename,
              MPS& psi, MPO& H,
              int& step)
{
    ifstream ifs = open_file (filename);
    itensor::read (ifs, psi);
    itensor::read (ifs, H);
    itensor::read (ifs, step);
}

vector<Real> order_by_bias (const WireSystem& system, Real bias, Real threshold=0.01)
{
    // Set order
    Real en0 = 0.5 * abs(bias);
    vector<Real> order;
    for(int i = 0; i < system.orbs().size(); i++)
    {
        auto [seg, ind, en] = system.orbs().at(i);
        if (seg == "S")
        {
            order.push_back (0.);
        }
        else if (abs(en) < threshold)
        {
            order.push_back (en);
        }
        else if (en < 0.)
        {
            Real d = abs (en + en0);
            order.push_back (-d);
        }
        else
        {
            Real d = abs (en - en0);
            order.push_back (d);
        }
    }
    return order;
}

int main(int argc, char* argv[])
{
    timer.start();
    string infile = argv[1];
    InputGroup input (infile,"basic");
    auto bias_leadL   = input.getReal("bias_leadL");
    auto bias_leadR   = input.getReal("bias_leadR");
    auto bias_device  = input.getReal("bias_device");

    auto SubCorrN   = input.getInt("SubCorrN",-1);
    auto corr_cutoff = input.getReal("corr_cutoff");

    auto input_dir       = input.getString("input_dir");
    auto psi_file       = input.getString("psi_file");
    auto H_file         = input.getString("H_file");

    auto dt            = input.getReal("dt");
    auto time_steps    = input.getInt("time_steps");
    auto NumCenter     = input.getInt("NumCenter");
    auto Truncate      = input.getYesNo("Truncate");
    auto mixNumCenter  = input.getYesNo("mixNumCenter",false);
    auto sweeps        = Read_sweeps (infile);

    auto globExpanN         = input.getInt("globExpanN",std::numeric_limits<int>::max());
    auto globExpanCutoff    = input.getReal("globExpanCutoff",1e-8);
    auto globExpanDim       = input.getReal("globExpanDim",10);
    auto globExpanMethod    = input.getString("globExpanMethod","DensityMatrix");  // Can be also "Fit"
    auto globExpanKrylovDim = input.getInt("globExpanKrylovDim",3);
    auto globExpanHpsiCutoff = input.getReal("globExpanHpsiCutoff",1e-8);
    auto globExpanHpsiDim   = input.getInt("globExpanHpsiDim",10);

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

    auto reorder        = input.getYesNo("reorder",false);
    auto cutoff_reorder = input.getReal("cutoff_reorder",1e-12);

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
        // Reorder basis
        if (reorder)
        {
            Real bias = bias_leadR - bias_leadL;
            auto new_order = order_by_bias (*system, bias);
            reorder_basis (psi, *system, new_order, {"Cutoff",cutoff_reorder,"verbose",true});
            cout << "dim after reorder = " << maxLinkDim (psi) << endl;
        }
        // Site indices
        auto site_inds = get_site_inds (psi);
        sites = Fermion (site_inds);
        // Hamiltonian
        system->add_mu ("L", bias_leadL);
        system->add_mu ("R", bias_leadR);
        system->add_mu ("S", bias_device);
        auto ampo = system->ampo (sites);
        H = toMPO (ampo);
    }
    else
    {
        readAll (read_dir+"/"+read_file, psi, H, step);
        auto site_inds = get_site_inds (psi);
        sites = Fermion (site_inds);
    }
    system->print_orbs();
    cout << "H MPO dim = " << maxLinkDim(H) << endl;
    int idevL = system->idevL(),
        idevR = system->idevR();
    cout << "device site = " << idevL << " " << idevR << endl;

    // Args parameters
    Args args_tdvp  = {"Quiet",true,"NumCenter",NumCenter,"DoNormalize",true,"Truncate",Truncate,
                       "UseSVD",UseSVD,"SVDmethod",SVDmethod,"WriteDim",WriteDim,"mixNumCenter",mixNumCenter};

    // Observer
    auto obs = MyObserver (sites, psi);
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
    if (SubCorrN < 0 or SubCorrN > length(psi))
    {
        SubCorrN = length(psi);
    }
    else if (SubCorrN < idevR-idevL+1)
    {
        cout << "sub-correlation size is too small: " << SubCorrN << endl;
        cout << "device size = " << idevL-idevR+1 << endl;
        throw;
    }
    int l = (N - SubCorrN) / 2;
    int ibeg = l+1,
        iend = l+SubCorrN;
    auto sub_corr = SubCorr (system, ibeg, iend);

    // Time evolution
    cout << "Start time evolution" << endl;
    cout << sweeps << endl;
    psi.position(1);
    Real en, err;
    Args args_tdvp_expansion = {"Cutoff",globExpanCutoff,"MaxDim",globExpanDim,"Method",globExpanMethod,
                                "KrylovOrd",globExpanKrylovDim, "DoNormalize",true, "Quiet",true};
    for(int i = 0; i < time_steps; i++)
    {
        cout << "step = " << step << endl;

        // Global space expansion
        if (globExpanCutoff != 0. and i < globExpanN)
        {
            timer["glob expan"].start();
            // Global subspace expansion
            std::vector<Real> cutoffK (globExpanKrylovDim-1, globExpanHpsiCutoff);
            std::vector<int>  dimK (globExpanKrylovDim-1, globExpanHpsiDim);
            addBasis (psi, H, cutoffK, dimK, args_tdvp_expansion);
            timer["glob expan"].stop();
        }
        // Time evolution
        timer["tdvp"].start();
        tdvp (psi, H, (0.1+1_i)*dt, sweeps, obs, args_tdvp);
        timer["tdvp"].stop();
        auto d1 = maxLinkDim(psi);

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
        // Measure currents by MPO
        for(int j : spec_links)
        {
            timer["current mps"].start();
            auto J = get_current (JMPOs.at(j), psi);
            timer["current mps"].stop();
            cout << "\t*current spec " << j << " " << j+1 << " = " << J << endl;
        }

        if (write)
        {
            timer["write"].start();
            writeAll (write_dir+"/"+write_file, psi, H, step);
            timer["write"].stop();
        }
        step++;
    }
    timer.print();

    return 0;
}
