#include <iomanip>
#include "itensor/all.h"
#include "Timer.h"
Timers timer;
#include "ReadInput.h"
#include "IUtility.h"
#include "MyObserver.h"
#include "MPSUtility.h"
#include "SingleParticle.h"
#include "MixedBasis.h"
#include "SystemStruct.h"
#include "ContainerUtility.h"
#include "Corr.h"
#include "TDVPObserver.h"
#include "MeaCurrent.h"
#include "tdvp.h"
#include "basisextension.h"
using namespace vectool;
using namespace itensor;
using namespace std;
using namespace iutility;

MPO Make_NMPO (const MixedBasis& sites)
{
    AutoMPO ampo (sites);
    for(int i = 1; i <= length(sites); i++)
    {
        ampo += 1.0,"N",i;
    }
    return toMPO (ampo);
}

template <typename SiteType>
MPS make_initstate (const SiteType& sites, int Np)
{
    int N = length(sites);
    int Npar = Np;
    InitState init (sites);
    for(int i = 1; i <= N; i++)
    {
        string state;
        if (i % 2 == 0 and Np > 0)
        {
            state = "Occ";
            Np--;
        }
        else
            state = "Emp";
        init.set (i, state);
        cout << i << ": " << state << endl;
    }
    if (Np > 0)
    {
        for(int i = 1; i <= N; i += 2)
            if (Np > 0)
            {
                init.set (i,"Occ");
                Np--;
            }
    }
    auto psi = MPS (init);

    auto Nmpo = Make_NMPO (sites);
    auto Ntot = inner (psi,Nmpo,psi);
    if (Ntot != Npar)
    {
        cout << "particle number not match:" << Ntot << " " << Npar << endl;
        throw;
    }
    return psi;
}

template <typename SiteType>
MPS make_initstate_charge_site (const SiteType& sites, int Np, const WireSystem& sys)
{
    int charge_site = sys.to_glob ("C",1);
    int N = length(sites);
    int Npar = Np;
    InitState init (sites);
    int Np_device = 0;
    for(int i = 1; i <= N; i++)
    {
        string state;
        if (i != charge_site)
        {
            if (i % 2 == 0 and Np > 0)
            {
                state = "Occ";
                Np--;
                // count number of particle in device
                auto [seg,ind] = sys.to_loc (i);
                if (seg == "S")
                    Np_device++;
            }
            else
                state = "Emp";
            init.set (i, state);
            cout << i << ": " << state << endl;
        }
    }
    if (Np > 0)
    {
        for(int i = 1; i <= N; i += 2)
            if (Np > 0 and i != charge_site)
            {
                init.set (i,"Occ");
                Np--;
                cout << i << ": Occ" << endl;

                // count number of particle in device
                auto [seg,ind] = sys.to_loc (i);
                if (seg == "S")
                    Np_device++;
            }
    }
    // Set charge site
    init.set (charge_site, str(Np_device));
    cout << charge_site << ": " << Np_device << endl;
    auto psi = MPS (init);

    return psi;
}

template <typename T>
bool ordered (const vector<T>& ns, T resolution=1e-4)
{
    for(int i = 1; i < ns.size(); i++)
    {
        T n1 = ns.at(i-1),
          n2 = ns.at(i);
        if (abs(n1-n2) > resolution and n2 > n1)
            return false;
    }
    return true;
}

Real den (const SiteSet& sites, const MPS& psi, int i)
{
    AutoMPO ampo (sites);
    ampo += 1.,"N",i;
    return inner(psi,ampo,psi);
}

template <typename MPSType>
ITensor print_wf (const MPSType& psi)
{
    ITensor pp (1.);
    vector<Index> iis;
    for(int i = 1; i <= length(psi); i++)
    {
        pp *= psi(i);
        auto is = findIndex (psi(i), "Site,0");
        iis.push_back (is);
        if constexpr (is_same_v <MPO, MPSType>)
            iis.push_back (prime(is));
    }
    pp.permute (iis);
    PrintData(pp);
    return pp;
}

template <typename SiteType, typename NumType>
Mat<NumType> exact_corr2 (const MPS& psi, const WireSystem& sys)
{
    int N = length (psi);
    auto sites = SiteType (siteInds(psi));
    auto corr = Mat<NumType> (N,N);

    auto apply_op = [&sites] (MPS& phi, string op, int j)
    {
        phi.ref(j) *= sites.op(op,j);
        phi.ref(j).noPrime("Site");
    };

    int ic = sys.to_glob ("C",1);
    for(int i = 1; i <= N; i++)
        for(int j = i; j <= N; j++)
        {
            auto [seg1, iseg1] = sys.to_loc(i);
            auto [seg2, iseg2] = sys.to_loc(j);
            if (seg1 == "C" or seg2 == "C")
                continue;

            auto phi = psi;
            apply_op (phi, "C", j);
            apply_op (phi, "Cdag", i);
            if (seg1 != "S" and seg2 == "S")
                apply_op (phi, "A", ic);
            else if (seg1 == "S" and seg2 != "S")
                apply_op (phi, "Adag", ic);

            if constexpr (is_same_v <NumType, Real>)
                corr(i-1,j-1) = inner (psi, phi);
            else
                corr(i-1,j-1) = innerC (psi, phi);
            if (i != j)
                corr(j-1,i-1) = corr(i-1,j-1);
        }
    return corr;
}

void writeAll (const string& filename,
               const MPS& psi, const MPO& H,
               const WireSystem& system,
               int step)
{
    ofstream ofs (filename);
    itensor::write (ofs, psi);
    itensor::write (ofs, H);
    itensor::write (ofs, step);
    system.write (ofs);
}

void readAll (const string& filename,
              MPS& psi, MPO& H,
              WireSystem& system,
              int& step)
{
    ifstream ifs = open_file (filename);
    itensor::read (ifs, psi);
    itensor::read (ifs, H);
    itensor::read (ifs, step);
    system.read (ifs);
}

int main(int argc, char* argv[])
{
    string infile = argv[1];
    InputGroup input (infile,"basic");

    auto L_lead   = input.getInt("L_lead");
    auto L_device   = input.getInt("L_device");
    auto t_lead     = input.getReal("t_lead");
    auto t_device   = input.getReal("t_device");
    auto t_contactL = input.getReal("t_contactL");
    auto t_contactR = input.getReal("t_contactR");
    auto mu_leadL   = input.getReal("mu_leadL");
    auto mu_leadR   = input.getReal("mu_leadR");
    auto mu_device  = input.getReal("mu_device");
    auto Delta      = input.getReal("Delta");
    auto Ec         = input.getReal("Ec");
    auto Ng         = input.getReal("Ng");
    auto mu_biasL   = input.getReal("mu_biasL");
    auto mu_biasS   = input.getReal("mu_biasS");
    auto mu_biasR   = input.getReal("mu_biasR");
    auto damp_decay_length = input.getInt("damp_decay_length",10000000);

    auto ConserveN   = input.getYesNo("ConserveN",false);
    auto ConserveNs  = input.getYesNo("ConserveNs",true);
    auto conserveQN  = input.getYesNo("conserveQN",true);

    auto dt            = input.getReal("dt");
    auto time_steps    = input.getInt("time_steps");
    auto NumCenter     = input.getInt("NumCenter");
    auto Truncate      = input.getYesNo("Truncate");
    auto mixNumCenter  = input.getYesNo("mixNumCenter",false);
    auto globExpanN          = input.getInt("globExpanN",std::numeric_limits<int>::max());
    auto globExpanItv        = input.getInt("globExpanItv",1);
    auto globExpanCutoff     = input.getReal("globExpanCutoff",1e-8);
    auto globExpanKrylovDim  = input.getInt("globExpanKrylovDim",3);
    auto globExpanHpsiCutoff = input.getReal("globExpanHpsiCutoff",1e-8);
    auto globExpanHpsiMaxDim = input.getInt("globExpanHpsiMaxDim",300);
    auto globExpanMethod     = input.getString("globExpanMethod","DensityMatrix");
    auto sweeps_t      = Read_sweeps (infile, "time_evolution_sweeps");
    auto sweeps_dmrg   = Read_sweeps (infile, "DMRG_sweeps");

    auto UseSVD        = input.getYesNo("UseSVD",true);
    auto SVDmethod     = input.getString("SVDMethod","gesdd");  // can be also "ITensor"
    auto WriteDim      = input.getInt("WriteDim");

    auto write         = input.getYesNo("write",false);
    auto write_dir     = input.getString("write_dir",".");
    auto write_file    = input.getString("write_file","");
    auto read          = input.getYesNo("read",false);
    auto read_dir      = input.getString("read_dir",".");
    auto read_file     = input.getString("read_file","");

    MPS psi;
    MPO H;
    auto system = WireSystem();
    int step = 1;
    auto sites = MixedBasis();

    if (!read)
    {
        // Define basis
        Real damp_fac = exp(-1./damp_decay_length);
        cout << "H left lead" << endl;
        auto H_leadL = Hamilt_k (L_lead, t_lead, mu_leadL, damp_fac, true, true);
        cout << "H right lead" << endl;
        auto H_leadR = Hamilt_k (L_lead, t_lead, mu_leadR, damp_fac, false, true);
        cout << "H dev" << endl;
        auto H_dev   = Hamilt_k (L_device, t_device, mu_device, damp_fac, true, true);
        auto H_zero  = Matrix(1,1);

        // WireSystem
        system.add_chain ("L",H_leadL);
        system.add_chain ("R",H_leadR);
        system.add_chain ("S",H_dev);
        system.add_chain ("C",H_zero);
        //system.sort_basis (sort_by_energy_S_middle (system.part("S"), {system.part("L"), system.part("R")}));
        system.sort_basis (sort_by_energy_S_middle_charging (system.parts().at("S"), system.parts().at("C"), {system.parts().at("L"), system.parts().at("R")}));
        system.add_hopping ("L","S",-1,1,t_contactL);
        system.add_hopping ("R","S",1,-1,t_contactR);
        system.set_Ec (Ec, Ng);
        system.add_SC ("S", Delta);
        cout << "device site = " << system.idevL() << " " << system.idevR() << endl;

        // Define SiteSet
        int N = system.N();
        int charge_site = system.to_glob ("C",1);
        // Find sites in S
        vector<int> S_sites;
        for(int i = 0; i < system.orbs().size(); i++)
        {
            if (get<0>(system.orbs().at(i)) == "S")
                S_sites.push_back (i+1);
        }
        sites = MixedBasis (N, S_sites, charge_site,
                            {"MaxOcc",L_device,"ConserveN",ConserveN,"ConserveNs",ConserveNs,"conserveQN",conserveQN});
        cout << "charge site = " << charge_site << endl;

        // Get ground state
        {
            auto psi0 = get_ground_state (system, sites, mu_biasL, mu_biasS, mu_biasR, 0);
            auto psi1 = get_ground_state (system, sites, mu_biasL, mu_biasS, mu_biasR, 1);
            psi0.position(1);
            psi1.position(1);

            auto ampo = get_ampo (system, sites, false);
            auto H0 = toMPO (ampo);
            // DMRG
            MyObserver<decltype(sites)> myobs0 (sites, psi0);
            MyObserver<decltype(sites)> myobs1 (sites, psi1);
            Real en0 = dmrg (psi0, H0, sweeps_dmrg, myobs0, {"WriteDim",WriteDim});
            Real en1 = dmrg (psi1, H0, sweeps_dmrg, myobs1, {"WriteDim",WriteDim});
            psi = (en0 < en1 ? psi0 : psi1);
            cout << "Initial state dim = " << maxLinkDim(psi) << endl;
            // Print density
            int iC = system.to_glob ("C",1);
            Real nC = den (sites, psi, iC);
            cout << "Initial Nc = " << nC << endl;

            Real Np = 0.;
            for(int i = 1; i <= length(psi); i++)
            {
                if (i != iC)
                    Np += den (sites, psi, i);
            }
            cout << "Initial Np = " << Np << endl;
        }

        // Make Hamiltonian MPO
        auto ampo = get_ampo (system, sites);
        H = toMPO (ampo);
        cout << "MPO dim = " << maxLinkDim(H) << endl;
    }
    else
    {
        readAll (read_dir+"/"+read_file, psi, H, system, step);
        sites = MixedBasis (siteInds(psi));
    }
    // ======================= Time evolution ========================

    // Observer
    auto obs = TDVPObserver (sites, psi);
    // Current MPO
    int lenL = system.parts().at("L").L();
    int lenS = system.parts().at("S").L();
    vector<int> spec_links = {lenL, lenL+lenS};
    //for(int i = 1; i < system.N_phys(); i++)
    //    spec_links.push_back (i);
    int N = system.N();
    vector<MPO> JMPOs (N);
    for(int i : spec_links)
    {
        auto mpo = get_current_mpo (sites, system, i);
        JMPOs.at(i) = mpo;
    }
    // Correlation
    /*if (SubCorrN < 0 or SubCorrN > length(psi))
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
    auto sub_corr = SubCorr (system, ibeg, iend);*/

    // Time evolution
    cout << "Start time evolution" << endl;
    cout << sweeps_t << endl;
    psi.position(1);
    Real en, err;
    Args args_tdvp_expansion = {"Cutoff",globExpanCutoff, "Method","DensityMatrix",
                                "KrylovOrd",globExpanKrylovDim, "DoNormalize",true, "Quiet",true};
    Args args_tdvp  = {"Quiet",true,"NumCenter",NumCenter,"DoNormalize",true,"Truncate",Truncate,
                       "UseSVD",UseSVD,"SVDmethod",SVDmethod,"WriteDim",WriteDim,"mixNumCenter",mixNumCenter};
    LocalMPO PH (H, args_tdvp);
    while (step <= time_steps)
    {
        cout << "step = " << step << endl;

        // Subspace expansion
        if (maxLinkDim(psi) < sweeps_t.mindim(1) or (step < globExpanN and (step-1) % globExpanItv == 0))
        {
            timer["glob expan"].start();
            addBasis (psi, H, globExpanHpsiCutoff, globExpanHpsiMaxDim, args_tdvp_expansion);
            PH.reset();
            timer["glob expan"].stop();
        }
        // Time evolution
        timer["tdvp"].start();
        //tdvp (psi, H, 1_i*dt, sweeps, obs, args_tdvp);
        TDVPWorker (psi, PH, 1_i*dt, sweeps_t, obs, args_tdvp);
        timer["tdvp"].stop();
        auto d1 = maxLinkDim(psi);

        // Measure currents by correlations
        /*timer["current corr"].start();
        sub_corr.measure<MixedBasis> (psi, {"Cutoff",corr_cutoff});
        timer["current corr"].stop();
        for(int j = 1; j < system.N_phys(); j++)
        {
            timer["current corr"].start();
            auto J = sub_corr.get_current (j);
            timer["current corr"].stop();
            cout << "\t*current " << j << " " << j+1 << " = " << J << endl;
        }*/
        // Measure currents by MPO
        for(int j : spec_links)
        {
            timer["current mps"].start();
            auto J = get_current (JMPOs.at(j), psi);
            timer["current mps"].stop();
            cout << "\t*current spec " << j << " " << j+1 << " = " << J << endl;
        }
        step++;
        if (write)
        {
            timer["write"].start();
            writeAll (write_dir+"/"+write_file, psi, H, system, step);
            timer["write"].stop();
        }
    }
    timer.print();
    return 0;
}
