#include <iomanip>
#include "itensor/all.h"
#include "ReadInput.h"
#include "IUtility.h"
#include "MyObserver.h"
#include "MPSUtility.h"
#include "SingleParticle.h"
#include "MixedBasis.h"
#include "SystemStruct.h"
#include "ContainerUtility.h"
using namespace vectool;
using namespace itensor;
using namespace std;
using namespace iutility;

template <typename SitesType>
MPO current_correlation (const SitesType& sites, int i)
{
    int N = length(sites);
    AutoMPO ampo (sites);
    if constexpr (is_same_v <SitesType, Fermion>)
    {
        ampo += -1_i,"Cdag",i,"C",i+1;
        ampo +=  1_i,"Cdag",i+1,"C",i;
        ampo +=  1_i,"Cdag",N-i,"C",N-i+1;
        ampo += -1_i,"Cdag",N-i+1,"C",N-i;
    }
    return toMPO (ampo);
}

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
        if (i % 2 == 0 && Np-- > 0)
            state = "Occ";
        else
            state = "Emp";
        init.set (i, state);
        cout << i << ": " << state << endl;
    }
    if (Np > 0)
    {
        for(int i = 1; i <= N; i += 2)
            if (Np-- > 0)
                init.set (i,"Occ");
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

int main(int argc, char* argv[])
{
    string infile = argv[1];
    InputGroup input (infile,"basic");

    auto L_lead   = input.getInt("L_lead");
    auto L_device   = input.getInt("L_device");
    auto Np         = input.getInt("N_fermions",-1);
    if (Np == -1) Np = L_lead + L_device/2;
    auto t_lead     = input.getReal("t_lead");
    auto t_device   = input.getReal("t_device");
    auto t_contactL = input.getReal("t_contactL");
    auto t_contactR = input.getReal("t_contactR");
    auto mu_leadL   = input.getReal("mu_leadL");
    auto mu_leadR   = input.getReal("mu_leadR");
    auto mu_device  = input.getReal("mu_device");
    auto V_device   = input.getReal("V_device");
    auto damp_decay_length = input.getReal("damp_decay_length");

    auto do_write   = input.getYesNo("write_to_file");
    auto out_dir    = input.getString("outdir",".");
    auto out_minm   = input.getInt("out_minm",0);
    auto out_Hamilt = input.getString("out_Hamilt","H");
    auto ConserveNf = input.getYesNo("ConserveNf",false);
    auto WriteDim   = input.getInt("WriteDim",-1);
    auto verbose    = input.getYesNo("verbose",false);
    auto sweeps     = Read_sweeps (infile);

    auto reorderN   = input.getInt("reorderN",0);
    auto cutoff_reorder = input.getReal("cutoff_reorder",1e-12);

    // Define basis
    Real damp_fac = exp(-1./damp_decay_length);
    cout << "H left lead" << endl;
    auto H_leadL = Hamilt_k (L_lead, t_lead, mu_leadL, damp_fac, true, true);
    cout << "H right lead" << endl;
    auto H_leadR = Hamilt_k (L_lead, t_lead, mu_leadR, damp_fac, false, true);
    cout << "H dev" << endl;
    auto H_dev   = Hamilt_k (L_device, t_device, mu_device, 1., true, true);
    auto system = WireSystem (H_leadL, H_dev, H_leadR);
    system.attach_leads (t_contactL, t_contactR);
    if (do_write)
        system.write (out_dir+"/"+out_Hamilt);
    system.print_orbs();

    auto sites = Fermion (system.N(), {"ConserveNf",ConserveNf});
   // auto sites = MixedBasis (system.N(), system.N()/2, {"MaxOcc",10,"ConserveQNs",ConserveQNs});
    auto ampo = system.ampo (sites);
    auto H = toMPO (ampo);
    cout << "MPO dim = " << maxLinkDim(H) << endl;
    cout << "device site = " << system.idevL() << " " << system.idevR() << endl;

    // Initialze MPS
    MPS psi;
    if (true)//(ConserveNf)
        psi = make_initstate (sites, Np);
    else
        psi = randomMPS (sites, 10);
    psi.position(1);


    // DMRG
    MyObserver myobs (sites, psi, {"Write",do_write,"out_dir",out_dir,"out_minm",out_minm});
    auto en = dmrg (psi, H, sweeps, myobs, {"WriteDim",WriteDim});
    auto d1 = maxLinkDim(psi);

    // Reorder basis
    {
        auto ns = myobs.ns();
        for(int i = 0; i < reorderN; i++)
        {
            reorder_basis (psi, system, ns, {"Cutoff",cutoff_reorder,"reverse",true});
            auto d2 = maxLinkDim(psi);
            cout << "original/rotated dim = " << d1 << " " << d2 << endl;

            cout << "Redo DMRG" << endl;
            ampo = system.ampo (sites);
            H = toMPO (ampo);
            auto en2 = dmrg (psi, H, sweeps, myobs, {"WriteDim",WriteDim});
            cout << "old/new energy = " << en << " " << en2 << endl;

            ns = myobs.ns();
            if (ordered (ns)) break;
        }
    }
    return 0;
}
