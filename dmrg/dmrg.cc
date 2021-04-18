#include "itensor/all.h"
#include "ReadInput.h"
#include "IUtility.h"
#include "MyObserver.h"
#include "MPSUtility.h"
#include "SingleParticle.h"
#include "MixedBasis.h"
using namespace itensor;
using namespace std;
using namespace myutility;

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

MPS make_initstate (const Fermion& sites, int Np)
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

/*inline ostream& operator<< (ostream& os, const Vector& v)
{
    for(size_t i=0; i < v.size(); i++)
        os << v(i) << " ";
    return os;
}*/

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
    auto mu_lead    = input.getReal("mu_lead");
    auto mu_device  = input.getReal("mu_device");
    auto V_device   = input.getReal("V_device");
    auto damp_decay_length = input.getReal("damp_decay_length");

    auto do_write   = input.getYesNo("write_to_file");
    auto out_dir    = input.getString("outdir",".");
    auto out_minm   = input.getInt("out_minm",0);
    auto out_Hamilt = input.getString("out_Hamilt","H");
    auto ConserveQNs = input.getYesNo("ConserveQNs",false);
    auto WriteDim   = input.getInt("WriteDim",-1);
    auto verbose    = input.getYesNo("verbose",false);
    auto sweeps     = Read_sweeps (infile);

    // Define basis
    Real damp_fac = exp(-1./damp_decay_length);
    cout << "H left lead" << endl;
    auto H_leadL = Hamilt_k (L_lead, t_lead, mu_lead, damp_fac, true, true);
    cout << "H right lead" << endl;
    auto H_leadR = Hamilt_k (L_lead, t_lead, mu_lead, damp_fac, false, true);
    cout << "H dev" << endl;
    auto H_dev   = Hamilt_k (L_device, t_device, mu_device, 1., true, true);
    auto system = WireSystem (H_leadL, H_dev, H_leadR, verbose);
    system.attach_leads (t_contactL, t_contactR);
    if (do_write)
        system.write (out_dir+"/"+out_Hamilt);

    auto sites = Fermion (system.L(), {"ConserveQNs",ConserveQNs});
    auto ampo = system.ampo (sites);
    auto H = toMPO (ampo);
    cout << "MPO dim = " << maxLinkDim(H) << endl;
    cout << "device site = " << system.idevL() << " " << system.idevR() << endl;

    // Initialze MPS
    MPS psi;
    if (ConserveQNs)
        psi = make_initstate (sites, Np);
    else
        psi = randomMPS (sites, 10);
    psi.position(1);


    // DMRG
    MyObserver<Fermion> myobs (sites, psi, {"Write",do_write,"out_dir",out_dir,"out_minm",out_minm});
    dmrg (psi, H, sweeps, myobs, {"WriteDim",WriteDim});
    return 0;
}
