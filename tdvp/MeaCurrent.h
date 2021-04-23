#ifndef __MEACURRENT_H_CMC__
#define __MEACURRENT_H_CMC__
#include "MixedBasis.h"
#include "ContainerUtility.h"
#include "MPOGen.h"
#include "MPSUtility.h"
using namespace vectool;

// Current between i and i+1
tuple <string, string, vector<OpTerm2>>
get_current_ops (const WireSystem& sys, int i, Real cutoff=1e-18)
{
    auto [partL, i1] = get_loc (sys, i);
    auto [partR, i2] = get_loc (sys, i+1);
    auto ops12 = cross_hopping (sys.part (partL), sys.part (partR), i1, i2, 1., cutoff);
    return {partL, partR, ops12};
}

template <typename SiteType>
MPO get_current_mpo (const SiteType& sites, const WireSystem& sys, int i, Real cutoff=1e-18)
{
    auto [partL, partR, ops12] = get_current_ops (sys, i, cutoff);
    AutoMPO ampo (sites);
    sys.add_cross (ampo, partL, partR, ops12);
    auto mpo = toMPO (ampo);
    return mpo;
}

Real get_current (const MPO& JMPO, const MPS& psi)
{
    auto J = innerC (psi, JMPO, psi);
    return -2. * imag(J);
}

template <typename SiteType>
MPO get_current_mpo2 (const SiteType& sites, const WireSystem& sys, int i)
{
    auto [p1, i1] = get_loc (sys, i);
    auto [p2, i2] = get_loc (sys, i+1);
    auto u1 = conj (sys.part (p1).Ui(i1));
    auto u2 = sys.part (p2).Ui(i2);

    int N = length(sites);
    auto gen = MPOGen (N, 2);
    for(int k1 = 1; k1 <= u1.size(); k1++)
        for(int k2 = 1; k2 <= u2.size(); k2++)
        {
            int j1 = sys.to_glob (p1,k1);
            int j2 = sys.to_glob (p2,k2);
            if (j1 == j2)
            {
                auto Nop = sites.op("N",j1);
                auto id = Identity (Nop.inds());
                auto CC = id - 2.*Nop;
                gen.set_onsite (j1, CC);
            }
            else if (j1 < j2)
            {
                gen.set_beg (j1, u1(k1-1), "Cdag", 3);
                gen.set_end (j2, u2(k2-1), "C", 3);
            }
            else if (j1 > j2)
            {
                gen.set_beg (j2, -u2(k2-1), "C", 4);
                gen.set_end (j1, u1(k1-1), "Cdag", 4);
            }
        }
    for(int k = 1; k <= N; k++)
    {
        gen.set_long_range (k, 3);
        gen.set_long_range (k, 4);
    }
    auto mpo = gen.makeMPO (sites);
    return mpo;
}
#endif
