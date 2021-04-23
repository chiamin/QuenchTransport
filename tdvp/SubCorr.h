#ifndef __SUBCORR_H_CMC__
#define __SUBCORR_H_CMC__
#include <unordered_map>
#include "itensor/all.h"
#include "Corr.h"
#include "SystemStruct.h"
using namespace std;
using namespace itensor;

// Measure only correlations in a block for part of orbitals
// Keep track of the orbital indices in different orders
class SubCorr
{
    public:
        SubCorr (shared_ptr<const WireSystem> sys, int ibeg, int iend);

        void measure (const MPS& psi, Args args=Args::global());
        Real get_current (int site) const;

    private:
        shared_ptr<const WireSystem> _sys;
        CMatrix                      _subcorr;
        int                          _ibeg, _iend;  // Define the measured correlation block in the global k-basis

        // For each segment, the local k involved in the sub-correlation
        // _sub_local_ks.at(segment) = [k in sub-correlation]
        unordered_map<string,vector<int>> _sub_local_ks;

        // The correlations C_{s1,s2}(k1,k2) between two segments
        map<pair<string,string>,ITensor> _corr_segs;

        ITensor Corr_segs   (string s1, string s2);
        ITensor get_Uik     (string seg, int i_seg, Index ii) const;
        Cplx    get_corr    (int k1, int k2) const;        // (k1,k2) are global k indices
};

SubCorr :: SubCorr (shared_ptr<const WireSystem> sys, int ibeg, int iend)
: _sys (sys)
, _ibeg (ibeg)
, _iend (iend)
{
    mycheck (ibeg >= 1 && iend >= ibeg, "invalid range");
    cout << "SubCorr block: " << ibeg << " " << iend << endl;

    // Map the global k-indices to local (segment, ki)
    // Store them for each segment
    for(int i = ibeg; i <= iend; i++)
    {
        auto [seg, i_seg] = _sys->to_loc(i);
        _sub_local_ks[seg].push_back (i_seg);
    }
}

// Measure the correlation in the block between <ibeg> and <iend> in the global k-basis
inline void SubCorr :: measure (const MPS& psi, Args args)
{
    mycheck (_iend <= length(psi), "out of MPS size");
    args.add ("corr_ibeg",_ibeg);
    args.add ("corr_iend",_iend);
    _subcorr = Measure_corr <Fermion,Cplx> (psi, args);

    /*// Normalize
    auto Nmpo = Make_NMPO (siteInds(psi));
    auto Np = real(innerC(psi,Nmpo,psi));
    Real trace = 0.;
    for(int i = 0; i < nrows(_subcorr); i++)
        trace += real(_subcorr (i,i));
    cout << "trace = " << trace << " " << Np << endl;
    _subcorr *= (Np / trace);*/

    // Store the correlation between two segments
    _corr_segs.clear();
    vector<pair<string,string>> seg_pairs;
    seg_pairs.emplace_back ("L", "L");
    seg_pairs.emplace_back ("L", "S");
    seg_pairs.emplace_back ("S", "S");
    seg_pairs.emplace_back ("S", "R");
    seg_pairs.emplace_back ("R", "R");
    for(auto seg_pair : seg_pairs)
    {
        auto [s1, s2] = seg_pair;
        auto scorr = Corr_segs (s1, s2);
        _corr_segs.insert ({seg_pair, scorr});
    }
}

// Compute the current at real-space link (site, site+1) based on the sub-correlation
// The current between (i1,i2) = i * \sum_{k1,k2} [ (u*)(i1,k1) u(i2,k2) corr(k1,k2) - H.C.].
// Suppose i1 and i2 belong to seg1 and seg2 respectively;
// k1 and k2 sum over the basis states in seg1 and seg2.
Real SubCorr :: get_current (int site) const
{
    int i1 = site,
        i2 = site+1;

    // Get local (segment, i) from global site indices
    auto [seg1, i1_seg] = get_loc (*_sys, i1);
    auto [seg2, i2_seg] = get_loc (*_sys, i2);

    auto corr = _corr_segs.at(make_pair(seg1,seg2));
    auto u1 = get_Uik (seg1, i1_seg, corr.inds()(1));
    auto u2 = get_Uik (seg2, i2_seg, corr.inds()(2));
    auto Jr = eltC (conj(u1) * corr * u2);
    return -2. * imag(Jr);
}

// Collect the sub-correlation between segments <s1> and <s2>
ITensor SubCorr :: Corr_segs (string s1, string s2)
{
    // Find out the involved k in each segment
    auto ks1 = _sub_local_ks.at(s1);
    auto ks2 = _sub_local_ks.at(s2);

    int N1 = ks1.size();
    int N2 = ks2.size();
    Index ii1(N1), ii2(N2);
    ITensor scorr (ii1, ii2);

    // Collect the correlation
    // Note that the order is based on _sub_local_ks.at(seg)
    for(int i1 = 0; i1 < N1; i1++)
    {
        int k1 = ks1.at(i1);                                    // For k1 in segment 1,
        int k1_glob = _sys->to_glob (s1, k1);                   // find out the corresponding global index k1_glob
        for(int i2 = 0; i2 < N2; i2++)
        {
            int k2 = ks2.at(i2);                                // For k2 in segment 2,
            int k2_glob = _sys->to_glob (s2, k2);               // find out the corresponding global index k2_glob
            scorr.set (i1+1, i2+1, get_corr (k1_glob, k2_glob));    // Collect from the sub-correlation
        }
    }
    return scorr;
}

// Get the coeffients U(i,k),
// where c_i = \sum_k U(i,k) c_k
// <site> is a global site-index
ITensor SubCorr :: get_Uik (string seg, int i_seg, Index ii) const
{
    // Get U(i,k) for all k
    auto uk = _sys->part(seg).Ui(i_seg);
    // Get U(i,k) for k in sub-correlation
    auto ks = _sub_local_ks.at(seg);
    mycheck (dim(ii) == ks.size(), "dim not match");
    ITensor subu (ii);
    for(int i = 0; i < ks.size(); i++)
    {
        int ki = ks.at(i);
        subu.set (i+1, uk(ki-1));
    }
    return subu;
}

// (k1,k2) are global k indices
inline Cplx SubCorr :: get_corr (int k1, int k2) const
{
    mycheck (k1 >= _ibeg and k1 <= _iend and k2 >= _ibeg and k2 <= _iend, "request (k1,k2) is not measured in the sub-correlation");
    return _subcorr (k1-_ibeg, k2-_ibeg);
}

void mea_current (int ibeg, int iend, Real J_zero, const SubCorr& sub_corr, Real muL, Real muR)
{
    Real Vbias = muR - muL;
    bool backward = (iend < ibeg);
    int itv = (backward ? -1 : 1);
    for(int i = ibeg; ; i += itv)
    {
        auto J = sub_corr.get_current (i);
        cout << "\t*current " << i << " " << i+1 << " = " << J << endl;
        if (abs(J/Vbias) < J_zero) break;
        if (backward and i <= iend) break;
        if (!backward and i >= iend) break;
    }
}
#endif
