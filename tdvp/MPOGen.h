#ifndef __MPOGEN_H_CMC__
#define __MPOGEN_H_CMC__
#include "IUtility.h"
using namespace iutility;

struct OpInfo
{
    Real coef1=0.;
    Real coef2=0.;
    string op1="";
    string op2="";
    Real long_range_coef=0.;
};

class MPOGen
{
    public:
        MPOGen (int N, int opN=0);

        void set_opN (int N);
        void set_beg (int site, Real coef, string op, int n);
        void set_end (int site, Real coef, string op, int n);
        void set_onsite (int site, const ITensor& op);
        void set_long_range (int site, int n, Real coef=1.);
        void set_fermion (bool fermion) { _fermion = fermion; }

        int dim () const { return _Ws.at(1).size()-1; }
        MPO makeMPO (const SiteSet& sites) const;

    private:
        int _nbeg=1, _nend=2;
        vector<vector<OpInfo>>  _Ws; // _Ws[site][op] = OpInfo
        vector<ITensor>         _onsite_op;
        bool _fermion=true;
};

inline MPOGen :: MPOGen (int N, int opN)
: _Ws (N+1)
, _onsite_op (N+1)
{
    set_opN (opN);
}

inline void MPOGen :: set_opN (int N)
{
    for(int i = 1; i < _Ws.size(); i++)
    {
        _Ws.at(i).resize (N+3);
    }
}

inline void MPOGen :: set_beg (int site, Real coef, string op, int n)
{
    mycheck (n > 2, "invalid n");
    _Ws.at(site).at(n).coef1 = coef;
    _Ws.at(site).at(n).op1 = op;
}

inline void MPOGen :: set_end (int site, Real coef, string op, int n)
{
    mycheck (n > 2, "invalid n");
    _Ws.at(site).at(n).coef2 = coef;
    _Ws.at(site).at(n).op2 = op;
}

inline void MPOGen :: set_onsite (int site, const ITensor& op)
{
    _onsite_op.at(site) = op;
}

inline void MPOGen :: set_long_range (int site, int n, Real coef)
{
    mycheck (n > 2, "invalid n");
    _Ws.at(site).at(n).long_range_coef = coef;
}

MPO MPOGen :: makeMPO (const SiteSet& sites) const
{
    MPO mpo (sites);
    int N = length (sites);
    Index li (this->dim(), "Link");
    auto li0 = li;
    for(int i = 1; i <= N; i++)
    {
        // Make tensor
        Index ri = sim (li);
        Index si = sites(i);
        auto& W = mpo.ref(i);
        W = ITensor (dag(si), prime(si), li, dag(ri));

        // Initialize
        auto I = Identity (dag(si), prime(si));
        W += I * setElt(li=_nbeg) * setElt(ri=_nbeg);
        W += I * setElt(li=_nend) * setElt(ri=_nend);

        // Set operators
        auto const& ops = _Ws.at(i);
        for(int n = 3; n < ops.size(); n++)
        {
            auto const& info = ops.at(n);
            // Start operator
            if (info.coef1 != 0.)
            {
                auto op1 = info.coef1 * sites.op (info.op1, i);
                if (_fermion)
                {
                    op1 *= prime(sites.op("F",i));
                    op1.mapPrime (2,1,"Site");
                }
                W += op1 * setElt(li=_nbeg) * setElt(ri=n);
            }
            // End operator
            if (info.coef2 != 0.)
            {
                auto op2 = info.coef2 * sites.op (info.op2, i);
                W += op2 * setElt(li=n) * setElt(ri=_nend);
            }
            // Long-range coefficient
            if (info.long_range_coef != 0.)
            {
                auto op = (_fermion ? sites.op("F",i) : I);
                auto D = info.long_range_coef * op;
                W += D * setElt(li=n) * setElt(ri=n);
            }
        }
        // Onsite operator
        if (_onsite_op.at(i))
        {
            W += _onsite_op.at(i) * setElt(li=_nbeg) * setElt(ri=_nend);
        }

        li = ri;
    }
    mpo.ref(1) *= setElt(dag(li0)=_nbeg);
    mpo.ref(N) *= setElt(dag(li)=_nend);
    return mpo;
}
#endif
