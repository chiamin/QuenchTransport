#ifndef __MIXEDBASIS_H_CMC__
#define __MIXEDBASIS_H_CMC__
#include <unordered_map>
#include "GeneralUtility.h"
#include "StringUtility.h"
#include "ReadInput.h"
using namespace std;

template <typename T>
using OpTerm1T = tuple <T, string, int>;
template <typename T>
using OpTerm2T = tuple <T, string, int, string, int>;
using OpTerm1 = OpTerm1T<Real>;
using OpTerm2 = OpTerm2T<Real>;
using OpTerm1C = OpTerm1T<Cplx>;
using OpTerm2C = OpTerm2T<Cplx>;

enum Segment { SegL, SegR, SegS, SegNum };

namespace itensor
{
    template <typename T>
    void write (ostream& s, const vector<OpTerm1T<T>>& v)
    {
        auto size = v.size();
        itensor::write (s,size);
        for(auto const& [a1,a2,a3] : v)
        {
            itensor::write (s, a1);
            itensor::write (s, a2);
            itensor::write (s, a3);
        }
    }
    template <typename T>
    void read (istream& s, vector<OpTerm1T<T>>& v)
    {
        auto size = v.size();
        itensor::read (s,size);
        v.resize(size);
        for(auto& t : v)
        {
            T       a1;
            string  a2;
            int     a3;
            itensor::read (s, a1);
            itensor::read (s, a2);
            itensor::read (s, a3);
            t = make_tuple (a1, a2, a3);
        }
    }
    template <typename T>
    void write (ostream& s, const vector<OpTerm2T<T>>& v)
    {
        auto size = v.size();
        itensor::write (s,size);
        for(auto const& [a1,a2,a3,a4,a5] : v)
        {
            itensor::write (s, a1);
            itensor::write (s, a2);
            itensor::write (s, a3);
            itensor::write (s, a4);
            itensor::write (s, a5);
        }
    }
    template <typename T>
    void read (istream& s, vector<OpTerm2T<T>>& v)
    {
        auto size = v.size();
        itensor::read (s,size);
        v.resize(size);
        for(auto& t : v)
        {
            T       a1;
            string  a2, a4;
            int     a3, a5;
            itensor::read (s, a1);
            itensor::read (s, a2);
            itensor::read (s, a3);
            itensor::read (s, a4);
            itensor::read (s, a5);
            t = make_tuple (a1, a2, a3, a4, a5);
        }
    }
    void write (ostream& s, const pair<Segment,int>& p)
    {
        itensor::write (s, p.first);
        itensor::write (s, p.second);
    }
    void read (istream& s, pair<Segment,int>& p)
    {
        Segment a1;
        int    a2;
        itensor::read (s, a1);
        itensor::read (s, a2);
        p = make_pair (a1, a2);
    }
    void write (ostream& s, const vector<pair<Segment,int>>& m)
    {
        auto size = m.size();
        itensor::write (s,size);
        for (auto const& p : m)
        {
            itensor::write (s, p.first);
            itensor::write (s, p.second);
        }
    }
    void read (istream& s, vector<pair<Segment,int>>& m)
    {
        auto size = m.size();
        itensor::read (s,size);
        m.resize (size);
        for (int i = 0; i < size; i++)
        {
            pair<Segment,int> a;
            itensor::read (s, a.first);
            itensor::read (s, a.second);
            m.at(i) = (a);
        }
    }
}

class NonInterChain
{
    public:
        NonInterChain () {}
        NonInterChain (const Matrix& H0, Segment seg)
        : _seg (seg)
        {
            diagHermitian (H0, _Uik, _ens);
            for(int i = 0; i < _ens.size(); i++)
                _ops.emplace_back (_ens(i),"N",i+1);
        }

        void add_mu (Real mu)
        {
            for(int i = 0; i < _ens.size(); i++)
                _ops.emplace_back (-mu,"N",i+1);
        }

              Segment          seg ()      const { return _seg; }
        const Matrix&          Uik ()      const { return _Uik; }
              Vector           Ui  (int i) const { return Vector (row (_Uik, i-1)); }
        const Vector&          ens ()      const { return _ens; }
        const int              L   ()      const { return _ens.size(); }
        const vector<OpTerm1>& ops ()      const { return _ops; }

        void write (ostream& s) const
        {
            itensor::write(s,_Uik);
            itensor::write(s,_ens);
            itensor::write(s,_ops);
        }
        void read (istream& s)
        {
            itensor::read(s,_Uik);
            itensor::read(s,_ens);
            itensor::read(s,_ops);
        }

    private:
        Matrix _Uik;
        Vector _ens;
        vector<OpTerm1> _ops;
        Segment _seg;
};

// Hopping between two non-interacting chains
// Return: 1. Hoppings from chain1 to chain2
//         2.               chain2 to chain1
template <typename T>
vector<OpTerm2T<T>>
cross_hopping (const NonInterChain& chain1, const NonInterChain& chain2, int i1, int i2, T t, Real cutoff=1e-18)
{
    vector<OpTerm2T<T>> ops12;
    auto uik1 = conj (chain1.Ui(i1));
    auto uik2 = chain2.Ui(i2);
    for(int ki1 = 0; ki1 < uik1.size(); ki1++)
    {
        auto u1 = uik1 (ki1);
        for(int ki2 = 0; ki2 < uik2.size(); ki2++)
        {
            auto u2 = uik2 (ki2);
            auto coef = -t*u1*u2;
            if (abs(coef) > cutoff)
                ops12.emplace_back (coef,"Cdag",ki1+1,"C",ki2+1);     // Cdag_ki1 C_ki2
        }
    }
    return ops12;
}

template <typename T>
inline vector<OpTerm2T<T>> dagger (vector<OpTerm2T<T>> ops, T coef=1.)
{
    for(auto& op : ops)
    {
        auto [val, op1, i1, op2, i2] = op;
        get<0>(op) = coef * val;
        get<2>(op) = i2;
        get<4>(op) = i1;
    }
    return ops;
}

template <typename T>
inline void mult (vector<OpTerm2T<T>>& ops, T val)
{
    for(auto& op : ops)
        get<0>(op) *= val;
}

class WireSystem
{
    public:
        WireSystem () {}
        WireSystem (const Matrix& H0_L, const Matrix& H0_S, const Matrix& H0_R, bool verbose=false);

        void    attach_leads (Real tL, Real tR);
        void    add_mu (Segment part, Real mu);

        void    add_diag  (AutoMPO& ampo, Segment p, const NonInterChain& chain) const;
        template <typename T>
        void    add_cross (AutoMPO& ampo, Segment p1, Segment p2, const vector<OpTerm2T<T>>& hops) const;

        AutoMPO ampo (const Fermion& sites) const;

              int   L       () const { return _leadL.L() + _leadR.L() + _dev.L(); }
              int   idevL   () const { return _leadL.L()+1; }
              int   idevR   () const { return _leadL.L()+_dev.L(); }
        const auto& leadL   () const { return _leadL; }
        const auto& leadR   () const { return _leadR; }
        const auto& dev     () const { return _dev; }
        const auto& orb_map (int sec, int ki) const { return _orb_map.at(sec).at(ki); }
        const auto& glob_to_loc (int i) const { return _orb_map_to_local.at(i); }
        const auto& chain   (Segment p) const
        {
            if (p == SegL) return _leadL;
            else if (p == SegR) return _leadR;
            else if (p == SegS) return _dev;
            else throw;
        }

        void write (ostream& s) const;
        void read  (istream& s);
        void write (const string& fname) const;
        void read  (const string& fname);

    private:
        NonInterChain               _leadL, _leadR, _dev;
        vector<vector<int>>         _orb_map;
        vector<pair<Segment,int>>   _orb_map_to_local;  // ortical index -> {partition, ki}
        vector<OpTerm2>             _hop_LS, _hop_SL, _hop_RS, _hop_SR;  // contact hoppings
};

vector<tuple<Segment,int,Real>>
sort_by_energy
(std::initializer_list<NonInterChain> chains)
{
    using Info = tuple<Segment,int,Real>;    // partition, orbital index, energy
    vector<Info> orbs;

    for(auto const& chain : chains)
    {
        for(int i = chain.L(); i > 0; i--)
            orbs.emplace_back (chain.seg(), i, chain.ens()(i-1));
    }
    // Sort the orbitals based on the energies
    auto sort_func = [] (const Info& s1, const Info& s2)
    {
        return get<2>(s1) < get<2>(s2);
    };
    std::sort (orbs.begin(), orbs.end(), sort_func);
    return orbs;
}

vector<tuple<Segment,int,Real>>
sort_by_energy_S_middle
(const NonInterChain& chainS, std::initializer_list<NonInterChain> other_chains)
{
    auto orb_S = sort_by_energy ({chainS});
    auto orbs = sort_by_energy (other_chains);
    auto it = orbs.begin();
    for(; it != orbs.end(); it++)
    {
        auto [seg,i,en] = *it;
        if (en > 0.) break;
    }
    orbs.insert (it, orb_S.begin(), orb_S.end());
    return orbs;
}

WireSystem :: WireSystem (const Matrix& H0_L, const Matrix& H0_S, const Matrix& H0_R, bool verbose)
: _leadL (H0_L, SegL)
, _leadR (H0_R, SegR)
, _dev (H0_S, SegS)
, _orb_map (SegNum)
{
    // --- Find the order for the leads ---
    // Note: orbs[i]={segment,k,energy}, i is 0-index, k is 1-index
    //auto orbs = sort_by_energy ({_dev, _leadL, _leadR});
    auto orbs = sort_by_energy_S_middle (_dev, {_leadL, _leadR});

    // --- Dictionary from {segment, k index} to orbital index ---
    _orb_map[SegL].resize (_leadL.L()+1);
    _orb_map[SegR].resize (_leadR.L()+1);
    _orb_map[SegS].resize (_dev.L()+1);
    _orb_map_to_local.resize (orbs.size()+1);
    if (verbose)
        cout << "orbitals, segment, ki, energy" << endl;
    for(int i = 1; i <= orbs.size(); i++)
    {
        auto [seg, ki, en] = orbs.at(i-1);
        _orb_map.at(seg).at(ki) = i;
        _orb_map_to_local.at(i) = make_pair (seg, ki);
        if (verbose)
            cout << i << " " << seg << " " << ki << " " << en << endl;
    }
    // Note: _orb_map[seg,ki] = i,
    //       _orb_map_to_local[i] = {seg,ki},
    //       both ki and i are 1-index
}

void WireSystem :: attach_leads (Real tL, Real tR)
{
    _hop_LS = cross_hopping (_leadL, _dev, _leadL.L(), 1, tL);
    _hop_RS = cross_hopping (_leadR, _dev, 1, _dev.L(), tR);
    _hop_SL = dagger (_hop_LS);
    _hop_SR = dagger (_hop_RS);
}

void WireSystem :: add_mu (Segment part, Real mu)
{
    mycheck (part == SegL or part == SegR or part == SegS, "invalid keyword");
    NonInterChain *pt = 0;
    if (part == SegL)      pt = &_leadL;
    else if (part == SegR) pt = &_leadR;
    else                  pt = &_dev;
    pt->add_mu (mu);
}

void WireSystem :: add_diag (AutoMPO& ampo, Segment p, const NonInterChain& chain) const
{
    for(auto [coef, op, i] : chain.ops())
    {
        int j = this->orb_map (p,i);
        ampo += coef, op, j;
    }
}

template <typename T>
void WireSystem :: add_cross (AutoMPO& ampo, Segment p1, Segment p2, const vector<OpTerm2T<T>>& hops) const
{
    for(auto [coef, op1, i1, op2, i2] : hops)
    {
        int j1 = this->orb_map (p1,i1);
        int j2 = this->orb_map (p2,i2);
        ampo += coef, op1, j1, op2, j2;
    }
}

AutoMPO WireSystem :: ampo (const Fermion& sites) const
{
    mycheck (length(sites) == L(), "size not match");

    AutoMPO ampo (sites);
    // Diagonal terms
    add_diag (ampo, SegL, _leadL);
    add_diag (ampo, SegR, _leadR);
    add_diag (ampo, SegS, _dev);

    // Contact hoppings
    add_cross (ampo, SegL, SegS, _hop_LS);
    add_cross (ampo, SegS, SegL, _hop_SL);
    add_cross (ampo, SegR, SegS, _hop_RS);
    add_cross (ampo, SegS, SegR, _hop_SR);

    return ampo;
}

void WireSystem :: write (ostream& s) const
{
    itensor::write(s,_leadL);
    itensor::write(s,_leadR);
    itensor::write(s,_dev);
    itensor::write(s,_orb_map);
    itensor::write(s,_orb_map_to_local);
    itensor::write(s,_hop_LS);
    itensor::write(s,_hop_SL);
    itensor::write(s,_hop_RS);
    itensor::write(s,_hop_SR);
}

void WireSystem :: read (istream& s)
{
    itensor::read(s,_leadL);
    itensor::read(s,_leadR);
    itensor::read(s,_dev);
    itensor::read(s,_orb_map);
    itensor::read(s,_orb_map_to_local);
    itensor::read(s,_hop_LS);
    itensor::read(s,_hop_SL);
    itensor::read(s,_hop_RS);
    itensor::read(s,_hop_SR);
}

void WireSystem :: write (const string& fname) const
{
    ofstream ofs (fname);
    this->write (ofs);
}

void WireSystem :: read (const string& fname)
{
    ifstream ifs = open_file (fname);
    this->read (ifs);
}

tuple <Segment,int> get_loc (const WireSystem& sys, int i)
{
    Segment part;
    if (i < sys.idevL())
    {
        part = SegL;
    }
    else if (i >= sys.idevL() && i <= sys.idevR())
    {
        part = SegS;
        i -= sys.idevL()-1;
    }
    else
    {
        part = SegR;
        i -= sys.idevR();
    }
    return {part, i};
}
#endif
