#ifndef __SYSTEMSTRUCT_H_CMC__
#define __SYSTEMSTRUCT_H_CMC__
#include <unordered_map>
#include "GeneralUtility.h"
#include "StringUtility.h"
#include "ReadInput.h"
#include "ReadWriteFile.h"
#include "GateContainer.h"
#include "IUtility.h"
using namespace std;

template <typename T>
using OpTerm1T = tuple <T, string, int>;
template <typename T>
using OpTerm2T = tuple <T, string, int, string, int>;
using OpTerm1 = OpTerm1T<Real>;
using OpTerm2 = OpTerm2T<Real>;
using OpTerm1C = OpTerm1T<Cplx>;
using OpTerm2C = OpTerm2T<Cplx>;
using BasisInfo = tuple <string, int, Real>; // partition, orbital index, energy

class NonInterChain
{
    public:
        NonInterChain () {}
        NonInterChain (const Matrix& H0, string name)
        : _name (name)
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

        const string&          name ()      const { return _name; }
        const Matrix&          Uik  ()      const { return _Uik; }
              Vector           Ui   (int i) const { return Vector (row (_Uik, i-1)); }
        const Vector&          ens  ()      const { return _ens; }
        const int              L    ()      const { return _ens.size(); }
        const vector<OpTerm1>& ops  ()      const { return _ops; }

        void write (ostream& s) const
        {
            itensor::write(s,_Uik);
            itensor::write(s,_ens);
            iutility::write(s,_ops);
        }
        void read (istream& s)
        {
            itensor::read(s,_Uik);
            itensor::read(s,_ens);
            iutility::read(s,_ops);
        }

    private:
        Matrix _Uik;
        Vector _ens;
        vector<OpTerm1> _ops;
        string _name;
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
        WireSystem (const Matrix& H0_L, const Matrix& H0_S, const Matrix& H0_R);

        void    attach_leads (Real tL, Real tR);
        void    add_mu       (const string& part, Real mu);

        void    add_diag  (AutoMPO& ampo, const string& p, const NonInterChain& chain) const;
        template <typename T>
        void    add_cross (AutoMPO& ampo, const string& p1, const string& p2, const vector<OpTerm2T<T>>& hops) const;
        template <typename SiteType>
        AutoMPO ampo      (const SiteType& sites) const;

              int   N       () const;
              int   idevL   () const { return _parts.at("L").L()+1; }
              int   idevR   () const { return _parts.at("L").L()+_parts.at("S").L(); }

        const auto& to_glob (const string& seg, int ki) const { return _to_glob.at(seg).at(ki); }
        const auto& to_loc  (int i)                     const { return _to_local.at(i); }
        const auto& part    (const string& p)           const { return _parts.at(p); }
        const auto& orbs    ()                          const { return _orbs; }

        void print_orbs () const;

        void write (ostream& s) const;
        void read  (istream& s);
        void write (const string& fname) const;
        void read  (const string& fname);

        template <typename T>
        vector<int> reorder_basis (vector<T>& ns, bool reverse=false);

    private:
        using GlobOrbDict = unordered_map <string, vector<int>>;

        unordered_map <string, NonInterChain>   _parts;
        vector<BasisInfo>           _orbs;
        GlobOrbDict                 _to_glob;           // {partition, ki} -> ortical index
        vector<pair<string,int>>    _to_local;          // ortical index -> {partition, ki}
        vector<OpTerm2>             _hop_LS, _hop_SL, _hop_RS, _hop_SR;  // contact hoppings

        void update_order ();
};

// Input: arbitrary number of NonInterChain objects
// Combine all the basis states and sort by the energies
vector<BasisInfo> sort_by_energy (std::initializer_list<NonInterChain> chains)
{
    vector<BasisInfo> orbs;
    for(auto const& chain : chains)
    {
        for(int i = chain.L(); i > 0; i--)
            orbs.emplace_back (chain.name(), i, chain.ens()(i-1));
    }
    // Sort the orbitals based on the energies
    auto sort_func = [] (const BasisInfo& s1, const BasisInfo& s2)
    {
        return get<2>(s1) < get<2>(s2);
    };
    std::sort (orbs.begin(), orbs.end(), sort_func);
    return orbs;
}

// Sort all the basis states by energy; however put the states from <chainS> in zero energy of the other states
vector<BasisInfo>
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

// Reorder _orbs based on <ns>, from small to large
// Update  _to_glob and _to_local
// Return the positions of swap gates
template <typename T>
vector<int> WireSystem :: reorder_basis (vector<T>& ns, bool reverse)
{
    mycheck (ns.size() == _orbs.size(), "size not match");
    vector<int> swap_pos;
    for(int i = 1; i < ns.size(); i++)
        for(int j = i; j >= 1; j--)
        {
            bool do_swap = (ns.at(j) < ns.at(j-1));
            if (reverse)
                do_swap = !do_swap;
            if (do_swap)
            {
                std::swap (ns[j], ns[j-1]);
                std::swap (_orbs[j], _orbs[j-1]);
                swap_pos.push_back (j);
            }
        }
    update_order ();
    return swap_pos;
}

// Update _to_glob and _to_local based on _orbs
void WireSystem :: update_order ()
{
    // --- Dictionary from {segment, k index} to orbital index ---
    // Note: _to_glob[seg,ki] = i,
    //       _to_local[i] = {seg,ki},
    //       both ki and i are 1-index
    _to_glob.clear();
    for(auto const& [name, chain] : _parts)
        _to_glob[name].resize (chain.L()+1);
    _to_local.resize (_orbs.size()+1);
    for(int i = 1; i <= _orbs.size(); i++)
    {
        auto [seg, ki, en] = _orbs.at(i-1);
        _to_glob.at(seg).at(ki) = i;
        _to_local.at(i) = make_pair (seg, ki);
    }
}

inline void WireSystem :: print_orbs () const
{
    cout << "orbitals, segment, ki, energy" << endl;
    for(int i = 1; i <= _orbs.size(); i++)
    {
        auto [seg, ki, en] = _orbs.at(i-1);
        cout << i << " " << seg << " " << ki << " " << en << endl;
    }
}

WireSystem :: WireSystem (const Matrix& H0_L, const Matrix& H0_S, const Matrix& H0_R)
{
    _parts.emplace ("L", NonInterChain (H0_L, "L"));
    _parts.emplace ("R", NonInterChain (H0_R, "R"));
    _parts.emplace ("S", NonInterChain (H0_S, "S"));

    // --- Find the order for the leads ---
    // Note: orbs[i]={segment,k,energy}, i is 0-index, k is 1-index
    //_orbs = sort_by_energy ({_parts.at("S"), _parts.at("L"), _parts.at("R")});
    _orbs = sort_by_energy_S_middle (_parts.at("S"), {_parts.at("L"), _parts.at("R")});

    update_order ();
}

int WireSystem :: N () const
{
    int L = 0;
    for(auto const& [p, chain] : _parts)
        L += chain.L();
    return L;
}

void WireSystem :: attach_leads (Real tL, Real tR)
{
    _hop_LS = cross_hopping (_parts.at("L"), _parts.at("S"), _parts.at("L").L(), 1, tL);
    _hop_RS = cross_hopping (_parts.at("R"), _parts.at("S"), 1, _parts.at("S").L(), tR);
    _hop_SL = dagger (_hop_LS);
    _hop_SR = dagger (_hop_RS);
}

void WireSystem :: add_mu (const string& part, Real mu)
{
    _parts.at(part).add_mu (mu);
}

void WireSystem :: add_diag (AutoMPO& ampo, const string& p, const NonInterChain& chain) const
{
    for(auto [coef, op, i] : chain.ops())
    {
        int j = this->to_glob (p,i);
        ampo += coef, op, j;
    }
}

template <typename T>
void WireSystem :: add_cross (AutoMPO& ampo, const string& p1, const string& p2, const vector<OpTerm2T<T>>& hops) const
{
    for(auto [coef, op1, i1, op2, i2] : hops)
    {
        int j1 = this->to_glob (p1,i1);
        int j2 = this->to_glob (p2,i2);
        ampo += coef, op1, j1, op2, j2;
    }
}

template <typename SiteType>
AutoMPO WireSystem :: ampo (const SiteType& sites) const
{
    mycheck (length(sites) == this->N(), "size not match");

    AutoMPO ampo (sites);
    // Diagonal terms
    for(auto const& [p, chain] : _parts)
        add_diag (ampo, p, chain);

    // Contact hoppings
    add_cross (ampo, "L", "S", _hop_LS);
    add_cross (ampo, "S", "L", _hop_SL);
    add_cross (ampo, "R", "S", _hop_RS);
    add_cross (ampo, "S", "R", _hop_SR);

    return ampo;
}

void WireSystem :: write (ostream& s) const
{
    iutility::write(s,_parts);
    iutility::write(s,_orbs);
    iutility::write(s,_to_glob);
    iutility::write(s,_to_local);
    iutility::write(s,_hop_LS);
    iutility::write(s,_hop_SL);
    iutility::write(s,_hop_RS);
    iutility::write(s,_hop_SR);
}

void WireSystem :: read (istream& s)
{
    iutility::read(s,_parts);
    iutility::read(s,_orbs);
    iutility::read(s,_to_glob);
    iutility::read(s,_to_local);
    iutility::read(s,_hop_LS);
    iutility::read(s,_hop_SL);
    iutility::read(s,_hop_RS);
    iutility::read(s,_hop_SR);
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

tuple <string,int> get_loc (const WireSystem& sys, int i)
{
    string part;
    if (i < sys.idevL())
    {
        part = "L";
    }
    else if (i >= sys.idevL() && i <= sys.idevR())
    {
        part = "S";
        i -= sys.idevL()-1;
    }
    else
    {
        part = "R";
        i -= sys.idevR();
    }
    return {part, i};
}

template <typename T>
void reorder_basis (MPS& psi, WireSystem& system, vector<T>& ns, const Args& args=Args::global())
{
    bool reverse = args.getBool("reverse",false);
    bool verbose = args.getBool("verbose",false);
    auto sites = siteInds (psi);
    auto gates = GateContainer();
    gates.new_gate<2> ("swap", sites, iutility::SwapGate);
    vector<int> swap_pos = system.reorder_basis (ns, reverse);
    for(int i : swap_pos)
        gates.add ("swap",i);
    gates.apply (psi, args);
}
#endif
