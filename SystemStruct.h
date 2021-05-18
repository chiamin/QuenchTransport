#ifndef __SYSTEMSTRUCT_H_CMC__
#define __SYSTEMSTRUCT_H_CMC__
#include <unordered_map>
#include <map>
#include "GeneralUtility.h"
#include "StringUtility.h"
#include "ReadInput.h"
#include "ReadWriteFile.h"
#include "GateContainer.h"
#include "IUtility.h"
#include "NonInterChain.h"
#include "SortBasis.h"
#include "ContainerUtility.h"
#include "GeneralUtility.h"
using namespace std;

template <typename T>
vector<tuple<T,int,int>>
quadratic_term_coefs (const NonInterChain& chain1, const NonInterChain& chain2, int i1, int i2, bool dag1, bool dag2)
{
    vector<tuple<T,int,int>> ops12;
    auto uik1 = chain1.Ui(i1);
    auto uik2 = chain2.Ui(i2);
    for(int ki1 = 0; ki1 < uik1.size(); ki1++)
    {
        auto u1 = uik1 (ki1);
        if constexpr (is_same_v <decltype(u1), Cplx>)
            if (dag1)
                u1 = std::conj(u1);
        for(int ki2 = 0; ki2 < uik2.size(); ki2++)
        {
            auto u2 = uik2 (ki2);
            if constexpr (is_same_v <decltype(u2), Cplx>)
                if (dag2)
                    u2 = std::conj(u2);
            ops12.emplace_back (u1*u2,ki1+1,ki2+1);     // Cdag_ki1 C_ki2
        }
    }
    return ops12;
}

// Return: Term of coef * Cdag_1 C_2 where 1, 2 are in chain1 and chain2
template <typename T>
vector<OpTerm2T<T>>
CdagC_terms (const NonInterChain& chain1, const NonInterChain& chain2, int i1, int i2, T coef, Real cutoff=1e-18)
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
            auto coefi = coef*u1*u2;
            if (abs(coefi) > cutoff)
                ops12.emplace_back (coefi,"Cdag",ki1+1,"C",ki2+1);     // Cdag_ki1 C_ki2
        }
    }
    return ops12;
}

class WireSystem
{
    public:
        WireSystem () {}

        void  add_chain    (const string& name, const Matrix& H0)   { _parts.emplace (name, NonInterChain (H0, name)); }
        void  add_mu       (const string& part, Real mu)            { _parts.at(part).add_mu (mu); }
        void  add_hopping  (const string& p1, const string& p2, int i1, int i2, Real t)
                                                                    { _hops.emplace_back (p1, p2, i1, i2, -t); }
        void  set_Ec       (Real Ec, Real Ng)                       { _Ec = Ec; _Ng = Ng; }
        void  add_SC       (const string& p, Real Delta)            { _SC.emplace_back (p, Delta); }

        int   N       () const;
        int   N_phys  () const { return has_part("C") ? N()-1 : N(); }
        int   idevL   () const { return _parts.at("L").L()+1; }
        int   idevR   () const { return _parts.at("L").L()+_parts.at("S").L()+_parts.at("C").L(); }

        template <typename SortFunction, typename... FuncArgs>
        void sort_basis (SortFunction sort_func, FuncArgs... args);
        void sort_basis (const vector<BasisInfo>& orbs);


        const auto& to_glob    (const string& seg, int ki) const { return _to_glob.at(seg).at(ki); }
        const auto& to_loc     (int i)                     const { return _to_local.at(i); }
        const auto& parts      ()                          const { return _parts; }
        bool        has_part   (const string& p)           const { return _parts.count(p) != 0; }
        const auto& orbs       ()                          const { return _orbs; }
        const auto& hops       ()                          const { return _hops; }
        Real        Ec         ()                          const { return _Ec; }
        Real        Ng         ()                          const { return _Ng; }
        void        print_orbs ()                          const;
        const auto& SC         ()                          const { return _SC; }

        void write (ostream& s) const;
        void read  (istream& s);
        void write (const string& fname) const;
        void read  (const string& fname);

        template <typename T>
        vector<int> reorder_basis (vector<T>& ns, bool reverse=false);

    private:
        using GlobOrbDict = unordered_map <string, vector<int>>;

        unordered_map <string, NonInterChain>     _parts;
        vector<BasisInfo>                         _orbs;
        GlobOrbDict                               _to_glob;           // {partition, ki} -> ortical index
        vector<pair<string,int>>                  _to_local;          // ortical index -> {partition, ki}
        vector<tuple<string,string,int,int,Real>> _hops;
        Real                                      _Ec, _Ng;
        vector<pair<string,Real>>                 _SC;

        void update_order ();
};

template <typename SortFunction, typename... FuncArgs>
void WireSystem :: sort_basis (SortFunction sort_func, FuncArgs... args)
{
    _orbs = sort_func (args...);
    update_order ();
}

void WireSystem :: sort_basis (const vector<BasisInfo>& orbs)
{
    mycheck (orbs.size() == this->N(), "size not match");
    _orbs = orbs;
    update_order ();
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

template <typename NumType>
void add_CdagC (AutoMPO& ampo, const WireSystem& sys, const string& p1, const string& p2, int i1, int i2, NumType coef)
{
    const auto& chain1 = sys.parts().at(p1);
    const auto& chain2 = sys.parts().at(p2);
    if (i1 < 0) i1 += chain1.L()+1;
    if (i2 < 0) i2 += chain2.L()+1;
    auto terms = quadratic_term_coefs<Real> (chain1, chain2, i1, i2, true, false);

    // 
    bool has_charging = sys.has_part("C");  // charging energy
    int jc = (has_charging ? sys.to_glob ("C",1) : -1);
    string op_charge = "";
    if (has_charging)
    {
        if ((p1 == "L" and p2 == "S") or    // Cdag_L C_S
            (p1 == "R" and p2 == "S"))      // Cdag_R C_S
        {
            op_charge = "A";
        }
        else
        if ((p1 == "S" and p2 == "L") or    // Cdag_S C_L
            (p1 == "S" and p2 == "R"))      // Cdag_S C_R
        {
            op_charge = "Adag";
        }
    }
    // Hopping terms
    for(auto [c12, k1, k2] : terms)
    {
        int j1 = sys.to_glob (p1,k1);
        int j2 = sys.to_glob (p2,k2);
        // hopping
        if (op_charge != "")
        {
//cout << "* " << coef << " " << op1 << " " << j1 << " " << op_charge << " " << jc << " " << op2 << " " << j2 << endl;
            ampo += coef*c12, "Cdag", j1, op_charge, jc, "C", j2;
        }
        else
        {
            ampo += coef*c12, "Cdag", j1, "C", j2;
        }
    }
}

template <typename NumType>
void add_SC (AutoMPO& ampo, const WireSystem& sys, const string& p1, const string& p2, int i1, int i2, NumType coef)
{
    const auto& chain1 = sys.parts().at(p1);
    const auto& chain2 = sys.parts().at(p2);
    if (i1 < 0) i1 += chain1.L()+1;
    if (i2 < 0) i2 += chain2.L()+1;
    auto terms = quadratic_term_coefs<Real> (chain1, chain2, i1, i2, false, false);

    for(auto [c12, k1, k2] : terms)
    {
        int j1 = sys.to_glob (p1,k1);
        int j2 = sys.to_glob (p2,k2);
        auto c = coef * c12;
        auto cc = iutility::conjT (c);
        ampo += -c, "C", j1, "C", j2;
        ampo += -cc, "Cdag", j2, "Cdag", j1;
    }
}

void add_diag (AutoMPO& ampo, const WireSystem& sys, const string& p, const NonInterChain& chain)
{
    for(auto [coef, op, i] : chain.ops())
    {
        int j = sys.to_glob (p,i);
        ampo += coef, op, j;
    }
}

template <typename SiteType>
AutoMPO get_ampo (const WireSystem& sys, const SiteType& sites, bool hopping=true)
{
    mycheck (length(sites) == sys.N(), "size not match");

    AutoMPO ampo (sites);
    // Diagonal terms
    for(auto const& [p, chain] : sys.parts())
        add_diag (ampo, sys, p, chain);
    // Hopping terms
    if (hopping)
    {
        for(auto const& [p1,p2,i1,i2,coef] : sys.hops())
        {
            add_CdagC (ampo, sys, p1, p2, i1, i2, coef);
            add_CdagC (ampo, sys, p2, p1, i2, i1, coef);
        }
    }
    // charging energy
    if (sys.has_part("C"))
    {
        int jc = sys.to_glob ("C",1);
        ampo += sys.Ec(),"NSqr",jc;
        ampo += sys.Ec()*sys.Ng()*sys.Ng(),"I",jc;
        ampo += -2.*sys.Ec()*sys.Ng(),"N",jc;
        //ampo +=  0.5*sys.Ec(),"NSqr",jc;
        //ampo += -0.5*sys.Ec(),"N",jc;
    }
    // superconducting
    for(auto [p, Delta] : sys.SC())
    {
        auto const& chain = sys.parts().at(p);
        for(int i = 1; i < chain.L(); i++)
        {
            add_SC (ampo, sys, p, p, i, i+1, Delta);
        }
    }
    return ampo;
}

template <typename SiteType>
MPS get_ground_state (const WireSystem& sys, const SiteType& sites, Real muL=0., Real muS=0., Real muR=0., int NS=0)
{
    mycheck (length(sites) == sys.N(), "size not match");
    int Ns=0, Np=0;
    Real E = 0.;
    vector<string> state (sys.N()+1, "Emp");
    for(string p : {"L","S","R"})
    {
        auto const& chain = sys.parts().at(p);
        Real mu;
        if (p == "L")      mu = muL;
        else if (p == "R") mu = muR;
        else if (p == "S") mu = muS;
        for(int i = 1; i <= chain.L(); i++)
        {
            auto const& en = chain.ens()(i-1);
            int j = sys.to_glob (p, i);

            if (p == "S")
            {
                if (Ns < NS)
                {
                    state.at(j) = "Occ";
                    E += en;
                    Np++;
                    Ns++;
                }
            }
            else if (en-mu < 0.)
            {
                state.at(j) = "Occ";
                E += en;
                Np++;
            }
        }
    }
    // Capacity site
    {
        int j = sys.to_glob ("C",1);
        state.at(j) = str(Ns);
    }

    InitState init (sites);
    for(int i = 1; i <= sys.N(); i++)
        init.set (i, state.at(i));

    cout << "orbitals, segment, ki, energy, state" << endl;
    for(int i = 1; i <= sys.orbs().size(); i++)
    {
        auto [seg, ki, en] = sys.orbs().at(i-1);
        cout << i << " " << seg << " " << ki << " " << en << " " << state.at(i) << endl;
    }
    cout << "initial energy = " << E << endl;
    cout << "initial particle number = " << Np << endl;
    return MPS (init);
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

int WireSystem :: N () const
{
    int L = 0;
    for(auto const& [p, chain] : _parts)
        L += chain.L();
    return L;
}

void WireSystem :: write (ostream& s) const
{
    iutility::write(s,_parts);
    iutility::write(s,_orbs);
    iutility::write(s,_to_glob);
    iutility::write(s,_to_local);
    iutility::write(s,_hops);
    iutility::write(s,_Ec);
    iutility::write(s,_SC);
}

void WireSystem :: read (istream& s)
{
    iutility::read(s,_parts);
    iutility::read(s,_orbs);
    iutility::read(s,_to_glob);
    iutility::read(s,_to_local);
    iutility::read(s,_hops);
    iutility::read(s,_Ec);
    iutility::read(s,_SC);
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
    string part = "L";
    int L_L = sys.parts().at("L").L();
    int L_S = sys.parts().at("S").L();
    if (i > L_L)
    {
        part = "S";
        i -= L_L;
        if (i > L_S)
        {
            part = "R";
            i -= L_S;
        }
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

    // Make swap gate
    auto i1 = sites(1);
    auto i2 = sites(2);
    ITensor swap (dag(i1), dag(i2), prime(i1), prime(i2));
    swap.set (1,1,1,1,1.);
    swap.set (1,2,2,1,1.);
    swap.set (2,1,1,2,1.);
    swap.set (2,2,2,2,-1.);

    gates.new_gate ("swap", swap, i1, i2);
    vector<int> swap_pos = system.reorder_basis (ns, reverse);
    for(int i : swap_pos)
        gates.add ("swap",i);
    gates.apply (psi, args);
}
#endif
