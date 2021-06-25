#ifndef __SpecialMixedSiteSet_H_CMC__
#define __SpecialMixedSiteSet_H_CMC__
#include "itensor/all.h"
#include "ContainerUtility.h"
using namespace itensor;

class SpecialFermionSite
{
    Index s;
    public:

    SpecialFermionSite(Index I) : s(I) { }

    SpecialFermionSite(Args const& args = Args::global())
    {
        auto ts = TagSet("Site,Fermion");
        auto n = 1;
        if(args.defined("SiteNumber"))
        {
          n = args.getInt("SiteNumber");
          ts.addTags("n="+str(n));
        }
        auto SC_scatter = args.getBool("SC_scatter");
        auto in_scatter = args.getBool("in_scatter");
        if (SC_scatter)
        {
            if (in_scatter)
            {
                s = Index(QN({"Nf",0,-1},{"Ps",0,-2}),1,
                          QN({"Nf",0,-1},{"Ps",1,-2}),1,
                          Out,ts);
            }
            else
            {
                s = Index(QN({"Nf",0,-1},{"Ps",0,-2}),1,
                          QN({"Nf",1,-1},{"Ps",0,-2}),1,
                          Out,ts);
            }
        }
        else
        {
            if (in_scatter)
            {
                s = Index(QN({"Nf",0,-1},{"Ns",0,-1}),1,
                          QN({"Nf",1,-1},{"Ns",1,-1}),1,
                          Out,ts);
            }
            else
            {
                s = Index(QN({"Nf",0,-1},{"Ns",0,-1}),1,
                          QN({"Nf",1,-1},{"Ns",0,-1}),1,
                          Out,ts);
            }
        }
    }

    Index
    index() const { return s; }

    IndexVal
    state(std::string const& state)
        {
        if(state == "Emp" || state == "0") 
            {
            return s(1);
            }
        else 
        if(state == "Occ" || state == "1") 
            {
            return s(2);
            }
        else
            {
            throw ITError("State " + state + " not recognized");
            }
        return IndexVal{};
        }

	ITensor
	op(std::string const& opname,
	   Args const& args) const
        {
        auto sP = prime(s);

        auto Emp  = s(1);
        auto EmpP = sP(1);
        auto Occ  = s(2);
        auto OccP = sP(2);
         
        auto Op = ITensor(dag(s),sP);

        if(opname == "N" || opname == "n")
            {
            Op.set(Occ,OccP,1);
            }
        else
        if(opname == "C")
            {
            Op.set(Occ,EmpP,1);
            }
        else
        if(opname == "Cdag")
            {
            Op.set(Emp,OccP,1);
            }
        else
        if(opname == "A")
            {
            Op.set(Occ,EmpP,1);
            }
        else
        if(opname == "Adag")
            {
            Op.set(Emp,OccP,1);
            }
        else
        if(opname == "F" || opname == "FermiPhase")
            {
            Op.set(Emp,EmpP,1);
            Op.set(Occ,OccP,-1);
            }
        else
        if(opname == "projEmp")
            {
            Op.set(Emp,EmpP,1);
            }
        else
        if(opname == "projOcc")
            {
            Op.set(Occ,OccP,1); 
            }
        else
            {
            throw ITError("Operator \"" + opname + "\" name not recognized");
            }

        return Op;
        }
};

using SpecialFermion = BasicSiteSet<SpecialFermionSite>;
class SpecialBosonSite
{
    Index s;

    vector<int> _ns;

    public:

    int n (int i) const { return _ns.at(i-1); }

    SpecialBosonSite(Index I) : s(I)
    {
        int d = dim(I);
        int maxOcc = (d-1)/2;
        for(int n = -maxOcc; n <= maxOcc; n++)
            _ns.push_back (n);
    }

    SpecialBosonSite(Args const& args = Args::global())
    {
        auto SC_scatter = args.getBool("SC_scatter");

        auto tags = TagSet("Site,Boson");
        if (args.defined("SiteNumber"))
        {
            auto n = args.getInt("SiteNumber");
            tags.addTags("n="+str(n));
        }

        auto maxOcc = args.getInt("MaxOcc");
        for(int n = -maxOcc; n <= maxOcc; n++)
            _ns.push_back (n);

        auto qints = Index::qnstorage(_ns.size());
        for(int i = 0; i < _ns.size(); i++)
        {
            int n = _ns.at(i);
            if(SC_scatter)
            {
                int p = n % 2;
                qints[i] = QNInt(QN({"Nf",n,-1},{"Ps",p,-2}),1);
            }
            else
            {
                qints[i] = QNInt(QN({"Nf",0,-1},{"Ns",-n,-1}),1);
            }
        }
        s = Index(std::move(qints),Out,tags);
    }

    Index
    index() const { return s; }

    IndexVal
    state(std::string state)
    {
	    if(state == "Emp")
            state = "0";
        for(int i = 0; i < _ns.size(); i++)
        {
            if(state == str(_ns.at(i))) return s(i+1);
        }
        throw ITError("State " + state + " not recognized");
        return IndexVal{};
    }

	ITensor op (std::string const& opname, Args const& args) const
    {
        auto sP = prime(s);

        auto Op = ITensor(dag(s),sP);
        if(opname == "N" || opname == "n")
        {
            for(int i = 0; i < _ns.size(); i++)
            {
                int j = i+1;
                int n = _ns.at(i);
                Op.set(s=j,sP=j,n);
            }
        }
        else if(opname == "NSqr" || opname == "nSqr")
        {
            for(int i = 0; i < _ns.size(); i++)
            {
                int j = i+1;
                int n = _ns.at(i);
                Op.set(s=j,sP=j,n*n);
            }
        }
        else if(opname == "A" or opname == "C")
        {
            for(int i = 1; i < _ns.size(); i++)
            {
                Op.set(s=1+i,sP=i,1);
            }
        }
        else if(opname == "Adag" or opname == "Cdag")
        {
            for(int i = 1; i < _ns.size(); i++)
            {
                Op.set(s=i,sP=1+i,1);
            }
        }
        else if(opname == "I")
        {
            for(int i = 1; i <= _ns.size(); i++)
            {
                Op.set(s=i,sP=i,1);
            }
        }
        else
        {
            throw ITError("Operator \"" + opname + "\" name not recognized");
        }

        return Op;
    }
};

//--------

class MixedBasis : public SiteSet
{
    int _maxOcc;

    public:

    int maxOcc () const { return _maxOcc; }

    MixedBasis() {}

    MixedBasis (int N, int iL, int iR, int iC, Args const& args=Args::global())
    {
        _maxOcc = args.getInt("MaxOcc");
        auto sites = SiteStore(N);
        for(int j = 1; j <= N; ++j)
        {
            bool in_scatter = (j >= iL and j <= iR);
            if(j == iC) sites.set (j,SpecialBosonSite   ({args,"SiteNumber=",j}));
            else        sites.set (j,SpecialFermionSite ({args,"SiteNumber=",j,"in_scatter",in_scatter}));
        }
        SiteSet::init(std::move(sites));
    }

    MixedBasis (int N, const vector<int>& C_sites, int iC, Args const& args=Args::global())
    {
        _maxOcc = args.getInt("MaxOcc");
        auto sites = SiteStore(N);
        for(int j = 1; j <= N; ++j)
        {
            bool in_scatter = vectool::in_vector (C_sites, j);
            if(j == iC) sites.set (j,SpecialBosonSite   ({args,"SiteNumber=",j}));
            else        sites.set (j,SpecialFermionSite ({args,"SiteNumber=",j,"in_scatter",in_scatter}));
        }
        SiteSet::init(std::move(sites));
    }

    MixedBasis (IndexSet const& is)
    {
        int N = is.length();
        auto sites = SiteStore(N);
        for(auto j : range1(N))
        {
            auto ii = is(j);
            mycheck (hasTags(ii,"Boson") or hasTags(ii,"Fermion"), "unknown site index");
            if (hasTags(ii,"Boson"))
            {
                sites.set(j,SpecialBosonSite(ii));
                int d = dim(ii);
                _maxOcc = (d-1)/2;
            }
            else
                sites.set(j,SpecialFermionSite(ii));
        }
        SiteSet::init(std::move(sites));
    }
};
#endif
