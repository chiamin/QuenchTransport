#ifndef __SpecialMixedSiteSet_H_CMC__
#define __SpecialMixedSiteSet_H_CMC__
#include "itensor/all.h"
#include "SpecialFermion.h"
#include "ContainerUtility.h"
using namespace itensor;

class SpecialBosonSite
{
    Index s;

    vector<int> _ns;

    public:

    int n (int i) const { return _ns.at(i-1); }

    SpecialBosonSite(Index I) : s(I) { }

    SpecialBosonSite(Args const& args = Args::global())
    {
        auto conserveNb = args.getBool("ConserveN",true);
        auto conserveNs = args.getBool("ConserveNs",true);
        auto conserveQN = args.getBool("conserveQN",true);

        auto tags = TagSet("Site,Boson");
        if (args.defined("SiteNumber"))
        {
            auto n = args.getInt("SiteNumber");
            tags.addTags("n="+str(n));
        }

        auto maxOcc = args.getInt("MaxOcc",1);
        for(int n = -maxOcc; n <= maxOcc; n++)
            _ns.push_back (n);

        if(conserveNb)
        {
            auto qints = Index::qnstorage(_ns.size());
            for(int i = 0; i < _ns.size(); i++)
            {
                int n = _ns.at(i);
                qints[i] = QNInt(QN({"Nf",0,-1},{"Ns",-n,-1}),1);
            }
            s = Index(std::move(qints),Out,tags);
        }
        else
        {
            auto qints = Index::qnstorage(_ns.size());
            for(int i = 0; i < _ns.size(); i++)
            {
                int n = _ns.at(i);
                int p = n % 2;
                if (conserveNs)
                    qints[i] = QNInt(QN({"Pf",0,-2},{"Ns",-n,-1}),1);
                else if (conserveQN)
                    qints[i] = QNInt(QN({"Pf",0,-2},{"Ps",p,-2}),1);
                else
                    qints[i] = QNInt(QN({"Ps",p,-2}),1);
            }
            s = Index(std::move(qints),Out,tags);
        }
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
            bool special_qn = (j >= iL and j <= iR);
            if(j == iC) sites.set (j,SpecialBosonSite   ({args,"SiteNumber=",j}));
            else        sites.set (j,SpecialFermionSite ({args,"SiteNumber=",j,"special_qn",special_qn}));
        }
        SiteSet::init(std::move(sites));
    }

    MixedBasis (int N, const vector<int>& C_sites, int iC, Args const& args=Args::global())
    {
        _maxOcc = args.getInt("MaxOcc");
        auto sites = SiteStore(N);
        for(int j = 1; j <= N; ++j)
        {
            bool special_qn = vectool::in_vector (C_sites, j);
            if(j == iC) sites.set (j,SpecialBosonSite   ({args,"SiteNumber=",j}));
            else        sites.set (j,SpecialFermionSite ({args,"SiteNumber=",j,"special_qn",special_qn}));
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
