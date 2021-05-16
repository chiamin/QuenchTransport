#ifndef __SpecialMixedSiteSet_H_CMC__
#define __SpecialMixedSiteSet_H_CMC__
#include "itensor/all.h"
#include "SpecialFermion.h"
#include "ContainerUtility.h"
using namespace itensor;

class SpecialBosonSite
{
    Index s;
    public:

    SpecialBosonSite(Index I) : s(I) { }

    SpecialBosonSite(Args const& args = Args::global())
    {
        auto conserveNb = args.getBool("ConserveN",true);
        auto conserveNs = args.getBool("ConserveNs",true);

        auto tags = TagSet("Site,Boson");
        auto n = 1;
        if(args.defined("SiteNumber") )
        {
            n = args.getInt("SiteNumber");
            tags.addTags("n="+str(n));
        }

        auto maxOcc = args.getInt("MaxOcc",1);
        {
            if(conserveNb)
            {
                auto qints = Index::qnstorage(1+maxOcc);
                for(int n : range(1+maxOcc)) 
                {
                    qints[n] = QNInt(QN({"Nf",0,-1},{"Ns",-n,-1}),1);
                }
                s = Index(std::move(qints),Out,tags);
            }
            else
            {
                auto qints = Index::qnstorage(1+maxOcc);
                for(int n : range(1+maxOcc)) 
                {
                    int p = n % 2;
                    if (conserveNs)
                        qints[n] = QNInt(QN({"Pf",0,-2},{"Ns",-n,-1}),1);
                    else
                        qints[n] = QNInt(QN({"Pf",0,-2},{"Ps",p,-2}),1);
                }
                s = Index(std::move(qints),Out,tags);
            }
        }
    }

    Index
    index() const { return s; }

    IndexVal
    state(std::string const& state)
    {
        auto maxOcc = dim(s)-1;
        for(auto n : range(1+maxOcc))
        {
            if(state == str(n)) return s(1+n);
        }
	    if(state == "Emp")
        {
            return s(1);
        }
        throw ITError("State " + state + " not recognized");
        return IndexVal{};
    }

	ITensor op (std::string const& opname, Args const& args) const
    {
        auto sP = prime(s);
        auto maxOcc = dim(s)-1;

        auto Op = ITensor(dag(s),sP);

        if(opname == "N" || opname == "n")
        {
            for(auto n : range1(1+maxOcc))
            {
                Op.set(s=n,sP=n,n-1);
            }
        }
        else if(opname == "NSqr" || opname == "nSqr")
        {
            for(auto n : range1(1+maxOcc))
            {
                Op.set(s=n,sP=n,(n-1)*(n-1));
            }
        }
        else if(opname == "A" or opname == "C")
        {
            for(auto n : range1(maxOcc))
            {
                Op.set(s=1+n,sP=n,1);
            }
        }
        else if(opname == "Adag" or opname == "Cdag")
        {
            for(auto n : range1(maxOcc))
            {
                Op.set(s=n,sP=1+n,1);
            }
        }
        else if(opname == "I")
        {
            for(auto n : range1(1+maxOcc))
            {
                Op.set(s=n,sP=n,1);
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
    public:

    MixedBasis() {}

    MixedBasis (int N, int iL, int iR, int iC, Args const& args=Args::global())
    {
        auto sites = SiteStore(N);
        for(int j = 1; j <= N; ++j)
        {
            bool special_qn = (j >= iL and j <= iR);
            if(j == iC) sites.set (j,SpecialBosonSite   ({args,"SiteNumber=",j}));
            else        sites.set (j,SpecialFermionSite (special_qn, {args,"SiteNumber=",j}));
        }
        SiteSet::init(std::move(sites));
    }

    MixedBasis (int N, const vector<int>& C_sites, int iC, Args const& args=Args::global())
    {
        auto sites = SiteStore(N);
        for(int j = 1; j <= N; ++j)
        {
            bool special_qn = vectool::in_vector (C_sites, j);
            if(j == iC) sites.set (j,SpecialBosonSite   ({args,"SiteNumber=",j}));
            else        sites.set (j,SpecialFermionSite (special_qn, {args,"SiteNumber=",j}));
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
                sites.set(j,SpecialBosonSite(ii));
            else
                sites.set(j,SpecialFermionSite(ii));
        }
        SiteSet::init(std::move(sites));
    }
};
#endif
