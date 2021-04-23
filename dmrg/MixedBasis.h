#ifndef __SpecialMixedSiteSet_H_CMC__
#define __SpecialMixedSiteSet_H_CMC__
#include "itensor/all.h"

class SpecialBosonSite
{
    Index s;
    public:

    SpecialBosonSite(Index I) : s(I) { }

    SpecialBosonSite(Args const& args = Args::global())
    {
        auto conserveQNs = args.getBool("ConserveQNs",true);
        auto conserveNb = args.getBool("ConserveNb",conserveQNs);

        auto tags = TagSet("Site,Boson");
        auto n = 1;
        if(args.defined("SiteNumber") )
        {
            n = args.getInt("SiteNumber");
            tags.addTags("n="+str(n));
        }

        auto maxOcc = args.getInt("MaxOcc",1);
        if(conserveQNs)
        {
            if(conserveNb)
            {
                auto qints = Index::qnstorage(1+maxOcc);
                for(int n : range(1+maxOcc)) 
                {
                    qints[n] = QNInt(QN({"Nb",n}),1);
                }
                s = Index(std::move(qints),tags);
            }
            else
            {
                s = Index(QN(),1+maxOcc,tags);
            }
        }
        else
        {
            if(conserveNb) throw ITError("ConserveNb cannot be true when ConserveQNs=false");
            s = Index(1+maxOcc,tags);
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
else if(state == "Occ" || state == "1") 
{
    return s(2);
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
        else if(opname == "A" || opname == "C")
        {
            for(auto n : range1(maxOcc))
            {
                Op.set(s=1+n,sP=n,std::sqrt(n));
            }
        }
        else if(opname == "Adag" || opname == "Cdag")
        {
            for(auto n : range1(maxOcc))
            {
                Op.set(s=n,sP=1+n,std::sqrt(n));
            }
        }
        else if(opname == "F" || opname == "FermiPhase")
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

template<typename ASiteType, typename BSiteType>
class SpecialMixedSiteSet : public SiteSet
    {
    public:

    SpecialMixedSiteSet() { }

    SpecialMixedSiteSet (int N, int bsite, Args const& args=Args::global())
    {
        auto sites = SiteStore(N);
        for(int j = 1; j <= N; ++j)
        {
            if(j == bsite) sites.set(j,BSiteType({args,"SiteNumber=",j}));
            else           sites.set(j,ASiteType({args,"SiteNumber=",j}));
        }
        SiteSet::init(std::move(sites));
    }

    SpecialMixedSiteSet (IndexSet const& is, int bsite)
    {
        int N = is.length();
        auto sites = SiteStore(N);
        for(auto j : range1(N))
        {
            if(j == bsite) sites.set(j,BSiteType(is(j)));
            else           sites.set(j,ASiteType(is(j)));
        }
        SiteSet::init(std::move(sites));
    }
/*
    void read (std::istream& s)
    {
        int N = itensor::read<int>(s);
        if(N > 0)
        {
            auto store = SiteStore(N);
            for(int j = 1; j <= N; ++j) 
            {
                auto I = Index{};
                I.read(s);
                if(j == bsite) store.set(j,BSiteType(I));
                else           store.set(j,ASiteType(I));
            }
            init(std::move(store));
        }
    }
*/
};

using MixedBasis = SpecialMixedSiteSet <FermionSite, SpecialBosonSite>;

#endif
