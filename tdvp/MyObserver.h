#ifndef __MYOBSERVER_H_CMC__
#define __MYOBSERVER_H_CMC__
#include <iomanip>
#include <map>
#include "itensor/all.h"
#include "Entanglement.h"
#include "MixedBasis.h"
#include "Corr.h"
#include "ContainerUtility.h"
using namespace vectool;
using namespace iutility;

class MyObserver : public DMRGObserver
{
    public:
        MyObserver (const Fermion& sites, const MPS& psi, const Args& args=Args::global())
        : DMRGObserver (psi, args)
        , _sites (sites)
        , _ns (length(psi)+1,0.)
        , _Npar (0.)
        , _specs (length(psi))
        {
            _write = args.getBool ("Write",false);
            _out_dir = args.getString("out_dir",".");
        }

        void measure (const Args& args);

        Real Npar () const { return _Npar; }
        const Spectrum& spec (int i) const { return _specs.at(i); }

    private:
        bool                    _write;
        string                  _out_dir;	// empty string "" if not write
        Fermion                 _sites;

        // Observables
        vector<Real>        _ns;
        Real                _Npar;
        vector<Spectrum>    _specs;
};

inline Real Onsite_mea (const ITensor& A, const ITensor& op)
{
    ITensor re = A * op;
    re.noPrime ("Site");
    re *= dag(A);
    return toReal (re);
}

void MyObserver :: measure (const Args& args)
{
    DMRGObserver::measure (args);

    cout << scientific << setprecision(14);
    // Define your measurements below
    // Call psi() to access the MPS
    //
    auto N = length(psi());
    auto b = args.getInt("AtBond");
    auto sw = args.getInt("Sweep");
    auto ha = args.getInt("HalfSweep");
    auto energy = args.getReal("Energy",0);

    if (b != N)
        _specs.at(b) = spectrum();

    int oc = orthoCenter(psi());
    int nc = args.getInt("NumCenter");
    // measure during the second half of sweep
    if ((nc == 2 && oc == N) || ha == 2)
    {
        // Density
        ITensor n_op = _sites.op("N",oc);
        Real ni = Onsite_mea (psi().A(oc), n_op);
        cout << "\tn " << oc << " = " << ni << endl;

        // Entanglement entropy
        Real S = EntangEntropy (spectrum());
        cout << "\tentang entropy " << oc << " = " << S << endl;
    }

    // At the end of a sweep
    if (oc == 1 && ha == 2 && b == 1)
    {
        for(int i = 1; i < N; i++)
        {
            cout << "\tbond dim " << i << " = " << dim(rightLinkIndex (psi(), i)) << endl;
        }

        if (_write)
        {
            cout << "write MPS" << endl;
            writeToFile (_out_dir+"/psi.mps", psi());
        }
    }
}
#endif
