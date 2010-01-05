
#include <string>
#include <ostream>
#include <iostream>
#include <sstream>
#include <dai/factorgraph.h>
#include "STREAM/CBP/CBP.h"
#include "dai/properties.h"

using namespace std;
using namespace stream;
using namespace dai;

int main( int argc, char *argv[] ) {

	string fgString;
	fgString = string("2\n\n") +
			string("2\n1 3\n2 2\n4\n0 0.4\n1 1.4\n2 2.0\n3 1.2\n\n") +
			string("2\n5 3\n2 2\n4\n0 0.4\n1 1.4\n2 2.0\n3 1.2\n");

	istringstream iss (fgString,istringstream::in);

	CFactorGraph cfg;
	iss >> cfg;

	PropertySet opts;
	opts.Set("maxiter",size_t(1000));
	opts.Set("tol",double(1e-9));
	opts.Set("updates",string("PARALL"));
	opts.Set("logdomain",false);
	opts.Set("compressionAlg", string("Position"));
	opts.Set("verbose", size_t(3));
	CBP cbp(cfg, opts);
	cbp.init();
	cbp.run();

	for (size_t i=0; i<cbp.nrVars(); i++) {
		cout << "superVar: " << i << ": ";
		dai::operator <<(cout, cbp.clusterV(i)) << endl;
	}

}
