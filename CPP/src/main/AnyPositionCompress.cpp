/*
 * AnyPositionCompress.cpp
 *
 *  Created on: Aug 14, 2009
 *      Author: stream
 */

#include "STREAM/CBP/Compress.h"
#include <stdio.h>

namespace stream {

using namespace std;
using namespace dai;

size_t AnyPositionCompress::facColor(size_t i) {
	for (size_t j=0; j< _cfg.nbF(i).size(); j++) {
		_facInbox[i][j] = _varSigs[_cfg.nbF(i,_cfg.getSigma(i)[j])];
	}
	_facInbox[i][_facInbox[i].size() - 1] = _facSigs[i];
	sort(_facInbox[i].begin(), _facInbox[i].end());

	return hashVector(_facInbox[i]);
}


size_t AnyPositionCompress::varColor(size_t i) {
	foreach (const dai::BipartiteGraph::Neighbor tmpFac, _cfg.nbV(i)) {
		vector<size_t> hashAndPos(2);

		hashAndPos[0] = _facSigs[tmpFac];
		hashAndPos[1] = 0;
		_varInbox[i][tmpFac.iter] = hashVector(hashAndPos);
	}
	_varInbox[i][_varInbox[i].size() - 1] = _varSigs[i];
	sort(_varInbox[i].begin(), _varInbox[i].end());

	return hashVector(_varInbox[i]);
}

void AnyPositionCompress::initFacColors ()  {
	_facSigs = vector<Signature>(_cfg.nrFactors());
	stringstream buf;
	boost::hash<std::string> string_hash;
	map< Var, size_t > states;
	vector<size_t> iDims;
	vector<size_t> iSigma;
	for (size_t i=0; i<_cfg.nrFactors(); i++) {

		iDims.clear();
		for (vector<Var>::const_iterator iter = _cfg.factor(i).vars().begin(); iter != _cfg.factor(i).vars().end(); iter++ ) {
			iDims.push_back( iter->states() );
		}

		bool allTrue = false;
		for (size_t j=0; j<i; j++)  {
			// muss set _facColorVec[i] properly
			if (_cfg.factor(i).states() != _cfg.factor(j).states()) {
				continue;
			}

			iSigma.clear();
			for (size_t pos=0; pos<_cfg.factor(i).vars().size(); pos++) {
				iSigma.push_back(pos);
			}
			vector<vector<size_t> > perms;
			permute(iSigma, 0, iSigma.size(), perms);

			for (size_t k=0; k<perms.size(); k++) {
				allTrue = true;
				Permute perm (iDims, perms[k]);
				for (size_t l=0; l<_cfg.factor(i).states(); l++) {
					if (_cfg.factor(j)[l] != _cfg.factor(i)[perm.convertLinearIndex(l)]) {
						allTrue = false;
						break;
					}
				}

				if (allTrue) {
					break;
				}
			}
			if (allTrue) {
				_facSigs[i] = _facSigs[j];
				break;
			}
		}

		if (!allTrue) {
			vector<double> tmpVec(_cfg.factor(i).states());
			for (size_t j=0; j<_cfg.factor(i).states(); j++) {
				tmpVec[j] = _cfg.factor(i)[j];
			}

			buf.str("");
			for (size_t j=0; j<_cfg.factor(i).states(); j++) {
				buf << tmpVec[j] << "|";
			}
			_facSigs[i] = string_hash(buf.str());
		}
		_facInbox.push_back(std::vector<size_t>(_cfg.nbF(i).size()+ 1) );
	}
}

void AnyPositionCompress::print(const vector<size_t> & v, const int size) {

	for (int i = 0; i < size; i++) {
		printf("%4d", v[i] );
	}
    printf("\n");
}



void AnyPositionCompress::permute(vector<size_t> & v, const int start, const int n, vector<vector<size_t> > & perms) {
	if (start == n-1) {
		perms.push_back(vector<size_t>(v));
	} else {
		for (int i = start; i < n; i++) {
			int tmp = v[i];

			v[i] = v[start];
			v[start] = tmp;
			permute(v, start+1, n, perms);
			v[start] = v[i];
			v[i] = tmp;
		}
	}
}

std::vector<std::map<size_t, int> > AnyPositionCompress::createCounts(size_t &gndFactor, VarSet &superVarSet) {
	// create zero entries for each position
	map<long, map<size_t, int> > countMap;
	foreach (const dai::BipartiteGraph::Neighbor &tmpVar, _cfg.nbF(gndFactor)) {
		Var liftedVar = _varRepr[_varColorVec[tmpVar]];
		size_t pos = find(_cfg.factor(gndFactor).sigma().begin(), _cfg.factor(gndFactor).sigma().end(), tmpVar.iter) - _cfg.factor(gndFactor).sigma().begin();
		countMap[liftedVar.label()][pos] = 0;
	}

	vector<map<size_t, int> > counts;
	for (vector<Var>::const_iterator iter = superVarSet.begin(); iter < superVarSet.end(); iter++) {
		size_t count = 0;
		foreach(const dai::BipartiteGraph::Neighbor gndFac, _cfg.nbV(_cfg.findVar(*iter))) {
			if (_facRepr[_facColorVec[gndFac]] == gndFactor) {
				count++;
			}
		}

		map<size_t,int>::iterator posIter = countMap[iter->label()].begin();
		posIter->second = count;

//		for (map<size_t, int>::iterator posIter=countMap[iter->label()].begin();posIter != countMap[iter->label()].end();posIter++) {
//			posIter->second = count / countMap[iter->label()].size();
//		}
		counts.push_back(countMap[iter->label()]);
	}

	return counts;
}

}
