/*
 * AnyPositionCnfCompress.cpp
 *
 *  Created on: Oct 20, 2009
 *      Author: stream
 */

#include "STREAM/CBP/Compress.h"

namespace stream {

using namespace std;
using namespace dai;

void AnyPositionCnfCompress::initFacColors() {
	_facSigs = vector<Signature>(_cfg.nrFactors());
	_zeroStates = vector<size_t>(_cfg.nrFactors());
	boost::hash<vector<double> > vector_hash;
	vector<size_t> dims;

	for (size_t i=0; i<_cfg.nrFactors(); i++) {
		CompressInterface::rearrangeCnfFactor(_cfg.factor(i));

		dims.clear();
		dims.reserve(_cfg.factor(i).vars().size());
		for (vector<Var>::const_iterator iter = _cfg.factor(i).vars().begin(); iter != _cfg.factor(i).vars().end(); iter++ ) {
			dims.push_back( iter->states() );
		}

		Permute perm (dims, _cfg.factor(i).sigma());
		vector<double> potential(_cfg.factor(i).states());
		for (size_t j=0; j<_cfg.factor(i).states(); j++) {
			potential[perm.convertLinearIndex(j)] = _cfg.factor(i)[j];
		}

		for (size_t j=0; j<_cfg.factor(i).states(); j++) {
			if (_cfg.factor(i)[j] == 1) {
				_zeroStates[i] = perm.convertLinearIndex(j);
			}
		}

		_facSigs[i] = vector_hash(potential);
		_facInbox.push_back(std::vector<size_t>(_cfg.nbF(i).size() + 1));
	}
}

size_t AnyPositionCnfCompress::facColor(size_t i) {
	for (size_t j=0; j< _cfg.nbF(i).size(); j++) {
		_facInbox[i][j] = _varSigs[_cfg.nbF(i,_cfg.getSigma(i)[j])];
	}
	_facInbox[i][_facInbox[i].size() - 1] = _facSigs[i];

	double res = log(_cfg.factor(i).states() - _zeroStates[i]) /  log(2);
	size_t nrPosLiterals = size_t(res);
	sort(_facInbox[i].begin(), _facInbox[i].begin() + nrPosLiterals);
	sort(_facInbox[i].begin() + nrPosLiterals, _facInbox[i].end() - 1);

	return hashVector(_facInbox[i]);
}

std::vector<std::map<size_t, int> > AnyPositionCnfCompress::createCounts(size_t &gndFactor, VarSet &superVarSet) {
	// create zero entries for each position
	map<long, map<size_t, int> > countMap;
	foreach (const dai::BipartiteGraph::Neighbor &tmpVar, _cfg.nbF(gndFactor)) {
		Var liftedVar = _varRepr[_varColorVec[tmpVar]];
		size_t pos = find(_cfg.factor(gndFactor).sigma().begin(), _cfg.factor(gndFactor).sigma().end(), tmpVar.iter) - _cfg.factor(gndFactor).sigma().begin();
		countMap[liftedVar.label()][pos] = 0;
	}

	vector<map<size_t, int> > counts;
	size_t posCount;
	size_t negCount;
	for (vector<Var>::const_iterator iter = superVarSet.begin(); iter < superVarSet.end(); iter++) {
		posCount = 0;
		negCount = 0;

		foreach(const dai::BipartiteGraph::Neighbor tmpFac, _cfg.nbV(_cfg.findVar(*iter))) {
			if (_facRepr[_facColorVec[tmpFac]] == gndFactor) {
				size_t pos = find(_cfg.factor(tmpFac).sigma().begin(), _cfg.factor(tmpFac).sigma().end(), tmpFac.dual) - _cfg.factor(tmpFac).sigma().begin();
				double res = log(_cfg.factor(tmpFac).states() - _zeroStates[tmpFac]) /  log(2);
				size_t nrPosLiterals = size_t(res);
				bool sign = (pos < nrPosLiterals);
				if (sign) {
					posCount++;
				} else {
					negCount++;
				}
			}
		}
		map<size_t, int>::iterator posIter=countMap[iter->label()].begin();

		if (posCount > 0) {
			posIter->second = posCount;
		}

		if (negCount > 0) {
			for (size_t j=0;j<posCount; j++, posIter++) {}
			posIter->second = negCount;
		}

		counts.push_back(countMap[iter->label()]);
	}

	return counts;
}
}
