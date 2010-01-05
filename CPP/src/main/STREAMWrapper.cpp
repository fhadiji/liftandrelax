/*
    Copyright (C) 2009
    Babak Ahmadi [babak dot ahmadi at iais dot fraunhofer dot de]
    Fabian Hadiji [fabian dot hadiji at iais dot fraunhofer dot de]
    Kristian Kersting (coordination) [kristian dot kersting at iais dot fraunhofer dot de]

    STREAM Project at
        Fraunhofer IAIS, Sankt Augustin, Germany, and
        KDML, Unversity of Bonn, Germany

    This file is part of libSTREAM.

    libSTREAM is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    libSTREAM is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, see <http://www.gnu.org/licenses/>.
*/
/*
 * libSTREAMWrapper.cpp
 *
 *  Created on: Jul 13, 2009
 *      Author: Fabian Hadiji
 */

#include <boost/python/class.hpp>
#include <boost/python/module.hpp>
#include <boost/python/operators.hpp>
#include <boost/python/def.hpp>
#include <boost/python/pure_virtual.hpp>
#include <boost/python/errors.hpp>
#include <boost/python/wrapper.hpp>
#include <boost/python/call.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <boost/python/suite/indexing/map_indexing_suite.hpp>
#include <boost/python/docstring_options.hpp>
#include <boost/python/enum.hpp>
#include <iostream>
#include <string>
#include <fstream>
#include <dai/smallset.h>
#include <dai/varset.h>
#include <dai/var.h>
#include <dai/bp.h>
#include <dai/gibbs.h>
#include <STREAM/CBP/CBP.h>
#include <STREAM/CBP/CFactorGraph.h>

using namespace dai;
using namespace std;
using namespace stream;
using namespace boost::python;

void (CBP::*CBPInit1)() = &CBP::init;

void (BP::*BPInit1)() = &BP::init;

void (Gibbs::*GibbsInit1)() = &Gibbs::init;

size_t (Var::*const_label)() const = &Var::label;
size_t (Var::*const_states)() const = &Var::states;

const VarSet& (Factor::*const_vars)() const = &Factor::vars;

const Factor & (FactorGraph::*const_factor)(size_t) const = &FactorGraph::factor;
void (FactorGraph::*clamp_idx)(size_t, size_t, bool) = &FactorGraph::clamp;
const FactorGraph::Neighbors & (FactorGraph::*varNeighbors)(size_t) const = &FactorGraph::nbV;

const CFactor & (CFactorGraph::*const_cfactor)(size_t) const = &CFactorGraph::factor;

PropertySet& propertySet_set (PropertySet &ps, string key, string val) {
	return ps.Set(key,val);
}

void propertySet_read (PropertySet &ps, const char* filename) {
	ifstream propertyFile(filename);
	propertyFile >> ps;
	propertyFile.close();
}

void factor_setitem(Factor& factor, int index, double value)
{
	if (index >= 0 && index < (int)factor.states()) {
		factor[index] = value;
	}
	else {
		PyErr_SetString(PyExc_IndexError, "index out of range");
		boost::python::throw_error_already_set();
	}
}

double factor_getitem(Factor &factor, int index)
{
	if (index >= 0 && index < (int)factor.states()) {
		return factor[index];
	}
	else {
		PyErr_SetString(PyExc_IndexError, "index out of range");
		boost::python::throw_error_already_set();
		return -1;
	}
}

void cfactor_setitem(CFactor& factor, int index, double value)
{
	if (index >= 0 && index < (int)factor.states()) {
		factor[index] = value;
	}
	else {
		PyErr_SetString(PyExc_IndexError, "index out of range");
		boost::python::throw_error_already_set();
	}
}

double cfactor_getitem(CFactor &factor, int index)
{
	if (index >= 0 && index < (int)factor.states()) {
		return factor[index];
	}
	else {
		PyErr_SetString(PyExc_IndexError, "index out of range");
		boost::python::throw_error_already_set();
		return -1;
	}
}

Factor vector_factor_getitem(vector<Factor> &factors, int index)
{
	if (index >= 0 && index < (int)factors.size()) {
		return factors[index];
	}
	else {
		PyErr_SetString(PyExc_IndexError, "index out of range");
		boost::python::throw_error_already_set();
		return -1;
	}
}

size_t stateOfVar (State &state, const Var &var) {
	return state(var);
}

void incrState(State &state) {
	state++;
}

string fgToDotString (FactorGraph &fg) {
	stringstream out;
	fg.printDot(out);
	return out.str();
}

class CompressInterfaceWrap : public CompressInterface, wrapper<CompressInterface> {
	string name() {
		return this->get_override("name")();
	}
};

BOOST_PYTHON_MODULE(libSTREAMWrapper)
{

	docstring_options doc_options(DEMO_DOCSTRING_SHOW_ALL);

	def("calcState", calcState);
	def("calcLinearState", calcLinearState);

	/* std containers */
	class_<map<long, int> >("map_long_int")
		.def(map_indexing_suite<map<long, int> >())
		;

	class_<map<Var, size_t> >("map_var_sizet")
		.def(str(self))
		.def(map_indexing_suite<map<Var, size_t> >())
		;

	class_<vector<size_t> >("vector_sizet", init<size_t>())
		.def(vector_indexing_suite<vector<size_t> >())
		.def("push_back", &vector<size_t>::push_back)
		;

	class_<Var>("Var", init<long, size_t>())
		.def(str(self))
		.def("label", const_label)
		.def("states", const_states)
		;

	class_<vector<Var> >("vector_var")
		.def(str(self))
		.def(vector_indexing_suite<std::vector<Var> >())
		.def("push_back", &vector<Var>::push_back)
		.def("size", &vector<Var>::size)
		;

	class_<SmallSet<Var> >("smallSet_var")
		;

	class_<VarSet, bases<SmallSet<Var> > >("VarSet")
		.def(init<SmallSet<Var> >())
		.def(self |=  Var())
		.def(self & self)
		.def(self | self)
		.def("__iter__", boost::python::iterator<VarSet, return_internal_reference<> >())
		.def(str(self))
		.def("size", &VarSet::size)
		;

	class_<vector<Factor> >("vector_factor")
		.def(str(self))
		.def("__iter__", boost::python::iterator<vector<Factor>, return_internal_reference<> >())
		.def("__getitem__", &vector_factor_getitem)
		.def("push_back", &vector<Factor>::push_back)
		.def("size", &vector<Factor>::size)
		;

	class_<Factor>("Factor", init<VarSet>())
		.def(init<VarSet, double>())
		.def("__getitem__", &factor_getitem)
		.def("__setitem__", &factor_setitem)
		.def(str(self))
		.def("vars", const_vars, return_internal_reference<>())
		.def("states", &Factor::states)
		;

	class_<FactorGraph>("FactorGraph")
		.def(init<vector<Factor> >())
		.def(str(self))
		.def("ReadFromFile", &FactorGraph::ReadFromFile)
		.def("nrEdges", &FactorGraph::nrEdges)
		.def("nrFactors", &FactorGraph::nrFactors)
		.def("nrVars", &FactorGraph::nrVars)
		.def("factor", const_factor, return_internal_reference<>())
		.def("factors", &FactorGraph::factors, return_internal_reference<>())
		.def("var", &FactorGraph::var, return_internal_reference<>())
		.def("vars", &FactorGraph::vars, return_internal_reference<>())
		.def("findVar", &FactorGraph::findVar)
		.def("nbV", varNeighbors, return_internal_reference<>())
		.def("clamp", clamp_idx)
		.def("fgToDotString",&fgToDotString)
		;

	class_<BipartiteGraph::Neighbor>("Neighbor")
		.def_readonly("dual", &BipartiteGraph::Neighbor::dual)
		.def_readonly("iter", &BipartiteGraph::Neighbor::iter)
		;

	class_<BipartiteGraph::Neighbors>("Neighbors")
		.def(vector_indexing_suite<BipartiteGraph::Neighbors>())
		;

	class_<BP, bases<FactorGraph> >("BP", init<FactorGraph, PropertySet>())
		.def("init", BPInit1)
		.def("run", &BP::run)
		.def("beliefV", &BP::beliefV)
		.def("beliefF", &BP::beliefF)
		.def("getIterations", &BP::Iterations)
		;

	class_<Gibbs>("Gibbs", init<FactorGraph, PropertySet>())
		.def("init", GibbsInit1)
		.def("run", &Gibbs::run)
		.def("beliefV", &Gibbs::beliefV)
		.def("beliefF", &Gibbs::beliefF)
		.def("getIterations", &Gibbs::Iterations)
		;


	class_<PropertySet>("PropertySet")
		.def(str(self))
		.def("Set", propertySet_set, return_internal_reference<>())
		.def("read", propertySet_read)
		;

	class_<State>("State", init<VarSet>())
		.def("stateOfVar", stateOfVar)
		.def("incr", incrState)
		;

	class_<vector<CFactor> >("vector_cfactor")
		.def(vector_indexing_suite<std::vector<CFactor> >())
		.def("push_back", &vector<CFactor>::push_back)
		;

	class_<CFactor>("CFactor", init<VarSet>())
		.def(init<VarSet, double>())
		.def("__getitem__", &cfactor_getitem)
		.def("__setitem__", &cfactor_setitem)
		.def(str(self))
		.def("vars", &CFactor::vars, return_internal_reference<>())
		.def("states", &CFactor::states)
		;

    class_<CFactorGraph>("CFactorGraph")
		.def(init<vector<CFactor> >())
		.def(init<FactorGraph>())
		.def(str(self))
		.def("readFromFile", &CFactorGraph::readFromFile)
		.def("nrEdges", &CFactorGraph::nrEdges)
		.def("nrFactors", &CFactorGraph::nrFactors)
		.def("nrVars", &CFactorGraph::nrVars)
		.def("factor", const_cfactor, return_internal_reference<>())
		.def("factors", &CFactorGraph::factors, return_internal_reference<>())
		.def("var", &CFactorGraph::var, return_internal_reference<>())
		.def("vars", &CFactorGraph::vars, return_internal_reference<>())
		.def("findVar", &CFactorGraph::findVar)
		.def("clamp", &CFactorGraph::clamp)
		/* CFactorGraph specific functions */
		.def("setSigma", &CFactorGraph::setSigma)
		.def("createCleanedGraph", &CFactorGraph::createCleanedGraph)
		.staticmethod("createCleanedGraph")
		;

	class_<CBP, bases<CFactorGraph> >("CBP", init<CFactorGraph, PropertySet>())
		.def(init<FactorGraph, PropertySet>())
		.def("init", CBPInit1)
		.def("run", &CBP::run)
		.def("maxDiff",&CBP::maxDiff)
		.def("beliefV", &CBP::beliefV)
		.def("beliefF", &CBP::beliefF)
		.def("getIterations", &CBP::Iterations)
		/* CFactorGraph specific functions*/
		.def("reprV", &CBP::reprV)
		.def("createVarMapping", &CBP::createVarMapping)
		;

	class_<CompressInterface>("CompressInterface", no_init)
		.def("getVarColorVec", &CompressInterface::getVarColorVec)
	        ;

	class_<PositionCompress, bases<CompressInterface> >("PositionCompress")
		.def("setCfg", &PositionCompress::setCfg)
		.def("init", &PositionCompress::init)
		.def("iterate", &PositionCompress::iterate)
		.def("createCFactorGraph", &PositionCompress::createCFactorGraph)
		.def("hasConverged", &PositionCompress::hasConverged)
		;


}
