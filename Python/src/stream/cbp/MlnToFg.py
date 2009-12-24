from libSTREAMWrapper import *
from MLN import *
from FOL import *
import math
import operator
from os import remove

class MlnToFg:
	
	perms = []
	cfg = None

	def truthTableForGndFormula (self, gf, weight):
		
		bstr = lambda n: (n > 0 and str(n & 1) + bstr(n >> 1).rstrip('0')) or '0'		    
		atoms = gf.getGroundAtoms()  
		atoms.sort(key=operator.attrgetter("idx")) 
		worldValues = {}
		truthtable = []
		comb = 2 ** len(atoms)
		for actual in range(comb):
			truthValues = bstr(actual)
			truthValues = truthValues.ljust(len(atoms), '0')

			worldValues.clear()
			for i, ga in enumerate(atoms):
				worldValues[ga.idx] = bool(int(truthValues[i]))
			if gf.isTrue(worldValues):
				truthtable.append(math.exp(weight))
			else:
				truthtable.append(1)

		return truthtable

	def mlnToFactorGraph (self, mln, cnf=False):
		vars = [];
		allVarSet = set()

		for i in range(len(mln.gndAtomsByIdx)):
			vars.append(Var(i, 2))
			allVarSet.add(i)

		facs = vector_factor();
		
		for gf in mln.gndFormulas:
			
			if cnf:
				gndCnf = gf.toCNF()
				for atom in gndCnf.getGroundAtoms():
					if atom.idx in allVarSet:
						allVarSet.remove(atom.idx)
					
				if len(gndCnf.getGroundAtoms()) > 0:
					
					if isinstance(gndCnf, Conjunction):
						for disjunction in gndCnf.children:
							facs.push_back(self.formualToFactor(disjunction,mln.formulas[gf.idxFormula].weight))
					else:
						facs.push_back(self.formualToFactor(gndCnf,mln.formulas[gf.idxFormula].weight)) 
			else:
				atoms = gf.getGroundAtoms();
				varSet = VarSet()
				idxSet = [];        
				for atom in atoms:
					idxSet.append(atom.idx)     
					varSet |= vars[atom.idx]
					if atom.idx in allVarSet:
						allVarSet.remove(atom.idx)  

				# create permuatation for restoring original ordering
				# e.g. a formula with ground atoms 15, 6, 22 corresponds to 1,0,2                              
				sortedIdxSet = list(idxSet);
				sortedIdxSet.sort();
				perm = vector_sizet(len(idxSet))
				for i, idx in enumerate(idxSet):
					perm[i] = sortedIdxSet.index(idx)
		
				factor = Factor(varSet);
	
				truthtable = self.truthTableForGndFormula(gf, mln.formulas[gf.idxFormula].weight)	       
		
				for i, row in enumerate(truthtable):
					factor[i] = row

				self.perms.append(perm)
				
				facs.push_back(factor)

		# hack to have all ground atoms as variables in the factorgraph
		for var in allVarSet:
			varSet = VarSet()
			varSet |= vars[var]
			factor = Factor(varSet)
			factor[0] = 1;
			factor[1] = 1;
			facs.push_back(factor)

		fg = FactorGraph(facs);
		assert len(mln.gndAtoms) == fg.nrVars();
		if not cnf:
			assert len(mln.gndFormulas) == fg.nrFactors();
		return fg

	def formualToFactor(self, gf, weight):
		atoms = gf.getGroundAtoms();
		varSet = VarSet()
		idxSet = [];
		for atom in atoms:
			idxSet.append(atom.idx)
			varSet |= Var(atom.idx,2)

		# create permuatation for restoring original ordering
		# e.g. a formula with ground atoms 15, 6, 22 corresponds to 1,0,2                              
		sortedIdxSet = list(idxSet);
		sortedIdxSet.sort();
		perm = vector_sizet(len(idxSet))
		for i, idx in enumerate(idxSet):
			perm[i] = sortedIdxSet.index(idx)
		
		factor = Factor(varSet);
		
		truthtable = self.truthTableForGndFormula(gf, weight)	       
		
		for i, row in enumerate(truthtable):
			factor[i] = row		

		return factor;

	def writeMlnToFgFile(self, mln, outfilename):	
		f = open(outfilename, 'w')		
		bstr = lambda n: (n>0 and bstr(n>>1).lstrip('0')+str(n&1)) or '0'
		
		f.write("%d\n\n" % len(mln.gndFormulas));
		for gf in mln.gndFormulas:
			currentWeight = mln.formulas[gf.idxFormula].weight
			atoms = gf.getGroundAtoms()
			l = len(atoms)
			f.write("%d\n" % l)
			for ga in atoms:
				f.write('%d ' % ga.idx)
			f.write('\n')
			for i in range(l):
				f.write("2 ")
			f.write("\n%d\n" % 2**l)
			worldValues = {}
			for conf in range(2**l):
				truthValues = bstr(conf).rjust(l, '0')
				for i, ga in enumerate(atoms):
					worldValues[ga.idx] = bool(int(truthValues[l-i-1]))
				wt = exp(currentWeight) if gf.isTrue(worldValues) else 1.0
				f.write('%d %f\n' % (conf, wt))
			f.write('\n')
		f.close()

	def mlnFileToFgFile(self, mlnfilename, outfilename, eviFile=None):
		print "Transforming MLN to grounded FG (in libDAI format with ordered variables)..."
		
		mln = MLN(mlnfilename, verbose=False)

		if eviFile == None:
			f = open('temp.db', 'w')
			f.close()
			evidence = evidence2conjunction(mln.combineDB('temp.db', verbose=True))
			remove('temp.db')	
		else:
			evidence = evidence2conjunction(mln.combineDB(eviFile, verbose=True))

		self.mlnToFgFile(mln, outfilename)
		print "...done"

	def mlnFileToFactorGraph(self, file, eviFile, cnf=False):
		print "Create FactorGraph from MLN file"
		mln = MLN(file)
		
		if eviFile == None:
			f = open('temp.db', 'w')
			f.close()
			evidence = evidence2conjunction(mln.combineDB('temp.db', verbose=True))
			remove('temp.db')	
		else:
			evidence = evidence2conjunction(mln.combineDB(eviFile, verbose=True))
		
		fg = self.mlnToFactorGraph(mln,cnf);
		
		# set evidence
		for i, atom in enumerate(mln.gndAtomsByIdx):
			if mln.evidence[atom] == True:
				fg.clamp(i, 1, False);
			elif mln.evidence[atom] == False:
				fg.clamp(i, 0, False);
		
		return fg;

	def mlnFileToCFactorGraph(self, file, eviFile):
		print "Create CFactorGraph from MLN file"
		mln = MLN(file)
		
		if eviFile == None:
			f = open('temp.db', 'w')
			f.close()
			evidence = evidence2conjunction(mln.combineDB('temp.db', verbose=True))
			remove('temp.db')	
		else:
			evidence = evidence2conjunction(mln.combineDB(eviFile, verbose=True))
						 		              
		return mlnToCFactorGraph(mln);
	
	def mlnToCFactorGraph (self, mln, cnf=True):
		print "Create CFactorGraph from MLN"
		fg = self.mlnToFactorGraph(mln, cnf);
		
		self.cfg = CFactorGraph(fg)
		
		for i, perm in enumerate(self.perms):
			self.cfg.setSigma(i, perm);    
		for i, atom in enumerate(mln.gndAtomsByIdx):        
			if mln.evidence[atom] == True:
				self.cfg.clamp(Var(i,2), 1, False);
			elif mln.evidence[atom] == False:
				self.cfg.clamp(Var(i,2), 0, False);
		return self.cfg;

	def cFactorGraphToFgFile(self, cfg, file):	
		print "Write CFactorGraph to FG file"	
		open(file, 'w').write(cfg.__str__());	

		
