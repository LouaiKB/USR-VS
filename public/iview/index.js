/*!
 * iview is an interactive WebGL visualizer for protein-ligand complex. iview is based on GLmol, three.js, zlib.js and jQuery.
 * http://github.com/HongjianLi/istar
 * Copyright (c) 2012-2014 Chinese University of Hong Kong
 * License: Apache License 2.0
 * Hongjian Li, Kwong-Sak Leung, Takanori Nakane and Man-Hon Wong.
 * iview: an interactive WebGL visualizer for protein-ligand complex.
 * BMC Bioinformatics, 15(1):56, 2014.
 *
 * GLmol
 * https://github.com/biochem-fan/GLmol
 * Copyright (c) 2011-2012 biochem_fan
 * License: dual license of MIT or LGPL3
 *
 * three.js
 * https://github.com/mrdoob/three.js
 * Copyright (c) 2010-2012 three.js Authors. All rights reserved.
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 *
 * zlib.js
 * https://github.com/imaya/zlib.js
 * Copyright (c) 2012 imaya
 * License: MIT License
 *
 * jQuery
 * http://jquery.org
 * Copyright (c) 2011 John Resig
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to
 * the following conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
 * LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 * OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 * WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

$(function () {
	var vdwRadii = { // Hu, S.Z.; Zhou, Z.H.; Tsai, K.R. Acta Phys.-Chim. Sin., 2003, 19:1073.
		 H: 1.08,
		HE: 1.34,
		LI: 1.75,
		BE: 2.05,
		 B: 1.47,
		 C: 1.49,
		 N: 1.41,
		 O: 1.40,
		 F: 1.39,
		NE: 1.68,
		NA: 1.84,
		MG: 2.05,
		AL: 2.11,
		SI: 2.07,
		 P: 1.92,
		 S: 1.82,
		CL: 1.83,
		AR: 1.93,
		 K: 2.05,
		CA: 2.21,
		SC: 2.16,
		TI: 1.87,
		 V: 1.79,
		CR: 1.89,
		MN: 1.97,
		FE: 1.94,
		CO: 1.92,
		NI: 1.84,
		CU: 1.86,
		ZN: 2.10,
		GA: 2.08,
		GE: 2.15,
		AS: 2.06,
		SE: 1.93,
		BR: 1.98,
		KR: 2.12,
		RB: 2.16,
		SR: 2.24,
		 Y: 2.19,
		ZR: 1.86,
		NB: 2.07,
		MO: 2.09,
		TC: 2.09,
		RU: 2.07,
		RH: 1.95,
		PD: 2.02,
		AG: 2.03,
		CD: 2.30,
		IN: 2.36,
		SN: 2.33,
		SB: 2.25,
		TE: 2.23,
		 I: 2.23,
		XE: 2.21,
		CS: 2.22,
		BA: 2.51,
		LA: 2.40,
		CE: 2.35,
		PR: 2.39,
		ND: 2.29,
		PM: 2.36,
		SM: 2.29,
		EU: 2.33,
		GD: 2.37,
		TB: 2.21,
		DY: 2.29,
		HO: 2.16,
		ER: 2.35,
		TM: 2.27,
		YB: 2.42,
		LU: 2.21,
		HF: 2.12,
		TA: 2.17,
		 W: 2.10,
		RE: 2.17,
		OS: 2.16,
		IR: 2.02,
		PT: 2.09,
		AU: 2.17,
		HG: 2.09,
		TL: 2.35,
		PB: 2.32,
		BI: 2.43,
		PO: 2.29,
		AT: 2.36,
		RN: 2.43,
		FR: 2.56,
		RA: 2.43,
		AC: 2.60,
		TH: 2.37,
		PA: 2.43,
		 U: 2.40,
		NP: 2.21,
		PU: 2.56,
		AM: 2.56,
		CM: 2.56,
		BK: 2.56,
		CF: 2.56,
		ES: 2.56,
		FM: 2.56,
	};
	var covalentRadii = { // http://en.wikipedia.org/wiki/Covalent_radius
		 H: 0.31,
		HE: 0.28,
		LI: 1.28,
		BE: 0.96,
		 B: 0.84,
		 C: 0.76,
		 N: 0.71,
		 O: 0.66,
		 F: 0.57,
		NE: 0.58,
		NA: 1.66,
		MG: 1.41,
		AL: 1.21,
		SI: 1.11,
		 P: 1.07,
		 S: 1.05,
		CL: 1.02,
		AR: 1.06,
		 K: 2.03,
		CA: 1.76,
		SC: 1.70,
		TI: 1.60,
		 V: 1.53,
		CR: 1.39,
		MN: 1.39,
		FE: 1.32,
		CO: 1.26,
		NI: 1.24,
		CU: 1.32,
		ZN: 1.22,
		GA: 1.22,
		GE: 1.20,
		AS: 1.19,
		SE: 1.20,
		BR: 1.20,
		KR: 1.16,
		RB: 2.20,
		SR: 1.95,
		 Y: 1.90,
		ZR: 1.75,
		NB: 1.64,
		MO: 1.54,
		TC: 1.47,
		RU: 1.46,
		RH: 1.42,
		PD: 1.39,
		AG: 1.45,
		CD: 1.44,
		IN: 1.42,
		SN: 1.39,
		SB: 1.39,
		TE: 1.38,
		 I: 1.39,
		XE: 1.40,
		CS: 2.44,
		BA: 2.15,
		LA: 2.07,
		CE: 2.04,
		PR: 2.03,
		ND: 2.01,
		PM: 1.99,
		SM: 1.98,
		EU: 1.98,
		GD: 1.96,
		TB: 1.94,
		DY: 1.92,
		HO: 1.92,
		ER: 1.89,
		TM: 1.90,
		YB: 1.87,
		LU: 1.87,
		HF: 1.75,
		TA: 1.70,
		 W: 1.62,
		RE: 1.51,
		OS: 1.44,
		IR: 1.41,
		PT: 1.36,
		AU: 1.36,
		HG: 1.32,
		TL: 1.45,
		PB: 1.46,
		BI: 1.48,
		PO: 1.40,
		AT: 1.50,
		RN: 1.50,
		FR: 2.60,
		RA: 2.21,
		AC: 2.15,
		TH: 2.06,
		PA: 2.00,
		 U: 1.96,
		NP: 1.90,
		PU: 1.87,
		AM: 1.80,
		CM: 1.69,
	};
	var atomColors = { // http://jmol.sourceforge.net/jscolors
		 H: new THREE.Color(0xFFFFFF),
		HE: new THREE.Color(0xD9FFFF),
		LI: new THREE.Color(0xCC80FF),
		BE: new THREE.Color(0xC2FF00),
		 B: new THREE.Color(0xFFB5B5),
		 C: new THREE.Color(0x909090),
		 N: new THREE.Color(0x3050F8),
		 O: new THREE.Color(0xFF0D0D),
		 F: new THREE.Color(0x90E050),
		NE: new THREE.Color(0xB3E3F5),
		NA: new THREE.Color(0xAB5CF2),
		MG: new THREE.Color(0x8AFF00),
		AL: new THREE.Color(0xBFA6A6),
		SI: new THREE.Color(0xF0C8A0),
		 P: new THREE.Color(0xFF8000),
		 S: new THREE.Color(0xFFFF30),
		CL: new THREE.Color(0x1FF01F),
		AR: new THREE.Color(0x80D1E3),
		 K: new THREE.Color(0x8F40D4),
		CA: new THREE.Color(0x3DFF00),
		SC: new THREE.Color(0xE6E6E6),
		TI: new THREE.Color(0xBFC2C7),
		 V: new THREE.Color(0xA6A6AB),
		CR: new THREE.Color(0x8A99C7),
		MN: new THREE.Color(0x9C7AC7),
		FE: new THREE.Color(0xE06633),
		CO: new THREE.Color(0xF090A0),
		NI: new THREE.Color(0x50D050),
		CU: new THREE.Color(0xC88033),
		ZN: new THREE.Color(0x7D80B0),
		GA: new THREE.Color(0xC28F8F),
		GE: new THREE.Color(0x668F8F),
		AS: new THREE.Color(0xBD80E3),
		SE: new THREE.Color(0xFFA100),
		BR: new THREE.Color(0xA62929),
		KR: new THREE.Color(0x5CB8D1),
		RB: new THREE.Color(0x702EB0),
		SR: new THREE.Color(0x00FF00),
		 Y: new THREE.Color(0x94FFFF),
		ZR: new THREE.Color(0x94E0E0),
		NB: new THREE.Color(0x73C2C9),
		MO: new THREE.Color(0x54B5B5),
		TC: new THREE.Color(0x3B9E9E),
		RU: new THREE.Color(0x248F8F),
		RH: new THREE.Color(0x0A7D8C),
		PD: new THREE.Color(0x006985),
		AG: new THREE.Color(0xC0C0C0),
		CD: new THREE.Color(0xFFD98F),
		IN: new THREE.Color(0xA67573),
		SN: new THREE.Color(0x668080),
		SB: new THREE.Color(0x9E63B5),
		TE: new THREE.Color(0xD47A00),
		 I: new THREE.Color(0x940094),
		XE: new THREE.Color(0x429EB0),
		CS: new THREE.Color(0x57178F),
		BA: new THREE.Color(0x00C900),
		LA: new THREE.Color(0x70D4FF),
		CE: new THREE.Color(0xFFFFC7),
		PR: new THREE.Color(0xD9FFC7),
		ND: new THREE.Color(0xC7FFC7),
		PM: new THREE.Color(0xA3FFC7),
		SM: new THREE.Color(0x8FFFC7),
		EU: new THREE.Color(0x61FFC7),
		GD: new THREE.Color(0x45FFC7),
		TB: new THREE.Color(0x30FFC7),
		DY: new THREE.Color(0x1FFFC7),
		HO: new THREE.Color(0x00FF9C),
		ER: new THREE.Color(0x00E675),
		TM: new THREE.Color(0x00D452),
		YB: new THREE.Color(0x00BF38),
		LU: new THREE.Color(0x00AB24),
		HF: new THREE.Color(0x4DC2FF),
		TA: new THREE.Color(0x4DA6FF),
		 W: new THREE.Color(0x2194D6),
		RE: new THREE.Color(0x267DAB),
		OS: new THREE.Color(0x266696),
		IR: new THREE.Color(0x175487),
		PT: new THREE.Color(0xD0D0E0),
		AU: new THREE.Color(0xFFD123),
		HG: new THREE.Color(0xB8B8D0),
		TL: new THREE.Color(0xA6544D),
		PB: new THREE.Color(0x575961),
		BI: new THREE.Color(0x9E4FB5),
		PO: new THREE.Color(0xAB5C00),
		AT: new THREE.Color(0x754F45),
		RN: new THREE.Color(0x428296),
		FR: new THREE.Color(0x420066),
		RA: new THREE.Color(0x007D00),
		AC: new THREE.Color(0x70ABFA),
		TH: new THREE.Color(0x00BAFF),
		PA: new THREE.Color(0x00A1FF),
		 U: new THREE.Color(0x008FFF),
		NP: new THREE.Color(0x0080FF),
		PU: new THREE.Color(0x006BFF),
		AM: new THREE.Color(0x545CF2),
		CM: new THREE.Color(0x785CE3),
		BK: new THREE.Color(0x8A4FE3),
		CF: new THREE.Color(0xA136D4),
		ES: new THREE.Color(0xB31FD4),
		FM: new THREE.Color(0xB31FBA),
	};
	var defaultAtomColor = new THREE.Color(0xCCCCCC);
	var defaultBoxColor  = new THREE.Color(0x1FF01F);
	var defaultBondColor = new THREE.Color(0x2194D6);
	var defaultBackgroundColor = new THREE.Color(0x000000);
	var sphereGeometry = new THREE.SphereGeometry(1, 64, 64);
	var cylinderGeometry = new THREE.CylinderGeometry(1, 1, 1, 64, 1);
	var sphereRadius = 1.5;
	var cylinderRadius = 0.3;
	var linewidth = 2;
	var catalogs = {
		'ACB Blocks': 'http://www.acbblocks.com',
		'Acorn PharmaTech': 'http://www.acornpharmatech.com',
		'Acros Organics': 'http://www.acros.be',
		'Active BioPharma': 'http://www.activebiopharma.com',
		'Adesis': 'http://www.adesisinc.com',
		'AF ChemPharm': 'http://www.afchempharm.co.uk',
		'AK Scientific': 'http://www.aksci.com',
		'Aldrich CPR': 'http://www.sigmaaldrich.com',
		'Alfa-Aesar': 'http://www.alfa.com',
		'Amadis Chemical': 'http://www.amadischem.com',
		'Ambinter': 'http://www.ambinter.com',
		'Ambinter Natural Products': 'http://www.ambinter.com',
		'American Custom Chemicals Corp.': 'http://www.acccorporation.com',
		'Amidohydrolase AH-EFI': 'http://www.enzymefunction.org',
		'Amino acid derivatives (EFI)': 'http://www.enzymefunction.org',
		'AmpC ligands': 'http://shoichetlab.compbio.ucsf.edu',
		'AmpC non-binders': 'http://shoichetlab.compbio.ucsf.edu',
		'AnalytiCon Discovery Natural Derivatives': 'http://www.ac-discovery.com',
		'AnalytiCon Discovery NP': 'http://www.ac-discovery.com',
		'AnalytiCon Discovery NP BB': 'http://www.ac-discovery.com',
		'Angene Building Blocks': 'http://an-gene.com',
		'Anward': 'http://www.anward.com',
		'Apeiron Synthesis': 'http://www.apeiron-synthesis.com',
		'Apexmol Building Blocks': 'http://www.apexmol.com',
		'APIChem': 'http://www.apichemistry.com',
		'Apollo Scientific': 'http://www.apolloscientific.co.uk',
		'Apollo Scientific Bioactives': 'http://www.apolloscientific.co.uk',
		'Ark Pharm Building Blocks': 'http://www.arkpharminc.com',
		'Aromsyn': 'http://www.aromsyn.com',
		'Aronis': 'http://www.aronis.ru',
		'Aronis (Make on Request)': 'http://www.aronis.ru',
		'Aronis BB Make-on-demand': 'http://www.aronis.ru',
		'Aronis BuildingBlocks': 'http://www.aronis.ru',
		'Asinex': 'http://www.asinex.com',
		'Asinex Building Blocks': 'http://www.asinex.com',
		'AsisChem': 'http://www.asischem.com',
		'AsisChem Building Blocks': 'http://www.asischem.com',
		'Bachem': 'http://www.bachem.com',
		'Beijing Advanced Technology': 'http://www.chemkingdom.com',
		'BePharm Building Blocks': 'http://www.bepharm.com',
		'BindingDB.org': 'http://www.bindingdb.org',
		'BioBlocks': 'http://www.bioblocks.com',
		'BioSynth': 'http://www.biosynth.ch',
		'Bitter DB': 'http://bitterdb.agri.huji.ac.il/bitterdb/',
		'Bosche Scientific': 'http://www.boschesci.com',
		'BroadPharm': 'http://www.broadpharm.com',
		'Capot Chemical': 'http://www.capotchem.com',
		'Cayman Chemical': 'http://www.caymanchem.com',
		'CCP W191G binders': 'http://shoichetlab.compbio.ucsf.edu',
		'CCP W191G non-binders': 'http://shoichetlab.compbio.ucsf.edu',
		'ChEBI': 'http://www.ebi.ac.uk/chebi/',
		'ChEMBL DrugStore': 'http://www.ebi.ac.uk',
		'ChEMBL12': 'http://www.ebi.ac.uk',
		'ChEMBL12 10uM': 'http://www.ebi.ac.uk',
		'ChEMBL13': 'http://www.ebi.ac.uk',
		'ChEMBL14': 'http://www.ebi.ac.uk',
		'ChEMBL15': 'http://www.ebi.ac.uk',
		'Chembo Pharma': 'http://www.chembopharma.com',
		'ChemBridge': 'http://www.chembridge.com',
		'ChemBridge BuildingBlocks': 'http://www.chembridge.com',
		'ChemDB': 'http://cdb.ics.uci.edu',
		'ChemDiv': 'http://www.chemdiv.com',
		'ChemDiv BuildingBlocks': 'http://www.chemdiv.com',
		'ChemFuture PharmTech': 'http://www.chemfuture.com',
		'Chemical Block': 'http://www.chemblock.com',
		'Chemical Block BB': 'http://www.chemblock.com',
		'Chemik Building Blocks': 'http://www.chemik.com',
		'Chemivate': 'http://www.chemivate.com',
		'ChemMol': 'http://www.chemmol.com',
		'ChiralBlock BioScience BB': 'http://www.chiralblock.com',
		'CiVentiChem': 'http://www.cvchem.com',
		'Collaborative Drug Discovery': 'http://www.collaborativedrug.com',
		'Combi-Blocks': 'http://www.combi-blocks.com',
		'CombiUgi': 'http://usefulchem.wikispaces.com/combiugi',
		'DiPeptides (Hao)': '',
		'DrugBank-approved': 'http://www.drugbank.ca',
		'DrugBank-experimental': 'http://www.drugbank.ca',
		'DrugBank-nutriceuticals': 'http://www.drugbank.ca',
		'DrugBank-Street Drugs': 'http://drugbank.ca',
		'DrugBank-withdrawn': 'http://www.drugbank.ca',
		'E. coli Metabolome Database': 'http://www.ecmdb.ca',
		'EDASA Scientific': 'http://www.edasascientific.com',
		'eMolecules': 'http://www.emolecules.com',
		'Enamine': 'http://www.enamine.net',
		'Enamine BB Make on Demand': 'http://www.enamine.net',
		'Enamine Building Blocks': 'http://www.enamine.net',
		'Enamine-REAL': 'http://www.enamine.net',
		'EndoTherm': 'http://www.endotherm-lsm.com',
		'Enolase EN-EFI': 'http://www.enzymefunction.org',
		'Enolase via KEGG (EFI)': 'http://www.enzymefunction.org',
		'EvoBlocks': 'http://www.evoblocks.com',
		'FDA-approved drugs (via DSSTOX)': 'http://www.epa.gov/nheerl/dsstox/',
		'FineTech': 'http://www.finetechnology-ind.com',
		'Florida Heterocyclic Compounds': 'http://www.ark.chem.ufl.edu',
		'Fluorochem': 'http://www.fluorochem.co.uk',
		'Focus Synthesis BB': 'http://focussynthesis.com',
		'Focus Synthesis BB Make-on-Demand': 'http://focussynthesis.com',
		'Fragmenta': 'http://www.fragmenta.com',
		'Frinton': 'http://frinton.com',
		'Frontier Scientific Services': 'http://www.frontierssi.com',
		'Frontier Scientific Services BB': 'http://www.frontierssi.com',
		'Georganics': 'http://georganics.sk',
		'Glutathione Transferrase GST-EFI': 'http://www.enzymefunction.org',
		'Haloacid dehalogenase HAD-EFI': 'http://www.enzymefunction.org',
		'Herbal Ingredients In-Vivo Metabolism': 'http://58.40.126.120/him/',
		'Herbal Ingredients Targets': 'http://lifecenter.sgst.cn/hit',
		'Human Metabolome Database': 'http://www.hmdb.ca',
		'IBM Patent Data': 'http://www-01.ibm.com/software/data/industry/life-sciences.html',
		'IBScreen': 'http://www.ibscreen.com',
		'IBScreen Bioactives': 'http://www.ibscreen.com',
		'IBScreen BuildingBlocks': 'http://www.ibscreen.com',
		'IBScreen NP': 'http://www.ibscreen.com',
		'Indofine': 'http://www.indofinechemical.com',
		'Indofine Natural Products': 'http://www.indofinechemical.com',
		'Infarmatik': 'http://www.infarmatik.com',
		'Infarmatik (make-on-demand)': 'http://www.infarmatik.com',
		'Inhibitor2': 'http://www.inhibitor2.com',
		'Innovapharm': 'http://www.innovapharm.com.ua',
		'Innovapharm BB Make on Demand': 'http://www.innovapharm.com.ua',
		'Innovapharm Building Blocks': 'http://www.innovapharm.com.ua',
		'Innovapharm Make-on-Demand': 'http://www.innovapharm.com.ua',
		'Isoprenoid synthase IS-EFI': 'http://www.enzymefunction.org',
		'Isoprenoids': 'http://www.isoprenoids.com',
		'IUPHAR Database': 'http://www.iuphar-db.org',
		'J&K Chemical': 'http://www.jk-scientific.com',
		'KaironKem': 'http://www.kaironkem.com',
		'KEGG via PubChem': 'http://pubchem.ncbi.nlm.nih.gov',
		'Key Organics Building Blocks': 'http://www.keyorganics.ltd.uk',
		'KeyOrganics': 'http://www.keyorganics.ltd.uk',
		'KeyOrganics Bioactives': 'http://www.keyorganics.ltd.uk',
		'L99A binders': 'http://shoichetlab.compbio.ucsf.edu',
		'L99A non-binders': 'http://shoichetlab.compbio.ucsf.edu',
		'L99A/M102Q binders': 'http://shoichetlab.compbio.ucsf.edu',
		'L99A/M102Q non-binders': 'http://shoichetlab.compbio.ucsf.edu',
		'Labotest': 'http://www.labotest.com',
		'Labotest Building Blocks': 'http://www.labotest.com',
		'Life Chemicals': 'http://www.lifechemicals.com',
		'Life Chemicals (Virtual)': 'http://www.lifechemicals.com',
		'Life Chemicals BB Make-on-Demand': 'http://www.lifechemicals.com',
		'Life Chemicals Building Blocks': 'http://www.lifechemicals.com',
		'Matrix Scientific': 'http://www.matrixscientific.com',
		'Maybridge': 'http://www.maybridge.com',
		'Maybridge Building Blocks': 'http://www.maybridge.com',
		'Maybridge Hit Finder': 'http://www.maybridge.com',
		'Mcule': 'http://www.mcule.com',
		'MicroCombiChem': 'http://www.microcombichem.com',
		'MicroCombiChem BB': 'http://www.microcombichem.com',
		'MicroCombiChem BB Make-on-demand': 'http://www.microcombichem.com',
		'MicroCombiChem Make-on-demand': 'http://www.microcombichem.com',
		'MicroSource Pharmakon': 'http://www.msdicovery.com',
		'MicroSource Spectrum': 'http://www.msdicovery.com',
		'MicroSource US Drugs': 'http://www.msdicovery.com',
		'MicroSource World Drugs': 'http://www.msdicovery.com',
		'Molcan': 'http://www.molcan.com',
		'MolMall (formerly Molecular Diversity Preservation International)': 'http://www.molmall.net',
		'Molport': 'http://www.molport.com',
		'Molport BB': 'http://www.molport.com',
		'Nagase': 'http://www.nagase-nam.com',
		'NCI Diversity 3': 'http://dtp.nci.nih.gov',
		'NCI Plated 2007': 'http://dtp.nci.nih.gov',
		'NIH Clinical Collection': 'http://www.nihclinicalcollection.com',
		'NIH Clinical Collection via PubChem': 'http://www.nihclinicalcollection.com',
		'Novochemy Building Blocks': 'http://www.novochemy.com',
		'NPACT Database': 'http://crdd.osdd.net/raghava/npact/',
		'NPC (NCGC Pharma)': 'http://tripod.nih.gov/npc',
		'Nubbe Natural Products': 'http://nubbe.iq.unesp.br',
		'Oakwood Chemical': 'http://www.oakwoodchemical.com',
		'Otava': 'http://www.otavachemicals.com',
		'Otava Premium BB': 'http://www.otavachemicals.com',
		'PBMR Labs': 'http://pbmr.com.ua',
		'PBMR Labs Building Blocks': 'http://pbmr.com.ua',
		'PDSP via PubChem': '',
		'Peakdale': 'http://www.peakdale.com',
		'PepTech': 'http://www.peptechcorp.com',
		'Pharmeks': 'http://www.pharmeks.com',
		'Phosphate sugars (EFI)': 'http://www.enzymefunction.org',
		'PKChem': 'http://www.pkchem.ru',
		'PKChem Building Blocks': 'http://www.pkchem.ru',
		'Prestwick Chemical': 'http://www.prestwickchemical.com',
		'Princeton BioMolecular BuildingBlocks': 'http://www.princetonbio.com',
		'Princeton BioMolecular Research': 'http://www.princetonbio.com',
		'Princeton NP': 'http://www.princetonbio.com',
		'Prous via PubChem': '',
		'ProVence': 'http://www.provetech.com',
		'ProVence Building Blocks': 'http://www.provetech.com',
		'PubChem': 'http://pubchem.ncbi.nlm.nih.gov',
		'Rare Chemicals': 'http://www.rarechem.de',
		'Ryan Scientific BB': 'http://www.ryansci.com',
		'Scientific Exchange': 'http://www.htscompounds.com',
		'Scientific Exchange (make on demand)': 'http://www.htscompounds.com',
		'Scientific Exchange Building Blocks': 'http://www.htscompounds.com',
		'Selleck BioChemicals': 'http://www.selleckbio.com',
		'Selleck BioChemicals NP': 'http://www.selleckbio.com',
		'Selleck Chemicals': 'http://www.selleckchem.com',
		'Sequoia Research Products': 'http://www.seqchem.com',
		'ShangHai Biochempartner': 'http://www.biochempartner.com',
		'Shanghai Sinofluoro Scientific': 'http://www.sinofluoro.com',
		'Sigma Aldrich (Building Blocks)': 'http://www.sigmaaldrich.com',
		'SMDC Asinex (building blocks)': '',
		'SMDC CDiv Carboxamide': '',
		'SMDC CDiv Diverse': '',
		'SMDC CDiv Kinase': '',
		'SMDC ChBr Diverse': '',
		'SMDC ChBr Premium': '',
		'SMDC Iconix': '',
		'SMDC Life and Maybridge (building blocks)': '',
		'SMDC Life Kinase': '',
		'SMDC MicroSource': '',
		'SMDC Pharmakon': '',
		'Specs': 'http://www.specs.net',
		'Specs Building Blocks': 'http://www.specs.net',
		'Specs Natural Products': 'http://www.specs.net',
		'Sphinx': 'http://www.sphinxscientificlab.com',
		'Sphinx Make-on-demand': 'http://www.sphinxscientificlab.com',
		'Squarix': 'http://www.squarix.de',
		'StreptomeDB': 'http://www.pharmaceutical-bioinformatics.de/streptomedb/',
		'Synergy Scientific BB': 'http://www.synergy-scientific.com',
		'SynQuest Building Blocks': 'http://www.synquestlabs.com',
		'Synthon-Lab': 'http://www.synthon-lab.com',
		'Synthonix Building Blocks': 'http://www.synthonix.com',
		'TCI': 'http://www.tci.co.uk',
		'TCM Database @ Taiwan': 'http://tcm.cmu.edu.tw',
		'Tetrahedron Building Blocks': 'http://www.thsci.com',
		'TimTec': 'http://www.timtec.com',
		'TimTec BB Make on Demand': 'http://www.timtec.net',
		'TimTec Building Blocks': 'http://www.timtec.com',
		'TimTec Make-on-Demand': 'http://www.timtec.net',
		'TimTec Natural Derivatives': 'http://www.timtec.net',
		'Toronto Research Chemicals': 'http://www.trc-canada.com',
		'Toslab': 'http://www.toslab.com',
		'Toslab Building Blocks': 'http://www.toslab.com',
		'Tractus': 'http://www.tractuschem.com',
		'TTD via PubChem': '',
		'Tyger Building Blocks': 'http://www.tygersci.com',
		'Ubichem': 'http://www.ubichem.com',
		'UEFS Natural Products': 'http://www.uefs.br',
		'UMBBD': 'http://umbbd.msi.umn.edu',
		'UORSY': 'http://www.ukrorgsynth.com',
		'UORSY BB Make-on-demand': 'http://www.ukrorgsynth.com',
		'UORSY Make-on-demand': 'http://www.ukrorgsynth.com',
		'Vitas-M': 'http://www.vitasmlab.com',
		'Vitas-M BB': 'http://www.vitasmlab.com',
		'Wisdom Chemicals': 'http://www.wisdompharma.com',
		'Yu Chen 1': '',
		'Yu Chen 2': '',
		'Yu Chen 3': '',
		'Zelinsky Institute': 'http://zelinsky.com',
		'Zelinsky Institute Building Blocks': 'http://zelinsky.com',
		'Zelinsky Institute Make on Demand': 'http://zelinsky.com',
		'ZereneX': 'http://www.zerenex-molecular.com',
		'ZereneX Building Blocks': 'http://www.zerenex-molecular.com',
		'Zylexa Pharma': 'http://www.zylexa-pharma.com',
		'Zylexa Pharma BB': 'http://www.zylexa-pharma.com',
	};

	var canvas = $('canvas');
	canvas.widthInv  = 1 / canvas.width();
	canvas.heightInv = 1 / canvas.height();
	var renderer = new THREE.WebGLRenderer({
		canvas: canvas.get(0),
		antialias: true,
	});
	renderer.setSize(canvas.width(), canvas.height());
	renderer.setClearColor(defaultBackgroundColor);
	var scene = new THREE.Scene();
	var directionalLight = new THREE.DirectionalLight(0xFFFFFF, 1.2);
	directionalLight.position.set(0.2, 0.2, -1).normalize();
	var ambientLight = new THREE.AmbientLight(0x202020);
	var rot = new THREE.Object3D();
	var mdl = new THREE.Object3D();
	rot.add(mdl);
	scene.add(directionalLight);
	scene.add(ambientLight);
	scene.add(rot);
	scene.fog = new THREE.Fog(defaultBackgroundColor, 100, 200);
	var camera = new THREE.PerspectiveCamera(20, canvas.width() / canvas.height(), 1, 800);
	camera.position.set(0, 0, -150);
	camera.lookAt(new THREE.Vector3(0, 0, 0));
	var sn, sf;
	var labelVertexShader = '\
uniform float width, height;\n\
varying vec2 vUv;\n\
void main()\n\
{\n\
	mat4 mv = modelViewMatrix;\n\
	mv[0][0] = mv[1][1] = mv[2][2] = 1.0;\n\
	mv[0][1] = mv[0][2] = mv[1][0] = mv[1][2] = mv[2][0] =  mv[2][1] = 0.0;\n\
	mat4 mat = projectionMatrix * mv;\n\
	vUv = uv;\n\
	float aspect = projectionMatrix[1][1] / projectionMatrix[0][0];\n\
	gl_Position = mat * vec4(position, 1.0);\n\
	gl_Position /= gl_Position.w;\n\
	gl_Position += vec4(uv.x * width * 1e-3, uv.y * height * aspect * 1e-3, 0.0, 0.0);\n\
	gl_Position.z = -0.9;\n\
}';
	var labelFragmentShader = '\
uniform sampler2D map;\n\
varying vec2 vUv;\n\
void main()\n\
{\n\
	gl_FragColor = texture2D(map, vec2(vUv.x, 1.0 - vUv.y));\n\
	if (gl_FragColor.a < 0.5) discard;\n\
}';
	var labelGeo = new THREE.Geometry();
	for (var i = 0; i < 6; ++i) {
		labelGeo.vertices.push(new THREE.Vector3(0, 0, 0));
	}
	labelGeo.faces.push(new THREE.Face3(0, 1, 2));
	labelGeo.faces.push(new THREE.Face3(0, 2, 3));
	labelGeo.faceVertexUvs[0].push([new THREE.Vector2(0, 0), new THREE.Vector2(1, 1), new THREE.Vector2(0, 1)]);
	labelGeo.faceVertexUvs[0].push([new THREE.Vector2(0, 0), new THREE.Vector2(1, 0), new THREE.Vector2(1, 1)]);

	var createSphere = function (atom, defaultRadius, forceDefault, scale) {
		var mesh = new THREE.Mesh(sphereGeometry, new THREE.MeshLambertMaterial({ color: atom.color }));
		mesh.scale.x = mesh.scale.y = mesh.scale.z = forceDefault ? defaultRadius : (vdwRadii[atom.elem] || defaultRadius) * (scale ? scale : 1);
		mesh.position.copy(atom.coord);
		return mesh;
	};
	var createCylinder = function (p0, p1, radius, color) {
		var mesh = new THREE.Mesh(cylinderGeometry, new THREE.MeshLambertMaterial({ color: color }));
		mesh.position.copy(p0).add(p1).multiplyScalar(0.5);
		mesh.matrixAutoUpdate = false;
		mesh.lookAt(p0);
		mesh.updateMatrix();
		mesh.matrix.multiply(new THREE.Matrix4().makeScale(radius, radius, p0.distanceTo(p1))).multiply(new THREE.Matrix4().makeRotationX(Math.PI * 0.5));
		return mesh;
	};
	var createLabel = function (text, size, color) {
		var canvas = document.createElement('canvas');
		canvas.style.backgroundColor = 'rgba(0, 0, 0, 0.0)';
		var ctx = canvas.getContext('2d');
		ctx.font = size + 'px Arial';
		canvas.width = ctx.measureText(text).width;
		canvas.height = size;
		ctx.font = size + 'px Arial';
		ctx.fillStyle = color;
		ctx.fillText(text, 0, size);
		var tex = new THREE.Texture(canvas);
		tex.flipY = false;
		tex.needsUpdate = true;
		return new THREE.Mesh(labelGeo, new THREE.ShaderMaterial({
			vertexShader: labelVertexShader,
			fragmentShader: labelFragmentShader,
			uniforms: {
				map: {
					type: 't',
					value: tex,
				},
				width: {
					type: 'f',
					value: tex.image.width,
				},
				height: {
					type: 'f',
					value: tex.image.height,
				},
			},
		}));
	};
	var createRepresentationSub = function (atoms, f0, f01) {
		var ged = new THREE.Geometry();
		for (var i in atoms) {
			var atom0 = atoms[i];
			f0 && f0(atom0);
			for (var j in atom0.bonds) {
				var atom1 = atom0.bonds[j];
				if (atom1.serial < atom0.serial) continue;
				if (atom1.chain === atom0.chain && ((atom1.resi === atom0.resi) || (atom0.name === 'C' && atom1.name === 'N') || (atom0.name === 'O3\'' && atom1.name === 'P'))) {
					f01 && f01(atom0, atom1);
				} else {
					ged.vertices.push(atom0.coord);
					ged.vertices.push(atom1.coord);
				}
			}
		}
		ged.computeLineDistances();
		return new THREE.Line(ged, new THREE.LineDashedMaterial({ linewidth: linewidth, color: defaultBondColor, dashSize: 0.25, gapSize: 0.125 }), THREE.LinePieces);
	};
	var createSphereRepresentation = function (atoms) {
		var obj = new THREE.Object3D();
		obj.add(createRepresentationSub(atoms, function (atom0) {
			obj.add(createSphere(atom0, sphereRadius));
		}));
		return obj;
	};
	var createStickRepresentation = function (atoms, atomR, bondR) {
		var obj = new THREE.Object3D();
		obj.add(createRepresentationSub(atoms, function (atom0) {
			obj.add(createSphere(atom0, atomR, true));
		}, function (atom0, atom1) {
			if (atom0.color === atom1.color) {
				obj.add(createCylinder(atom0.coord, atom1.coord, bondR, atom0.color));
			} else {
				var mp = atom0.coord.clone().add(atom1.coord).multiplyScalar(0.5);
				obj.add(createCylinder(atom0.coord, mp, bondR, atom0.color));
				obj.add(createCylinder(atom1.coord, mp, bondR, atom1.color));
			}
		}));
		return obj;
	};
	var createLineRepresentation = function (atoms) {
		var obj = new THREE.Object3D();
		var geo = new THREE.Geometry();
		obj.add(new THREE.Line(geo, new THREE.LineBasicMaterial({ linewidth: linewidth, vertexColors: true }), THREE.LinePieces));
		obj.add(createRepresentationSub(atoms, function (atom0) {
			if (atom0.solvent) {
				obj.add(createSphere(atom0, sphereRadius, false, 0.2));
			}
		}, function (atom0, atom1) {
			if (atom0.color === atom1.color) {
				geo.vertices.push(atom0.coord);
				geo.vertices.push(atom1.coord);
				geo.colors.push(atom0.color);
				geo.colors.push(atom1.color);
			} else {
				var mp = atom0.coord.clone().add(atom1.coord).multiplyScalar(0.5);
				geo.vertices.push(atom0.coord);
				geo.vertices.push(mp);
				geo.vertices.push(atom1.coord);
				geo.vertices.push(mp);
				geo.colors.push(atom0.color);
				geo.colors.push(atom0.color);
				geo.colors.push(atom1.color);
				geo.colors.push(atom1.color);
			}
		}));
		return obj;
	};
	var createLabelRepresentation = function (atoms) {
		var obj = new THREE.Object3D();
		for (var i in atoms) {
			var atom = atoms[i];
			var bb = createLabel(atom.name === 'CA' ? atom.chain + ':' + atom.resn + ':' + atom.resi : atom.name, 24, '#dddddd');
			bb.position.copy(atom.coord);
			obj.add(bb);
		}
		return obj;
	};
	var refreshLigand = function(ligand) {
		if (ligand.representations === undefined) {
			var atoms = ligand.atoms;
			var labels = ligand.labels = {};
			ligand.representations = {
				label: createLabelRepresentation(labels),
			};
		}
		mdl.add(ligand.representations.label);
		mdl.position.set(0, 0, 0);
		rot.quaternion.set(1, 0, 0, 0);
		var data = $('#data');
		$('span', data).each(function() {
			var $this = $(this);
			$this.text(ligand[$this.attr('id')]);
		});
		$('#id', data).parent().attr('href', '//zinc.docking.org/substance/' + ligand.id);
/*		$('#suppliers', data).html(ligand.suppliers.map(function(supplier) {
			var link = catalogs[supplier];
			return '<li><a' + (link === undefined || link.length === 0 ? '' : ' href="' + link + '"') + '>' + supplier + '</a></li>';
		}).join(''));*/
	};
	var render = function () {
		var center = rot.position.z - camera.position.z;
		if (center < 1) center = 1;
		camera.near = center + sn;
		if (camera.near < 1) camera.near = 1;
		camera.far = center + sf;
		if (camera.near + 1 > camera.far) camera.far = camera.near + 1;
		camera.updateProjectionMatrix();
		scene.fog.near = camera.near + 0.4 * (camera.far - camera.near);
		scene.fog.far = camera.far;
		renderer.render(scene, camera);
	};
	var hasCovalentBond = function (atom1, atom2) {
		var r = covalentRadii[atom1.elem] + covalentRadii[atom2.elem];
		return atom1.coord.distanceToSquared(atom2.coord) < 1.3 * r * r;
	};
	var path = '/usr/jobs/' + location.search.substr(1) + '/';
	$('#downloads a').each(function () {
		var t = $(this);
		t.attr('href', path + t.text());
	});
	var ligand;
	$.ajax({
		url: path + 'ligands.sdf.gz',
		mimeType: 'application/octet-stream; charset=x-user-defined',
	}).done(function (lsrcz) {
		var gunzipWorker = new Worker('/gunzip.js');
		gunzipWorker.addEventListener('message', function (e) {
			var ligands = [];
			var lines = e.data.split('\n');
			for (var offset = 0, l = lines.length - 1; offset < l;)
			{
				var ligand = {
					atoms: {},
					refresh: function() {
						var r = ligand.representations[ligand.active];
						if (r === undefined) {
							switch (ligand.active) {
								case 'line':
									r = createLineRepresentation(ligand.atoms);
									break;
								case 'stick':
									r = createStickRepresentation(ligand.atoms, cylinderRadius, cylinderRadius);
									break;
								case 'ball & stick':
									r = createStickRepresentation(ligand.atoms, cylinderRadius, cylinderRadius * 0.5);
									break;
								case 'sphere':
									r = createSphereRepresentation(ligand.atoms);
									break;
							}
							ligand.representations[ligand.active] = r;
						}
						mdl.add(r);
					},
					id: lines[offset],
				}, atoms = ligand.atoms;
				offset += 3;
				var atomCount = parseInt(lines[offset].substr(0, 3));
				var bondCount = parseInt(lines[offset].substr(3, 3));
				for (var i = 1; i <= atomCount; ++i) {
					var line = lines[++offset];
					var atom = {
						serial: i,
						coord: new THREE.Vector3(parseFloat(line.substr( 0, 10)), parseFloat(line.substr(10, 10)), parseFloat(line.substr(20, 10))),
						elem: line.substr(31, 2).replace(/ /g, '').toUpperCase(),
						bonds: [],
					};
					if (atom.elem === 'H') continue;
					atom.color = atomColors[atom.elem] || defaultAtomColor;
					atoms[atom.serial] = atom;
				}
				for (var i = 1; i <= bondCount; ++i) {
					var line = lines[++offset];
					var atom0 = atoms[parseInt(line.substr(0, 3))];
					if (atom0 === undefined) continue;
					var atom1 = atoms[parseInt(line.substr(3, 3))];
					if (atom1 === undefined) continue;
					atom0.bonds.push(atom1);
					atom1.bonds.push(atom0);
				}
				while (lines[offset++] !== "$$$$");
				ligands.push(ligand);
			}
			$('#nligands').text(ligands.length);
			var ids = $('#ids');
			ids.html(ligands.map(function(ligand) {
				return '<label class="btn btn-primary"><input type="radio">' + ligand.id + '</label>';
			}).join(''));
			$(':first', ids).addClass('active');
			$('> .btn', ids).click(function(e) {
				mdl.remove(ligand.representations.label);
				mdl.remove(ligand.representations[ligand.active]);
				ligands.forEach(function(l) {
					if (l.id.toString() === $(e.target).text().trim()) {
						ligand = l;
					}
				});
				refreshLigand(ligand);
				ligand.active = $('#ligand .active').text().trim();
				ligand.refresh();
				render();
			});
			refreshLigand(ligand = ligands[0]);
			ligand.active = $('#ligand .active').text().trim();
			ligand.refresh();
			$('#ligand').click(function (e) {
				mdl.remove(ligand.representations[ligand.active]);
				ligand.active = $(e.target).text().trim();
				ligand.refresh();
				render();
			});
			var lmin = new THREE.Vector3( 9999, 9999, 9999);
			var lmax = new THREE.Vector3(-9999,-9999,-9999);
			atoms = ligand.atoms;
			for (var i in atoms) {
				var atom = atoms[i];
				var coord = atom.coord;
				lmin.min(coord);
				lmax.max(coord);
			}
			var maxD = lmax.distanceTo(lmin) + 4;
			sn = -maxD;
			sf =  maxD;
			rot.position.z = maxD * 0.35 / Math.tan(Math.PI / 180.0 * 10) - 140;
			render();
		});
		gunzipWorker.postMessage(lsrcz);
	});
	var dg, wh, cx, cy, cq, cz, cp, cn, cf;
	canvas.bind('contextmenu', function (e) {
		e.preventDefault();
	});
	canvas.bind('mouseup touchend', function (e) {
		dg = false;
	});
	canvas.bind('mousedown touchstart', function (e) {
		e.preventDefault();
		var x = e.pageX;
		var y = e.pageY;
		if (e.originalEvent.targetTouches && e.originalEvent.targetTouches[0]) {
			x = e.originalEvent.targetTouches[0].pageX;
			y = e.originalEvent.targetTouches[0].pageY;
		}
		dg = true;
		wh = e.which;
		cx = x;
		cy = y;
		cq = rot.quaternion.clone();
		cz = rot.position.z;
		cp = mdl.position.clone();
		cn = sn;
		cf = sf;
	});
	canvas.bind('mousemove touchmove', function (e) {
		e.preventDefault();
		if (!dg) return;
		var x = e.pageX;
		var y = e.pageY;
		if (e.originalEvent.targetTouches && e.originalEvent.targetTouches[0]) {
			x = e.originalEvent.targetTouches[0].pageX;
			y = e.originalEvent.targetTouches[0].pageY;
		}
		var dx = (x - cx) * canvas.widthInv;
		var dy = (y - cy) * canvas.heightInv;
		if (!dx && !dy) return;
		if (e.ctrlKey && e.shiftKey) { // Slab
			sn = cn + dx * 100;
			sf = cf + dy * 100;
		} else if (e.ctrlKey || wh == 3) { // Translate
			var scaleFactor = Math.max((rot.position.z - camera.position.z) * 0.85, 20);
			mdl.position.copy(cp).add(new THREE.Vector3(-dx * scaleFactor, -dy * scaleFactor, 0).applyQuaternion(rot.quaternion.clone().inverse().normalize()));
		} else if (e.shiftKey || wh == 2) { // Zoom
			var scaleFactor = Math.max((rot.position.z - camera.position.z) * 0.85, 80);
			rot.position.z = cz - dy * scaleFactor;
		} else { // Rotate
			var r = Math.sqrt(dx * dx + dy * dy);
			var rs = Math.sin(r * Math.PI) / r;
			rot.quaternion.set(1, 0, 0, 0).multiply(new THREE.Quaternion(Math.cos(r * Math.PI), 0, rs * dx, rs * dy)).multiply(cq);
		}
		render();
	});
	canvas.bind('mousewheel', function (e) {
		e.preventDefault();
		rot.position.z -= e.originalEvent.wheelDelta * 0.025;
		render();
	});
	canvas.bind('DOMMouseScroll', function (e) {
		e.preventDefault();
		rot.position.z += e.originalEvent.detail;
		render();
	});
	$('#exportCanvas').click(function (e) {
		render();
		window.open(renderer.domElement.toDataURL('image/png'));
	});
});
