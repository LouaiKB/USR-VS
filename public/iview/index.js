/*!
 * iview is an interactive WebGL visualizer for protein-molecule complex. iview is based on GLmol, three.js and jQuery.
 * http://github.com/HongjianLi/istar
 * Copyright (c) 2012-2014 Chinese University of Hong Kong
 * License: Apache License 2.0
 * Hongjian Li, Kwong-Sak Leung, Takanori Nakane and Man-Hon Wong.
 * iview: an interactive WebGL visualizer for protein-molecule complex.
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

	// Get the latest job status
	var status = $('#status');
	var jobid = location.search.substr(1);
	var tick = function () {
		$.get('../job', { id: jobid }, function (job) {
			job.usrF = ['USR', 'USRCAT'][job.usr];
			['submitted', 'started', 'completed'].forEach(function (key) {
				if (job[key] === undefined) return;
				job[key] = new Date(job[key]);
				job[key+'F'] = $.format.date(job[key], 'yyyy/MM/dd HH:mm:ss.SSS');
			});
			job.info = job.completed ? (job.error ? ["", "Failed to parse the query file. Please choose a file in SDF format."][parseInt(job.error)] : 'Completed ' + job.nqueries + ' ' + (job.nqueries == 1 ? 'query' : 'queries') + ' in ' + (runtime=((job['completed']-job['started'])*0.001)).toFixed(3) + ' seconds.<br>Screening speed was ' + (93.903333*parseInt(job.nqueries)/runtime).toFixed(0) + ' million 3D conformers per second.') : (job.started ? 'Execution in progress <img src="loading.gif" style="width: 16px; height: 16px;">' : 'Queued for execution');
			$('span', status).each(function (d) {
				var t = $(this);
				var c = job[t.attr('id')];
				if (t.html() !== c) {
					t.html(c).hide().fadeIn('slow');
				}
			});
			var path = '../jobs/' + jobid + '/';
			$('#filename', status).parent().attr('href', path + 'query.sdf');
			if (!job.completed) {
				setTimeout(tick, 1000);
				return;
			}
			if (job.error) return;
			$('#results').removeClass('hidden');

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
				'AmpC molecules': 'http://shoichetlab.compbio.ucsf.edu',
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
				'Zelinsky Institute': 'http://zelinsky.com',
				'Zelinsky Institute Building Blocks': 'http://zelinsky.com',
				'Zelinsky Institute Make on Demand': 'http://zelinsky.com',
				'ZereneX': 'http://www.zerenex-molecular.com',
				'ZereneX Building Blocks': 'http://www.zerenex-molecular.com',
				'Zylexa Pharma': 'http://www.zylexa-pharma.com',
				'Zylexa Pharma BB': 'http://www.zylexa-pharma.com',
			};
			var atomColors = { // http://jmol.sourceforge.net/jscolors
				 H: new THREE.Color(0xFFFFFF),
				 C: new THREE.Color(0x909090),
				 N: new THREE.Color(0x3050F8),
				 O: new THREE.Color(0xFF0D0D),
				 F: new THREE.Color(0x90E050),
				NA: new THREE.Color(0xAB5CF2),
				MG: new THREE.Color(0x8AFF00),
				 P: new THREE.Color(0xFF8000),
				 S: new THREE.Color(0xFFFF30),
				CL: new THREE.Color(0x1FF01F),
				 K: new THREE.Color(0x8F40D4),
				CA: new THREE.Color(0x3DFF00),
				MN: new THREE.Color(0x9C7AC7),
				FE: new THREE.Color(0xE06633),
				CO: new THREE.Color(0xF090A0),
				NI: new THREE.Color(0x50D050),
				CU: new THREE.Color(0xC88033),
				ZN: new THREE.Color(0x7D80B0),
				AS: new THREE.Color(0xBD80E3),
				SE: new THREE.Color(0xFFA100),
				BR: new THREE.Color(0xA62929),
				SR: new THREE.Color(0x00FF00),
				MO: new THREE.Color(0x54B5B5),
				CD: new THREE.Color(0xFFD98F),
				 I: new THREE.Color(0x940094),
				CS: new THREE.Color(0x57178F),
				HG: new THREE.Color(0xB8B8D0),
				 U: new THREE.Color(0x008FFF),
			};
			var defaultAtomColor = new THREE.Color(0xCCCCCC);
			var defaultBackgroundColor = new THREE.Color(0x000000);
			var sphereGeometry = new THREE.SphereGeometry(1, 64, 64);
			var cylinderGeometry = new THREE.CylinderGeometry(1, 1, 1, 64, 1);
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
			var createSphere = function (atom, radius) {
				var mesh = new THREE.Mesh(sphereGeometry, new THREE.MeshLambertMaterial({ color: atom.color }));
				mesh.scale.x = mesh.scale.y = mesh.scale.z = radius;
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
				canvas.width = size;
				canvas.height = size;
				var ctx = canvas.getContext('2d');
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
			var createStickRepresentation = function (atoms, atomR, bondR) {
				var obj = new THREE.Object3D();
				for (var i in atoms) {
					var atom0 = atoms[i];
					obj.add(createSphere(atom0, atomR, false, 0.4));
					for (var j in atom0.bonds) {
						var atom1 = atom0.bonds[j];
						if (atom1.serial < atom0.serial) continue;
						if (atom0.color === atom1.color) {
							obj.add(createCylinder(atom0.coord, atom1.coord, bondR, atom0.color));
						} else {
							var mp = atom0.coord.clone().add(atom1.coord).multiplyScalar(0.5);
							obj.add(createCylinder(atom0.coord, mp, bondR, atom0.color));
							obj.add(createCylinder(atom1.coord, mp, bondR, atom1.color));
						}
					}
				}
				return obj;
			};
			var createLabelRepresentation = function (atoms) {
				var obj = new THREE.Object3D();
				for (var i in atoms) {
					var atom = atoms[i];
					var bb = createLabel(atom.elem, 64, '#EEEEEE');
					bb.position.copy(atom.coord);
					obj.add(bb);
				}
				return obj;
			};
			var iview = function (canvas) {
				this.canvas = $(canvas);
				this.canvas.height(this.canvas.width());
				this.canvas.widthInv  = 1 / this.canvas.width();
				this.canvas.heightInv = 1 / this.canvas.height();
				this.renderer = new THREE.WebGLRenderer({
					canvas: this.canvas.get(0),
					antialias: true,
				});
				this.renderer.setSize(this.canvas.width(), this.canvas.height());
				this.renderer.setClearColor(defaultBackgroundColor);
				this.directionalLight = new THREE.DirectionalLight(0xFFFFFF, 1.2);
				this.directionalLight.position.set(0.2, 0.2, -1).normalize();
				this.ambientLight = new THREE.AmbientLight(0x202020);
				this.mdl = new THREE.Object3D();
				this.rot = new THREE.Object3D();
				this.rot.add(this.mdl);
				this.scene = new THREE.Scene();
				this.scene.add(this.directionalLight);
				this.scene.add(this.ambientLight);
				this.scene.add(this.rot);
				this.scene.fog = new THREE.Fog(defaultBackgroundColor, 100, 200);
				this.camera = new THREE.PerspectiveCamera(20, this.canvas.width() / this.canvas.height(), 1, 800);
				this.camera.position.set(0, 0, -150);
				this.camera.lookAt(new THREE.Vector3(0, 0, 0));
				var me = this;
				this.canvas.bind('contextmenu', function (e) {
					e.preventDefault();
				});
				this.canvas.bind('mouseup touchend', function (e) {
					me.dg = false;
				});
				this.canvas.bind('mousedown touchstart', function (e) {
					e.preventDefault();
					var x = e.pageX;
					var y = e.pageY;
					if (e.originalEvent.targetTouches && e.originalEvent.targetTouches[0]) {
						x = e.originalEvent.targetTouches[0].pageX;
						y = e.originalEvent.targetTouches[0].pageY;
					}
					me.dg = true;
					me.wh = e.which;
					me.cx = x;
					me.cy = y;
					me.cq = me.rot.quaternion.clone();
					me.cz = me.rot.position.z;
					me.cp = me.mdl.position.clone();
					me.cn = me.sn;
					me.cf = me.sf;
				});
				this.canvas.bind('mousemove touchmove', function (e) {
					e.preventDefault();
					if (!me.dg) return;
					var x = e.pageX;
					var y = e.pageY;
					if (e.originalEvent.targetTouches && e.originalEvent.targetTouches[0]) {
						x = e.originalEvent.targetTouches[0].pageX;
						y = e.originalEvent.targetTouches[0].pageY;
					}
					var dx = (x - me.cx) * me.canvas.widthInv;
					var dy = (y - me.cy) * me.canvas.heightInv;
					if (!dx && !dy) return;
					if (e.ctrlKey && e.shiftKey) { // Slab
						me.sn = me.cn + dx * 100;
						me.sf = me.cf + dy * 100;
					} else if (e.ctrlKey || me.wh == 3) { // Translate
						var scaleFactor = Math.max((me.rot.position.z - me.camera.position.z) * 0.85, 20);
						me.mdl.position.copy(me.cp).add(new THREE.Vector3(-dx * scaleFactor, -dy * scaleFactor, 0).applyQuaternion(me.rot.quaternion.clone().inverse().normalize()));
					} else if (e.shiftKey || me.wh == 2) { // Zoom
						var scaleFactor = Math.max((me.rot.position.z - me.camera.position.z) * 0.85, 80);
						me.rot.position.z = me.cz - dy * scaleFactor;
					} else { // Rotate
						var r = Math.sqrt(dx * dx + dy * dy);
						var rs = Math.sin(r * Math.PI) / r;
						me.rot.quaternion.set(1, 0, 0, 0).multiply(new THREE.Quaternion(Math.cos(r * Math.PI), 0, rs * dx, rs * dy)).multiply(me.cq);
					}
					me.render();
				});
				this.canvas.bind('mousewheel', function (e) {
					e.preventDefault();
					me.rot.position.z -= e.originalEvent.wheelDelta * 0.025;
					me.render();
				});
				this.canvas.bind('DOMMouseScroll', function (e) {
					e.preventDefault();
					me.rot.position.z += e.originalEvent.detail;
					me.render();
				});
			};
			iview.prototype = {
				constructor: iview,
				reset: function (molecule) {
					var maxD = molecule.maxD;
					if (maxD === undefined) {
						var cmin = new THREE.Vector3( 10000, 10000, 10000);
						var cmax = new THREE.Vector3(-10000,-10000,-10000);
						var csum = new THREE.Vector3();
						var atoms = molecule.atoms;
						for (var i in atoms) {
							var coord = atoms[i].coord;
							csum.add(coord);
							cmin.min(coord);
							cmax.max(coord);
						}
						molecule.maxD = maxD = cmax.distanceTo(cmin) + 4;
						molecule.ctrV = csum.clone().multiplyScalar(-1 / molecule.nha);
					}
					this.sn = -maxD;
					this.sf =  maxD;
					this.mdl.position.copy(molecule.ctrV);
					this.rot.position.z = maxD * 0.35 / Math.tan(Math.PI / 180.0 * 10) - 140;
				},
				render: function () {
					var center = this.rot.position.z - this.camera.position.z;
					if (center < 1) center = 1;
					this.camera.near = center + this.sn;
					if (this.camera.near < 1) this.camera.near = 1;
					this.camera.far = center + this.sf;
					if (this.camera.near + 1 > this.camera.far) this.camera.far = this.camera.near + 1;
					this.camera.updateProjectionMatrix();
					this.scene.fog.near = this.camera.near + 0.4 * (this.camera.far - this.camera.near);
					this.scene.fog.far = this.camera.far;
					this.renderer.render(this.scene, this.camera);
				},
				exportCanvas: function () {
					this.render();
					window.open(this.renderer.domElement.toDataURL('image/png'));
				},
			};
			var parseSDF = function (src) {
				var molecules = [];
				for (var lines = src.split(/\r?\n/), l = lines.length - 1, offset = 0; offset < l;) {
					var molecule = {
						atoms: {},
						id: lines[offset],
					}, atoms = molecule.atoms;
					offset += 3;
					var atomCount = parseInt(lines[offset].substr(0, 3));
					var bondCount = parseInt(lines[offset].substr(3, 3));
					for (var i = 1; i <= atomCount; ++i) {
						var line = lines[++offset];
						var atom = {
							serial: i,
							coord: new THREE.Vector3(parseFloat(line.substr(0, 10)), parseFloat(line.substr(10, 10)), parseFloat(line.substr(20, 10))),
							elem: line.substr(31, 2).replace(/ /g, '').toUpperCase(),
							bonds: [],
						};
						if (atom.elem === 'H') continue;
						atom.color = atomColors[atom.elem] || defaultAtomColor;
						atoms[atom.serial] = atom;
					}
					molecule.nha = Object.keys(atoms).length;
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
					molecules.push(molecule);
				}
				return molecules;
			};
			var refreshMolecule = function (molecule, iv) {
				if (molecule.representations === undefined) {
					molecule.representations = {
						stick: createStickRepresentation(molecule.atoms, 0.3, 0.3),
						label: createLabelRepresentation(molecule.atoms),
					};
				}
				iv.mdl.children = [];
				iv.mdl.add(molecule.representations.stick);
				iv.mdl.add(molecule.representations.label);
				iv.reset(molecule);
				iv.render();
			};

			var iviews = $('canvas').map(function (index) {
				var iv = new iview(this);
				$('#exportButton' + index).click(function (e) {
					iv.exportCanvas();
				});
				return iv;
			});
			$.ajax({
				url: path + 'query.sdf',
			}).done(function (qsdf) {
				var qmolecules = parseSDF(qsdf).slice(0, 1);
				if (qmolecules.length !== job.nqueries) throw Error("qmolecules.length !== job.nqueries");
				$('#qids_label').text(qmolecules.length + ' query molecule' + (qmolecules.length == 1 ? '' : 's'));
				var qindex;
				var refreshQuery = function (qidx) {
					refreshMolecule(qmolecules[qindex = qidx], iviews[0]);
					$('#downloads a').each(function () {
						var t = $(this);
						t.attr('href', path + qindex + '/' + t.text());
					});
					$.ajax({
						url: path + qindex + '/hits.sdf',
					}).done(function (hsdf) {
						var hmolecules = parseSDF(hsdf);
						if (hmolecules.length !== 100) throw Error("hmolecules.length !== 100");
						$.ajax({
							url: path + qindex + '/hits.csv',
						}).done(function (hcsv) {
							var logs = hcsv.split(/\r?\n/).slice(1, 101);
							if (logs.length !== hmolecules.length) throw Error("logs.length !== hmolecules.length");
							var propNames = [ 'usr_score', 'usrcat_score', 'mwt', 'lgp', 'ads', 'pds', 'hbd', 'hba', 'psa', 'chg', 'nrb', 'smiles', 'suppliers' ];
							$.each(hmolecules, function (i, molecule) {
								var properties = logs[i].split(',');
								if (molecule.id !== properties[0]) throw Error("molecule.id !== properties[0]");
								$.each(propNames, function (j, propName) {
									molecule[propName] = properties[1+j];
								});
								molecule.suppliers = molecule.suppliers.split(' | ').slice(1);
								molecule.nsuppliers = molecule.suppliers.length;
							});
							$('#hids_label').text(hmolecules.length + ' hit molecules sorted by ' + job.usrF + ' score');
							var hindex;
							var refreshHit = function (hidx) {
								var molecule = hmolecules[hindex = hidx];
								refreshMolecule(molecule, iviews[1]);
								var output = $('#output');
								$('span', output).each(function () {
									var t = $(this);
									t.text(molecule[t.attr('id')]);
								});
								$('#id', output).parent().attr('href', '//zinc.docking.org/substance/' + molecule.id);
								$('#suppliers', output).html(molecule.suppliers.map(function (supplier) {
									var link = catalogs[supplier];
									return '<li><a' + (link === undefined ? '' : ' href="' + link + '"') + '>' + supplier + '</a></li>';
								}).join(''));
							};
							var hids = $('#hids');
							hids.html(hmolecules.map(function (molecule, index) {
								return '<label class="btn btn-primary"><input type="radio">' + index + '</label>';
							}).join(''));
							$('> .btn', hids).click(function (e) {
								var hidx = $(e.target).text();
								if (hidx == hindex) return;
								refreshHit(hidx);
							});
							$(':first', hids).addClass('active');
							refreshHit(0);
						});
					});
				};
				var qids = $('#qids');
				qids.html(Array.apply(0, Array(qmolecules.length)).map(function (value, index) {
					return '<label class="btn btn-primary"><input type="radio">' + index + '</label>';
				}).join(''));
				$('> .btn', qids).click(function (e) {
					var qidx = $(e.target).text();
					if (qidx == qindex) return;
					refreshQuery(qidx);
				});
				$(':first', qids).addClass('active');
				refreshQuery(0);
			});
		});
	};
	tick();
});
