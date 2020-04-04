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

$(() => {

	// Get the latest job status
	var status = $('#status');
	var jobid = location.search.substr(1);
	var tick = () => {
		$.get('../job', { id: jobid }, (job) => {
			job.usrF = ['USR', 'USRCAT'][job.usr];
			['submitted', 'started', 'completed'].forEach((key) => {
				if (job[key] === undefined) return;
				job[key] = new Date(job[key]);
				job[key+'F'] = $.format.date(job[key], 'yyyy/MM/dd HH:mm:ss.SSS');
			});
			job.info = job.completed ? (job.error ? ["", "Failed to parse the query file. Please choose a file in SDF format."][parseInt(job.error)] : 'Completed ' + job.nqueries + ' ' + (job.nqueries == 1 ? 'query' : 'queries') + ' in ' + (runtime=((job['completed']-job['started'])*0.001)).toFixed(3) + ' seconds.<br>Screening speed was ' + (job.numConformers*0.001*parseInt(job.nqueries)/runtime).toFixed(0) + 'K 3D conformers per second.') : (job.started ? 'Execution in progress <img src="loading.gif" style="width: 16px; height: 16px;">' : 'Queued for execution');
			$('span', status).each(function (d) { // 'this' binding is used.
				var t = $(this);
				var c = job[t.attr('id')];
				if (t.html() !== c) {
					t.html(c).hide().fadeIn('slow');
				}
			});
			const path = '../jobs/' + jobid + '/';
			$('#filename', status).parent().attr('href', path + 'query.sdf');
			if (!job.completed) {
				setTimeout(tick, 1000);
				return;
			}
			if (job.error) return;
			$('#results').removeClass('hidden');
			const gaugeScreeningSpeed = echarts.init(document.getElementById('gaugeScreeningSpeed'));
			gaugeScreeningSpeed.setOption({
				series: [{
					name: 'Screening speed',
					type: 'gauge',
					min: 0,
					max: 60,
					splitNumber: 12,
//					radius: '100%',
					axisLine: {
						lineStyle: {
							color: [[1/6, '#91c7ae'], [5/6, '#63869e'], [1, '#c23531']],
							width: 8,
						},
					},
					axisTick: {
						length: 15,
					},
					splitLine: {
						length: 20,
						lineStyle: {
							color: 'auto',
						},
					},
					title: {
						fontWeight: 'bold',
					},
					detail: {
						formatter: (value) => {
							return value.toFixed(2);
						},
						fontWeight: 'bold',
					},
					data: [{
						value: job.numConformers*(10**-6)*parseInt(job.nqueries)/runtime,
						name: 'Million\nconformers / s',
					}],
				}],
			});

			const atomColors = { // http://jmol.sourceforge.net/jscolors
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
			const defaultAtomColor = new THREE.Color(0xCCCCCC);
			const defaultBackgroundColor = new THREE.Color(0x000000);
			const sphereGeometry = new THREE.SphereBufferGeometry(1, 64, 64);
			const cylinderGeometry = new THREE.CylinderBufferGeometry(1, 1, 1, 64, 1);
			const labelVertexShader = '\
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
			const labelFragmentShader = '\
uniform sampler2D map;\n\
varying vec2 vUv;\n\
void main()\n\
{\n\
	gl_FragColor = texture2D(map, vec2(vUv.x, 1.0 - vUv.y));\n\
	if (gl_FragColor.a < 0.5) discard;\n\
}';
			const labelGeo = new THREE.Geometry();
			for (var i = 0; i < 6; ++i) {
				labelGeo.vertices.push(new THREE.Vector3(0, 0, 0));
			}
			labelGeo.faces.push(new THREE.Face3(0, 1, 2));
			labelGeo.faces.push(new THREE.Face3(0, 2, 3));
			labelGeo.faceVertexUvs[0].push([new THREE.Vector2(0, 0), new THREE.Vector2(1, 1), new THREE.Vector2(0, 1)]);
			labelGeo.faceVertexUvs[0].push([new THREE.Vector2(0, 0), new THREE.Vector2(1, 0), new THREE.Vector2(1, 1)]);
			var fontLoader = new THREE.FontLoader();
			var helvetiker_regular;
			fontLoader.load('helvetiker_regular.typeface.json', (font) => {
				helvetiker_regular = font;
			});
			var createSphere = (atom, radius) => {
				var mesh = new THREE.Mesh(sphereGeometry, new THREE.MeshLambertMaterial({ color: atom.color }));
				mesh.scale.x = mesh.scale.y = mesh.scale.z = radius;
				mesh.position.copy(atom.coord);
				return mesh;
			};
			var createCylinder = (p0, p1, radius, color) => {
				var mesh = new THREE.Mesh(cylinderGeometry, new THREE.MeshLambertMaterial({ color: color }));
				mesh.position.copy(p0).add(p1).multiplyScalar(0.5);
				mesh.lookAt(p0);
				mesh.updateMatrix();
				mesh.matrixAutoUpdate = false;
				mesh.matrix.multiply(new THREE.Matrix4().makeScale(radius, radius, p0.distanceTo(p1))).multiply(new THREE.Matrix4().makeRotationX(Math.PI * 0.5));
				return mesh;
			};
			var createLabel = (text, size, color) => {
				var mesh = new THREE.Mesh(new THREE.TextBufferGeometry(text, {
					font: helvetiker_regular,
					size: 0.6,
					height: 0.1,
				}), new THREE.MeshBasicMaterial({ color: color }));
//				return mesh;
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
			var createStickRepresentation = (atoms, atomR, bondR) => {
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
			var createLabelRepresentation = (atoms) => {
				var obj = new THREE.Object3D();
				for (var i in atoms) {
					var atom = atoms[i];
					if (atom.elem === 'C') continue;
					var bb = createLabel(atom.elem, 64, '#EEEEEE');
					bb.position.copy(atom.coord.clone().add(new THREE.Vector3(0.2, 0.2, 0.2)));
					obj.add(bb);
				}
				return obj;
			};
			var iview = function (canvas) { // 'this' binding is used.
				this.canvas = $(canvas);
				this.canvas.height(this.canvas.width());
				this.canvas.widthInv  = 1 / this.canvas.width();
				this.canvas.heightInv = 1 / this.canvas.height();
				this.renderer = new THREE.WebGLRenderer({
					canvas: this.canvas.get(0),
//					context: canvas.getContext('webgl2'),
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
				this.canvas.bind('contextmenu', (e) => {
					e.preventDefault();
				});
				this.canvas.bind('mouseup touchend', (e) => {
					me.dg = false;
				});
				this.canvas.bind('mousedown touchstart', (e) => {
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
				this.canvas.bind('mousemove touchmove', (e) => {
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
				this.canvas.bind('mousewheel', (e) => {
					e.preventDefault();
					me.rot.position.z -= e.originalEvent.wheelDelta * 0.025;
					me.render();
				});
				this.canvas.bind('DOMMouseScroll', (e) => {
					e.preventDefault();
					me.rot.position.z += e.originalEvent.detail;
					me.render();
				});
			};
			iview.prototype = {
				constructor: iview,
				reset: function (molecule) { // 'this' binding is used.
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
				render: function () { // 'this' binding is used.
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
				exportCanvas: function () { // 'this' binding is used.
					this.render();
					window.open(this.renderer.domElement.toDataURL('image/png'));
				},
			};
			var parseSDF = (src) => {
				var molecules = [];
				for (var lines = src.split(/\r\n|\n|\r/), l = lines.length - 1, offset = 0; offset < l;) {
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
					for (var line = lines[offset++]; line !== '$$$$'; line = lines[offset++]) {
						if (line[0] === '>') {
							const prop = line.split('<')[1].split('>')[0];
							molecule[prop] = lines[offset++];
							if (prop.startsWith('num')) {
								molecule[prop] = parseInt(molecule[prop]);
							} else if (['exactMW', 'clogP', 'tPSA'].includes(prop)) {
								molecule[prop] = parseFloat(molecule[prop]);
							}
						}
					}
					molecules.push(molecule);
				}
				return molecules;
			};
			var refreshMolecule = (molecule, iv) => {
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

			var iviews = $('canvas.three').map(function (index) { // 'this' binding is used.
				var iv = new iview(this);
				$('#exportButton' + index).click(function (e) {
					iv.exportCanvas();
				});
				return iv;
			});
			const smilesDrawer = new SmilesDrawer.Drawer({
				width: 540,
				height: 540,
			});
			$.ajax({
				url: path + 'query.sdf',
			}).done((qsdf) => {
				var qmolecules = parseSDF(qsdf).slice(0, 1); // Take out only the first query molecule and discard subsequent query molecules, if any.
				if (qmolecules.length !== job.nqueries) throw Error("qmolecules.length !== job.nqueries");
				$('#qids_label').text(qmolecules.length + ' query molecule' + (qmolecules.length == 1 ? '' : 's'));
				var qindex;
				var refreshQuery = (qidx) => {
					const qmolecule = qmolecules[qindex = qidx];
					refreshMolecule(qmolecule, iviews[0]);
					var qpath = path + qindex + '/';
					var output = $('#output');
					$('#downloads a', output).each(function () { // 'this' binding is used.
						var t = $(this);
						t.attr('href', qpath + t.text());
					});
					var qproperties = $('#qproperties');
					$('span', qproperties).each(function () { // 'this' binding is used.
						var t = $(this);
						var prop = t.attr('id');
						var idx = ['tPSA', 'exactMW', 'clogP'].indexOf(prop);
						if (idx === -1) {
							t.text(qmolecule[prop]);
						} else {
							t.text(qmolecule[prop].toFixed(2 + idx)); // Display tPSA with 2 digits. Display exactMW with 3 digits. Display clogP with 4 digits.
						}
					});
					if (!qmolecule['canonicalSmilesTree']) {
						SmilesDrawer.parse(qmolecule['canonicalSMILES'], (qtree) => { // SmilesDrawer.parse() is a static function.
							qmolecule['canonicalSmilesTree'] = qtree;
		                }, (err) => {
							// TODO: noty()
						});
					}
					smilesDrawer.draw(qmolecule['canonicalSmilesTree'], 'qdrawer', 'dark');
					$.ajax({
						url: qpath + 'hits.sdf',
					}).done((hsdf) => {
						var hmolecules = parseSDF(hsdf);
						if (hmolecules.length !== 100) throw Error("hmolecules.length !== 100");
						$.ajax({
							url: qpath + 'hits.csv',
						}).done((hcsv) => {
							var logs = hcsv.split(/\r?\n/).slice(1, 101);
							if (logs.length !== hmolecules.length) throw Error("logs.length !== hmolecules.length");
							var propNames = [ 'usr_score', 'usrcat_score', 'tanimoto_score', 'canonicalSMILES', 'molFormula', 'numAtoms', 'numHBD', 'numHBA', 'numRotatableBonds', 'numRings', 'exactMW', 'tPSA', 'clogP' ];
							$.each(hmolecules, (i, molecule) => {
								var properties = logs[i].split(',');
								if (molecule.id !== properties[0]) throw Error("molecule.id !== properties[0]");
								$.each(propNames, (j, propName) => {
									molecule[propName] = properties[1+j];
								});
							});
							$('#hids_label').text(hmolecules.length + ' hit molecules sorted by ' + job.usrF + ' score');
							var hindex;
							var refreshHit = (hidx) => {
								var molecule = hmolecules[hindex = hidx];
								refreshMolecule(molecule, iviews[1]);
								$('span', output).each(function () { // 'this' binding is used.
									var t = $(this);
									t.text(molecule[t.attr('id')]);
								});
								var hproperties = $('#hproperties');
								$('span', hproperties).each(function () { // 'this' binding is used.
									var t = $(this);
									t.text(molecule[t.attr('id')]);
								});
								$('#id', hproperties).parent().attr('href', '//zinc.docking.org/substance/' + molecule.id);
								if (!molecule['canonicalSmilesTree']) {
									SmilesDrawer.parse(molecule['canonicalSMILES'], (htree) => { // SmilesDrawer.parse() is a static function.
										molecule['canonicalSmilesTree'] = htree;
					                }, (err) => {
										// TODO: noty()
									});
								}
								smilesDrawer.draw(molecule['canonicalSmilesTree'], 'hdrawer', 'dark');
							};
							var hids = $('#hids');
							hids.html(hmolecules.map((molecule, index) => {
								return `<button type="button" class="btn btn-primary">${index}</button>`;
							}).join(''));
							$('> button', hids).click((e) => {
								var hidx = $(e.target).text();
								if (hidx == hindex) return;
								$('> button.active', hids).removeClass('active');
								refreshHit(hidx);
							});
							$(':first', hids).addClass('active');
							refreshHit(0);
						});
					});
				};
				var qids = $('#qids');
				qids.html(Array.apply(0, Array(qmolecules.length)).map((value, index) => {
					return `<button type="button" class="btn btn-primary">${index}</button>`;
				}).join(''));
				$('> .btn', qids).click((e) => {
					var qidx = $(e.target).text();
					if (qidx == qindex) return;
					$('> button.active', qids).removeClass('active');
					refreshQuery(qidx);
				});
				$(':first', qids).addClass('active');
				refreshQuery(0);
			});
		});
	};
	tick();
});
