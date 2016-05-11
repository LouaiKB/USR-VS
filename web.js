#!/usr/bin/env node
var cluster = require('cluster');
if (cluster.isMaster) {
	var numWorkerProcesses = 4;
	console.log('Forking %d worker processes', numWorkerProcesses);
	for (var i = 0; i < numWorkerProcesses; i++) {
		cluster.fork();
	}
	cluster.on('death', function(worker) {
		console.error('Worker process %d died. Restarting...', worker.pid);
		cluster.fork();
	});
} else {
	var mongodb = require('mongodb');
	new mongodb.MongoClient.connect('mongodb://' + process.argv[2] + '/istar', function(err, db) {
		if (err) throw err;
		db.authenticate(process.argv[3], process.argv[4], function(err, authenticated) {
			if (err) throw err;
			var usr = db.collection('usr2');
			var express = require('express');
			var bodyParser = require('body-parser');
			var errorHandler = require('errorhandler');
			var app = express();
			app.use(bodyParser.urlencoded({ limit: '60kb', extended: false }));
			app.use(errorHandler({ dumpExceptions: true, showStack: true }));
			app.use(express.static(__dirname + '/public'));
			var validator = require('./public/validator');
			var fs = require('fs');
			var cp = require('child_process');
			app.route('/job').get(function(req, res) {
				var v = new validator(req.query);
				if (v
					.field('id').message('must be a valid job ID').objectid().copy()
					.failed()) {
					res.json(v.err);
					return;
				};
				usr.find({_id: new mongodb.ObjectID(v.res.id)}, {
					fields: {
						'_id': 0,
						'filename': 1,
						'usr': 1,
						'submitted': 1,
						'started': 1,
						'completed': 1,
						'nqueries': 1,
						'error': 1,
					},
				}).limit(1).next(function(err, doc) {
					if (err) throw err;
					res.json(doc);
				});
			}).post(function(req, res) {
				var v = new validator(req.body);
				if (v
					.field('filename').message('must be provided, at most 20 characters').length(1, 20).xss().copy()
					.field('query').message('must be provided, at most 50KB').length(1, 50000)
					.field('usr').message('must be 0 or 1').int(0).min(0).max(1).copy()
					.failed()) {
					res.json(v.err);
					return;
				}
				var validatesdf = cp.spawn(__dirname + '/bin/validatesdf');
				validatesdf.on('close', function (code, signal) {
					if (code) {
						res.json(code);
						return;
					} else if (signal) {
						res.json(signal);
						return;
					}
					v.res.version = 1;
					v.res.submitted = new Date();
					v.res._id = new mongodb.ObjectID();
					var dir = __dirname + '/public/jobs/' + v.res._id;
					fs.mkdir(dir, function (err) {
						if (err) throw err;
						fs.writeFile(dir + '/query.sdf', req.body['query'], function(err) {
							if (err) throw err;
							usr.insert(v.res, { w: 0 });
							res.json(v.res._id);
						});
					});
				});
				validatesdf.stdin.write(req.body['query']);
				validatesdf.stdin.end();
			});
			app.route('/embed').post(function(req, res) {
				var v = new validator(req.body);
				if (v
					.field('smiles').message('must be provided, at most 1KB').length(1, 1000)
					.failed()) {
					res.json(v.err);
					return;
				}
				var mktemp = cp.spawn('mktemp', ['-p', __dirname + '/public/help/tmp', '--suffix', '.svg', 'XXXXXXXXXX']);
				var mktemp_out = new Buffer(0);
				mktemp.stdout.on('data', function (data) {
					mktemp_out = Buffer.concat([mktemp_out, data]);
				});
				mktemp.on('close', function (code, signal) {
					if (code) {
						res.json({
							code: code,
						});
					} else if (signal) {
						res.json({
							signal: signal,
						});
					} else {
						var tmp = mktemp_out.slice(0, -1).toString();
						var embed = cp.spawn(__dirname + '/bin/embed', [tmp]);
						var embed_out = new Buffer(0);
						embed.stdout.on('data', function (data) {
							embed_out = Buffer.concat([embed_out, data]);
						});
						embed.on('close', function (code, signal) {
							if (code) {
								res.json({
									code: code,
								});
							} else if (signal) {
								res.json({
									signal: signal,
								});
							} else {
								res.json({
									tmp: tmp.substr(tmp.length - 14),
									embed_out: embed_out.toString(),
								});
							}
						});
						embed.stdin.write(req.body['smiles']);
						embed.stdin.end();
					}
				});
			});
			var http_port = 4000;
			app.listen(http_port);
			console.log('Worker %d listening on HTTP port %d in %s mode', process.pid, http_port, app.settings.env);
		});
	});
}
