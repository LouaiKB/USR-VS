#!/usr/bin/env node
var cluster = require('cluster');
if (cluster.isMaster) {
	var numWorkerProcesses = 4;
	console.log('Forking %d worker processes', numWorkerProcesses);
	for (var i = 0; i < numWorkerProcesses; i++) {
		cluster.fork();
	}
	cluster.on('death', (worker) => {
		console.error('Worker process %d died. Restarting...', worker.pid);
		cluster.fork();
	});
} else {
	const mongodb = require('mongodb');
	mongodb.MongoClient.connect('mongodb://localhost:27017', { useNewUrlParser: true, useUnifiedTopology: true }).then((mongoClient) => { // poolSize is 5 by default
		const jstar = mongoClient.db('jstar');
		const usr = jstar.collection('usr2');
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
		app.route('/job').get((req, res) => {
			var v = new validator(req.query);
			if (v
				.field('id').message('must be a valid job ID').objectid().copy()
				.failed()) {
				res.json(v.err);
				return;
			};
			usr.findOne({_id: new mongodb.ObjectID(v.res.id)}, {
				projection: {
					'_id': 0,
					'filename': 1,
					'usr': 1,
					'submitted': 1,
					'started': 1,
					'completed': 1,
					'nqueries': 1,
					'numConformers': 1,
				},
			}).then((doc) => {
				res.json(doc);
			});
		}).post((req, res) => {
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
			var validatesdf_out = new Buffer(0);
			validatesdf.stdout.on('data', (data) => {
				validatesdf_out = Buffer.concat([validatesdf_out, data]);
			});
			validatesdf.on('close', (code, signal) => {
				if (code) {
					res.json(code);
					return;
				} else if (signal) {
					res.json(signal);
					return;
				}
				v.res.submitted = new Date();
				v.res._id = new mongodb.ObjectID();
				var dir = __dirname + '/public/jobs/' + v.res._id;
				fs.mkdir(dir, (err) => {
					if (err) throw err;
					fs.writeFile(dir + '/query.sdf', validatesdf_out.toString(), (err) => {
						if (err) throw err;
						usr.insertOne(v.res, { w: 0 });
						res.json(v.res._id);
					});
				});
			});
			validatesdf.stdin.write(req.body['query']);
			validatesdf.stdin.end();
		});
		app.route('/embed').post((req, res) => {
			var v = new validator(req.body);
			if (v
				.field('smiles').message('must be provided, at most 1KB').length(1, 1000)
				.failed()) {
				res.json(v.err);
				return;
			}
			var embed = cp.spawn(__dirname + '/bin/embed');
			var embed_out = new Buffer(0);
			embed.stdout.on('data', (data) => {
				embed_out = Buffer.concat([embed_out, data]);
			});
			embed.on('close', (code, signal) => {
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
						embed_out: embed_out.toString(),
					});
				}
			});
			embed.stdin.write(req.body['smiles']);
			embed.stdin.end();
		});
		var http_port = 4001;
		app.listen(http_port);
		console.log('Worker %d listening on HTTP port %d in %s mode', process.pid, http_port, app.settings.env);
	});
}
