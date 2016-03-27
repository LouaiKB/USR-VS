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
					.field('query').message('must be provided, at most 100KB').length(1, 102400)
					.field('usr').message('must be 0 or 1').int(0).min(0).max(1).copy()
					.failed()) {
					res.json(v.err);
					return;
				}
				v.res.version = 1;
				v.res.submitted = new Date();
				v.res._id = new mongodb.ObjectID();
				var dir = __dirname + '/public/jobs/' + v.res._id;
				fs.mkdir(dir, function (err) {
					if (err) throw err;
					fs.writeFile(dir + '/query.sdf', req.body['query'].split(/\r\n|\n|\r/).map(function (line) {
						if (line[5] == '.' && line[15] == '.' && line[25] == '.' && line[31] == ' ') {
							return line.substr(0, 31) + line.substr(33);
						} else {
							return line;
						}
					}).join('\n'), function(err) {
						if (err) throw err;
						usr.insert(v.res, { w: 0 });
						res.json(v.res._id);
					});
				});
			});
			var http_port = 4000;
			app.listen(http_port);
			console.log('Worker %d listening on HTTP port %d in %s mode', process.pid, http_port, app.settings.env);
		});
	});
}
