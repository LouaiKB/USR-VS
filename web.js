#!/usr/bin/env node
var fs = require('fs'),
	cluster = require('cluster');
if (cluster.isMaster) {
	// Fork worker processes with cluster
	var numCPUs = require('os').cpus().length;
	console.log('Forking %d worker processes', numCPUs);
	for (var i = 0; i < numCPUs; i++) {
		cluster.fork();
	}
	cluster.on('death', function(worker) {
		console.error('Worker process %d died. Restarting...', worker.pid);
		cluster.fork();
	});
} else {
	// Connect to MongoDB
	var mongodb = require('mongodb');
	new mongodb.MongoClient.connect('mongodb://' + process.argv[2] + '/istar', function(err, db) {
		if (err) throw err;
		db.authenticate(process.argv[3], process.argv[4], function(err, authenticated) {
			if (err) throw err;
			var usr = db.collection('usr2');
			// Configure express server
			var express = require('express');
			var bodyParser = require('body-parser');
			var favicon = require('serve-favicon');
			var errorHandler = require('errorhandler');
			var app = express();
			app.use(bodyParser.urlencoded({ limit: '100kb', extended: false }));
			app.use(errorHandler({ dumpExceptions: true, showStack: true }));
			var env = process.env.NODE_ENV || 'development';
			if (env == 'development') {
				app.use(express.static(__dirname + '/public'));
				app.use(favicon(__dirname + '/public/favicon.ico'));
			} else if (env == 'production') {
				var oneDay = 1000 * 60 * 60 * 24;
				var oneYear = oneDay * 365.25;
				app.use(express.static(__dirname + '/public', { maxAge: oneDay }));
				app.use(favicon(__dirname + '/public/favicon.ico', { maxAge: oneYear }));
			};
			// Define helper variables and functions
			var validator = require('./public/validator');
			app.route('/jobs').get(function(req, res) {
				var v = new validator(req.query);
				if (v
					.field('skip').message('must be a non-negative integer').int(0).min(0).copy()
					.field('count').message('must be a non-negative integer').int(0).min(0).copy()
					.failed() || v
					.range('skip', 'count')
					.failed()) {
					res.json(v.err);
					return;
				};
				usr.count(function(err, count) {
					if (err) throw err;
					if (v
						.field('count').message('must be no greater than ' + count).max(count)
						.failed()) {
						res.json(v.err);
						return;
					}
					usr.find({}, {
						fields: v.res.count == count ? {
							'_id': 0,
							'started': 1,
							'done': 1,
						} : {
							'description': 1,
							'submitted': 1,
							'started': 1,
							'done': 1,
						},
						sort: {'submitted': 1},
						skip: v.res.skip,
						limit: count - v.res.skip
					}).toArray(function(err, docs) {
						if (err) throw err;
						res.json(docs);
					});
				});
			}).post(function(req, res) {
				var v = new validator(req.body);
				if (v
					.field('email').message('must be valid').email().copy()
					.field('description').message('must be provided, at most 20 characters').length(1, 20).xss().copy()
					.field('ligand').message('must be provided and must not exceed 100KB').length(1, 102400)
					.failed()) {
					res.json(v.err);
					return;
				}
				v.res.submitted = new Date();
				v.res._id = new mongodb.ObjectID();
				var dir = process.cwd() + '/public/jobs/' + v.res._id;
				fs.mkdir(dir, function (err) {
					if (err) throw err;
					fs.writeFile(dir + '/ligand.sdf', req.body['ligand'], function(err) {
						if (err) throw err;
						usr.insert(v.res, { w: 0 });
						res.json({});
					});
				});
			});
			// Start listening
			var http_port = 4000;
			app.listen(http_port);
			console.log('Worker %d listening on HTTP port %d in %s mode', process.pid, http_port, app.settings.env);
		});
	});
}
