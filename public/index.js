$(function() {

	// Initialize pager
	var pager = $('#pager');
	pager.pager('init', [ 'Description', 'Submitted', 'Status', 'Result' ], function(job) {
		var status, result = '<a href="iview/?' + job._id + '"><img src="iview/logo.png" alt="iview"></a>';
		if (!job.started) {
			status = 'Queued for execution';
		} else if (!job.done) {
			status = 'Execution in progress';
		} else {
			status = 'Done ' + $.format.date(new Date(job.done), 'yyyy/MM/dd HH:mm:ss');
			result += '<a href="jobs/' + job._id + '/log.csv.gz"><img src="/excel.png" alt="log.csv.gz"></a><a href="jobs/' + job._id + '/hits.sdf.gz"><img src="/molecule.png" alt="hits.sdf.gz"></a>';
		}
		return [
			job.description,
			$.format.date(new Date(job.submitted), 'yyyy/MM/dd HH:mm:ss'),
			status,
			result
		];
	});

	// Refresh the table of jobs and its pager every second
	var jobs = [], skip = 0;
	var tick = function() {
		$.get('jobs', { skip: skip, count: jobs.length }, function(res) {
			if (res.length) {
				for (var i = skip; i < jobs.length; ++i) {
					var job = res[i - skip];
					jobs[i].started = job.started;
					jobs[i].done = job.done;
				}
				pager.pager('refresh', skip, jobs.length, 2, 5, false);
				if (res.length > jobs.length - skip) {
					var len = jobs.length;
					jobs = jobs.concat(res.slice(jobs.length - skip));
					pager.pager('source', jobs);
					pager.pager('refresh', len, jobs.length, 0, 5, true);
				}
				for (; skip < jobs.length && jobs[skip].done; ++skip);
			}
			setTimeout(tick, 1000);
		});
	};
	tick();

	// Load query ligand locally
	var query;
	$('input[type="file"]').change(function() {
		var file = this.files[0];
		if (file === undefined) return;
		$('#description').val(file.name);
		var reader = new FileReader();
		reader.onload = function () {
			query = reader.result;
		};
		reader.readAsText(file);
	});

	// Initialize tooltips
	$('.form-group a').tooltip();

	// Process submission
	var submit = $('#submit');
	submit.click(function() {
		// Hide tooltips
		$('.form-group a').tooltip('hide');
		// Do client side validation
		var job = {
			query: query,
			description: $('#description').val(),
			email: $('#email').val(),
		};
		var v = new validator(job);
		if (v
			.field('query').message('must be provided and must not exceed 100KB').length(1, 102400)
			.field('description').message('must be provided, at most 20 characters').length(1, 20)
			.field('email').message('must be valid').email()
			.failed()) {
			var keys = Object.keys(v.err);
			keys.forEach(function(key) {
				$('#' + key + '_label').tooltip('show');
			});
			$('#' + keys[0]).focus();
			return;
		}
		job.query = query;
		// Disable the submit button for a while
		submit.prop('disabled', true);
		// Post a new job with server side validation
		$.post('jobs', job, function(res) {
			var keys = Object.keys(res);
			// If server side validation fails, show the tooltips
			if (keys.length) {
				keys.forEach(function(key) {
					$('#' + key + '_label').tooltip('show');
				});
				$('#' + keys[0]).focus();
			} else {
				$('html, body').animate({ scrollTop: pager.offset().top });
			}
		}, 'json').always(function() {
			submit.prop('disabled', false);
		});
	});

	// Apply accordion to tutorials
	$('.ui-accordion').accordion({
		collapsible: true,
		active: false,
		heightStyle: 'content',
		activate: function(event, ui) {
			$('img', this).trigger('expand');
		}
	});
	$('.ui-accordion img').lazyload({
		event: 'expand',
		effect: 'fadeIn',
	});
});
