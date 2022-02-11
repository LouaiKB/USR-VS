$(function() {
	$('.form-group a').tooltip();
	var query_label = $('#query_label');
	var submit = $('#submit');
	submit.click(function() {
		var file = $('input[type="file"]')[0].files[0];
		if (file === undefined || file.size > 50000) {
			query_label.tooltip('show');
			return;
		}
		var reader = new FileReader();
		reader.onload = function () {
			$.post('job', {
				query: reader.result,
				filename: file.name.substr(0, 20),
				usr: $('input[name="usrRadioOptions"]:checked').val(),
			}, function(res) {
				if (res.length !== 24) {
					query_label.tooltip('show');
					return;
				}
				location.assign('iview/?' + res);
			});
		};
		reader.readAsText(file);
	});
});
