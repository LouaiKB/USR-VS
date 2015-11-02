$(function() {
	$('.form-group a').tooltip();
	var submit = $('#submit');
	submit.click(function() {
		var file = $('input[type="file"]')[0].files[0];
		if (file === undefined || file.size > 102400) {
			$('#query_label').tooltip('show');
			return;
		}
		var reader = new FileReader();
		reader.onload = function () {
			$.post('job', {
				query: reader.result,
				filename: file.name.substr(0, 20),
				usr: $('input[name="usrRadioOptions"]:checked').val(),
			}, function(res) {
				location.replace('iview/?' + res);
			});
		};
		reader.readAsText(file);
	});
});
