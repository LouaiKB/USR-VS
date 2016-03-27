$(function() {
	var embed_in = $('#embed_in');
	var embed_out = $('#embed_out');
	embed_out.hide();
	var submit = $('#submit');
	submit.click(function() {
		$.post('/embed', {
			smiles: embed_in.val(),
		}, function(res) {
			if (typeof res === 'object') {
				embed_out.attr('rows', 1);
				embed_out.text(JSON.stringify(res));
			} else {
				embed_out.attr('rows', res.split(/\n/).length);
				embed_out.text(res);
			}
			embed_out.show();
		});
	});
});
