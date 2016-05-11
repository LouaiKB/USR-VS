$(function() {
	var embed_in = $('#embed_in');
	var embed_out = $('#embed_out');
	var embed_svg = $('#embed_svg');
	embed_out.hide();
	embed_svg.hide();
	var submit = $('#submit');
	submit.click(function() {
		submit.prop('disabled', true);
		$.post('/embed', {
			smiles: embed_in.val(),
		}, function(res) {
			if (res.tmp) {
				embed_out.attr('rows', res.embed_out.split(/\n/).length);
				embed_out.text(res.embed_out);
				embed_svg.attr('src', 'tmp/' + res.tmp);
				embed_svg.show();
			} else {
				embed_svg.hide();
				embed_out.attr('rows', 1);
				embed_out.text(JSON.stringify(res));
			}
			embed_out.show();
		}).always(function() {
			submit.prop('disabled', false);
		});
	});
});
