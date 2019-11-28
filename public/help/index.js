$(() => {
	var embed_in = $('#embed_in');
	var embed_out = $('#embed_out');
	embed_out.hide();
	const smilesDrawer = new SmilesDrawer.Drawer({
		width: 540,
		height: 540,
	});
	var submit = $('#submit');
	submit.click(() => {
		submit.prop('disabled', true);
		$.post('/embed', {
			smiles: embed_in.val(),
		}, (res) => {
			if (res.embed_out) { // embed exited normally.
				embed_out.attr('rows', res.embed_out.split(/\n/).length);
				embed_out.text(res.embed_out);
				SmilesDrawer.parse(embed_in.val(), (tree) => { // SmilesDrawer.parse() is a static function.
					smilesDrawer.draw(tree, 'drawer');
				}, (err) => {
					// TODO: noty()
				});
			} else { // embed did not run normally and exited with either a code or a signal.
				embed_out.attr('rows', 1);
				embed_out.text(JSON.stringify(res));
			}
			embed_out.show();
		}).always(() => {
			submit.prop('disabled', false);
		});
	});
});
