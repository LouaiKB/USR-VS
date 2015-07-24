(function($) {

	var tr, tds_empty, tbody_trs, frstpage, lastpage, prevpage, nextpage, whatpage;
	var num_pages, page = 0;
	var methods = {
		init: function(thead, tr_arg) {
			// Set up tr function.
			var cols = thead.length, i;
			tds_empty = new Array(cols);
			for (i = 0; i < cols; ++i) tds_empty[i] = '&nbsp;';
			tr = function(record) {
				if (!record) return tds_empty;
				return tr_arg(record);
			};
			// Set up table
			var tr_tag = '<tr>';
			for (i = 0; i < cols; ++i) tr_tag += '<td>&nbsp;</td>';
			tr_tag += '</tr>';
			var tbody = '';
			for (i = 0; i < 8; ++i) tbody += tr_tag;
			this.html('<div class="table-responsive"><table class="table table-hover"><thead><tr><th>' + thead.join('</th><th>') + '</th></tr></thead><tbody>' + tbody + '</tbody></table></div><div class="text-center"><ul class="pagination"><li id="frstpage"><a>&laquo;&laquo;</a></li><li id="prevpage"><a>&laquo;</a></li><li class="form-group"><input type="text" id="whatpage" class="form-control"></li><li id="nextpage"><a>&raquo;</a></li><li id="lastpage"><a>&raquo;&raquo;</a></li></ul></div>');
			tbody_trs = $('tbody tr', this);
			// Set up pager events.
			frstpage = $('#frstpage', this);
			lastpage = $('#lastpage', this);
			prevpage = $('#prevpage', this);
			nextpage = $('#nextpage', this);
			whatpage = $('#whatpage', this);
			frstpage.click(function() {
				if (frstpage.hasClass('disabled')) return;
				page = 1;
				render();
			});
			lastpage.click(function() {
				if (lastpage.hasClass('disabled')) return;
				page = num_pages;
				render();
			});
			prevpage.click(function() {
				if (prevpage.hasClass('disabled')) return;
				page -= 1;
				render();
			});
			nextpage.click(function() {
				if (nextpage.hasClass('disabled')) return;
				page += 1;
				render();
			});
			whatpage.change(function() {
				var p = parseInt(whatpage.val());
				if ((p < 1) || (p > num_pages) || (p === page)) return;
				page = p;
				render();
			});
			return this;
		}, source: function(records_arg) {
			records = records_arg;
			num_pages = records.length ? ((records.length + 7) >> 3) : 1;
			page = num_pages;
			render();
			return this;
		}, refresh: function(row_lb, row_ub, col_lb, col_ub, fade) {
			var offset = 8 * (page - 1);
			tbody_trs.each(function(row) {
				row += offset;
				if (!(row_lb <= row && row < row_ub)) return;
				var tds = tr(records[row]);
				$('td', this).each(function(col) {
					if (!(col_lb <= col && col < col_ub)) return;
					var t = $(this);
					if (t.html() !== tds[col]) t.html(tds[col]);
					if (fade) t.hide().fadeIn('slow');
				});
			});
			return this;
		}
	};

	$.fn.pager = function(method) {
		return methods[method].apply(this, Array.prototype.slice.call(arguments, 1));
	};

	// Render the table and its pager
	function render() {
		var offset = 8 * (page - 1);
		tbody_trs.each(function(i) {
			var tds = tr(records[offset + i]);
			$('td', this).each(function(j) {
				$(this).html(tds[j]);
			});
		});
		if (page === 1) {
			frstpage.addClass('disabled');
			prevpage.addClass('disabled');
		} else {
			frstpage.removeClass('disabled');
			prevpage.removeClass('disabled');
		}
		if (page === num_pages) {
			lastpage.addClass('disabled');
			nextpage.addClass('disabled');
		} else {
			lastpage.removeClass('disabled');
			nextpage.removeClass('disabled');
		}
		whatpage.val(page);
	}

})(jQuery);
