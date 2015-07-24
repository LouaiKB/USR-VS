importScripts('gunzip.min.js');

self.addEventListener('message', function (e) {
	var srcz = e.data;
	var srczr = new Uint8Array(srcz.length);
	for (var i = 0, l = srcz.length; i < l; ++i) {
		srczr[i] = srcz.charCodeAt(i);
	}
	var srcr = new Zlib.Gunzip(srczr).decompress();
	var src = '';
	for (var i = 0, l = srcr.length; i < l; ++i) {
		src += String.fromCharCode(srcr[i]);
	}
    self.postMessage(src);
});
