function _load_odex() {
	var exports = {};
	function onLoad() {
		//console.log(this.responseText);
		eval(this.responseText);
		console.log("loaded");
	};
	var oReq = new XMLHttpRequest();
	oReq.addEventListener("load", onLoad);
	oReq.open("GET", "https://raw.githubusercontent.com/littleredcomputer/odex-js/master/src/odex.js");
	oReq.send();
	return exports;
}

var odex = _load_odex();
console.log("loading");