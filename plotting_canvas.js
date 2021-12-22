var PlottingCanvas = function(canvas_id, min_x, min_y, max_x, max_y) {
	this._canvas_id = canvas_id;
	this.min_x = min_x;
	this.min_y = min_y;
	this.max_x = max_x;
	this.max_y = max_y;
	this._canvas = document.getElementById(this._canvas_id);
	this._ctx = this._canvas.getContext("2d");
	this._height = this._canvas.height;
	this._width = this._canvas.width;
	this._x_scale = this._width / (max_x - min_x);
	this._y_scale = this._height / (max_y - min_y);
	this._x_shift = (max_x + min_x) / 2;
	this._y_shift = (max_y + min_y) / 2;
	this._transform = function(x, y) {
		return {
			x: this._width / 2 + (x - this._x_shift) * this._x_scale,
			y: this._height / 2 - (y - this._y_shift) * this._y_scale
		};
	};
	this._vector_transform = function(x, y) {
		return {
			x: x * this._x_scale,
			y: y * this._y_scale
		};
	};
	this.set_color = function(color) {
		this._ctx.fillStyle = color;
		this._ctx.strokeStyle = color;
	};
	this.beginPath = function() {
		this._ctx.beginPath();
	};
	this.closePath = function() {
		this._ctx.closePath();
	};
	this.fill = function() {
		this._ctx.fill();
	};
	this.stroke = function() {
		this._ctx.stroke();
	}
	this.moveTo = function(x, y) {
		let t = this._transform(x, y);
		this._ctx.moveTo(t.x, t.y);
	};
	this.lineTo = function(x, y) {
		let t = this._transform(x, y);
		this._ctx.lineTo(t.x, t.y);
	};
	this.addDot = function(x, y) {
		let t = this._transform(x, y);
		this._ctx.fillRect(t.x, t.y, 20, 20)
	};
};

var meshgrid = function(x_arr, y_arr) {
	var x_len = x_arr.length;
	var y_len = y_arr.length;
	X = Array(y_len).fill([...x_arr]);
	Y = Array(y_len);
	for (let i = 0; i < y_len; i++) {
		Y[i] = Array(x_len);
		Y[i].fill(y_arr[i]);
	}
	return [X, Y];
};

var linspace = function(start, stop, num) {
	var step = (stop - start) / (num - 1);
	var lin = Array(num);
	var cur = start;
	lin[0] = cur;
	for (let i = 1; i < num; i++) {
		cur += step;
		lin[i] = cur;
	}
	return lin;
};