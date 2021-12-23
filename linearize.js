function ODESystem2d(dot_x, dot_y) {
	// \dot x = f(x, y)
	// \dot y = g(x, y)
	this._f = nerdamer(dot_x);
	this._g = nerdamer(dot_y);
	this._f_comp = this._f.buildFunction(["x", "y"]);
	this._g_comp = this._g.buildFunction(["x", "y"]);
	this._jacobi_tl = nerdamer.diff(this._f, "x");
	this._jacobi_tr = nerdamer.diff(this._f, "y");
	this._jacobi_bl = nerdamer.diff(this._g, "x");
	this._jacobi_br = nerdamer.diff(this._g, "y");
	this._jacobi_tl_comp = this._jacobi_tl.buildFunction(["x", "y"]);
	this._jacobi_tr_comp = this._jacobi_tr.buildFunction(["x", "y"]);
	this._jacobi_bl_comp = this._jacobi_bl.buildFunction(["x", "y"]);
	this._jacobi_br_comp = this._jacobi_br.buildFunction(["x", "y"]);

	this.linearized = function(at_x, at_y) {
		var lin_f = this._jacobi_tl.evaluate({ x: at_x, y: at_y }).multiply(nerdamer("x").subtract(nerdamer(at_x))).add(this._jacobi_tr.evaluate({ x: at_x, y: at_y }).multiply(nerdamer("y").subtract(nerdamer(at_y))));
		var lin_g = this._jacobi_bl.evaluate({ x: at_x, y: at_y }).multiply(nerdamer("x").subtract(nerdamer(at_x))).add(this._jacobi_br.evaluate({ x: at_x, y: at_y }).multiply(nerdamer("y").subtract(nerdamer(at_y))));
		return new ODESystem2d(lin_f, lin_g);
	};

	this.focus_type = function(at_x, at_y) {
		var eigval = nerdamer.solve(this._jacobi_tl.subtract(nerdamer("eig")).multiply(this._jacobi_br.subtract(nerdamer("eig"))).subtract(this._jacobi_tr.multiply(this._jacobi_bl)).evaluate({ x: at_x, y: at_y }), "eig");
		if (nerdamer.imagpart(nerdamer.vecget(eigval, 0)).eq("0")) {
			if (nerdamer.size(eigval).eq("1")) {
				if (this._jacobi_tr.evaluate({ x: at_x, y: at_y }).eq("0") && this._jacobi_bl.evaluate({ x: at_x, y: at_y }).eq("0")) {
					return "proper node"; // дикритический узел
				} else {
					return "degenerate node"; // вырожденный узел
				}
			} else {
				if (nerdamer.sign(nerdamer.vecget(eigval, 0)).eq(nerdamer.sign(nerdamer.vecget(eigval, 1)))) {
					return "node";
				} else {
					return "saddle";
				}
			}
		} else {
			return "center/focus";
			/* if (nerdamer.realpart(nerdamer.vecget(eigval, 0)).eq("0")) {
				return "center (or focus for nonlinear)";
			} else {
				return "focus (or center for nonlinear)";
			} */
		}
	};
}