function _det2x2(a11, a12, a21, a22) {
    return a11 * a22 - a12 * a21;
}

function _kramer2x2(a11, a12, a21, a22, b1, b2) {
    var mdet = _det2x2(a11, a12, a21, a22);
    if (Math.abs(mdet) < 1e-4) {
        return null;
    }
    return [_det2x2(b1, a12, b2, a22) / mdet, _det2x2(a11, b1, a21, b2) / mdet];
}

function ODESystem2d(dot_x, dot_y, root_threshold=1e-3, newton_max_steps=10) {
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
    this._root_threshold = root_threshold;
    this._newton_max_steps = newton_max_steps;

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

    this.right_part = function(x, y) {
        return [this._f_comp(x, y), this._g_comp(x, y)];
    };

    this._equilibrium_point = function() {
        var threshold = this._root_threshold;
        try {
            nerdamer.set('SOLUTIONS_AS_OBJECT', true);
            var solution = nerdamer.solveEquations([this._f.toString(), this._g.toString()]);
            // if (this._f.evaluate({x: solution.x, y: solution.y}).eq("0") && this._g.evaluate({x: solution.x, y: solution.y}).eq("0")) {
            if (nerdamer.abs(this._f.evaluate({x: solution.x, y: solution.y})).lt(threshold) && nerdamer.abs(this._g.evaluate({x: solution.x, y: solution.y})).lt(threshold)) {
                return solution;
            } else {
                return null;
            }
        } catch (err) {
            return null;
        }
    };

    this._newtons_method = function(x0, y0) {
        var x = x0;
        var y = y0;
        var iter = 0;
        var res;
        while(Math.abs(this._f_comp(x, y)) > this._root_threshold && Math.abs(this._g_comp(x, y)) > this._root_threshold && iter < this._newton_max_steps) {
            res = _kramer2x2(this._jacobi_tl_comp(x, y), this._jacobi_tr_comp(x, y), this._jacobi_bl_comp(x, y), this._jacobi_br_comp(x, y),
                -this._f_comp(x, y) + this._jacobi_tl_comp(x, y) * x + this._jacobi_tr_comp(x, y) * y,
                -this._g_comp(x, y) + this._jacobi_bl_comp(x, y) * x + this._jacobi_br_comp(x, y) * y);
            if (res == null) {
                return null;
            }
            x = res[0];
            y = res[1];
            iter++;
        }
        if (Math.abs(this._f_comp(x, y)) <= this._root_threshold && Math.abs(this._g_comp(x, y)) <= this._root_threshold) {
            return res;
        } else {
            return null;
        }
    };

    this.equilibrium_points = function(x_arr, y_arr) {
        var p_list = [];
        var res;
        /* var res = this._equilibrium_point();
        if (res != null) {
            p_list.push({x: res[0], y: res[1]});
        } */
        //console.log(p_list[0]);
        for (var i = 0; i < x_arr.length; i++) {
            for (var j = 0; j < x_arr.length; j++) {
                res = this._newtons_method(x_arr[i], y_arr[j]);
                if (res != null) {
                    var is_root_new = true;
                    for (var k = 0; k < p_list.length; k++) {
                        if ((res[0] - p_list[k].x) * (res[0] - p_list[k].x) +  (res[1] - p_list[k].y) * (res[1] - p_list[k].y) < this._root_threshold) {
                            is_root_new = false;
                            break;
                        }
                    }
                    if (is_root_new) {
                        p_list.push({x: res[0], y: res[1]});
                    }
                }
            }
        }
        return p_list;
    };
}