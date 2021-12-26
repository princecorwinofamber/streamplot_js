var streamplotjs = {};

streamplotjs._load_odex = function() {
    var exports = {};
    function onLoad() {
        eval(this.responseText);
    };
    var oReq = new XMLHttpRequest();
    oReq.addEventListener("load", onLoad);
    oReq.open("GET", "https://raw.githubusercontent.com/littleredcomputer/odex-js/master/src/odex.js");
    oReq.send();
    return exports;
}

var odex = streamplotjs._load_odex();

streamplotjs._det2x2 = function(a11, a12, a21, a22) {
    return a11 * a22 - a12 * a21;
}

streamplotjs._kramer2x2 = function(a11, a12, a21, a22, b1, b2) {
    var mdet = streamplotjs._det2x2(a11, a12, a21, a22);
    if (Math.abs(mdet) < 1e-4) {
        return null;
    }
    return [streamplotjs._det2x2(b1, a12, b2, a22) / mdet, streamplotjs._det2x2(a11, b1, a21, b2) / mdet];
}

streamplotjs.ODESystem2d = function(dot_x, dot_y, root_threshold=1e-3, newton_max_steps=10) {
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
        return new streamplotjs.ODESystem2d(lin_f, lin_g);
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
        while((Math.abs(this._f_comp(x, y)) > this._root_threshold || Math.abs(this._g_comp(x, y)) > this._root_threshold) && iter < this._newton_max_steps) {
            res = streamplotjs._kramer2x2(this._jacobi_tl_comp(x, y), this._jacobi_tr_comp(x, y), this._jacobi_bl_comp(x, y), this._jacobi_br_comp(x, y),
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

streamplotjs.PlottingCanvas = function(canvas_id, min_x, min_y, max_x, max_y) {
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
    this.get_color = function() {
        return {fill: this._ctx.fillStyle, stroke: this._ctx.strokeStyle};
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
        //this._ctx.fillRect(t.x, t.y, 20, 20);
        this._ctx.beginPath();
        this._ctx.arc(t.x, t.y, 10, 0, 2 * Math.PI, false);
        this._ctx.fill();
    };
    this.clear = function() {
        this._ctx.clearRect(0, 0, this._canvas.width, this._canvas.height);
    }

    this.drawArrow = function(curx, cury, dx, dy) {
        this._ctx.beginPath();
        let t = this._transform(curx, cury);
        curx = t.x;
        cury = t.y;
        this._ctx.moveTo(curx, cury);
        this._ctx.lineTo(curx-dx*0.5-dy*0.1, cury-dy*0.5+dx*0.1);
        this._ctx.lineTo(curx-dx*0.5+dy*0.1, cury-dy*0.5-dx*0.1);
        this._ctx.closePath();
        this._ctx.fill();
    }

    this.update = function(min_x, min_y, max_x, max_y) {
        this.min_x = min_x;
        this.min_y = min_y;
        this.max_x = max_x;
        this.max_y = max_y;
        this._x_scale = this._width / (max_x - min_x);
        this._y_scale = this._height / (max_y - min_y);
        this._x_shift = (max_x + min_x) / 2;
        this._y_shift = (max_y + min_y) / 2;
    }
};

streamplotjs.meshgrid = function(x_arr, y_arr) {
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

streamplotjs.linspace = function(start, stop, num) {
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

streamplotjs.default_palette = ["#40db40", "#db8040", "#cc40db", "#2b2b28", "#24c793", "#909196", "#150999", "#7b7a80", "#a10306", "#7e9635"];

streamplotjs.Streamplot = function(f, startx=-3, starty=-1.7, endx=3, endy=1.7, stepx=0.5, stepy=0.5, canvasName="canvas", color="RoyalBlue", iter_steps=1000, arrowSize=30, lr=0.001, zoomstep=2, trajectoryColor="#edcb42", startTime = 0, endTime=1000, root_search_discretization_x=25, root_search_discretization_y=3, palette=streamplotjs.default_palette) {
    this.startx = startx;
    this.starty = starty;
    this.endx = endx;
    this.endy = endy;
    this.stepx = stepx;
    this.stepy = stepy;
    this.canvasName = canvasName;
    this.f = f;
    this.iter_steps = iter_steps;
    this.arrowSize = arrowSize;
    this.lr = lr;
    this.ctx = new streamplotjs.PlottingCanvas(canvasName, startx, starty, endx, endy);
    this.canvas = this.ctx._canvas;
    this.mx=0;
    this.my=0;
    this.redraw = false;
    this.zoom_step = zoomstep;
    this.trajectoryPoints = [];
    this.trajectoryColor = trajectoryColor;
    this.mode = true;
    this.solver = false;
    this.startTime = startTime;
    this.endTime = endTime;
    this.root_search_discretization_x = root_search_discretization_x;
    this.root_search_discretization_y = root_search_discretization_y;
    this.palette = palette;
    this.equilibrium_points_output = null;
    
    this.setColor = function (color) {
        this.color = color;
        this.ctx.set_color(this.color);
    }

    this.setColor(color);

    this.streamplot = function () {
        for (var x = this.startx-this.startx%this.stepx-this.stepx; x < this.endx+this.stepx; x += this.stepx) {
            for (var y = this.starty-this.starty%this.stepy-this.stepy; y < this.endy+this.stepy; y += this.stepy) {
                this.ctx.beginPath();
                curx = x;
                cury = y;
                this.ctx.moveTo(curx, cury);
                for (var i = 0; i < iter_steps; i++) {
                    aaa = this.f.right_part(curx, cury);
                    curx += aaa[0]*lr;
                    cury += aaa[1]*lr;
                    this.ctx.lineTo(curx, cury);
                }
                this.ctx.stroke();
                var dx = aaa[0]*lr;
                var dy = aaa[1]*lr;
                var l = Math.sqrt(dx*dx + dy*dy);
                dx = arrowSize*dx/l;
                dy = arrowSize*dy/l;
                this.ctx.drawArrow(curx, cury, dx, -dy);
            }
        }
    }

    this.trajectoryPlot = function () {
        var start = true;
        this.ctx.set_color(this.trajectoryColor);
        this.ctx._ctx.lineWidth = 3;
        this.ctx.beginPath();
        for (i in this.trajectoryPoints) {
            point = this.trajectoryPoints[i];
            if (start) {
                start = false;
                this.ctx.moveTo(point.x, point.y);
            } else {
                this.ctx.lineTo(point.x, point.y);
            }
        }
        this.ctx.stroke();
        this.ctx.set_color(this.color);
        this.ctx._ctx.lineWidth = 1;
    }

    this.equlibriumPlot = function() {
        if (this.equilibrium_points_output != null) {
            document.getElementById(this.equilibrium_points_output).innerHTML = "";
        }
        var last_color = this.color;
        var cur_color_ind = 0;
        var equilibrium_points = this.f.equilibrium_points(streamplotjs.linspace(this.startx, this.endx, this.root_search_discretization_x), streamplotjs.linspace(this.starty, this.endy, this.root_search_discretization_y));
        for (var equilibrium_point_ind = 0; equilibrium_point_ind < equilibrium_points.length; equilibrium_point_ind++) {
            var point = equilibrium_points[equilibrium_point_ind];
            if ((point.x < this.startx) || (point.x > this.endx) || (point.y < this.starty) || (point.y > this.endy)) {
                continue;
            }
            /* var new_color = null;
            if (cur_color_ind < this.palette.length) {
                new_color = this.palette[cur_color_ind];
            } else {
                new_color = "#123456";
            } */
            var new_color = "hsl(".concat(Math.round(360*(point.y-this.starty)/(this.endy-this.starty)).toString()).concat(", 100%, ").concat(Math.round(80*(point.x-this.startx)/(this.endx-this.startx)).toString()).concat("%)");
            //console.log(new_color);
            cur_color_ind++;
            this.setColor(new_color);
            this.ctx.addDot(point.x, point.y);
            if (this.equilibrium_points_output != null) {
                var span = document.getElementById(this.equilibrium_points_output).appendChild(document.createElement("span"));
                span.setAttribute("style", "margin-top: 0; color:".concat(new_color));
                span.textContent = this.f.focus_type(point.x, point.y);
                document.getElementById(this.equilibrium_points_output).appendChild(document.createElement("br"));
            }
        };
        this.setColor(last_color);
    };

    this.draw = function () {
        this.ctx.clear();
        this.streamplot();
        this.trajectoryPlot();
        this.equlibriumPlot();
    }

    this.switchMode = function() {
        this.mode = !this.mode;
    }

    this.addLongTrajectory = function(posx, posy) {
        if (this.solver === false) {
            this.solver = new odex.Solver(2);
            this.solver.maxStepSize = this.lr;
            this.solver.maxSteps = 30000;
        }
        this.trajectoryPoints = [];
        this.solver.solve((t, x) => {return this.f.right_part(x[0], x[1])}, this.startTime, [posx, posy], this.endTime, (n, x0, x1, X) => { this.trajectoryPoints.push({x:X[0], y:X[1]}) });
        this.draw();
    };

    this.canvas.addEventListener('mousedown', e => {
        this.mx = e.offsetX;
        this.my = e.offsetY;
        if (this.mode) {
            this.redraw = true;
        } else {
            var posx = (this.mx/this.canvas.width)*(this.endx-this.startx)+this.startx;
            var posy = (1-this.my/this.canvas.height)*(this.endy-this.starty)+this.starty;
            this.addLongTrajectory(posx, posy);
        }
    });

    this.canvas.addEventListener('mousemove', e => {
        var dx = e.offsetX-this.mx;
        var dy = e.offsetY-this.my;
        if (this.redraw === true) {
            dx = (this.endx-this.startx)*dx/this.canvas.width;
            dy = (this.endy-this.starty)*dy/this.canvas.height;
            this.startx -= dx;
            this.endx -= dx;
            this.starty += dy;
            this.endy += dy;
            this.mx = e.offsetX;
            this.my = e.offsetY;
            this.ctx.update(this.startx, this.starty, this.endx, this.endy);
            this.draw();
        }
    });

    window.addEventListener('mouseup', e => {
        this.redraw = false;
    });

    this.canvas.addEventListener('wheel', e => {
        var width = this.endx - this.startx;
        var height = this.endy - this.starty;
        var ratio = height/width;
        var scale = width/this.canvas.width;
        scale *= 1+e.deltaY*0.0005;
        width = this.canvas.width*scale;
        height = width*ratio;
        var centerx = (this.endx + this.startx)/2;
        var centery = (this.endy + this.starty)/2;
        this.startx = centerx - width/2;
        this.endx = centerx + width/2;
        this.starty = centery - height/2;
        this.endy = centery + height/2;

        while (this.stepx * this.zoom_step < width/10) {
            this.stepx *= this.zoom_step;
            this.stepy *= this.zoom_step;
        }
        while (this.stepx > width*this.zoom_step/10) {
            this.stepx /= this.zoom_step;
            this.stepy /= this.zoom_step;
        }

        this.ctx.update(this.startx, this.starty, this.endx, this.endy);
        this.draw();
    });

    this.set_equilibrium_points_output = function(div_id) {
        var old_equilibrium_points_output = this.equilibrium_points_output;
        this.equilibrium_points_output = div_id;
        if (old_equilibrium_points_output != null) {
            document.getElementById(old_equilibrium_points_output).innerHTML = "";
        }
    };
}

