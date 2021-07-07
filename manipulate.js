class ManipulateDisplay {
    constructor(canvas, extRenderCallback) {
        this.extRenderCallback = extRenderCallback;
        this.canvas = canvas;
        this.ctx = canvas.getContext('2d');
        canvas.style.cursor = 'pointer';

        this.width = canvas.width;
        this.height = canvas.height;
        this.handleTol = 15;
        this.scrollAmt = 0.1;

        // initial state
        this.sys = null;
        this.svd = null;
        this.cx = 0;
        this.cy = 0;
        this.z = 2;

        this.handles = [];

        this.mouseState = {
            // clicking
            mousedown: false,
            xdown: -1,          // coordinates of click
            ydown: -1,
            savex: 0.0,         // saved coordinates as needed
            savey: 0.0,
            // handle
            onhandle: false,    // whether mouse is over a handle
            handleclick: false, // whether handle was clicked
            handleid: -1,       // ID of clicked handle
        };
        var me = this;
        // mouse interaction functions
        canvas.onmousemove = function(e) {
            var xx = e.offsetX;
            var yy = e.offsetY;
            if (me.mouseState.mousedown) {
                // dragging
                let hid = me.mouseState.handleid;
                if (me.mouseState.handleclick) {
                    if (!me.isSecondHandle(hid) || !e.ctrlKey) {
                        let xy = me.screenToWorld(xx, yy);
                        me.changeHandleLoc(hid, xy[0], xy[1], me.pressedKeys(e));
                        me.updateFromHandle(hid);
                        me.update();
                        me.extRenderCallback(me.sys);
                    }
                }
                else {
                    let mxy = me.screenToWorld(me.mouseState.xdown, me.mouseState.ydown);
                    let xy = me.screenToWorld(xx, yy);
                    let dx = mxy[0] - xy[0];
                    let dy = mxy[1] - xy[1];
                    me.cx = me.mouseState.savex + dx;
                    me.cy = me.mouseState.savey + dy;
                    me.update();
                }
            }
            else if (me.sys !== null) {
                // check for nearby handles
                var mindist = me.width + me.height;
                var minloc = -1;
                for (let i = 0; i < me.handles.length; i++) {
                    let ss = me.handles[i].loc;
                    let xy = me.worldToScreen(ss[0], ss[1]);
                    let d = Math.hypot(xx - xy[0], yy - xy[1]);
                    if (d < mindist) {
                        mindist = d;
                        minloc = i;
                    }
                }
                if (mindist <= me.handleTol && me.mouseState.handleid != minloc) {
                    me.mouseState.onhandle = true;
                    me.mouseState.handleid = minloc;
                    me.update();
                }
                else if (mindist > me.handleTol && me.mouseState.onhandle) {
                    me.mouseState.onhandle = false;
                    me.mouseState.handleid = -1;
                    me.update();
                }
            }
        };
        canvas.onmousedown = function(e) {
            me.canvas.style.cursor = 'grabbing';
            me.mouseState.mousedown = true;
            var xx = e.offsetX;
            var yy = e.offsetY;
            me.mouseState.xdown = xx;
            me.mouseState.ydown = yy
            if (me.mouseState.onhandle) {
                me.mouseState.handleclick = true;
            }
            else {
                me.mouseState.savex = me.cx;
                me.mouseState.savey = me.cy;
            }
        };
        canvas.onmouseup = function(e) {
            if (me.mouseState.handleclick && me.isSecondHandle(me.mouseState.handleid) && e.ctrlKey) {
                me.changeHandleLoc(me.mouseState.handleid, 0, 0, me.pressedKeys(e));
                me.updateFromHandle(me.mouseState.handleid);
                me.update();
                me.extRenderCallback(me.sys);
            }
            me.mouseState.handleclick = false;
            me.mouseState.mousedown = false;
            me.mouseState.xdown = -1;
            me.mouseState.ydown = -1;
            me.mouseState.savex = 0;
            me.mouseState.savey = 0;
            me.canvas.style.cursor = 'pointer';
        };
        canvas.onwheel = function(e) {
            e.preventDefault();
            var delta = Math.sign(e.deltaY) * me.scrollAmt;
            me.z += delta;
            me.z = me.z < 0.5 ? 0.5 : me.z;
            me.update();
        }
    }

    pressedKeys(mouseE) {
        return {ctrl: mouseE.ctrlKey, shift: mouseE.shiftKey, alt: mouseE.altKey};
    }

    isSecondHandle(handleId) {
        return this.handles[handleId].param.includes('2');
    }

    changeHandleLoc(handleId, x, y, keys) {
        // x and y are in world coordinates
        var hd = this.handles[handleId];
        var loc = hd.loc;
        if (hd.param === 'T') {
            // translate all handles belonging to this system
            let delta = vecSub([x,y], loc);
            for (let i = 0; i < this.handles.length; i++)
                if (this.handles[i].sysindex === hd.sysindex)
                    this.handles[i].loc = vecAdd(this.handles[i].loc, delta);
        } else if (this.isSecondHandle(handleId)) {
            if (keys.ctrl) {
                // flip (180 degree rotation) this handle around its base
                let ref = this.handles[hd.ref];
                hd.loc = vecAddScaled(ref.loc, vecSub(loc, ref.loc), -1);
            } else {
                let isU = hd.param[0] === 'U';
                let base = this.handles[hd.ref];
                let hv = vecSub(hd.loc, base.loc);
                let other = this.handles[handleId + (isU ? 2 : -2)];
                let ov = vecSub(other.loc, base.loc);
                let sigother = vecNorm(ov);
                let max = isU ? 1 : sigother;
                let min = isU ? sigother : 0.01;
                let mousev = vecSub([x, y], base.loc);
                let d = vecInner(mousev, hv) / vecNorm(hv);
                d = d > max ? max : d;
                d = d < min ? min : d;
                hd.loc = vecAddScaled(base.loc, hv, d / vecNorm(hv));
                other = this.handles[handleId - 1];
                ov = vecSub(other.loc, base.loc);
                other.loc = vecAddScaled(base.loc, ov, d / vecNorm(ov));
            }
        } else {
            // rotate handles around their base
            let ref = this.handles[hd.ref];
            let xy = vecSub(loc, ref.loc);
            let mxy = vecSub([x, y], ref.loc);
            let theta = angleBetween(xy, mxy);
            for (let i = 0; i < this.handles.length; i++) {
                let hh = this.handles[i];
                if (hh.sysindex === hd.sysindex && hh.param !== 'T') {
                    if (!keys.ctrl && hh.param[0] != hd.param[0]) continue;
                    xy = vecSub(hh.loc, ref.loc);
                    xy = vecRotate(xy, theta);
                    xy = vecAdd(xy, ref.loc);
                    this.handles[i].loc = xy;
                }
            }
        }
    }

    updateFromHandle(handleId) {
        var h = this.handles[handleId];
        if (h.param === 'T') {
            this.sys[h.sysindex][0][2] = h.loc[0];
            this.sys[h.sysindex][1][2] = h.loc[1];
        } else {
            let svd = this.svd[h.sysindex];
            let base = this.handles[h.ref];
            for (let i = 0; i < this.handles.length; i++) {
                let hh = this.handles[i];
                if (hh.sysindex === h.sysindex && hh.param !== 'T') {
                    let p = hh.param[0];
                    let k = parseInt(hh.param[1], 10) - 1;
                    let v = vecSub(hh.loc, base.loc);
                    let d = vecNorm(v);
                    if (hh.param === 'U1') {
                        svd.S[0] = d;
                    } else if (hh.param === 'V1') {
                        svd.S[1] = d;
                    }
                    v = vecScale(v, 1 / d);
                    svd[p][0][k] = v[0];
                    svd[p][1][k] = v[1];
                }
            }
            let m = isvd(svd);
            for (let i = 0; i < 2; i++)
                for (let j = 0; j < 2; j++)
                    this.sys[h.sysindex][i][j] = m[i][j];
        }
    }

    worldToScreen(x, y) {
        var scx = this.width / 2;
        var scy = this.height / 2;
        var xs = (x - this.cx) * scx / this.z + scx;
        var ys = -(y - this.cy) * scy / this.z + scy;
        return [xs, ys];
    }

    screenToWorld(x, y) {
        var scx = this.width / 2;
        var scy = this.height / 2;
        var xw = (x - scx) * this.z / scx + this.cx;
        var yw = -(y - scy) * this.z / scy + this.cy;
        return [xw, yw];
    }

    drawLine(x1, y1, x2, y2, color=null) {
        var ctx = this.ctx;
        color = color === null ? 'black' : color;
        ctx.strokeStyle = color;
        ctx.fillStyle = color;
        ctx.beginPath();
        ctx.moveTo(x1, y1);
        ctx.lineTo(x2, y2);
        ctx.stroke();
    }

    drawCircle(x, y, r, fill=false, color=null) {
        var ctx = this.ctx;
        color = color === null ? 'black' : color;
        ctx.beginPath();
        ctx.arc(x, y, r, 0, Math.PI * 2);
        if (fill) {
            ctx.fillStyle = color;
            ctx.fill();
        } else {
            ctx.strokeStyle = color;
            ctx.stroke();
        }
    }

    drawArrow(x1, y1, x2, y2, a=7, fill=true, color=null) {
        var ctx = this.ctx;
        color = color === null ? 'black' : color;
        let r = Math.sqrt(3) / 6 * a;
        let x = x2 - x1;
        let y = y2 - y1;
        let d = Math.hypot(x, y);

        let xy = [x / d, y / d];
        let xyt = [-xy[1], xy[0]];
        let sp = [x1, y1];

        let p1 = vecAddScaled(sp, xy, 2*r);
        let p2 = vecAddScaled(sp, xy, -r)
        let p3 = vecAddScaled(p2, xyt, a/2);
        p2 = vecAddScaled(p2, xyt, -a/2);

        ctx.beginPath();
        ctx.moveTo(p1[0], p1[1]);
        ctx.lineTo(p2[0], p2[1]);
        ctx.lineTo(p3[0], p3[1]);
        ctx.lineTo(p1[0], p1[1]);
        if (fill) {
            ctx.fillStyle = color;
            ctx.fill();
        } else {
            ctx.strokeStyle = color;
            ctx.stroke();
        }
    }

    drawCoordinates() {
        var ctx = this.ctx;
        var fsize = 9;
        ctx.font = fsize.toString() + 'px Courier';
        ctx.fillStyle = '#F2F2F2';
        ctx.strokeStyle = 'white';
        ctx.lineWidth = 2;

        var w = this.width;
        var h = this.height;
        ctx.fillRect(0, 0, w, h);
        ctx.fillStyle = 'black';

        var ng = 8;
        for (let i = 1; i < ng; i++) {
            let u = w * i / ng;
            let v = h * i / ng;

            // horizontal
            ctx.textBaseline = 'middle';
            ctx.textAlign = 'left';
            ctx.beginPath();
            ctx.moveTo(0, v);
            ctx.lineTo(w, v);
            ctx.stroke();
            let loc = this.cy - (i - ng / 2) * 2 * this.z / ng;
            let loc_s = (loc > 0 ? '+' : loc == 0 ? ' ' : '') + loc.toFixed(2);
            ctx.fillText(loc_s, 2, v);

            // vertical
            ctx.textBaseline = 'top';
            ctx.textAlign = 'center';
            ctx.beginPath();
            ctx.moveTo(u, 0);
            ctx.lineTo(u, h);
            ctx.stroke();
            loc = this.cx + (i - ng / 2) * 2 * this.z / ng;
            loc_s = (loc > 0 ? '+' : loc == 0 ? ' ' : '') + loc.toFixed(2);
            ctx.fillText(loc_s, u, 2);
        }
        let orig = this.worldToScreen(0, 0);
        this.drawCircle(orig[0], orig[1], 3, true, '#555555');
    }

    drawHandles() {
        var ctx = this.ctx;
        var sys = this.sys;
        var w = this.width;
        var h = this.height;

        for (let i = 0; i < this.handles.length; i++) {
            ctx.lineWidth = 1;
            var hd = this.handles[i];
            var sc = this.worldToScreen(hd.loc[0], hd.loc[1]);

            let color;
            switch (hd.param[0]) {
                case 'T': color = 'black'; break;
                case 'U': color = 'red'; break;
                case 'V': color = 'blue'; break;
            }

            // draw
            var bs = null;
            if (hd.ref >= 0) {
                bs = this.handles[hd.ref].loc;
                bs = this.worldToScreen(bs[0], bs[1]);
                this.drawLine(bs[0], bs[1], sc[0], sc[1], color);
            }
            let isSecond = this.isSecondHandle(i) ? true : false; 
            let isHover = i === this.mouseState.handleid;
            if (isSecond) {
                if (isHover) {
                    this.drawArrow(sc[0], sc[1], bs[0], bs[1], 15, true, color);
                } else {
                    this.drawCircle(sc[0], sc[1], 4, true, color);
                }
            } else {
                ctx.lineWidth = isHover ? 3 : 1;
                this.drawCircle(sc[0], sc[1], 7, false, color);
            }
        }
    }

    addHandle(sysindex, param, loc, ref=-1) {
        this.handles.push({
            sysindex: sysindex, // index of system that the handle belongs to
            param: param,       // parameter of system that the handle manipulates
            loc: loc,           // location (world coords) of handle
            ref: ref            // index of another handle that a line should be drawn to
        });
    }

    update(sys=null, force=false) {
        if (sys === null)
            sys = this.sys;
        else if (force || this.sys !== sys) {
            this.sys = sys;
            this.svd = [];
            this.handles = [];
            let hPerS = 5;
            for (let i = 0; i < this.sys.length; i++) {
                let s = svd2x2(sys[i]);
                this.svd.push(s);
                let ref = hPerS*(i + 1) - 1;
                let v = getColumn(this.sys[i], 2);
                this.addHandle(i, 'U1', vecAddScaled(v, getColumn(s.U, 0), s.S[0]), ref);
                this.addHandle(i, 'U2', vecAddScaled(v, getColumn(s.U, 1), s.S[0]), ref);
                this.addHandle(i, 'V1', vecAddScaled(v, getColumn(s.V, 0), s.S[1]), ref);
                this.addHandle(i, 'V2', vecAddScaled(v, getColumn(s.V, 1), s.S[1]), ref);
                this.addHandle(i, 'T', v);
            }
        }

        this.drawCoordinates();
        this.drawHandles();
    }
}

// ------------------------------------------
// Linear Algebra functions (in 2 dimensions)
// ------------------------------------------
function angleBetween(v1, v2) {
    // return the angle between two vectors
    let dp = v1[0] * v2[0] + v1[1] * v2[1];
    let n1 = Math.hypot(v1[0], v1[1]);
    let n2 = Math.hypot(v2[0], v2[1]);
    let ct = dp / (n1 * n2);
    ct = ct > 1 ? 1 : ct;
    ct = ct < -1 ? -1 : ct;
    let theta = Math.acos(ct);
    let s = Math.sign(v1[0]*v2[1]-v1[1]*v2[0]);
    return s*theta;
}

function vecRotate(v, theta) {
    // rotate vector v by angle theta
    let c = Math.cos(theta);
    let s = Math.sin(theta);
    let q = [v[0]*c - v[1]*s, v[0]*s + v[1]*c];
    return q;
}

function vecScale(v, a) {
    // scale vector v by scalar a
    return [a * v[0], a * v[1]];
}

function vecAdd(v1, v2) {
    // add two vectors
    return [v1[0] + v2[0], v1[1] + v2[1]];
}

function vecAddScaled(v1, v2, a) {
    // v1 + a*v2
    return [v1[0] + a * v2[0], v1[1] + a * v2[1]];
}

function vecSub(v1, v2) {
    // add two vectors
    return [v1[0] - v2[0], v1[1] - v2[1]];
}

function vecNorm(v) {
    return Math.hypot(v[0], v[1]);
}

function vecInner(v1, v2) {
    // length of v1 in the direction of v2
    return v1[0] * v2[0] + v1[1] * v2[1];
}

function getColumn(m, i) {
    // return the ith column of matrix m as a vector
    return [m[0][i], m[1][i]];
}

function isvd(svd) {
    let U = svd.U;
    let S = svd.S;
    let V = svd.V;
    let m = [[
        S[0] * U[0][0] * V[0][0] + S[1] * U[0][1] * V[1][0],
        S[0] * U[0][0] * V[0][1] + S[1] * U[0][1] * V[1][1],
        ],[
        S[0] * U[1][0] * V[0][0] + S[1] * U[1][1] * V[1][0],
        S[0] * U[1][0] * V[0][1] + S[1] * U[1][1] * V[1][1],
    ]];
    return m;
}

function svd2x2(A) {
    // adapted from https://scicomp.stackexchange.com/a/11646
    var u, s, v;
    var a = A[0][0], b = A[0][1], c = A[1][0], d = A[1][1];
    const EPSILON = 1E-6;
    // If it is diagonal, SVD is trivial
    if (Math.abs(b - c) < EPSILON && Math.abs(b) < EPSILON) {
        v = [[a < 0 ? -1 : 1, 0], [0, d < 0 ? -1 : 1]];
        s = [Math.abs(a), Math.abs(d)];
        u = [[1.0, 0.0], [0.0, 1.0]];
    }
    // Otherwise, we need to compute A^T*A
    else {
        let j = a * a + b * b;
        let k = c * c + d * d;
        let u_c = a * c + b * d;
        // Check to see if A^T*A is diagonal
        if (Math.abs(u_c) < EPSILON) {
            let s1 = Math.sqrt(j);
            let s2 = Math.abs(j-k) < EPSILON ? s1 : Math.sqrt(k);
            s = [s1, s2];
            u = [[1.0, 0.0], [0.0, 1.0]];
            v = [[a / s1, b / s1], [c / s2, d / s2]];
        }
        // Otherwise, solve quadratic for eigenvalues
        else {
            let jmk = j - k;
            let jpk = j + k;
            let root = Math.sqrt(jmk*jmk + 4*u_c*u_c);
            let eig = (jpk+root)/2;
            let s1 = Math.sqrt(eig);
            let s2 = Math.abs(root) < EPSILON ? s1 : Math.sqrt((jpk - root)/2);
            s = [s1, s2];
            // Use eigenvectors of A^T*A as U
            let u_s = eig - j;
            let le = Math.sqrt(u_s*u_s + u_c*u_c);
            u_c /= le;
            u_s /= le;
            u = [[u_c, -u_s], [u_s, u_c]];
            // Compute V matrix as AU/s
            v = [
                [(a * u_c + c * u_s) / s1,
                 (b * u_c + d * u_s) / s1],
                [(c * u_c - a * u_s) / s2,
                 (d * u_c - b * u_s) / s2]
            ];
        }
    }
    return {U: u, S: s, V: v};
}