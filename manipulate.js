class ManipulateDisplay {
    constructor(canvas) {
        this.canvas = canvas;
        this.ctx = canvas.getContext('2d');

        this.width = canvas.width;
        this.height = canvas.height;

        // initial state
        this.cx = 0;
        this.cy = 0;
        this.z = 2;
        this.sys = null;

        // set up interactive
        var me = this;
        canvas.onmousemove = function(e) {
            if (me.sys !== null) {
                var xx = e.offsetX;
                var yy = e.offsetY;
                // var xy = me.screenToWorld(xx, yy);
                var mindist = me.width + me.height;
                var minloc = -1;
                for (let i = 0; i < me.sys.length; i++) {
                    let ss = me.sys[i];
                    let xy = me.worldToScreen(ss[0][2], ss[1][2])
                    let d = Math.hypot(xx - xy[0], yy - xy[1]);
                    if (d < mindist) {
                        mindist = d;
                        minloc = i;
                    }
                }
                if (mindist < 10) {
                    console.log('Near to sys' + minloc);
                }
            }
        };
    }

    worldToScreen(x, y) {
        var scx = this.width / 2;
        var scy = this.height / 2;
        var xs = (x + this.cx) * scx / this.z + scx;
        var ys = -(y + this.cy) * scy / this.z + scy;
        return [xs, ys];
    }

    screenToWorld(x, y) {
        var scx = this.width / 2;
        var scy = this.height / 2;
        var xw = (x - scx) * this.z / scx - this.cx;
        var yw = (y - scy) * this.z / scy - this.cy;
        return [xw, yw];
    }

    drawVec(x1, y1, x2, y2, color=null) {
        var ctx = this.ctx;
        color = color === null ? 'black' : color;
        ctx.strokeStyle = color;
        ctx.beginPath();
        ctx.moveTo(x1, y1);
        ctx.lineTo(x2, y2);
        ctx.stroke();

        ctx.beginPath();
        ctx.arc(x2, y2, 5, 0, Math.PI * 2);
        ctx.stroke();
    }

    drawOrigin(x, y, color=null) {
        var ctx = this.ctx;
        color = color === null ? 'black' : color;
        ctx.strokeStyle = color;
        ctx.beginPath();
        ctx.arc(x, y, 7, 0, Math.PI * 2)
        ctx.stroke();
    }

    drawCoordinates() {
        var ctx = this.ctx;
        ctx.fillStyle = '#F5F5F5';
        ctx.strokeStyle = 'white';
        var fsize = 10;
        ctx.font = fsize.toString() + 'px Courier';

        var w = this.width;
        var h = this.height;
        ctx.fillRect(0, 0, w, h);
        ctx.fillStyle = 'black';

        var ng = 8;
        for (let i = 1; i < ng; i++) {
            let u = w * i / ng;
            let v = h * i / ng;

            // horizontal
            ctx.beginPath();
            ctx.moveTo(0, v);
            ctx.lineTo(w, v);
            ctx.stroke();
            let loc = this.cy + (i - ng / 2) * 2 * this.z / ng;
            let loc_s = (loc > 0 ? '+' : loc == 0 ? ' ' : '') + loc.toFixed(2);
            ctx.fillText(loc_s, 2, v + 0.25*fsize);

            // vertical
            ctx.beginPath();
            ctx.moveTo(u, 0);
            ctx.lineTo(u, h);
            ctx.stroke();
            loc = this.cx + (i - ng / 2) * 2 * this.z / ng;
            loc_s = (loc > 0 ? '+' : loc == 0 ? ' ' : '') + loc.toFixed(2);
            ctx.fillText(loc_s, u - 2*fsize, 2+fsize*(1 + (i+1) % 2));
        }
    }

    showSystemDecomp(sys) {
        this.sys = sys;
        var ctx = this.ctx;

        this.drawCoordinates();

        var w = this.width;
        var h = this.height;

        ctx.moveTo(w, 0);
        ctx.lineTo(w, h);
        ctx.stroke();

        for (let i = 0; i < sys.length; i++) {
            var svd = svd2x2(sys[i]);
            var s1 = svd.S[0];
            var s2 = svd.S[1];
            var bx = sys[i][0][2];
            var by = sys[i][1][2];
            var bw = this.worldToScreen(bx, by, w, h);

            // draw
            if (bw[0] >= 0 && bw[0] < w && bw[1] >= 0 && bw[1] < h) {
                // left singular vectors (U)
                var xy = this.worldToScreen(bx + svd.U[0][0] * s1, by + svd.U[1][0] * s2);
                this.drawVec(bw[0], bw[1], xy[0], xy[1], 'red');
                xy = this.worldToScreen(bx + svd.U[0][1] * s1, by + svd.U[1][1] * s2);
                this.drawVec(bw[0], bw[1], xy[0], xy[1], 'red');
                // right singular vectors (V^T)
                xy = this.worldToScreen(bx + svd.V[0][0], by + svd.V[1][0], w, h);
                this.drawVec(bw[0], bw[1], xy[0], xy[1], 'blue');
                xy = this.worldToScreen(bx + svd.V[0][1], by + svd.V[1][1], w, h);
                this.drawVec(bw[0], bw[1], xy[0], xy[1], 'blue');
                this.drawOrigin(bw[0], bw[1]);
            }
                
        }
    }

}

function svd2x2(A) {
    var a = A[0][0];
    var b = A[0][1];
    var c = A[1][0];
    var d = A[1][1];

    var x, y, z, c1, c2, s1, s2, d1, d2;

    // Calculate RQ decomposition of A
    if (c == 0) {
        x = a;
        y = b;
        z = d;
        c2 = 1;
        s2 = 0;
    }
    else {
        let maxden = Math.max(Math.abs(c), Math.abs(d));

        let rcmaxden = 1/maxden;
        c = c * rcmaxden;
        d = d * rcmaxden;

        let den = 1 / Math.sqrt(c*c + d*d);

        let numx = (-b*c + a*d);
        let numy = (a*c + b*d);
        x = numx * den;
        y = numy * den;
        z = maxden/den;

        s2 = -c * den;
        c2 = d * den;
    }

    // Calculate tangent of rotation on R[x,y0,z] to diagonalize R^T*R
    let scaler = 1 / Math.max(Math.abs(x), Math.abs(y));
    let x_ = x*scaler;
    let y_ = y*scaler;
    let z_ = z*scaler;
    let numer = ((z_-x_)*(z_+x_)) + y_*y_;
    let gamma = x_*y_;
    gamma = numer == 0 ? 1 : gamma;
    var t = 0;
    if (gamma != 0) {
        let zeta = numer / gamma;
        t = 2 * (zeta >= 0 ? 1 : -1) / (Math.abs(zeta) + Math.sqrt(zeta*zeta+4));
    }

    // Calculate sines and cosines
    c1 = 1 / Math.sqrt(1 + t*t);
    s1 = c1*t;

    // Calculate U*S = R*R(c1,s1)
    var usa = c1*x - s1*y ;
    var usb = s1*x + c1*y;
    var usc = -s1*z;
    var usd = c1*z;

    // Update V = R(c1,s1)^T*Q
    t = c1*c2 + s1*s2;
    s2 = c2*s1 - c1*s2;
    c2 = t;

    // Separate U and S
    d1 = Math.hypot(usa, usc);
    d2 = Math.hypot(usb, usd);
    var dmax = Math.max(d1, d2);
    var usmax1 = d2 > d1 ? usd : usa;
    var usmax2 = d2 > d1 ? usb : -usc;

    var signd1 = x*z >= 0 ? 1 : -1;
    dmax *= d2 > d1 ? signd1 : 1;
    d2 *= signd1;
    var rcpdmax = 1 / dmax;

    c1 = dmax != 0 ? usmax1 * rcpdmax : 1;
    s1 = dmax != 0 ? usmax2 * rcpdmax : 0;

    U = [[c1, s1], [-s1, c1]];
    S = [d1, d2];
    V = [[c2, -s2], [s2, c2]];

    return {"U": U, "S": S, "V": V};
}