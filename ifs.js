function randomSystem(n=null, bval=1) {
    if (n === null)
        n = [2, 5];
    if (Array.isArray(n)) {
        n = Math.floor(Math.random() * (n[1] - n[0])) + n[0];   
    }

    var matrices = [];
    for (let i=0; i<n; i++) {
        var sv1 = (Math.random() + 1) * 0.5;
        var sv2 = Math.random() * (sv1 - 0.1) + 0.1;
        var b1 = (Math.random() - 0.5) * 2 * bval;
        var b2 = (Math.random() - 0.5) * 2 * bval;

        var uv = [];
        for (let j=0; j < 2; j++) {
            x1 = Math.sign(Math.random() - 0.5);
            x2 = Math.sign(Math.random() - 0.5);

            var angle = Math.PI * (Math.random() * 2 - 1);
            var s = Math.sin(angle);
            var c = Math.cos(angle);

            var rmat = [[c * x1, s * x1], [-s * x2, c * x2]];
            uv.push(rmat);
        }
        let u = uv[0];
        let v = uv[1];
        
        var m = [
            [u[0][0]*sv1*v[0][0] + u[0][1]*sv2*v[1][0], u[0][0]*sv1*v[0][1] + u[0][1]*sv2*v[1][1], b1],
            [u[1][0]*sv1*v[0][0] + u[1][1]*sv2*v[1][0], u[1][0]*sv1*v[0][1] + u[1][1]*sv2*v[1][1], b2],
        ];

        matrices.push(m);
    }
    return matrices;
}

function iterateIFS(ws, n_iter, v0=null) {
    var dets = [];
    var ps = [];
    var total_sum = 0;
    for (let i = 0; i < ws.length; i++) {
        let det = ws[i][0][0] * ws[i][1][1] - ws[i][0][1] * ws[i][1][0];
        dets.push(det);
        ps.push(Math.abs(det));
        total_sum += ps[i];
    }
    ps[0] /= total_sum;
    for (let i = 1; i < ws.length; i++) {
        ps[i] = ps[i-1] + ps[i] / total_sum;
    }

    var x, y;
    if (v0 === null) {
        let s = 1 / (1 + dets[0] - ws[0][0][0] - ws[0][1][1]);
        x = s * ((1 - ws[0][1][1]) * ws[0][0][2] + ws[0][0][1] * ws[0][1][2]);
        y = s * ((1 - ws[0][0][0]) * ws[0][1][2] + ws[0][1][0] * ws[0][0][2]);
    } else {
        x = v0[0];
        y = v0[0];
    }

    var coords = [];

    var w, xt;
    for (let i = 0; i < n_iter; i++) {
        r = Math.random();
        var k;
        for (k = 0; k < ps.length; k++) {
            if (r < ps[k]) break;
        }
        w = ws[k];
        xt = x;
        x = w[0][0] * xt + w[0][1] * y + w[0][2];
        y = w[1][0] * xt + w[1][1] * y + w[1][2];
        coords.push([x, y]);
    }
    return coords;
}

function minmax(points) {
    var mms = [points[0][0], points[0][1], points[0][0], points[0][1]];
    for (let i = 1; i < points.length; i++) {
        let x = points[i][0];
        let y = points[i][1];
        if (x < mms[0]) mms[0] = x;
        if (y < mms[1]) mms[1] = y;
        if (x > mms[2]) mms[2] = x;
        if (y > mms[3]) mms[3] = y;
    }
    return mms;
}