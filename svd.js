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