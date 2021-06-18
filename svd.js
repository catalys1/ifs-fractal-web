// function svd2x2(mat) {
//     var a = mat[0][0];
//     var b = mat[0][1];
//     var c = mat[1][0];
//     var d = mat[1][1];

//     var aa = a*a;
//     var bb = b*b;
//     var cc = c*c;
//     var dd = d*d;

//     // U
//     var theta = 0.5 * Math.atan2(2*a*c + 2*b*d, aa + bb - cc - dd);
//     var ctheta = Math.cos(theta);
//     var stheta = Math.sin(theta);
//     var U = [[ctheta, -stheta], [stheta, ctheta]];

//     // SIGMA
//     var s1 = aa + bb + cc + dd;
//     var s2 = Math.sqrt(Math.pow(aa + bb - cc - dd, 2) + 4*Math.pow(a*c + b*d, 2));
//     var sig1 = Math.sqrt((s1 + s2) / 2);
//     var sig2 = Math.sqrt((s1 - s2) / 2);
//     var S = [sig1, sig2];

//     // V
//     var phi = 0.5 * Math.atan2(2*a*b + 2*c*d, aa - bb + cc - dd);
//     var cphi = Math.cos(phi);
//     var sphi = Math.sin(phi);
//     var s11 = Math.sign((a*ctheta + c*stheta)*cphi + (b*ctheta + d*stheta)*sphi);
//     var s22 = Math.sign((a*stheta - c*ctheta)*sphi + (-b*stheta + d*ctheta)*cphi);
//     var V = [[s11*cphi, -s22*sphi], [s11*sphi, s22*cphi]];

//     return {'U': U, 'S': S, 'V': V};
// }

// function svd2x2(mat) {
//     var y1 = mat[1][0] + mat[0][1];
//     var x1 = mat[0][0] - mat[1][1];
//     var y2 = mat[1][0] - mat[0][1];
//     var x2 = mat[0][0] + mat[1][1];

//     var h1 = Math.hypot(y1, x1);
//     var h2 = Math.hypot(y2, x2);

//     var t1 = x1 / h1;
//     var t2 = x2 / h2;

//     var cc = Math.sqrt((1 + t1) * (1 + t2));
//     var ss = Math.sqrt((1 - t1) * (1 - t2));
//     var cs = Math.sqrt((1 + t1) * (1 - t2));
//     var sc = Math.sqrt((1 - t1) * (1 + t2));

//     var c1 = (cc - ss) / 2;
//     var s1 = (sc + cs) / 2;
//     var U = [[c1, -s1], [s1, c1]];

//     var S = [(h1 + h2) / 2, (h1 - h2) / 2];

//     if (h1 != h2) {
        
//     }
// }

function svd2x2(mat) {
    // s[0] = (sqrt(pow(a[0] - a[3], 2) + pow(a[1] + a[2], 2)) + sqrt(pow(a[0] + a[3], 2) + pow(a[1] - a[2], 2))) / 2;
    // s[1] = fabs(s[0] - sqrt(pow(a[0] - a[3], 2) + pow(a[1] + a[2], 2)));
    // v[2] = (s[0] > s[1]) ? sin((atan2(2 * (a[0] * a[1] + a[2] * a[3]), a[0] * a[0] - a[1] * a[1] + a[2] * a[2] - a[3] * a[3])) / 2) : 0;
    // v[0] = sqrt(1 - v[2] * v[2]);
    // v[1] = -v[2];
    // v[3] = v[0];
    // u[0] = (s[0] != 0) ? (a[0] * v[0] + a[1] * v[2]) / s[0] : 1;
    // u[2] = (s[0] != 0) ? (a[2] * v[0] + a[3] * v[2]) / s[0] : 0;
    // u[1] = (s[1] != 0) ? (a[0] * v[1] + a[1] * v[3]) / s[1] : -u[2];
    // u[3] = (s[1] != 0) ? (a[2] * v[1] + a[3] * v[3]) / s[1] : u[0];
    var a = mat[0][0];
    var b = mat[0][1];
    var c = mat[1][0];
    var d = mat[1][1];

    var aa = a*a;
    var bb = b*b;
    var cc = c*c;
    var dd = d*d;

    var sig1 = (Math.hypot(a - d, b + c) + Math.hypot(a + d, b - c)) / 2;
    var sig2 = Math.abs(sig1 - Math.hypot(a - d, b + c));
    var S = [sig1, sig2];

    var v10 = sig1 > sig2 ? Math.sin((Math.atan2(2 * (a * b + c * d), a*a - b*b + c*c - d*d)) / 2) : 0;
    var v00 = Math.sqrt(1 - v10*v10);
    var V = [[v00, -v10], [v10, v00]];

    var u00 = (sig1 != 0) ? (a * v00 + b * v10) / sig1 : 1;
    var u10 = (sig1 != 0) ? (c * v00 + d * v10) / sig1 : 0;
    var u01 = (sig2 != 0) ? (a * -v10 + b * v00) / sig2 : -u10;
    var u11 = (sig2 != 0) ? (c * -v10 + d * v00) / sig2 : u00;
    var U = [[u00, u01], [u10, u11]];

    return {'U': U, 'S': S, 'V': V};
}

function rq2x2(mat) {
    var a = mat[0][0];
    var b = mat[0][1];
    var c = mat[1][0];
    var d = mat[1][1];

    if (c == 0) {
        return [a, b, d, 1, 0];
    }
    var maxden = Math.max(Math.abs(c), Math.abs(d));
    c = c / maxden;
    d = d / maxden;
    var den = 1 / Math.hypot(c, d);
    var numx = (-b*c + a*d);
    var numy = (a*c + b*d);
    return [numx * den, numy * den, maxden / den, d * den, -c * den];
}
// template <class T>
// void Rq2x2Helper(const Matrix<T, 2, 2>& A, T& x, T& y, T& z, T& c2, T& s2) {
//     T a = A(0, 0);
//     T b = A(0, 1);
//     T c = A(1, 0);
//     T d = A(1, 1);

//     if (c == 0) {
//         x = a;
//         y = b;
//         z = d;
//         c2 = 1;
//         s2 = 0;
//         return;
//     }
//     T maxden = std::max(abs(c), abs(d));

//     T rcmaxden = 1/maxden;
//     c *= rcmaxden;
//     d *= rcmaxden;

//     T den = 1/sqrt(c*c + d*d);

//     T numx = (-b*c + a*d);
//     T numy = (a*c + b*d);
//     x = numx * den;
//     y = numy * den;
//     z = maxden/den;

//     c2 = d * den;
//     s2 = -c * den;
// }

function svd2x2(mat) {
    var rq = rq2x2(mat);
    var x = rq[0];
    var y = rq[1];
    var z = rq[2];
    var c2 = rq[3];
    var s2 = rq[4];

    var scaler = 1 / Math.max(Math.abs(x), Math.abs(y));
    var x_ = x * scaler;
    var y_ = y * scaler;
    var z_ = z * scaler;
    var numer = ((z_ - x_) * (z_ + x_)) + y_*y_;
    var gamma = numer == 0 ? 1 : x_ * y_;
    var zeta = numer / gamma;
    
    var t = 2 * (zeta >= 0 ? 1 : -1) / (Math.abs(zeta) + Math.sqrt(zeta*zeta + 4));
    var c1 = 1 / Math.sqrt(1 + t*t);
    var s1 = c1 * t;

    var usa = c1*x - s1*y; 
    var usb = s1*x + c1*y;
    var usc = -s1*z;
    var usd = c1*z;

    t = c1*c2 + s1*s2;
    s2 = c2*s1 - c1*s2;
    c2 = t;

    var d1 = Math.hypot(usa, usc);
    var d2 = Math.hypot(usb, usd);
    var dmax = Math.max(d1, d2);
    var usmax1 = d2 > d1 ? usd : usa;
    var usmax2 = d2 > d1 ? usb : -usc;

    var signd1 = (x*z >= 0 ? 1 : -1);
    dmax *= d2 > d1 ? signd1 : 1;
    d2 *= signd1;
    var rcpdmax = 1 / dmax;

    c1 = dmax != 0 ? usmax1 * rcpdmax : 1;
    s1 = dmax != 0 ? usmax2 * rcpdmax : 0;
}
// template <class T>
// void Svd2x2Helper(const Matrix<T, 2, 2>& A, T& c1, T& s1, T& c2, T& s2, T& d1, T& d2) {
//     // Calculate RQ decomposition of A
//     T x, y, z;
//     Rq2x2Helper(A, x, y, z, c2, s2);

//     // Calculate tangent of rotation on R[x,y;0,z] to diagonalize R^T*R
//     T scaler = T(1)/std::max(abs(x), abs(y));
//     T x_ = x*scaler, y_ = y*scaler, z_ = z*scaler;
//     T numer = ((z_-x_)*(z_+x_)) + y_*y_;
//     T gamma = x_*y_;
//     gamma = numer == 0 ? 1 : gamma;
//     T zeta = numer/gamma;

//     T t = 2*impl::sign_nonzero(zeta)/(abs(zeta) + sqrt(zeta*zeta+4));

//     // Calculate sines and cosines
//     c1 = T(1) / sqrt(T(1) + t*t);
//     s1 = c1*t;

//     // Calculate U*S = R*R(c1,s1)
//     T usa = c1*x - s1*y; 
//     T usb = s1*x + c1*y;
//     T usc = -s1*z;
//     T usd = c1*z;

//     // Update V = R(c1,s1)^T*Q
//     t = c1*c2 + s1*s2;
//     s2 = c2*s1 - c1*s2;
//     c2 = t;

//     // Separate U and S
//     d1 = std::hypot(usa, usc);
//     d2 = std::hypot(usb, usd);
//     T dmax = std::max(d1, d2);
//     T usmax1 = d2 > d1 ? usd : usa;
//     T usmax2 = d2 > d1 ? usb : -usc;

//     T signd1 = impl::sign_nonzero(x*z);
//     dmax *= d2 > d1 ? signd1 : 1;
//     d2 *= signd1;
//     T rcpdmax = 1/dmax;

//     c1 = dmax != T(0) ? usmax1 * rcpdmax : T(1);
//     s1 = dmax != T(0) ? usmax2 * rcpdmax : T(0);
// }