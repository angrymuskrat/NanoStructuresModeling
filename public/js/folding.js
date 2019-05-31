const EPS = 1e-7;
const defaultRoundAmount = 7;

function equals(a, b) {
    return Math.abs(a - b) < EPS;
}

function sign(a) {
    return isZero(a) ? 0 : (a < 0 ? -1 : 1);
}


function isZero(a) {
    return equals(a, 0);
}

function round(a, amount = defaultRoundAmount) {
    let factor = Math.pow(10, amount);
    return Math.round(a * factor) / factor;
}

function Point(x, y, z) {
    this.x = x;
    this.y = y;
    this.z = z;

    this.add = (a) => new Point(this.x + a.x, this.y + a.y, this.z + a.z);
    this.increment = (a) => new Point(this.x - a.x, this.y - a.y, this.z - a.z);
    this.toString = () => `(${this.x}, ${this.y}, ${this.z})`;

    this.toPolarPoint = () => new PolarPoint(Math.sqrt(sqr(this.x) + sqr(this.y)), Math.atan2(this.y, this.x), z);
}
function PolarPoint(r, fi, h) {
    this.r = r;
    this.fi = fi;
    this.h = h;

    this.toString = () => `(${this.r}, ${this.fi}, ${this.h})`;
    this.toPoint = () => new Point(r * Math.cos(this.fi), r * Math.sin(this.fi), h);

    this.increment = (p) => new PolarPoint(this.r - p.r, this.fi - p.fi, this.h - p.h);
}

//{x, y, z}
Point.round = (pointLike) => ['x', 'y', 'z'].map(name => pointLike[name] = round(pointLike[name]));


function Vector(x, y, z) {
    Point.call(this, x, y, z);

    this.length = () => Math.sqrt(this.x * this.x + this.y * this.y + this.z * this.z);
    this.multiplication = (val) => new Vector(this.x * val, this.y * val, this.z * val);
    this.normalize = () => this.multiplication(1 / this.length())
    this.dot = (v) => (this.x * v.x + this.y * v.y + this.z * v.z) / (this.length() * v.length());
    this.cross = (v) => new Vector(
        this.y * v.z - this.z * v.y, -this.x * v.z + this.z * v.x, this.x * v.y - this.y * v.x)
        / (this.length() * v.length());

    this.remove = (x, y, z) => {
        this.x = x;
        this.y = y;
        this.z = z;
    }
}
Vector.fromPointLike = ({x, y, z}) => new Vector(x, y, z);

Vector.makeFromPoints =  (a, b) => {
    let tmp = b.increment(a);
    return new Vector(tmp.x, tmp.y, tmp.z);
};

let newPointError = (mes) => new Error('Невозможно вычислить точку. ' + mes || '');
let sqr = a => a * a;

//init newton iteration
//return function which return subtrahend for newton iteration x(k+1) by x(k)
function initNewtonIteration(p1, p2, l1, l2, ri, g, r0) {
    return function ({ r, fi, h }) {
        let matrix = [], fx = {};
        let cos1 = Math.cos(fi - p1.fi), sin1 = Math.sin(fi - p1.fi);
        let cos2 = Math.cos(fi - p2.fi), sin2 = Math.sin(fi - p2.fi);

        //set yakobi matrix
        matrix.push([-1, ri, g]);
        matrix.push([2 * r - 2 * p1.r * cos1, 2 * r * p1.r * sin1, -2 * (-h + p1.h)]);
        matrix.push([2 * r - 2 * p2.r * cos2, 2 * r * p2.r * sin2, -2 * (-h + p2.h)]);

        //calc value of function
        fx.r = -r + r0 + fi * ri + h * g;
        fx.fi = -sqr(l1) + sqr(p1.h - h) + sqr(r) + sqr(p1.r) - 2 * r * p1.r * Math.cos(p1.fi - fi);
        fx.h  = -sqr(l2) + sqr(p2.h - h) + sqr(r) + sqr(p2.r) - 2 * r * p2.r * Math.cos(p2.fi - fi);
        return inverseMatrixAndMul(matrix, fx);
    }
}

//inverse matrix 3x3 and multiply on vector 3x1
function inverseMatrixAndMul(m, {r, fi, h}) {
    let matrix = [], factor = 0, ans = {};
    factor += m[0][2] * m[1][1] * m[2][0] - m[0][1] * m[1][2] * m[2][0] - m[0][2] * m[1][0] * m[2][1];
    factor += m[0][0] * m[1][2] * m[2][1] + m[0][1] * m[1][0] * m[2][2] - m[0][0] * m[1][1] * m[2][2];
    matrix.push([m[1][2] * m[2][1] - m[1][1] * m[2][2], m[0][1] * m[2][2] - m[0][2] * m[2][1], m[0][2] * m[1][1] - m[0][1] * m[1][2]]);
    matrix.push([m[1][0] * m[2][2] - m[1][2] * m[2][0], m[0][2] * m[2][0] - m[0][0] * m[2][2], m[0][0] * m[1][2] - m[0][2] * m[1][0]]);
    matrix.push([m[1][1] * m[2][0] - m[1][0] * m[2][1], m[0][0] * m[2][1] - m[0][1] * m[2][0], m[0][1] * m[1][0] - m[0][0] * m[1][1]]);
    ans.r = (matrix[0][0] * r + matrix[0][1] * fi + matrix[0][2] * h) / factor;
    ans.fi = (matrix[1][0] * r + matrix[1][1] * fi + matrix[1][2] * h) / factor;
    ans.h = (matrix[2][0] * r + matrix[2][1] * fi + matrix[2][2] * h) / factor;
    return ans;
}

let t = 0, operation = 0;

function check(p1, p2) {
    return Math.abs(p1.r - p2.r) < EPS && Math.abs(p2.fi - p1.fi) < EPS && Math.abs(p1.h - p2.h) < EPS;
}

//solve system of equations
function solve(p1, p2, l1, l2, ri, g, r0, p0) {
    if (operation === 69)
        console.log(1);
    let last = p0 || (new PolarPoint((p1.r + p2.r) + 2, (p1.fi + p2.fi) / 2 + 0.12, (p1.fi + p2.fi) / 2));
    let iteration = initNewtonIteration(p1, p2, l1, l2, ri, g, r0);
    let p = last.increment(iteration(last)), ind = 0;
    t++;
    while (!check(p, last)) {
        last = p;
        p = p.increment(iteration(p));
        t++;
        ind++;
        if (ind > 50) {
            console.log(p1.toString(), p2.toString());
            console.log(p.toString());
            console.log('problem ' + operation);
            return p;
        }
    }
    operation++;
    return p;
}

function findAndSumVectorPolar(p0, p1, p2) {
    let v = Vector.makeFromPoints(p0.toPoint(), p2.toPoint());
    //return p1.toPoint().add(v).toPolarPoint();
    return new PolarPoint(p2.r - p0.r + p1.r, p2.fi - p0.fi + p1.fi, p2.h - p0.h + p1.h);
}

function folding(r0, ri, n, m, alpha, beta, a) {
    let b = Math.sqrt(2) * a;
    let h = a * Math.cos(alpha);
    let r = r0 + a * Math.sin(alpha);
    let points = [[], []];
    points[0].push(new PolarPoint(r0, 0, 0));
    points[1].push(new PolarPoint(r, 0, h));
    let p = points[1][0].toPoint();
    p.y -= a;
    points[1][-1] = p.toPolarPoint();

    let riPi = ri / (2 * Math.PI);
    let g = Math.tan(alpha);

    for (let i = 1; i <= n; i++) {
        points[0].push(solve(points[0][i - 1], points[1][i - 1], a, b, riPi, g, r0,
            findAndSumVectorPolar(points[1][i - 2], points[0][i - 1], points[1][i - 1])));
        points[1].push(solve(points[1][i - 1], points[0][i], a, a, riPi, g, r0,
            findAndSumVectorPolar(points[0][i - 1], points[1][i - 1], points[0][i])));
    }

    for (let j = 2; j <= m; j++) {
        points.push([]);
        h = j * a * Math.cos(alpha);
        r = r0 + j * a * Math.sin(alpha);
        points[j].push(new PolarPoint(r, 0, h));

        p = points[j][0].toPoint();
        p.y -= a;
        points[1][-1] = p.toPolarPoint();
        for (let i = 1; i <= n; i++) {
            points[j].push(solve(points[j][i - 1], points[j - 1][i], a, a, riPi, g, r0,
                findAndSumVectorPolar(points[j - 1][i - 1], points[j][i - 1], points[j - 1][i])));
        }
    }
    console.log('newton iteration on one point: ', t / operation);
    points.map((elem, ind) => points[ind] = elem.map(point => point.toPoint()));
    return points;
}

/*
//let r0 = 10, ri = 3, alpha = Math.PI / 30, n = 36, m = 0, beta = 0;
let r0 = 10, ri = 1, alpha = Math.PI / 300, n = 200, m = 3, beta = 0;
folding(r0, ri, n, m, alpha, beta)//.map(point => console.log(point.toString()));
console.log(t / 2 / n);*/
/*
let p0 = new PolarPoint(12, 0.4, 3);
let p1 = new PolarPoint(5, 0.1, 8);
let p2 = new PolarPoint(1, 2.4, -6);
let p4 = new PolarPoint(p2.r - p0.r + p1.r, p2.fi - p0.fi + p1.fi, p2.h - p0.h + p1.h);

console.log(findAndSumVectorPolar(p0, p1, p2).toString());
console.log(p4.toString());*/




/*function folding(p1, p2, p3, fi1, fi2, amount) {
    let makeMyPoint = ({ x, y, z }) => new Point(x, y, z);
    let points = [makeMyPoint(p1), makeMyPoint(p2), makeMyPoint(p3)];

    for (let i = 3; i <= amount; i++) {
        let fi = i % 2 ? fi1 : fi2;

        points.push(foldingOnePoint(points[i - 3], points[i - 2], points[i - 1], fi, i % 2 === 1));
        fi1 += fi1 / (amount * 8);
        fi2 += fi2 / (amount * 9);
    }
    points.map(point => Point.round(point, 3));
    return points;
}*/
/*function foldingOnePoint(p1, p2, p3, fi, isDiagonal) {
    // solve this system
    //1) cos(fi) = (v31, v24) / (|v31| * |v24|) = x * v31.x + y * v31.y + z * v31.z
    //2) cos(v24^v23) = t = (v23, v24) / (|v23| * |v24|) = x * v23.x + y * v23.y + z * v23.z
    //3) 1 = x*x + y*y + z*z

    let v23 = Vector.makeFromPoints(p2, p3);
    let v31 = Vector.makeFromPoints(p3, p1);
    let t = isDiagonal ? v31.length() / v23.length() : 0;
    let length = v31.length();
    let cosFi = Math.cos(fi), sinFi = Math.sin(fi);
    v31 = v31.normalize();
    v23 = v23.normalize();

    let ans = [];
    let addAns = ({x, y, z}) => ans.push({x, y, z});

    // (4) = v23.z * z = t - v23.x * x - v23.y * y - (2) in other form
    // (1) * v23.z and insert  right part of (4) instead v23.z * z  =>  k = a * x + b * y
    let k = v23.z * cosFi - v31.z * t;
    let a = v23.z * v31.x - v31.z * v23.x;
    let b = v23.z * v31.y - v31.z * v23.y;

    if (isZero(a) || isZero(b)) {
        //cos(fi) - v31.y * y = v31.x * x + v31.z * z (1)'
        //t - v23.y * y = v23.x * x + v23.z * z (2)'
        //1 - y^2 = x^2 + z^2 (3)'
        //y = k / b
        //p.s a = 0 => b != 0 and one of parameter before x, z isn't zero
        if (isZero(a) && isZero(b))
            throw newPointError('Нулевые a и b');
        let [found, name] = isZero(a) ? [k/b, 'y'] : [k/a, 'x'];
        k = cosFi - v31[name] * found;
        let l = t - v23[name] * found;
        let [first, second] = name === 'y' ? ['x', 'z'] : ['y', 'z'];

        let findRoots = (v, root, free) => {
            let decision = {};
            decision[name] = found;
            let other = root === first ? second : first;
            let disc = sqr(v[root] * free) - (sqr(v[root]) + sqr(v[other]))
                * ((sqr(found) - 1) * sqr(v[other]) + sqr(free));

            if (!isZero(disc) && disc < 0)
                throw newPointError('Отрицательный дискриминант');

            let r1 = (Math.sqrt(disc) + v[root]) / (sqr(v[root]) + sqr(v[other]));
            let r2 = (-Math.sqrt(disc) + v[root]) / (sqr(v[root]) + sqr(v[other]));
            decision[root] = r1;
            decision[other] = (free - v[root] * decision[root]) / v[other];
            addAns(decision);
            decision[root] = r2;
            decision[other] = (free - v[root] * decision[root]) / v[other];
            addAns(decision);
        };
        if (!isZero(v31[first]) || !isZero(v23[first])) {
            if (!isZero(v31[first]))
                findRoots(v31, second, k);
            else
                findRoots(v23, second, l);
        }
        else {
            if (!isZero(v31[second]))
                findRoots(v31, first, k);
            else {
                findRoots(v23, first, l);
            }
        }
    }
    else {
        //right part of (4) to (3) * v23.z^2
        //l = c * x^2 + d * y^2 - 2 * v23.x * x - 2 * v23.y * y + 2 * v23.x * v23.y * x * y (5)
        let l = v23.z * v23.z - t * t;
        let c = v23.x * v23.x + v23.z * v23.z;
        let d = v23.y * v23.y + v23.z * v23.z;

        //a*x = k - b * y (6)
        //right part of (6) to (5) * a^2
        //f*y^2 + e*y + j = 0
        let f = c * b * b + a * a * d - 2 * v23.x * v23.y * a * b;
        let e = (v23.x * v23.y * a * k + v23.x * a * b * t - v23.y * a * a * t - c * k * b);
        let j = (c * k * k - a * a * l - 2 * v23.x * a * k * t);

        let findRoots = (y) => {
            let decision = { y };
            decision.x = (k - b * y) / a;
            decision.z = Math.sqrt(1 - sqr(decision.x) - sqr(y)) || 0;
            addAns(decision);
            decision.z = -decision.z;
            addAns(decision);
        };

        if (isZero(f)) {
            if (isZero(e)) {
                throw newPointError();
            }
            findRoots(j / e);
        }
        else {
            let disc = e * e - f * j;
            if (!isZero(d) && d < 0)
                throw newPointError('Отрицательный дискреминант');
            disc = Math.sqrt(disc) || 0;
            let r1 = (disc - e) / f;
            let r2 = (-disc - e) / f;
            findRoots(r1);
            findRoots(r2);
        }
    }
    //check one root in init system
    let check = ({ x, y, z }, p) => {
        let main = x * v31.x + y * v31.y + z * v31.z;
        let second = x * v23.x + y * v23.y + z * v23.z;
        let third = sqr(x) + sqr(y) + sqr(z);
        let angleSign = planeChecker(v23, v31, p1, { x: x + p.x, y: y + p.y, z: z + p.z }) * (isDiagonal ? -1 : 1);
        return equalsSign(angleSign, sinFi) && equals(third, 1) && equals(cosFi, main) && equals(t, second);
    };

    //console.log('next');
    for (let i in ans) {
        if (check(ans[i], p2)) {
            //Point.round(ans[i]);
            return p2.add(Vector.fromPointLike(ans[i]).multiplication(length));
        }
    }
    throw newPointError('Уравнение не имеет решений');
}*/