const EPS = 1e-4;
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
    this.decrement = (a) => new Point(this.x - a.x, this.y - a.y, this.z - a.z);
    this.toString = () => `(${this.x}, ${this.y}, ${this.z})`;

    this.toPolarPoint = () => new PolarPoint(Math.sqrt(sqr(this.x) + sqr(this.y)), Math.atan2(this.y, this.x), z);
}
function PolarPoint(r, fi, h) {
    this.r = r;
    this.fi = fi;
    this.h = h;

    this.toString = () => `(${this.r}, ${this.fi}, ${this.h})`;
    this.toPoint = () => new Point(r * Math.cos(this.fi), r * Math.sin(this.fi), h);

    this.decrement = (p) => new PolarPoint(this.r - p.r, this.fi - p.fi, this.h - p.h);
    this.equals = (p) => equals(p.r, this.r) && equals(p.fi, this.fi) && equals(p.h, this.h);
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
    let tmp = b.decrement(a);
    return new Vector(tmp.x, tmp.y, tmp.z);
};

let sqr = a => a * a;

//init newton iteration
//return function which return subtrahend for newton iteration x(k+1) by x(k)
function initNewtonIteration(p1, p2, l1, l2, ri, g, r0) {
    return function ({ r, fi, h }) {
        let matrix = [], fx = [];
        let cos1 = Math.cos(fi - p1.fi), sin1 = Math.sin(fi - p1.fi);
        let cos2 = Math.cos(fi - p2.fi), sin2 = Math.sin(fi - p2.fi);

        //set yakobi matrix
        matrix.push([-1, ri, g]);
        matrix.push([2 * r - 2 * p1.r * cos1, 2 * r * p1.r * sin1, -2 * (-h + p1.h)]);
        matrix.push([2 * r - 2 * p2.r * cos2, 2 * r * p2.r * sin2, -2 * (-h + p2.h)]);

        //calc value of function
        fx[0] = -r + r0 + fi * ri + h * g;
        fx[1] = -sqr(l1) + sqr(p1.h - h) + sqr(r) + sqr(p1.r) - 2 * r * p1.r * Math.cos(p1.fi - fi);
        fx[2]  = -sqr(l2) + sqr(p2.h - h) + sqr(r) + sqr(p2.r) - 2 * r * p2.r * Math.cos(p2.fi - fi);
        return inverseMatrixAndMul(matrix, fx);
    }
}

//inverse matrix 3x3 and multiply on vector 3x1
function inverseMatrixAndMul(m, fx) {
    let matrix = [], factor = 0, ans = {};
    factor += m[0][2] * m[1][1] * m[2][0] - m[0][1] * m[1][2] * m[2][0] - m[0][2] * m[1][0] * m[2][1];
    factor += m[0][0] * m[1][2] * m[2][1] + m[0][1] * m[1][0] * m[2][2] - m[0][0] * m[1][1] * m[2][2];
    matrix.push([m[1][2] * m[2][1] - m[1][1] * m[2][2], m[0][1] * m[2][2] - m[0][2] * m[2][1], m[0][2] * m[1][1] - m[0][1] * m[1][2]]);
    matrix.push([m[1][0] * m[2][2] - m[1][2] * m[2][0], m[0][2] * m[2][0] - m[0][0] * m[2][2], m[0][0] * m[1][2] - m[0][2] * m[1][0]]);
    matrix.push([m[1][1] * m[2][0] - m[1][0] * m[2][1], m[0][0] * m[2][1] - m[0][1] * m[2][0], m[0][1] * m[1][0] - m[0][0] * m[1][1]]);
    ans.r = (matrix[0][0] * fx[0] + matrix[0][1] * fx[1] + matrix[0][2] * fx[2]) / factor;
    ans.fi = (matrix[1][0] * fx[0] + matrix[1][1] * fx[1] + matrix[1][2] * fx[2]) / factor;
    ans.h = (matrix[2][0] * fx[0] + matrix[2][1] * fx[1] + matrix[2][2] * fx[2]) / factor;
    return ans;
}

let t = 0, operation = 0;

function check(p1, p2) {
    return Math.abs(p1.r - p2.r) < EPS && Math.abs(p2.fi - p1.fi) < EPS && Math.abs(p1.h - p2.h) < EPS;
}

//solve system of equations
function solve(p1, p2, l1, l2, ri, g, r0, p0) {
    let last = p0;
    let iteration = initNewtonIteration(p1, p2, l1, l2, ri, g, r0);
    let p = last.decrement(iteration(last));
    t++;
    while (!p.equals(last)) {
        last = p;
        p = p.decrement(iteration(p));
        t++;
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
    b = a;
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
    console.log('newton iteration on one point: ', t / operation, t, operation);
    points.map((elem, ind) => points[ind] = elem.map(point => point.toPoint()));
    return points;
}
