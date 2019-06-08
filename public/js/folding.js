const EPS = 1e-7;

function Point(x, y, z, isEmpty = false, isBimbo = false) {
    this.x = x;
    this.y = y;
    this.z = z;
    this.isBimbo = isBimbo;
    this.isEmpty = isEmpty;

    this.decrement = (a) => new Point(this.x - a.x, this.y - a.y, this.z - a.z);

    this.lengthTo = (p) => Math.sqrt(sqr(this.x - p.x) + sqr(this.y - p.y) + sqr(this.z - p.z));
    this.sqrLengthTo = (p) => sqr(this.x - p.x) + sqr(this.y - p.y) + sqr(this.z - p.z);

    this.toString = () => `(${this.x}, ${this.y}, ${this.z})`;

    this.toPolarPoint = () => new PolarPoint(Math.sqrt(sqr(this.x) + sqr(this.y)), Math.atan2(this.y, this.x), z, this.isEmpty, this.isBimbo);
}

function PolarPoint(r, fi, h, isEmpty = false, isBimbo = false) {
    this.r = r;
    this.fi = fi;
    this.h = h;
    this.isBimbo = isBimbo;
    this.isEmpty = isEmpty;

    this.toString = () => `(${this.r}, ${this.fi}, ${this.h})`;
    this.toPoint = () => new Point(r * Math.cos(this.fi), r * Math.sin(this.fi), h, this.isEmpty, this.isBimbo);

    this.decrement = (p) => new PolarPoint(this.r - p.r, this.fi - p.fi, this.h - p.h);
    this.equals = (p) => equals(p.r, this.r) && equals(p.fi, this.fi) && equals(p.h, this.h);
}

function CircleList() {
    function Elem(p0, p1, p2, prev, next) {
        this.prev = prev;
        this.next = next;
        this.p0 = p0;
        this.p1 = p1;
        this.p2 = p2;
    }
    let root = new Elem();

    let pointer = root;
    root.prev = root.next = root;

    this.push = (p0, p1, p2) => {
        let elem = new Elem(p0, p1, p2, root.prev, root);
        root.prev.next = elem;
        root.prev = elem;
    };

    this.shift = (p0, p1, p2) => {
        let elem = new Elem(p0, p1, p2, root, root.next);
        root.next.prev = elem;
        root.next = elem;
    };

    this.headValue = () => {
        if (pointer === root)
            return null;
        return pointer;
    };

    this.next = () => {
        pointer = pointer.next;
        return this.headValue();
    };

    this.deleteHead = () => {
        let val = this.headValue();
        if (val) {
            val.prev.next = val.next;
            val.next.prev = val.prev;
        }
        return val;
    }
}

function equals(a, b) {
    return Math.abs(a - b) < EPS;
}

function isZero(a) {
    return equals(a, 0);
}

function sqr(a) {
    return a * a;
}

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


function findAndSumVectorPolar(p0, p1, p2) {
    return new PolarPoint(p2.r - p0.r + p1.r, p2.fi - p0.fi + p1.fi, p2.h - p0.h + p1.h)
}

//solve system of equations
function solveEquations(p0, p1, p2, l1, l2, ri, g, r0) {
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

function folding(r0, ri, n, m, alpha, beta, a) {
    let aHelper = Math.sqrt(3) / 2 * a;
    let bHelper = a / 2;
    let cHelper = 2 * sqr(aHelper) * (1 - Math.cos(beta));
    let h = aHelper * Math.cos(alpha);
    let r = r0 + aHelper * Math.sin(alpha);
    let riPi = ri / (2 * Math.PI);
    let g = Math.tan(alpha);
    let points = [[], []];
    points[0].push(new PolarPoint(r0, 0, 0));

    let helpPoint = new PolarPoint(r, 0, h);

    if (!isZero(beta))
        helpPoint = solveEquations(new PolarPoint(helpPoint.r, beta > 0 ? 0.1 : -0.1, helpPoint.h),
            helpPoint, points[0][0], cHelper, aHelper, riPi, g, r0);

    let helpPointLeft = new PolarPoint(helpPoint.r, helpPoint.fi + 0.1, helpPoint.h);
    let helpPointRight = new PolarPoint(helpPoint.r, helpPoint.fi - 0.1, helpPoint.h);

    let solve = (p0, p1, p2) =>
        solveEquations(findAndSumVectorPolar(p0, p1, p2), p1, p2, a, a, riPi, g, r0);

    points[1].push(solveEquations(helpPointRight, points[0][0], helpPoint, a, bHelper, riPi, g, r0));
    points[1].push(solveEquations(helpPointLeft, points[0][0], helpPoint, a, bHelper, riPi, g, r0));
    points[0].push(solve(points[1][0], points[0][0], points[1][1], a));

    let p = points[1][0].toPoint();
    p.y -= a;


    for (let i = 2; i <= n; i++) {
        points[1].push(solve(points[0][i - 2], points[0][i - 1], points[1][i - 1]));
        points[0].push(solve(points[1][i - 1], points[1][i], points[0][i - 1]));
    }

    for (let j = 2; j <= m; j++) {
        points.push([[], []]);
        let [ti, ni] = j % 2 === 0 ? [0, 1] : [1, 0];

        points[j][ti] = solve(points[j - 2][ti], points[j - 1][ni], points[j - 1][ti]);
        points[j][ni] = solve(points[j - 1][ti], points[j - 1][ni], points[j][ti]);

        for (let i = 2; i <= n; i++) {
            if (i % 2 || i === n)
                points[j].push(solve(points[j - 1][i - 1], points[j - 1][i], points[j][i - 1]));
            else {
                if (j % 2)
                    points[j].push(solve(points[j - 2][i], points[j - 1][i], points[j - 1][i - 1]));
                else
                    points[j].push(solve(points[j - 2][i], points[j - 1][i], points[j - 1][i + 1]));
            }
        }
    }
    console.log('newton iteration on one point: ', t / operation, t, operation);
    return points;
}

function getCubsCoordinate(point, size) {
    let x = Math.floor(point.x / size);
    let y = Math.floor(point.y / size);
    let z = Math.floor(point.z / size);
    return [x, y, z];
}

function addPointToCub(cubs, size, point, i, j) {
    let p = point.toPoint();
    let [x, y, z] = getCubsCoordinate(p, size);

    cubs[x] = cubs[x] === undefined ? [] : cubs[x];
    cubs[x][y] = cubs[x][y] === undefined ? [] : cubs[x][y];
    cubs[x][y][z] = cubs[x][y][z] === undefined ? [] : cubs[x][y][z];
    p.i = i;
    p.j = j;

    cubs[x][y][z].push(p);
}

function initCubs(points, size) {
    let cubs = [];
    points.map((array, i) => array.map((point, j) => addPointToCub(cubs, size, point, i, j)));
    return cubs;
}

function isNotEmpty(i, j) {
    return ((j + (i % 2 + 2) % 2 + 1) % 3 + 3) % 3 > 0;
}

function getNeighbours(cubs, size, args) {
    let [x, y, z] = getCubsCoordinate(args[0], size);
    let xmin = x - 1, xmax = x + 1, ymin = y - 1, ymax = y + 1, zmin = z - 1, zmax = z + 1;
    for (let i in args) {
        [x, y, z] = getCubsCoordinate(args[i], size);
        xmin = Math.min(xmin, x - 1);
        xmax = Math.max(xmax, x + 1);
        ymin = Math.min(ymin, y - 1);
        ymax = Math.max(ymax, y + 1);
        zmin = Math.min(zmin, z - 1);
        zmax = Math.max(zmax, z + 1);
    }
    let ans = [];
    for (let i = xmin; i <= xmax; i++) {
        for (let j = ymin; j <= ymax; j++) {
            for (let l = zmin; l <= zmax; l++) {
                if (!cubs[i]) continue;
                if (!cubs[i][j]) continue;
                if (!cubs[i][j][l]) continue;
                ans.push(...cubs[i][j][l]);
            }
        }
    }
    return ans;
}

function twoPointPotential(a, b, s1, s2, size, force = false) {
    if (!force && (a.isEmpty || b.isEmpty || a.isBimbo || b.isBimbo)) return 0;

    let sqrLen = a.sqrLengthTo(b);

    if (isZero(sqrLen) || sqrLen > size * size) return 0;

    let len6 = Math.pow(sqrLen, 3);
    return s1 / sqr(len6) - s2 / len6;
}

function getPLJ(cubs, size, args, s1, s2) {
    let points = args.map(point => point.toPoint());
    let near = getNeighbours(cubs, size, points);
    let potential = 0;
    for (let i in near) {
        for (let j in points) {
            potential += twoPointPotential(near[i], points[j], s1, s2, size);
        }
    }
    for (let i in points) {
        for (let j in points) {
            if ((i > j) && (points[i].isBimbo || points[j].isBimbo)) {
                potential += twoPointPotential(points[i], points[j], s1, s2, size, true);
            }
            if ((i < j) && (!points[i].isBimbo && !points[j].isBimbo)) {
                potential -= twoPointPotential(points[i], points[j], s1, s2, size, true);
            }
        }
    }
    //console.log(potential)
    return potential;
}

function makeRealPoint(cubs, size, elem, point, isEmpty) {
    point.isBimbo = false;
    point.isEmpty = isEmpty;
    if (!isEmpty) {
        addPointToCub(cubs, size, point, elem.i, elem.j);
    }
}

function buildHex(cubs, size, s1, s2, points, solve, energyMaximum = 0, ind0, ind1, ind2) {
    let getPoint = elem => points[elem.i][elem.j];
    let makePoint = (ind0, ind1, ind2) => {
        let i = ind1.i + ind2.i - ind0.i;
        let j = ind1.j + ind2.j - ind0.j;
        if (ind1.i === ind2.i)
            j = ind0.j;
        if (!points[i])
            points[i] = [];
        if (!points[i][j])
            points[i][j] = solve(getPoint(ind0), getPoint(ind1), getPoint(ind2));
        return {i, j};
    };
    let hex = [ind1, ind2];
    hex[-1] = makePoint(ind0, ind1, ind2);
    hex[5] = makePoint(hex[1], hex[0], hex[-1]);
    hex[2] = makePoint(hex[0], hex[1], hex[-1]);
    hex[3] = makePoint(hex[1], hex[2], hex[-1]);
    hex[4] = makePoint(hex[0], hex[5], hex[-1]);

    let isNew = false, isNegate = false, args = [];
    for (let i = 0; i < 6; i++) {
        if (getPoint(hex[i]).isBimbo)
            isNew = true;
        if (hex[i].i < 0 || hex[i].j < 0) {
            isNegate = true;
        }
        args.push(getPoint(hex[i]));
    }
    let energy = getPLJ(cubs, size, args, s1, s2);

    if (isNew && !isNegate && energy < energyMaximum) {
        for (let i = 0; i < 6; i++) {
            makeRealPoint(cubs, size, hex[i], getPoint(hex[i]), false);
        }
        makeRealPoint(cubs, size, hex[-1], getPoint(hex[-1]), true);
        return hex;
    }
    return null;
}
let tt = 0;
function buildup(cubs, size, points, r0, ri, alpha, a, amountHex, energyMaximum = 0, s1, s2) {
    let riPi = ri / (2 * Math.PI);
    let g = Math.tan(alpha);
    let solve = (p0, p1, p2) => {
        let point = solveEquations(findAndSumVectorPolar(p0, p1, p2), p1, p2, a, a, riPi, g, r0);
        point.isBimbo = true;
        point.isEmpty = true;
        return point;
    };
    let makeElem = (i, j) => ({i, j});
    let makeHex = (ind0, ind1, ind2) => buildHex(cubs, size, s1, s2, points, solve, energyMaximum, ind0, ind1, ind2);
    let len = points[0].length - 1, last, init;
    switch (len % 3) {
        case 0: last = len; init = 0; break;
        case 1: last = len - 1; init = 0; break;
        default: last = len; init = 1;
    }
    let list = new CircleList();
    for (let i = init; i < points.length - 1; i += 2) {
        list.push(makeElem(i, last - 1), makeElem(i, last), makeElem(i + 1, last - init));
    }
    last = points.length - 1;
    init = last % 2;
    for (let i = 0; i < points[last].length - 3; i++) {
        if (!points[last][i].isEmpty && !points[last - 1][i + 1 - init].isEmpty) {
            list.push(makeElem(last - 1, i - init), makeElem(last, i), makeElem(last - 1, i + 1 - init));
        }
    }
    let cur, hex, buildupHex = 0, wasChange = false;
    while(buildupHex < amountHex) {
        cur = list.next();
        if (!cur) {
            if (!wasChange)
                break;
            wasChange = false;
            continue;
        }

        hex = makeHex(cur.p0, cur.p1, cur.p2);
        if (hex) {
            wasChange = true;
            list.deleteHead();
            buildupHex++;
            list.shift(hex[-1], hex[3], hex[4]);
        }
    }
    console.log('newton iteration on one point: ', t / operation, t, operation);
    for (let i in points) {
        for (let j in points[i]) {
            points[i][j] = points[i][j].toPoint();
        }
    }
    return points;
}