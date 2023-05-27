function multilateralizationExecution(refPoints, coordsUnknownPoint){
    const distances = [];
    let A_pseudo_inv = 0;

    // Matrices para el sist. de ecuaciones
    let A = [];
    let b = [];

    // Matriz Identidad
    const I = createIdentityMatrix(refPoints[0].length);

    console.log("Puntos Ancla: " + refPoints)
    console.log("Coords punto desconocido: " + coordsUnknownPoint)

    calculateDistances(distances, coordsUnknownPoint);
    generateAEcs(A, refPoints);
    generateBEcs(b, refPoints, distances);

    A_pseudo_inv = math.pinv(A)
    console.log(A_pseudo_inv)
    const x_p = math.multiply(A_pseudo_inv, b)

    // SVD (Descomposición en valores singulares). U matriz unitaria. S matriz diagonal. VT matriz transpuesta conjugada.
    const { u, v, q } = SVDJS.SVD(A)
    let VT = math.transpose(v);
    
    // Vector x_h
    let x_h = VT[VT.length - 2];

    console.log("x_p =", x_p)
    console.log("x_h =", x_h)

    let secondGradeVars = secondGradeVariables(x_p, x_h)
    
    let solutionSecondGradeEc = solve_quadratic(secondGradeVars)
    
    if(solutionSecondGradeEc.length){
        const resultado1 = x_h.map(elemento => solutionSecondGradeEc[0] * elemento);
        const resultado2 = x_h.map(elemento => solutionSecondGradeEc[1] * elemento);
        let sol1 = math.add(x_p, resultado1);
        let sol2 = math.add(x_p, resultado2);
        calculateDifferences(sol1, sol2);
        
        sol1 = sol1.slice(1); // Me interesa quedarme con las coordenadas desconocidas, ya que las primeras lo son.
        sol2 = sol2.slice(1);
            
        const N1 = math.multiply(sol1, I);
        const N2 = math.multiply(sol2, I);
        
        console.log("N1 =", N1);
        console.log("N2 =", N2);
    }else{
        const resultado = x_h.map(elemento => solutionSecondGradeEc[0] * elemento);
        let sol = math.add(x_p, resultado);
        
        sol = sol.slice(1); // Me interesa quedarme con las coordenadas desconocidas, ya que las primeras lo son.
            
        const N = math.multiply(sol, I);
        
        console.log("N =", N);
    }

    // findSolution();
}

function createIdentityMatrix(n) {
    return Array.from({length: n}, (_, i) => {
        return Array.from({length: n}, (_, j) => i === j ? 1 : 0);
    });
}

// Función que calcule las distancias
function calculateDistances(distances, coordsUnknownPoint) {
    for (const coords of refPoints) {
        let value = 0;
        for (let i = 0; i < coordsUnknownPoint.length; i++) {
        value += (coords[i] - coordsUnknownPoint[i]) ** 2;
        }
        const s = Math.sqrt(value);
        distances.push(s);
    }
    console.log("Distancias: " + distances);
}
  
// Generar las ecuaciones de A en función de los puntos de referencia
function generateAEcs(A, refPoints) {
    let newRow = [];
    for (const coords of refPoints) {
    newRow = [1];
    for (const value of coords) {
        newRow.push(-2*value);
    }
    console.log(newRow);
    A.push(newRow);
    }
    console.log("Matriz A = ", A);
}

// Generar las ecuaciones de B en función de los puntos de referencia y distancias
function generateBEcs(b, refPoints, distances) {
  for (let i = 0; i < refPoints.length; i++) {
    const coords = refPoints[i];
    const dist = distances[i];
    let newRow = dist**2;
    for (const value of coords) {
      newRow -= value**2;
    }
    b.push(newRow);
  }
  console.log("Matriz b = ", b);
}

// Calculate EC 2º Grado values
function secondGradeVariables(x_p, x_h) {
    if (x_p.length !== x_h.length) {
        throw new Error("La longitud debe ser la misma");
    } else {
        A_segundoGrado = x_h.slice(1).reduce((sum, val) => sum + val ** 2, 0);

        B_segundoGrado = x_p.slice(1).reduce((sum, val, i) => sum + 2 * val * x_h[i + 1], 0) - x_h[0];

        C_segundoGrado = x_p.slice(1).reduce((sum, val) => sum + val ** 2, 0) - x_p[0];
    }
    return [A_segundoGrado, B_segundoGrado, C_segundoGrado]
}
  

// Ec 2º Grado
function solve_quadratic(secondGradeVariables) {
    let a = secondGradeVariables[0]
    let b = secondGradeVariables[1]
    let c = secondGradeVariables[2]

    let insideSqrt = (b ** 2) - (4 * a * c);

    if (insideSqrt < 0) {
        let realPart = -b / (2 * a);
        return realPart;
    } else {
        let solucion1 = (-b + math.sqrt(insideSqrt)) / (2*a);
        let solucion2 = (-b - math.sqrt(insideSqrt)) / (2*a);
        return [solucion1, solucion2]
    }
}

function calculateDifferences(x1, x2) {
    d1 = x1[0] - x1.slice(1).reduce((total, current) => total + current ** 2, 0);
    d2 = x2[0] - x2.slice(1).reduce((total, current) => total + current ** 2, 0);
}
  
function findSolution() {
    const X = N1;
    const Y = N2;
    const P = refPoints;

    let res1 = 0;
    let res2 = 0;
    let res = [];

    for (let i = 0; i < P.length; i++) {
        res1 += (norm(X, P[i]) ** 2 - distances[i]) ** 2;
        res2 += (norm(Y, P[i]) ** 2 - distances[i]) ** 2;
    }

    res.push(res1);
    res.push(res2);

    const indexMinValue = res.indexOf(Math.min(...res));

    console.log("minimice =", res[indexMinValue]);

    if (indexMinValue == 0) {
        console.log("La solución es N1");
    } else {
        console.log("La solución es N2");
    }
}

function norm(a, b) {
    let distance = 0;
    for (let i = 0; i < a.length; i++) {
        distance += (a[i] - b[i]) ** 2;
    }
    return Math.sqrt(distance);
}