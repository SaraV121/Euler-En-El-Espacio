
/* =====================================================
   EULER EN EL ESPACIO
   Métodos Numéricos
   Euler - Heun - Runge-Kutta 4 - Solución Analítica
===================================================== */


/* =====================================================
   1. REFERENCIAS A LOS ELEMENTOS HTML
===================================================== */

const masaInput = document.getElementById("masa");
const empujeInput = document.getElementById("empuje");
const gravedadInput = document.getElementById("gravedad");
const resistenciaInput = document.getElementById("resistencia");

const velocidadInicialInput = document.getElementById("velocidadInicial");
const alturaInicialInput = document.getElementById("alturaInicial");

const tiempoFinalInput = document.getElementById("tiempoFinal");
const pasoInput = document.getElementById("paso");

const btnSimular = document.getElementById("btnSimular");
const btnRestablecer = document.getElementById("btnRestablecer");

const mensaje = document.getElementById("mensaje");

const velocidadResultado = document.getElementById("velocidadResultado");
const alturaResultado = document.getElementById("alturaResultado");

const velocidadExacta = document.getElementById("velocidadExacta");
const alturaExacta = document.getElementById("alturaExacta");

const tablaResultados = document.getElementById("tablaResultados");
const tablaErrores = document.getElementById("tablaErrores");

const precisionMensaje = document.getElementById("precisionMensaje");


/* =====================================================
   2. VARIABLES GLOBALES
===================================================== */

let graficaVelocidad = null;
let graficaAltura = null;


/* =====================================================
   3. OBTENER Y VALIDAR LOS PARÁMETROS
===================================================== */

function obtenerParametros() {

    const parametros = {
        masa: Number(masaInput.value),
        empuje: Number(empujeInput.value),
        gravedad: Number(gravedadInput.value),
        resistencia: Number(resistenciaInput.value),

        velocidadInicial: Number(velocidadInicialInput.value),
        alturaInicial: Number(alturaInicialInput.value),

        tiempoFinal: Number(tiempoFinalInput.value),
        paso: Number(pasoInput.value)
    };

    const valores = Object.values(parametros);

    if (!valores.every(Number.isFinite)) {
        throw new Error("Todos los parámetros deben ser números válidos.");
    }

    if (parametros.masa <= 0) {
        throw new Error("La masa debe ser mayor que cero.");
    }

    if (parametros.empuje < 0) {
        throw new Error("El empuje no puede ser negativo.");
    }

    if (parametros.gravedad < 0) {
        throw new Error("La gravedad no puede ser negativa.");
    }

    if (parametros.resistencia < 0) {
        throw new Error("La resistencia no puede ser negativa.");
    }

    if (parametros.tiempoFinal <= 0) {
        throw new Error("El tiempo final debe ser mayor que cero.");
    }

    if (parametros.paso <= 0) {
        throw new Error("El tamaño de paso debe ser mayor que cero.");
    }

    if (parametros.tiempoFinal / parametros.paso > 10000) {
        throw new Error("El número de iteraciones es demasiado grande.");
    }

    return parametros;
}


/* =====================================================
   4. DEFINICIÓN DEL MODELO DIFERENCIAL
===================================================== */

// v' = T/m - g - (k/m)v
// y' = v

function derivadas(v, parametros) {

    const aceleracion =
        parametros.empuje / parametros.masa
        - parametros.gravedad
        - (parametros.resistencia / parametros.masa) * v;

    return {
        dv: aceleracion,
        dy: v
    };
}


/* =====================================================
   5. SOLUCIÓN ANALÍTICA
===================================================== */

function solucionExacta(t, parametros) {

    const A =
        parametros.empuje / parametros.masa
        - parametros.gravedad;

    const B =
        parametros.resistencia / parametros.masa;

    const v0 = parametros.velocidadInicial;
    const y0 = parametros.alturaInicial;

    let v, y;

    if (Math.abs(B) < 1e-12) {

        v = v0 + A * t;

        y = y0 + v0 * t + 0.5 * A * t * t;

    } else {

        const velocidadLimite = A / B;

        const exponencial = Math.exp(-B * t);

        v = velocidadLimite
            + (v0 - velocidadLimite) * exponencial;

        y = y0
            + velocidadLimite * t
            + (v0 - velocidadLimite) * (1 - exponencial) / B;

    }

    return { v, y };
}


/* =====================================================
   6. MÉTODO DE EULER
===================================================== */

function metodoEuler(v, y, h, parametros) {

    const d = derivadas(v, parametros);

    const nuevaVelocidad = v + h * d.dv;

    const nuevaAltura = y + h * d.dy;

    return {
        v: nuevaVelocidad,
        y: nuevaAltura
    };
}


/* =====================================================
   7. MÉTODO DE HEUN
===================================================== */

function metodoHeun(v, y, h, parametros) {

    // Primera pendiente

    const k1 = derivadas(v, parametros);

    // Predicción de Euler

    const vPredicho = v + h * k1.dv;

    const yPredicho = y + h * k1.dy;

    // Segunda pendiente

    const k2 = derivadas(vPredicho, parametros);

    // Corrección mediante el promedio

    const nuevaVelocidad =
        v + (h / 2) * (k1.dv + k2.dv);

    const nuevaAltura =
        y + (h / 2) * (k1.dy + k2.dy);

    return {
        v: nuevaVelocidad,
        y: nuevaAltura
    };
}


/* =====================================================
   8. MÉTODO RUNGE-KUTTA DE CUARTO ORDEN
===================================================== */

function metodoRK4(v, y, h, parametros) {

    // Pendiente 1

    const k1 = derivadas(v, parametros);

    // Pendiente 2

    const k2 = derivadas(
        v + (h / 2) * k1.dv,
        parametros
    );

    // Pendiente 3

    const k3 = derivadas(
        v + (h / 2) * k2.dv,
        parametros
    );

    // Pendiente 4

    const k4 = derivadas(
        v + h * k3.dv,
        parametros
    );

    // Cálculo de la nueva velocidad

    const nuevaVelocidad = v + (h / 6) * (
        k1.dv + 2 * k2.dv + 2 * k3.dv + k4.dv
    );

    // Cálculo de la nueva altura

    const nuevaAltura = y + (h / 6) * (
        k1.dy + 2 * k2.dy + 2 * k3.dy + k4.dy
    );

    return {
        v: nuevaVelocidad,
        y: nuevaAltura
    };
}


/* =====================================================
   9. EJECUTAR LOS CUATRO MÉTODOS
===================================================== */

function ejecutarSimulacion(parametros) {

    let vEuler = parametros.velocidadInicial;
    let yEuler = parametros.alturaInicial;

    let vHeun = parametros.velocidadInicial;
    let yHeun = parametros.alturaInicial;

    let vRK4 = parametros.velocidadInicial;
    let yRK4 = parametros.alturaInicial;

    let tiempo = 0;

    const resultados = [];

    // Registrar las condiciones iniciales

    resultados.push({
        tiempo: 0,

        eulerV: vEuler,
        eulerY: yEuler,

        heunV: vHeun,
        heunY: yHeun,

        rk4V: vRK4,
        rk4Y: yRK4,

        exactaV: parametros.velocidadInicial,
        exactaY: parametros.alturaInicial
    });

    // Realizar las iteraciones

    while (tiempo < parametros.tiempoFinal - 1e-12) {

        const h = Math.min(
            parametros.paso,
            parametros.tiempoFinal - tiempo
        );

        // Euler

        const resultadoEuler = metodoEuler(
            vEuler,
            yEuler,
            h,
            parametros
        );

        vEuler = resultadoEuler.v;
        yEuler = resultadoEuler.y;

        // Heun

        const resultadoHeun = metodoHeun(
            vHeun,
            yHeun,
            h,
            parametros
        );

        vHeun = resultadoHeun.v;
        yHeun = resultadoHeun.y;

        // Runge-Kutta de cuarto orden

        const resultadoRK4 = metodoRK4(
            vRK4,
            yRK4,
            h,
            parametros
        );

        vRK4 = resultadoRK4.v;
        yRK4 = resultadoRK4.y;

        // Actualizar tiempo

        tiempo = Math.min(
            tiempo + h,
            parametros.tiempoFinal
        );

        // Solución analítica

        const exacta = solucionExacta(tiempo, parametros);

        // Guardar resultados

        resultados.push({

            tiempo: tiempo,

            eulerV: vEuler,
            eulerY: yEuler,

            heunV: vHeun,
            heunY: yHeun,

            rk4V: vRK4,
            rk4Y: yRK4,

            exactaV: exacta.v,
            exactaY: exacta.y

        });

    }

    return resultados;
}


/* =====================================================
   10. MOSTRAR RESULTADOS PRINCIPALES
===================================================== */

function mostrarResultados(resultados) {

    const ultimo = resultados[resultados.length - 1];

    velocidadResultado.textContent =
        ultimo.eulerV.toFixed(4);

    alturaResultado.textContent =
        ultimo.eulerY.toFixed(4);

    velocidadExacta.textContent =
        ultimo.exactaV.toFixed(4);

    alturaExacta.textContent =
        ultimo.exactaY.toFixed(4);
}


/* =====================================================
   11. GENERAR TABLA DE ITERACIONES
===================================================== */

function generarTabla(resultados) {

    tablaResultados.innerHTML = "";

    resultados.forEach(r => {

        const fila = document.createElement("tr");

        const valores = [
            r.tiempo.toFixed(2),

            r.eulerV.toFixed(4),
            r.eulerY.toFixed(4),

            r.heunV.toFixed(4),
            r.heunY.toFixed(4),

            r.rk4V.toFixed(4),
            r.rk4Y.toFixed(4),

            r.exactaV.toFixed(4),
            r.exactaY.toFixed(4)
        ];

        valores.forEach(valor => {

            const celda = document.createElement("td");

            celda.textContent = valor;

            fila.appendChild(celda);

        });

        tablaResultados.appendChild(fila);

    });
}


/* =====================================================
   12. CALCULAR Y MOSTRAR ERRORES
===================================================== */

function calcularErrores(resultados) {

    const ultimo = resultados[resultados.length - 1];

    const metodos = [
        {
            nombre: "Euler",
            velocidad: ultimo.eulerV,
            altura: ultimo.eulerY
        },

        {
            nombre: "Heun",
            velocidad: ultimo.heunV,
            altura: ultimo.heunY
        },

        {
            nombre: "Runge-Kutta 4",
            velocidad: ultimo.rk4V,
            altura: ultimo.rk4Y
        }
    ];

    tablaErrores.innerHTML = "";

    const errores = [];

    metodos.forEach(metodo => {

        const errorV = Math.abs(
            ultimo.exactaV - metodo.velocidad
        );

        const errorY = Math.abs(
            ultimo.exactaY - metodo.altura
        );

        errores.push({
            nombre: metodo.nombre,
            velocidad: errorV,
            altura: errorY
        });

        const fila = document.createElement("tr");

        fila.innerHTML = `
            <td>${metodo.nombre}</td>
            <td>${errorV.toFixed(6)} m/s</td>
            <td>${errorY.toFixed(6)} m</td>
        `;

        tablaErrores.appendChild(fila);

    });

    return errores;
}


/* =====================================================
   13. EVALUAR LA PRECISIÓN
===================================================== */

function evaluarPrecision(errores) {

    const euler = errores.find(e => e.nombre === "Euler");

    if (euler.velocidad < 0.5 && euler.altura < 5) {

        precisionMensaje.textContent =
            "✓ EULER CUMPLE LOS DOS REQUISITOS DE PRECISIÓN.";

        precisionMensaje.style.color = "#34d399";

        precisionMensaje.style.borderColor = "#34d399";

    } else {

        precisionMensaje.textContent =
            "⚠ LA PRECISIÓN DE EULER TODAVÍA ES INSUFICIENTE. PRUEBA UN PASO MENOR.";

        precisionMensaje.style.color = "#fbbf24";

        precisionMensaje.style.borderColor = "#fbbf24";

    }

}


/* =====================================================
   14. CREAR GRÁFICA DE VELOCIDAD
===================================================== */

function crearGraficaVelocidad(resultados) {

    const ctx = document.getElementById("graficaVelocidad");

    if (graficaVelocidad) {
        graficaVelocidad.destroy();
    }

    graficaVelocidad = new Chart(ctx, {

        type: "line",

        data: {

            labels: resultados.map(r => r.tiempo.toFixed(2)),

            datasets: [

                {
                    label: "Euler",
                    data: resultados.map(r => r.eulerV),
                    borderColor: "#38bdf8",
                    backgroundColor: "#38bdf8",
                    borderWidth: 2,
                    pointRadius: 2,
                    tension: 0.2
                },

                {
                    label: "Heun",
                    data: resultados.map(r => r.heunV),
                    borderColor: "#fbbf24",
                    borderWidth: 2,
                    pointRadius: 1,
                    tension: 0.2
                },

                {
                    label: "RK4",
                    data: resultados.map(r => r.rk4V),
                    borderColor: "#34d399",
                    borderWidth: 2,
                    pointRadius: 1,
                    tension: 0.2
                },

                {
                    label: "Solución exacta",
                    data: resultados.map(r => r.exactaV),
                    borderColor: "#f472b6",
                    borderWidth: 2,
                    pointRadius: 0,
                    tension: 0.2
                }

            ]

        },

        options: {

            responsive: true,

            maintainAspectRatio: true,

            interaction: {
                intersect: false,
                mode: "index"
            },

            plugins: {

                legend: {
                    labels: {
                        color: "#cbd5e1"
                    }
                }

            },

            scales: {

                x: {
                    title: {
                        display: true,
                        text: "Tiempo (segundos)",
                        color: "#cbd5e1"
                    },

                    ticks: {
                        color: "#94a3b8"
                    },

                    grid: {
                        color: "#263b57"
                    }
                },

                y: {
                    title: {
                        display: true,
                        text: "Velocidad (m/s)",
                        color: "#cbd5e1"
                    },

                    ticks: {
                        color: "#94a3b8"
                    },

                    grid: {
                        color: "#263b57"
                    }
                }

            }

        }

    });

}


/* =====================================================
   15. CREAR GRÁFICA DE ALTURA
===================================================== */

function crearGraficaAltura(resultados) {

    const ctx = document.getElementById("graficaAltura");

    if (graficaAltura) {
        graficaAltura.destroy();
    }

    graficaAltura = new Chart(ctx, {

        type: "line",

        data: {

            labels: resultados.map(r => r.tiempo.toFixed(2)),

            datasets: [

                {
                    label: "Euler",
                    data: resultados.map(r => r.eulerY),
                    borderColor: "#38bdf8",
                    borderWidth: 2,
                    pointRadius: 2,
                    tension: 0.2
                },

                {
                    label: "Heun",
                    data: resultados.map(r => r.heunY),
                    borderColor: "#fbbf24",
                    borderWidth: 2,
                    pointRadius: 1,
                    tension: 0.2
                },

                {
                    label: "RK4",
                    data: resultados.map(r => r.rk4Y),
                    borderColor: "#34d399",
                    borderWidth: 2,
                    pointRadius: 1,
                    tension: 0.2
                },

                {
                    label: "Solución exacta",
                    data: resultados.map(r => r.exactaY),
                    borderColor: "#f472b6",
                    borderWidth: 2,
                    pointRadius: 0,
                    tension: 0.2
                }

            ]

        },

        options: {

            responsive: true,

            maintainAspectRatio: true,

            interaction: {
                intersect: false,
                mode: "index"
            },

            plugins: {

                legend: {
                    labels: {
                        color: "#cbd5e1"
                    }
                }

            },

            scales: {

                x: {
                    title: {
                        display: true,
                        text: "Tiempo (segundos)",
                        color: "#cbd5e1"
                    },

                    ticks: {
                        color: "#94a3b8"
                    },

                    grid: {
                        color: "#263b57"
                    }
                },

                y: {
                    title: {
                        display: true,
                        text: "Altura (metros)",
                        color: "#cbd5e1"
                    },

                    ticks: {
                        color: "#94a3b8"
                    },

                    grid: {
                        color: "#263b57"
                    }
                }

            }

        }

    });

}


/* =====================================================
   16. FUNCIÓN PRINCIPAL DE SIMULACIÓN
===================================================== */

function iniciarSimulacion() {

    try {

        // Obtener parámetros

        const parametros = obtenerParametros();

        // Ejecutar métodos

        const resultados = ejecutarSimulacion(parametros);

        // Mostrar resultados

        mostrarResultados(resultados);

        // Generar tabla

        generarTabla(resultados);

        // Calcular errores

        const errores = calcularErrores(resultados);

        // Evaluar precisión

        evaluarPrecision(errores);

        // Crear gráficas

        crearGraficaVelocidad(resultados);

        crearGraficaAltura(resultados);

        // Mensaje de confirmación

        mensaje.textContent =
            "✓ Simulación ejecutada correctamente. " +
            resultados.length + " puntos calculados.";

    } catch (error) {

        mensaje.textContent = "⚠ " + error.message;

        mensaje.style.color = "#f87171";

    }

}


/* =====================================================
   17. RESTABLECER PARÁMETROS
===================================================== */

function restablecerValores() {

    masaInput.value = 1000;
    empujeInput.value = 15000;
    gravedadInput.value = 9.81;
    resistenciaInput.value = 100;

    velocidadInicialInput.value = 0;
    alturaInicialInput.value = 0;

    tiempoFinalInput.value = 20;
    pasoInput.value = 2;

    mensaje.style.color = "#34d399";

    mensaje.textContent =
        "Parámetros restablecidos correctamente.";

    iniciarSimulacion();

}


/* =====================================================
   18. EVENTOS DE LOS BOTONES
===================================================== */

btnSimular.addEventListener("click", iniciarSimulacion);

btnRestablecer.addEventListener("click", restablecerValores);


/* =====================================================
   19. EJECUCIÓN INICIAL
===================================================== */

document.addEventListener("DOMContentLoaded", function() {

    iniciarSimulacion();

});
