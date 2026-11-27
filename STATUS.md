# Estado del trabajo autónomo hacia la alternativa open source a PIPEPHASE

Objetivo: prestaciones básicas de PIPEPHASE para redes de recolección de
petróleo, gas y agua en régimen permanente.  Fuente de verdad del plan:
`ROADMAP.md`.  Este archivo resume cada release en español para el
maintainer y recoge lo que espera su decisión.

## Pendiente del maintainer

1. **PyPI**: último paso, decisión del maintainer.
2. **Peneloux en el modelo de flujo composicional** (ROADMAP, limitaciones
   de `Well`): <1 % en el condensado del tutorial 09, 12 % en un crudo rico
   en n-C16.  Cambiaría resultados composicionales existentes; sería opt-in.
3. **Salto laminar/turbulento** (ROADMAP, etapa 5, limitaciones).  64/Re
   por debajo de Re 2300 y Colebrook por encima (0.0278 frente a 0.0475 a
   ε/D 3e-4).  Un tubo con Re ≈ 2300 puede no tener solución.  Nuevo dato
   de v1.28.0: ocurre también en un solve **isotérmico** (black-oil, red
   de dos lazos a 5, 5.5 y 6 °C; converge a 4.5 y 7 °C), así que no es un
   efecto del acoplamiento térmico sino de la discontinuidad misma.
   Opciones: (a) dejarlo con el diagnóstico; (b) interpolar entre Re 2300
   y 4000; (c) Churchill (1977), una ecuación para todo Re (cifras de
   memoria, a leer en la fuente).  (b) y (c) cambian todo tubo con Re
   entre 2300 y 4000.
4. **Viscosidad black-oil por debajo de 21.1 °C** (ROADMAP, etapa 5,
   limitaciones): (a) extrapolar y avisar (actual); (b) congelar en
   21.1 °C; (c) otra correlación con datos a baja temperatura; (d) rechazar
   el caso.  Impacto ≤ 1.4 % en los tutoriales salvo el 07 (31.8 % → 11.2 %).

## Releases

### v1.37.0 — Ítem 2b.3: inhibidores de hidratos, metanol y MEG

**Qué cambió.**  Módulo `inhibitor`: actividad del agua en soluciones de
metanol y MEG (0–50 % en masa en el agua libre).  Todas las funciones de
hidratos, `hydrate_partition` (el inhibidor se concentra al formarse
hidrato), los fluidos black-oil y composicional, el margen de hidratos de
ambos solvers y la GUI (diálogos de fluido y curva de la envolvente)
aceptan un inhibidor.  En black-oil, sin composición, la curva de Katz se
desplaza la menor de dos depresiones de referencia (metano sI y gas
C1/C2/C3 sII).  Ningún resultado existente cambia (es opcional).

**Comparado con fuente (mismo laboratorio).**  Metanol frente a muestras
nuevas del mismo laboratorio (Semenov et al. 2024, Data in Brief 53,
110138, leída en Europe PMC; mismo equipo y método que los datos del
ajuste, así que no es validación independiente), 5–50 %, 3–13 MPa: error
en la depresión ≤ 0.07 K al 5 %, 0.19 al 10 %, 0.39 al 20 %, 0.78 al 30 %,
0.98 al 40 %, 1.39 K al 50 %.  Signo: a 3 MPa la subestima (lado seguro);
por encima de ~7 MPa la sobrestima (lado inseguro, +0.3 a +1.0 K a
9.7–12.8 MPa), y la temperatura absoluta de hidrato, la que usa el margen,
sale hasta 1.89 K fría (50 %, 12.76 MPa).  La depresión medida es plana en
presión; la del modelo crece (30.2 K a 3 MPa, 32.7 a 13, 34.4 a 30): por
encima de 13 MPa no validado y probablemente más del lado inseguro.  Metano
puro 274–289 K: el modelo vdW-P sale 0.60–0.91 K frío.

**Solo consistencia.**  MEG: ajustado al 10/30/50 % y comprobado al 20 y
40 % del mismo conjunto (≤ 0.13 y ≤ 0.41 K).  No hay datos de MEG de otra
fuente leídos.

**No validado / supuesto.**  Otros gases: ajustado solo con metano (a
7 MPa y 30 % de metanol la depresión del modelo es 16.4 K para metano,
14.9 K para C1/C2/C3 y 20.9 K para CO2, sin datos); el
inhibidor no pasa al gas ni al aceite; sin sales; sin fase hielo (con
inhibidor, el agua se trata como líquida hasta −60 °C).  Los puntos de
hielo que implica la actividad del metanol salen 1.0–8.9 K altos al
20–50 % (por eso no se usan como límite).  El desplazamiento black-oil
ignora la dependencia del gas (1.45 K entre los dos gases de referencia al
30 % de metanol; se toma el menor, del lado seguro).

**Revisión antes del tag.**  `hydrate_partition` con inhibidor y sin agua
libre no daba hidrato (lado inseguro; corregido); el desplazamiento
black-oil bajo 3 MPa usaba el valor de 3 MPa (0.14–0.53 K de más; ahora la
presión del tubo); se corrigieron la etiqueta "validado", el signo del
error a alta presión, la tendencia con la presión y la dependencia del gas.

**Qué se probó y se descartó.**  Actividad a partir de los puntos de hielo
o de UNIFAC: sobrestimaba la depresión hasta 6.3 K al 50 % de metanol (lado
inseguro).  Ng & Robinson (1985), el criterio del ROADMAP, no fue
accesible; se usaron los datos de Semenov.

### v1.36.0 — Ítem 2b.2: el hidrato como fase (cantidad formada)

**Qué cambió.**  `hydrate_partition(names, zs, T, P)`: reparte una
alimentación con agua entre fluido e hidrato en equilibrio (igual
potencial químico del agua, mínimo de Gibbs con vdW-P).  Da el agua y los
huéspedes en el hidrato, el fluido restante (agua libre, vapor de agua,
cada hidrocarburo), el número de hidratación y qué detuvo el crecimiento:
equilibrio, agua agotada, formador agotado o ninguno.  Vale para gas–agua
y petróleo–gas–agua (hidrocarburos bifásicos; el heptano no entra).  Solo
API: ningún resultado de solver, GUI o tutorial cambia.

**Comparado, no validado.**  Número de hidratación del metano: Sum,
Burruss y Sloan (1997), 6.04 ± 0.03 (extracto del resumen vía buscador;
condiciones no leídas).  Modelo: 5.895 a 2 °C y 70 bar (−2.4 %, unas 5
veces su incertidumbre).  Discrepancia registrada.

**Validado por consistencia.**  La cantidad incipiente es cero en la curva
vdW-P y crece linealmente debajo (0.198 mol/K); balances de agua y de cada
formador cierran a 1e-12; el gas restante queda en su propia curva de
hidrato (o hidrato–vapor) a T, a 1e-5 K.

**No validado / supuesto.**  Una sola estructura (si la otra también es
estable al final se marca `other_structure_stable`; no se calcula el
equilibrio con dos hidratos: entonces el resultado no es el mínimo de
Gibbs y cerca del cambio de estructura la cantidad varía bruscamente,
0.300 → 0.160 mol entre 6.85 y 6.90 °C en un caso); un gas denso que la
EOS etiqueta líquido no lleva vapor de agua y la cantidad salta (+32 %
entre 130 y 140 bar en un caso); con muy poco gas (≤ 1e-5) tarda 13–27 s;
sin hielo (T ≥ 0 °C); solo composicional (black-oil no tiene composición
que repartir).

**Revisión antes del tag.**  Hidrocarburo sin vapor (propano o CO2
líquidos, gas denso) salía "equilibrium" con G < 0 (ahora "water"); con
poco gas la bisección dejaba 0.4 % del gas sin consumir (ahora se acota
por décadas); un fallo del flash ya no supone todo vapor en silencio.

### v1.35.0 — Ítem 2b.1: contenido de agua del gas e hidratos con gas no saturado

**Qué cambió.**  Módulo `water_content` (contenido de agua en ppm, mg/Sm³
o fracción de saturación; punto de rocío de agua).  Las funciones de
hidratos aceptan el contenido de agua: bajo saturación usan el equilibrio
hidrato–vapor.  En composicional, una corriente toda vapor con agua se
evalúa con su contenido de agua; con líquido, se mantiene la curva
saturada (conservador).  El cambio en VF = 0.999 es un salto: un gas con
60 ppm pasa de 0.42 °C a 17.03 °C con una traza de condensado.  Los
tutoriales cambian como máximo 7·10⁻¹⁰ K.  La curva de hidratos de la GUI
aplica la misma regla punto a punto (`hydrate_curve_at_water_content`).

**Revisión antes del tag.**  La GUI usaba el agua del feed sin mirar las
fases (ahora regla del solver); entradas fuera de rango (rocío fuera de
−50 … 373.85 °C, NaN, P ≤ 0, fracciones > 1, agua negativa) dan error en
vez de un número recortado; el test del tubo compara el valor exacto
(0.420 °C) con una llamada directa; 42.21 mol/Sm³ en el docstring.

**Validado contra fuente.**  Presión de saturación del agua frente a NIST
WebBook (IAPWS-95), 1–91 °C: ≤ 7.2·10⁻⁵; el punto de rocío la invierte a
0.002 K.

**Validado por consistencia.**  Gas saturado: la curva saturada exacta.  Gas
no saturado: la temperatura de hidrato queda entre su punto de rocío y la
curva saturada (70 bar: rocío 5 °C → hidrato 7.22 °C; saturada 17.03 °C).

**No validado / supuesto.**  El contenido de agua a presión usa vapor ideal
(Raoult) y subestima el agua real a alta presión; el equilibrio
hidrato–vapor no se ha comparado con datos medidos.  Sesgo conocido de la
actividad del agua y_w·P/P_sat: omite Poynting (+5.5 % en a_w a 70 bar y
10 °C, calculado) y el coeficiente de fugacidad del agua en el vapor
(+10–15 %, estimado): la temperatura de hidrato del gas no saturado sale
≈ 0.6–2 K más alta (estimación; del lado seguro).  Una fase densa sobre la
cricondenbara que el flash etiqueta líquida toma la curva saturada.

### v1.34.0 — Ítem 7.4: informe de perfiles por conexión

**Qué cambió.**  `profile_table(case, result)`: para cada conexión, sus
componentes en orden (los tramos de un tubo segmentado), con distancia y
elevación acumuladas, presión y temperatura en cada extremo, caída de
presión, caudal y, en tubos bifásicos, patrón, holdup, velocidad de mezcla
y relación erosional.  Sección/hoja "Profiles" en los informes CSV y Excel
(la GUI la incluye).

**Revisión independiente.**  Encontró que la GUI podía exportar el perfil
de una red distinta de la resuelta (la escena editada después del
cálculo) y que el agrupamiento fallaba con ids vacíos o con un prefijo
"link_" que no encadena.  Corregido: se guarda el caso de cada ejecución y
un caso que no corresponde se rechaza.

**Validado por consistencia.**  Los tramos encajan (nodo, distancia,
elevación, presión); las distancias suman la longitud, las elevaciones el
desnivel y las caídas la caída total.  No hay nada que validar contra
fuente: el informe reordena la solución.

### v1.33.0 — Ítem 7.3: factores de calibración y dimensionamiento de línea

**Qué cambió.**  Multiplicadores por tubo del factor de fricción (todas las
correlaciones) y del holdup de Beggs-Brill; `size_line` prueba diámetros
candidatos de un tubo contra caída de presión, velocidad, velocidad
erosional o presión mínima y elige el menor que cumple.  GUI: campos del
tubo y *Studies → Line Sizing…*.

**Validado analíticamente (medido).**  A caudal fijo en un tubo horizontal
la caída de presión escala exactamente con el multiplicador (≤ 5.9·10⁻⁹,
turbulento y laminar).

**Validado por consistencia.**  Factor 1: todos los tutoriales idénticos.
En red, más fricción da menos caudal, de forma monótona.  Observación
física: en un tubo horizontal con Beggs-Brill, más holdup baja la fricción
bifásica y da más caudal (+21 % con 1.3 en la red de prueba); en subida
añade gravedad.

**Revisión independiente.**  Física correcta: el multiplicador equivale a
escalar la longitud de fricción en todas las correlaciones y modos.  Fallos
de uso corregidos: el criterio erosional se saltaba en silencio en tubos
monofásicos (y el diálogo lo rellenaba), una línea segmentada se
dimensionaba por un solo tramo, los campos no aparecían en escenas
antiguas.  Matiz: ΔP ∝ c solo en incompresible; en gas y multifásico no.

**No validado.**  Frente a datos de campo (los factores son para ajustar a
ellos).

### v1.32.0 — Ítem 7.2: estudios de sensibilidad y análisis nodal

**Qué cambió.**  `angelica.studies`: variar uno o dos parámetros del caso
(presión o caudal de un contorno, temperatura, cualquier campo numérico de
un componente) y resolver cada punto; análisis nodal de un pozo (curva del
pozo frente a la presión que la red necesita en su cabeza, cruce en la rama
estable).  En la GUI, menú *Studies* con *Sensitivity…* y *Nodal
Analysis…* (gráficas).  Tutorial black-oil 12, Parte G.

**Validado contra fuente.**  Nada.

**Validado por consistencia (medido).**  El punto de operación del análisis
nodal es el de la red completa: tutorial 12, pozo 1 6.6933 frente a 6.6932
kg/s y 4.04543 frente a 4.04548 MPa, pozo 2 2.2189 frente a 2.2190; tutorial
composicional 09, 6.2268 = 6.2268 kg/s.  Un pozo muerto no tiene cruce.

**Revisión independiente.**  Encontró que mi primera versión del análisis
nodal daba puntos de operación equivocados sin avisar en cinco casos
(correlación PVT del caso sustituida por Standing, yacimiento compartido,
solver isotérmico, yacimiento sin condición térmica, cabeza de pozo
compartida: errores del 2 % al 127 %) y uno con solo una nota.  Corregido:
los cuatro primeros coinciden con la red a ≤ 3.8·10⁻⁵; la cabeza compartida
y la bisección sin terminar ya no dan punto.  También: un tubo de diámetro
negativo se aceptaba (corregido).

**No validado.**  Frente a un ejemplo publicado de análisis nodal.
**Supuesto / limitaciones.**  Pozos con gas lift no soportados en el
análisis nodal; la GUI de estudios trabaja en SI.

### v1.31.0 — Ítem 7.1: unidades de campo y caudales estándar

**Qué cambió.**  Caudales estándar por tubo (petróleo, gas, agua; Sm³/d a
60 °F y 1 atm, flash de una etapa) en black-oil y composicional, impresos y
exportados también en STB/d, MMscf/d y bbl/d; módulo `angelica.units` con
las definiciones exactas; sistema de unidades "Oil field" en la GUI (psia,
°F, ft, in, lb/s; las temperaturas de contorno, ambiente, consigna del
calentador y del canvas en °F).

**Validado contra fuente.**  Factores de conversión frente a NIST SP 811,
Apéndice B.8 (leído): todos dentro de media unidad de su 7.ª cifra.

**Validado por consistencia (medido).**  En black-oil los caudales estándar
del flujo que entra al separador del tutorial 14 son los del separador:
petróleo 705.49 Sm³/d en ambos, agua a 1.6·10⁻⁶; el gas que sale conserva el
de entrada a 0.11 % (la corriente de solo gas en black-oil lleva GOR 10⁶,
un poco de petróleo).

**Revisión independiente.**  Encontró dos errores físicos míos, corregidos:
el gas estándar de black-oil restaba el Rs que Standing aún da a 60 °F y
1 atm (−0.7 % en el tutorial 14, −13 % con GOR 10); ahora va por
definición, q_g = GOR·q_o, comprobado con la fórmula escrita aparte.  Y en
composicional faltaba el agua que el flash deja en el gas (24 % del agua
de un condensado); ahora cuenta toda el agua como líquido.

**No validado.**  Caudales estándar composicionales frente a un flash
publicado.  Supuesto: flash de una sola etapa a condiciones estándar, que no
es el camino de varias etapas de una planta (da otra cifra de líquido).

### v1.30.0 — Ítem 5.2: compresor

**Qué cambió.**
- `Compressor` en gas, black-oil y composicional: especificado por relación
  de presiones, subida de presión o potencia; eficiencia isentrópica o
  politrópica.  El trabajo isentrópico ∫v·dP se integra en ln P sobre la
  h y la v de cada fluido, así que no necesita una función de entropía y
  sirve igual para los tres.  Resultados en `SolveResult.compressors`
  (potencia, relación, temperaturas de descarga real e isentrópica), en la
  GUI (paleta y canvas) y en la Parte B del tutorial de gas 06.
- Corregido: los solvers de presión absoluta dejaban llegar un nodo a 0 Pa
  y la EOS dividía por cero; ahora el suelo es 1 kPa (los líquidos siguen
  en 0 Pa, usan presión manométrica).

**Validado contra fuente.**  Ninguna fuente leída (los ejemplos de GPSA no
están en el repositorio).

**Validado analíticamente o contra la EOS (medido).**
- Gas ideal con cp constante: trabajo y T₂ₛ = T₁·r^(R/cp) a ≤ 8.8·10⁻⁵.
- La comprobación pedida del límite isentrópico: el final de la
  integración cae sobre la isentrópica de la EOS, s(P₂,T₂ₛ) = s(P₁,T₁), con
  la entropía calculada aparte: T₂ₛ·|Δs| ≤ 1.3·10⁻⁴ del trabajo (PR
  metano) y ≤ 4.7·10⁻⁵ (gas natural composicional, entropía del flash de
  thermo).

**Validado solo por consistencia.**  En red: relación o subida cumplidas en
las presiones de los nodos a 10⁻⁴, potencia a 3·10⁻⁴; potencia = ṁ·Δh
recalculada con el fluido; masa coherente a 2·10⁻³; balance de energía ≤
1.5·10⁻³ %.  Barrido en gas:
relaciones 1.7–6 y potencias 0.1–10 MW convergen; relación ≤ 1.5 (45 bar
contra 50) da caudal nulo, como debe.

**No validado.**  Frente a un caso publicado de compresor; la compresión
en black-oil (no hay función de entropía con qué comprobarla).

**Medido frente a supuesto.**  Cifras medidas.  Supuesto: un compresor sin
mapa (sin curvas de altura y eficiencia frente a caudal, sin surge), una
etapa, sin interenfriamiento.

**Hallazgo del tutorial 06.**  La "bomba como compresor" que se usaba hasta
ahora (w = ΔP/ρ) da 92.4 kJ/kg donde la isentrópica del gas da 108.5
(+17.5 %) y el gas sale 6.5 K más caliente: para gas, conviene el
`Compressor`.

**Revisión independiente.**  Encontró fallos serios de mi primera versión,
todos corregidos: daba por convergida una especificación imposible (ambos
extremos en contorno: 5 370 kg/s y 553 MW) y no cumplía la relación con
una tobera pequeña, porque comparaba la subida impuesta en vez de las
presiones de los nodos; el suelo de 1 kPa pisaba un contorno de 0.5 kPa;
un nodo de caudal impuesto llegaba a −72 kPa y lanzaba excepción;
máquinas paradas se daban por no convergidas.  Además añadí el aviso de
alimentación no gaseosa y la exportación de resultados.

**Limitaciones.**  Relación 10 o 30 MW en la estación de prueba (descarga
≈ 280 °C, succión casi vacía) y potencias de 1–30 kW (menos del 1 % del
caudal de diseño): no convergido con diagnóstico.  Compresores en
paralelo con especificaciones distintas no convergen (iguales sí).

### v1.29.0 — Ítem 5.3: calentador/enfriador por temperatura de salida

**Qué cambió.**
- `HeatSource(outlet_temperature_c=…, max_duty_w=…)` en los cuatro modos
  térmicos: el solver calcula el duty, ṁ·(h_out − h_in) con la entalpía del
  fluido (calor latente incluido en composicional), calentando o enfriando
  según haga falta.  Por encima de `max_duty_w` trabaja al máximo y el
  resultado dice "duty limited".  Resultados en `SolveResult.heaters`, en
  el informe, CSV, Excel y en el canvas de la GUI.
- Tutorial black-oil 16 (estación de calentamiento a mitad de un oleoducto
  de 50 km) y Parte D del tutorial composicional 11 (precalentador del city
  gate a 45 °C con 600 kW como máximo).
- Corregidos por el camino (eran previos):
  - un calentador de duty fijo en gas no conservaba la energía (salida con
    cp constante): con 500 kW en la red de prueba, 7.1 % de error en metano
    y 5.9 % en el condensado (salida de la red 1.06 y 1.67 K demasiado
    caliente); 0.39 % en el tutorial 11.  Ahora ≤ 4.3·10⁻⁴ %.  Black-oil no
    cambia.  El tutorial 11 Parte B cambia hasta 0.24 K;
  - un calentador sin caída de presión bloqueaba una red de agua (246 kg/s
    por el calentador en una línea de 64);
  - una red de gas en reposo con un accesorio lanzaba una excepción
    (−310 °C en un nodo);
  - los solvers de gas y de líquido térmico no daban diagnóstico al no
    converger.

**Validado contra fuente.**  Duty del agua frente a entalpías de NIST
WebBook a 0.5 MPa: +0.048 % (20 → 60 °C) y −0.044 % (80 → 10 °C);
tolerancia 0.51 %, la dispersión del cp de NIST alrededor del cp constante
del modelo entre 5 y 80 °C.

**Revisión independiente.**  Encontró una regresión que yo había
introducido: mi primer arreglo del acoplamiento del calentador (suelo de
0.01 m/s) impedía converger a un calentador tras un pozo muerto.  Lo cambié
a 10⁻⁴ m/s y quedó como test.  También corrigió mis cifras: el error del
calentador en gas era 7.1 %, no solo el 0.39 % del tutorial 11; y la
cercanía del nodo a la consigna (0.0040 K) es de esta red, no una cota
general (≈ 10⁻³ × T en °C).

**Validado solo por consistencia (medido).**  En la red de dos lazos con un
calentador en cada modo: salida exactamente en la consigna; el nodo
siguiente a ≤ 0.0040 K; duty = ṁ·Δh a ≤ 7.4·10⁻⁶; balance de energía ≤
1.9·10⁻⁴ %; límite de duty, límite cero y flujo invertido correctos;
consignas de −60 a 300 °C convergen (salvo black-oil bajo 0 °C).

**No validado.**  El duty en gas, black-oil y composicional frente a un
calentador publicado: descansa en el modelo de entalpía de cada fluido.

**Medido frente a supuesto.**  Todas las cifras son medidas.  Se supone
que el calentador entrega todo su duty al fluido (sin rendimiento de
combustión ni aproximación del lado de proceso).

**Limitaciones registradas (sin decisión necesaria, salvo la última).**
- Red en reposo con un elemento de ley cuadrática (calentador, accesorio,
  choke): no se asienta (ruido de presión → caudal); no convergido con
  diagnóstico, sin excepción.  Remedio posible: ley lineal a muy baja
  velocidad — sería un cambio de modelo, no lo hice.
- Black-oil enfriado por debajo de ~0 °C: no converge (viscosidad
  extrapolada + salto laminar/turbulento), diagnosticado.  Depende de las
  decisiones 3 y 4 de arriba.

### v1.28.0 — Parte 1 del objetivo: lo que ya existe, comprobado con números

**Qué cambió.**
- Nueva opción isotérmica en gas, black-oil y composicional
  (`isothermal_temperature_c`; en la GUI, *Numerics → Isothermal at*).
  No existía: esos tres solvers resolvían siempre la energía.  Con ella,
  nodos, tubos, separadores y tubing de pozos quedan a esa temperatura; los
  calentadores y las condiciones térmicas a otra temperatura se ignoran con
  aviso.  Desactivada por defecto: ningún resultado existente cambia.
- Red de prueba con dos lazos independientes y dos entradas en los cinco
  modos (agua isotérmica y térmica, metano PR, black-oil con gas y agua,
  composicional de siete componentes con agua libre y condensado), en
  térmico y en isotérmico: 26 tests (más 19 del modo isotérmico).
- Corregido: las temperaturas fijadas (contornos y modo isotérmico) deben
  ser finitas y mayores que −273.15 °C.  Antes, NaN o −300 °C daban
  excepciones en gas y un black-oil "convergido" a −300 °C.
- Corregido: en el diálogo *Numerics* de la GUI la casilla "Max
  Temperature Iterations" estaba encima del combo del esquema de
  convección.

**Validado contra fuente.**  Nada nuevo.  La red mallada publicada sigue
siendo Hanoi (agua; resultados de EPANET, tutorial 07).  No encontré una
red mallada de gas ni multifásica publicada con resultados.

**Validado solo por consistencia (medido).**
- Balance de masa en cada nodo ≤ 1.07·10⁻⁵ (criterio del solver 10⁻³:
  casi automático si converge).
- La ley de lazos se cumple por construcción (una presión por nodo).  Lo
  comprobado es que cada ΔP es el de su correlación a su caudal: cada tubo
  resuelto solo, con el caudal y el estado aguas arriba de la red, da la
  presión aguas abajo de la red a ≤ 3.02·10⁻⁵ de su ΔP (tolerancia 10⁻⁴).
  La revisión mostró que mi test de "suma en el lazo" era circular; lo
  quité.
- Balance de energía de las corridas térmicas ≤ 2.77·10⁻⁴ % (cifra del
  ROADMAP 0.01 %).
- Margen de hidratos por tubo (black-oil y composicional) igual a la
  función independiente en el extremo más frío.
- Isotérmico: tubing de pozo a temperatura constante (con Ramey se
  enfría); fuera de rango (−40 a 200 °C, red en reposo, flujo invertido,
  pozos, gas lift, separadores, reguladores) todo converge salvo black-oil
  a 5–6 °C (salto laminar/turbulento, diagnosticado).  Ninguna excepción.
  Los avisos nombran lo que el modo isotérmico deja fuera: compresores,
  efecto Joule-Thomson en válvulas, intercambio de los pozos con la
  formación, contornos de gradiente.

**Revisión independiente (antes del tag).**  Encontró: la GUI seguía
exigiendo un contorno de temperatura fija en black-oil y composicional;
temperaturas no físicas aceptadas; equipos cuyo cambio de temperatura se
perdía sin aviso; el test circular; tests de rango que no exigían
convergencia.  Todo corregido.  También señaló algo previo, no corregido
aquí: el solver de gas no da diagnóstico cuando no converge (anotado en el
ROADMAP).

**No validado.**  Lazos de gas y multifásicos frente a un caso publicado;
la envolvente de mezclas frente a datos medidos (solo metano puro frente a
NIST); la curva de hidratos de mezclas naturales (ver `test_hydrate.py`).

**Medido frente a supuesto.**  Todas las cifras de arriba son medidas en
esta versión.  Se supone que la coherencia entre la red y el tubo aislado
basta para probar las leyes de Kirchhoff en los lazos: comprueba que cada
ΔP es el de su correlación a su caudal, no que la correlación sea correcta.
