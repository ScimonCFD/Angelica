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
