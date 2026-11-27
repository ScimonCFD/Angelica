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
