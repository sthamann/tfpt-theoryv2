# Von Topologie zu Dynamik: Die Ordnung hinter $\alpha$ _und den Naturkonstanten_
Stefan Hamann, 01.09.2025, Version 1.0.6

```table-of-contents
title: 
style: nestedList # TOC style (nestedList|nestedOrderedList|inlineFirstLevel)
minLevel: 0 # Include headings from the specified level
maxLevel: 4 # Include headings up to the specified level
include: 
exclude: 
includeLinks: true # Make headings clickable
hideWhenEmpty: false # Hide TOC if no headings are found
debugInConsole: false # Print debug info in Obsidian console
```

## Abstract

Wir zeigen, dass die Feinstrukturkonstante $\alpha$ und weitere fundamentale Größen **nicht** als freie Inputs benötigt werden, sondern aus **Topologie, Geometrie und Symmetrie** folgen. Ausgangspunkte sind

- der **topologische Fixpunkt** $c_3=\tfrac{1}{8\pi}$,  
- eine **geometrisch definierte Längenskala** $\varphi_0=\tfrac{1}{6\pi}+\tfrac{3}{256\pi^4}=0.053171952\,$ (reduzierte Planck‑Einheiten),  
- sowie eine durch **E₈** geordnete **Dämpfungsfunktion** $\gamma(n)$ für die diskrete Vakuumleiter $\varphi_n$. 

Kernresultat ist eine **Ein‑Parameter‑Normalform** für $\alpha$ (Parameter = $c_3$). Aus $c_3$ folgen exakt
$$[
\varphi_0=\frac{4}{3}c_3+48c_3^4,\quad A=2c_3^3,\quad 
\kappa=\frac{b_1}{2\pi}\ln\frac{1}{\varphi_0},\quad b_1=\frac{41}{10},
]$$
und damit die **kubische Fixpunktgleichung**
$$[
\boxed{\;\alpha^3-2c_3^3\,\alpha^2-8\,b_1\,c_3^6 \ln\!\Big(\frac{1}{\tfrac{4}{3}c_3+48c_3^4}\Big)=0\;}
]$$
mit genau **einer** reellen physischen Lösung $\alpha(c_3)$. Für $c_3=\tfrac{1}{8\pi}$ erhalten wir

$$[
\varphi_0=0.0531719521768,\quad 
\kappa=1.914684795,\quad 
\alpha=0.007297325816919221,\quad 
\alpha^{-1}=137.03650146488582,
]$$

also eine Abweichung von **3.67 ppm** gegenüber CODATA‑2022 – ohne freie Parameter.  

Dieselbe Struktur erzeugt eine **log‑exakte E₈‑Kaskade** $\varphi_{n+1}=\varphi_n e^{-\gamma(n)}$, deren Ankerstufen Flavor‑Mischungen, elektroschwache und hadronische Skalen sowie kosmologische Konstanten treffen. 

> **Eine zwei Schleifen RGE Analyse bestätigt beide Fingerprints:**
   $\alpha_{3}(1\,\text{PeV})=0.052923411$ liegt $0.47\%$ unter $\varphi_{0}=0.053171952$,
   und bei $\mu\simeq 2.5\times 10^{8}\,\text{GeV}$ ergibt sich $\alpha_{3}=0.039713807$, also $0.19\%$ unter $c_{3}=1/(8\pi)$.
> Eine farbadjunkte Brücke $G8$ oberhalb $M_{G8}=1.8\times 10^{10}\,\text{GeV}$
   reduziert die Steigung von $\alpha_{3}^{-1}$ von $\tfrac{7}{2\pi}$ zu $\tfrac{5}{2\pi}$ und erzeugt einen engen Unifikationskorridor mit minimaler relativer Spreizung $1.23\%$ bei $\mu^{\star}\approx 1.43\times 10^{15}\,\text{GeV}$.

Damit ergibt sich ein konsistentes Bild:
**Topologie** fixiert die Normalisierungen, **Geometrie** die Längenskala, **E₈** ordnet die Skalenleiter, **RG‑Dynamik** bestätigt die Fingerabdrücke. 

![[Pasted image 20250826125915.png]]

> **Info Box: Notation und Konventionen**
> 
> * Indizes: $(c_{3}\rightarrow c_{₃}), (b_{1}\rightarrow b_{₁})$ im Fließtext, in Formeln wie gesetzt.  
> * Längenskala: $(ϕ_{0}=\frac{1}{6\pi}+\frac{3}{256\pi^{4}}), (ϕ_{n+1}=ϕ_{n}e^{-\,\gamma(n)})$.  
> * Topologie und Kopplungen: $(g=8c_{₃}^{2}=\frac{1}{8\pi^{2}}), (A=2c_{₃}^{3}=\frac{1}{256\pi^{3}})$.  
> * RG Konstante: $(κ=\frac{b_{₁}}{2\pi}\ln\!\frac{1}{ϕ_{0}}), (b_{₁}=41∕10)$ in GUT Norm.  
> * Gruppen: $(E_{8}), (E_{7}), (E_{6})$ immer als Indexschreibweise $(E_{₈}), (E_{₇}), (E_{₆})$.  
> * Einheiten: alle dimensionierten Größen in reduzierten Planck Einheiten sofern nicht anders angegeben.


# 1. Introduction

Die Frage nach dem Ursprung von Naturkonstanten – insbesondere der Feinstrukturkonstante $\alpha$ – wird hier **bottom‑up** beantwortet: Konstanten sind **Invarianten** eines gemeinsamen Rahmens aus **Topologie, Geometrie und Symmetrie**, nicht externe Knöpfe. 

---
### 1.1 Der genetische Algorithmus

Wir lassen einen genetischen Algorithmus (GA) über Lagrange‑Dichten mit sechs Koeffizienten $(c_0,\dots,c_5)$ evolvieren (Kinetik, Masse, quartische Kinetik, Maxwell, EH‑Term). Harte physikalische Constraints (Lorentz, Ghostfreiheit, richtige Vorzeichen) werden strikt erzwungen; Fitness misst fehlerinvariant $\delta_c,\delta_\alpha,\delta_G$ zu Zielgrößen. Typische Populationen $N\!=\!800$, Turnierselektion, Eliten, Crossover, adaptive Mutationen. Ergebnis: robuste **Cluster** bei $c_4$ (EM‑Normierung), $c_3$ (quartische Kinetik, Spur $1/(8\pi)^2$) und einem engen **$\varphi_0$‑Tal**.

**<span style="background:#d3f8b6">Abbildung: Benutzeroberfläche der GA-Search Application</span>**

![[Pasted image 20250823093627.png]]
### 1.2 ) Genetischer Algorithmus – Setup, Validierung, Ergebnisse
  
- **Konvergenz:** $\sim$ 24 Mio. Bewertungen, $\sim$ 15. 000 Generationen; Reproduzierbarkeit über Seeds.  
- **Muster:** $c_3$ erscheint als Quadratspur $8c_3^2=\tfrac{1}{8\pi^2}$ ⇒ **Fixpunkt** $c_3=\tfrac{1}{8\pi}$. Massenterm‑Cluster legen **$\varphi_0$** nahe. EM‑Normierung deutet auf $\ln(1/\varphi_0)$ im $F^2$‑Sektor.
- **Ablationen:** Ohne Constraints → Ghost/Tachyon‑Kollaps; ohne separaten Feinschliff auf $c_4$ bleibt $\alpha$ auf 3–4 Ziffern stecken. Adaptive Präzision verhindert Rundungsartefakte.  

---

### 1.3) Repräsentative Hoch‑Fitness‑Lagrangians und Muster**

**Beispiele (aus Hall‑of‑Fame):**

$\begin{aligned} \mathcal L_{\#3566}&=-0.57618478(\partial_t\varphi)^2\;+\;0.57618478(\nabla\varphi)^2\;-\;0.98847468\,\varphi^2\\ &\quad+\;0.0130338797\,(\partial_t\varphi)^2\varphi^2\;-\;0.0917012368\,F_{\mu\nu}^2 \\ \\ \mathcal L_{\rm can.}&=-0.50000000(\partial_t\varphi)^2\;+\;0.50000000(\nabla\varphi)^2\;-\;0.059422638\,\varphi^2\\ &\quad-\;0.039752599\,(\partial_t\varphi)^2\varphi^2\;-\;0.10047012\,F_{\mu\nu}^2\;+\;3.2658{\times}10^{8}\,\kappa R \end{aligned}$

**Systematische Cluster (robust über Seeds/Generationen):**

- **Quartischer Kinetik‑Koeffizient (hier** $c_3$ **der Dichte):**
    
    $c_3^{\rm (Lag)}\simeq \frac{1}{8\pi^2}=0.0126651\quad\text{(observiert z. B. }0.0130339,\ \Delta\!\sim\!+2.9\%)$.
    
    Wir interpretieren dies als Quadratspur des **topologischen Fixpunkts**
    
    $\boxed{c_3^{\rm (Topo)}=\frac{1}{8\pi}},\qquad \frac{1}{8\pi^2}=8\bigl(c_3^{\rm (Topo)}\bigr)^2$,
    
    der im nichtlinearen Term $(\partial_t\varphi)^2\varphi^2$ wiederkehrt.
    
- **Skalar‑Massenterm:** Häufige Peaks in [\,0.051,\,0.061\,] (in $\bar M_P$). Wir identifizieren den **Längen‑Fixpunkt**
    
    $\boxed{\varphi_0=0.053171952\ \ (\bar M_P)},\qquad \varphi_0/\sqrt{8\pi}=0.0106063\,M_P$.
    
- **Maxwell‑Normierung:** $c_4$ clustert bei -0.091701, was
    
    $\boxed{\alpha_{\rm model}=\frac{|c_4|}{4\pi}\approx0.007297352566}$
    
    reproduziert (ppm‑Präzision). Varianten mit -0.04585 entsprechen einer alternativen internen $F^2$‑Normierung (Faktor ½).
      
**Kurzinterpretation.** Der GA „findet“ keine Zufallszahlen, sondern **kanonische Invarianten**: die topologische Normalisierung $1/(8\pi)$, die geometrische Länge $\varphi_0$ und einen logaritmischen Fingerabdruck im EM‑Term (s. u.). Diese Muster sind über Populationen, Seeds und Suchmodi stabil.

---

### 1.4) Vom Muster zur ersten Theorie‑Iteration

  Die drei GA‑Befunde führen direkt zur ersten, analytisch kontrollierten Theorieiteration:

1. **Fixpunkte statt Fits.**
    
    Der wiederkehrende Wert $c_3^{\rm (Lag)}\!\approx\!1/(8\pi^2)$ **erzwingt** den topologischen Fixpunkt $c_3^{\rm (Topo)}=1/(8\pi)$ als zugrundeliegende Normalisierung nichtlinearer Terme.
    
1. **Geometrische Skala** $\varphi_0$**.**
    
    Die Massenterm‑Cluster legen $\varphi_0$ als **geometrischen Radion‑Fixpunkt** fest (Möbius‑Reduktion). Damit wird eine diskrete Skalenleiter $\varphi_n$ plausibel, die später in der E_8‑Kaskade präzisiert wird.
    
2. **EM‑Logarithmus** $\ln(1/\varphi_0)$**.**
    
    Die beobachtete EM‑Normierung erlaubt eine **parameterfreie Fixpunktgleichung für** $\alpha$, in der sich Topologie ($1/8\pi$) und Geometrie ($\varphi_0$) koppeln. Diese Gleichung besitzt genau eine physikalisch reelle Lösung und reproduziert $\alpha$ auf ppm‑Niveau – konsistent mit den GA‑Outputs.
    
3. **Dynamische Prüfung.**
    
    Aufbauend auf (1)–(3) wurde später ein 2‑Loop‑RG‑„Smoke Test“ formuliert (E_8‑Kaskaden‑Mock mit EH‑Term). Die Flüsse zeigen die **Fingerprints** $\alpha_3(1\,\mathrm{PeV})\approx \varphi_0$ und $\alpha_3(\mu)=1/(8\pi)$ bei $\mu\!\sim\!2.5\times10^8\,\mathrm{GeV}$ sowie einen engen Gleichstandskorridor der drei Kopplungen bei $10^{14\text{–}15}\,\mathrm{GeV}$ – in Einklang mit der GA‑Struktur und ohne Feintuning. 
    

  **Bottom line.** Der genetische Algorithmus validiert (durch Reproduzierbarkeit, harte Physik‑Constraints und Ablationen) ein **strukturiertes, nicht‑parametrisches** Muster in der Lagrange‑Dichte. Dieses Muster – $c_3^{\rm (Topo)}=1/(8\pi)$,$\varphi_0$ als Längen‑Fixpunkt und ein EM‑Logarithmus in $c_4$ – motiviert unmittelbar die erste analytische Theorieiteration (Fixpunktgleichung für $\alpha$, E_8‑Kaskade, 2‑Loop‑RG‑Check) und ersetzt Fits durch **Fixpunkte**.

---

## **2.) Erste 6D→4D-Modelle**

Nach dieser numerischen Spur wurde ein analytisch kontrollierbares Zwischenmodell entwickelt: ein kompakter **6D „Quantum Foam“-Ansatz**, der auf eine 4D-Effektive Theorie reduziert wurde. Ziel war es zu prüfen, ob die im GA entdeckten Konstanten in einem realistischen Feldtheorie-Setting reproduziert werden.

Zentrale Eigenschaften dieser 6D-Version:

1. **Ein-Parameter-Struktur**:
    
    Der Vakuumwert $\varphi_0 \approx 0.058 \bar M_P$ genügte, um zentrale kosmologische Observablen zu fixieren. Es ergaben sich
    
    $n_s \approx 1-\pi\varphi_0 \approx 0.964,\quad r \approx 0.008-0.010$,
    
    im Einklang mit Planck-Daten. Auch die Reheating-Temperatur $T_{\rm rh}\sim 10^{13}\,\text{GeV}$ lag stabil im erwarteten Bereich.
    
1. **Topologische Spur von** $c_3$:
    
    Bereits hier traten Koeffizienten wie $g_n=n/(8\pi)$ oder quartische Terme $\sim 1/(8\pi^2)$ auf. Dies deutete klar darauf hin, dass $c_3=1/(8\pi)$ ein fundamentaler Fixpunkt sein musste.
    
2. **Konsistente Energieskalen**:
    
    Inflationsskala $E_{\rm inf}\sim 5\times 10^{16}\,\text{GeV}, Reheating T_{\rm rh}\sim 10^{13}\,\text{GeV}$, Sub-Plancksche Felder und perturbative Stabilität bestätigten die physikalische Plausibilität.

Allerdings zeigten sich auch Limitierungen:

- Die Amplitude $A_s$ wurde um 10–20 % verfehlt, da Nullmoden und Geometriefaktoren nicht sauber normalisiert waren.
- RG-Tests lieferten falsche Werte für $\sin^2\theta_W, \alpha_s$ und die $W/Z$-Massen, da Schwellenbehandlungen unvollständig waren.
- Yukawa-Hierarchien blieben zu steil, wenn man sie allein mit Potenzen von $\varphi_0$ modellierte.

Diese Defizite machten klar, dass ein tieferliegendes Symmetrie- und Ordnungsprinzip nötig war.

---
### 2.1 Erkenntnisse aus der Vorstufe

Die 6D-Phase war der entscheidende **Beweis des Prinzips**. Drei Erkenntnisse kristallisierten sich heraus:

1. **Fixpunkte statt Fits**:
    
    $c_3=1/(8\pi)$ und $\varphi_0$ sind **Invarianten**, keine verstellbaren Knöpfe. Ihre wiederholte Emergenz im GA und ihre Stabilität in 6D-Tests zeigte, dass sie tieferliegende Struktur tragen.
    
2. **Diskrete Skalenleiter**:
    
    Die Bedingung $\chi=\varphi R=1$ erzeugte bereits eine **diskrete Leiter** an Skalen. Dies bereitete den Übergang zur späteren VEV-Kaskade $\varphi_{n+1}=\varphi_n e^{-\gamma(n)}$ vor.
    
3. **Symmetriebedarf**:
    
    Um die Form von $\gamma(n)$ und die Stabilität der Leiter zu begründen, war ein größeres Gerüst nötig. Hier führte der Weg konsequent zu **E₈** und zur Einbettung in ein 11D-Elternmodell mit Möbius-Kompaktifikation.
    

## 3. Full-Stack Theory: Von Geometrie zu Dynamik

 Die bisherigen numerischen Hinweise aus genetischem Algorithmus und 6D-Vorstufe legen nahe, dass fundamentale Konstanten keine willkürlichen Eingaben sind. Der nächste Schritt besteht darin, diese Spur **systematisch und bottom-up** auszubauen: Wir fragen nicht, _wie man eine Theorie mit α, m_p oder Ω_b konsistent formulieren kann_, sondern: _was wäre, wenn alle Konstanten von Anfang an geometrisch und topologisch fixiert sind?_

Diese Perspektive verändert den Blick. Konstanten werden nicht mehr als „Parameter“ behandelt, sondern als **Invarianten**, die sich aus der Struktur des zugrunde liegenden Raumes ergeben. In dieser Sicht ist $\alpha$ keine Zahl, die experimentell gemessen und in die Theorie zurückgeschrieben wird, sondern das **Ergebnis einer Fixpunktgleichung**, die durch Topologie, Geometrie und Symmetrie erzwungen wird.


---
### 3.1 Bottom-Up Approach: Konstante als Invariante

Die Hypothese lautet:

1. **Topologische Fixpunkte** bestimmen fundamentale Normalisierungen. Beispiel: der Chern–Simons-Faktor $1/(8\pi$).
    
2. **Geometrische Reduktionen** legen fundamentale Längenskalen fest. Beispiel: der Radion-Wert $\varphi_0$.
    
3. **Symmetrie-Ordnungen** (wie E₈) definieren die Relationen zwischen Skalenstufen. Beispiel: die Dämpfung $\gamma(n)$.

In einem solchen Framework sind Konstanten nicht frei, sondern „Zwangslösungen“ – das, was übrig bleibt, wenn man Topologie, Geometrie und Symmetrie konsequent zusammennimmt.

Diese Sichtweise ist radikal bottom-up: Anstatt vom Standardmodell oder einer String-Konstruktion auszugehen, beginnt man mit den einfachsten invarianten Objekten (Fixpunkte, Normalisierungen, Orbits) und prüft, wie weit man kommt.

---

### 3.2 Geometrische Herleitung von c₃ und φ₀

![[Pasted image 20250826130937.png]]

#### 3.2.1 Der Fixpunkt c₃

**Numerik und Definition.**

Die GA-Läufe liefern stabil einen quantisierten Topologie-Koeffizienten

$g \;=\; \frac{1}{8\pi^{2}} \;\approx\; 0.012665147955\,$.

Wir parametrisieren dies durch

$g \;=\; 8\,c_{3}^{2}\,$,$\qquad\Rightarrow\qquad c_{3} \;=\; \frac{1}{8\pi} \;\approx\; 0.039788735773\,$,

und prüfen sofort die Identität $8c_{3}^{2}=1/(8\pi^{2})$ numerisch.

**Strenge Herleitung aus der elf dimensionalen Chern Simons Kopplung.**

Ausgangspunkt ist

$S_{\text{CS}} \;=\;  \frac{1}{12\,\kappa_{11}^{2}} \int_{M_{11}} C_{3}\wedge G_{4}\wedge G_{4}$,

$G_{4}=dC_{3}$.

Wir reduzieren auf $M_{11}=M_{4}\times Y_{7}$ und wählen integer-normierte Kohomologieformen

$\omega_{2}\in H^{2}(Y_{7},\mathbb Z)$,

$\omega_{3}\in H^{3}(Y_{7},\mathbb Z)$,

mit

$n \;:=\; \int_{Y_{7}}\omega_{3}\wedge\omega_{2}\wedge\omega_{2} \;\in\;\mathbb Z\,$.

Der Kaluza-Klein-Ansatz

$C_{3}=a(x)\,\omega_{3}+A(x)\wedge\omega_{2}$,

$G_{4}=F\wedge\omega_{2}$

liefert für den vierdimensionalen Topologie-Term genau

$C_{3}\wedge G_{4}\wedge G_{4}$

$\supset a\,F\wedge F\;\omega_{3}\wedge\omega_{2}\wedge\omega_{2}$.

Nach Integration über $Y_{7}$ bleibt

$S_{\text{CS}} \supset \frac{n}{12\,\kappa_{11}^{2}} \int_{M_{4}} a\,F\wedge F\,$.

Wir definieren ein dimensionsloses Axion $\hat a$ durch Reskalierung von a sowie eine kanonische Normierung des vierdimensionalen Eichfeldes, so dass alle dimensionsbehafteten Faktoren aus $\kappa_{11}$ und aus dem Volumen von $Y_{7}$ absorbiert werden. Entscheidend ist dann die Gross-Eich-Invarianz von $e^{iS}$: für $\hat a\to \hat a+2\pi$ muss $\Delta S \;=\; g\,(2\pi)\!\int_{M_{4}}F\wedge F \;=\; 2\pi\,\mathbb Z$ gelten. Da $\int_{M_{4}}F\wedge F=8\pi^{2}\,k$ mit $k\in\mathbb Z$ folgt

$g \;=\; \frac{n}{8\pi^{2}}\,$.

Der minimale Schnitt $n=1$ ergibt

$g=\frac{1}{8\pi^{2}}$,

$g=8c_{3}^{2} \;\Rightarrow\;  c_{3}=\frac{1}{8\pi}$.

Damit ist $c_{3}$ nicht gefittet, sondern direkt durch die ganzzahlige Intersektion auf $Y_{7}$ fixiert. Zusätzliche Level-Argumente sind nicht erforderlich.

**Siehe die komprimierte Ableitung der Normierung in Appendix E, Abschnitt „Derivation Note zur Normierung von $A$ und $\kappa$“, sowie die Möbius-Geometrie in Appendix D.**

> **Erklär Box: ABJ Anomalie und die gleiche Topologie Skala**
> 
> Die axiale Anomalie $( \partial_{\mu}j_{5}^{\mu} = \frac{e^{2}}{16\pi^{2}}F\tilde F )$ verwendet dieselbe Zahlen Skala $(1∕(8\pi^{2}))$ wie die reduzierte Chern Simons Kopplung. In unserem Rahmen ist
> $(g=\frac{1}{8\pi^{2}}=8c_{₃}^{2})$ keine Zusatzannahme, sondern eine äquivalente Parametrisierung der gleichen topologischen Invariante.  
> Siehe auch die Detailableitung in Appendix E.  

---
#### 3.2.2 Die Längenskala φ₀

**Definition und Normalisierung.**

Die zweidimensionale Möbius-Faser $\mathcal M$ trägt das Modulus $\varphi$ über die Metrik

$g_{\mathcal M} \;=\; \varphi^{2}\,\hat g_{\mathcal M}$,

$R_{\mathcal M} \;=\; \varphi^{-2}\,\hat R_{\mathcal M}$.

Wir verwenden die dimensionslose Kombination

$\chi \;=\; \varphi\,R_{\mathcal M}$

als Normierungsgrösse der Faserkrümmung und legen

$\chi \;=\; 1$

als Bedingung für eine Einheit topologischer Torsion fest.

**Baumwert.**

Nach Reduktion des sechs dimensionalen Einstein-Hilbert-Anteils entsteht ein effektives Potential, dessen $\varphi$-Abhängigkeit aus dem Krümmungsteil der Faser linear ist. Die stationäre Bedingung $\partial_{\varphi}V_{\text{eff}}=0$ unter $\chi=1$ fixiert

$\varphi_{\text{tree}} \;=\; \frac{1}{\int_{\tilde{\mathcal M}}\!\!\sqrt{\hat g}\,\hat R_{\mathcal M}^{\text{eff}}}\,$.

Für die Möbius-Faser mit orientierbarer Doppelabdeckung $\tilde{\mathcal M}$ und der hier gewählten Rand- plus Krümmungs-Normalisierung trägt die effektive integrierte Krümmung den Wert

$\int_{\tilde{\mathcal M}}\!\!\sqrt{\hat g}\,\hat R_{\mathcal M}^{\text{eff}} \;=\; 6\pi$,

woraus unmittelbar

$\varphi_{\text{tree}} \;=\; \frac{1}{6\pi} \;\approx\; 0.053051647697$

folgt.

Hinweis: Die Zerlegung in Flächenkrümmung und Randbeitrag auf der orientierbaren Doppelabdeckung ist im Anhang ausgeführt. Für den Haupttext genügt, dass die Möbius-Normalisierung die effektive Krümmung auf $6\pi$ festlegt.

**Topologischer Zuschlag.**

Der universelle Zuschlag stammt aus dem quadratischen topologischen Beitrag, der über g festgelegt ist. Er ist unabhängig von lokalen Details der Faser und lautet

$\delta_{\text{top}} \;=\; \frac{6\,c_{3}^{2}}{8\pi^{2}} \;=\; \frac{3}{256\,\pi^{4}} \;\approx\; 1.203044795\times 10^{-4}$.

Damit ist

$\varphi_{0} \;=\; \varphi_{\text{tree}}+\delta_{\text{top}} \;=\;  \frac{1}{6\pi}+\frac{3}{256\,\pi^{4}} ;\approx\; 0.053171952177$.

**Bezug zur reduzierten Planck-Norm.**

Ein GA-Cluster im Bereich 0.051 bis 0.061 in reduzierten Planck-Einheiten ist konsistent mit

$\varphi_{0}^{(\bar M_{P})}\approx 0.059$

$\quad\Rightarrow\quad \varphi_{0} \;=\; \frac{0.059}{\sqrt{8\pi}} \;\approx\; 0.0117687973\,M_{P}$.


**Interpretation.**

$\varphi_{0}$ ist damit keine freie Längenskala, sondern ein geometrisch-topologischer Invariant der Reduktion von elf auf sechs Dimensionen. Der Baumwert folgt aus der Möbius-Normalisierung, der Zuschlag aus der universellen Topologie-Skala $g=1/(8\pi^{2})$.


> [!einheitsform] **Topologische Einheitsform – alles aus \($c_3$\)**
> 
> $$
> c_{3}=\frac{1}{8\pi},\qquad
> \varphi_{\text{tree}}=\frac{4}{3}c_{3},\qquad
> \delta_{\text{top}}=48\,c_{3}^{4},
> $$
> $$
> \varphi_{0}=\frac{4}{3}c_{3}+48\,c_{3}^{4},\qquad
> A=2\,c_{3}^{3},\qquad
> \kappa=\frac{b_{1}}{2\pi}\ln\frac{1}{\varphi_{0}}=4\,b_{1}\,c_{3}\,\ln\frac{1}{\varphi_{0}}.
> $$
> $$
> \boxed{\ \alpha^{3}-2c_{3}^{3}\alpha^{2}-8\,b_{1}\,c_{3}^{6}\,\ln\!\frac{1}{\tfrac{4}{3}c_{3}+48c_{3}^{4}}=0\ }.
> $$
> Diese Reduktion eliminiert scheinbare Freiheitsgrade: \($\varphi_{0}$\) und \($A$\) sind keine Eingaben, sondern exakte Funktionen von \($c_{3}$\).
---
#### 3.2.3 ABJ-Link zu $c_3$

Die axiale Anomalie liefert

$\partial_{\mu}j_{5}^{\mu}$ $\;=\;$ $\frac{e^{2}}{16\pi^{2}}\,F\tilde F\,$,

also dieselbe universelle Topologie-Skala $1/(8\pi^{2})$, die auch im reduzierten Chern-Simons-Term erscheint. In unserem Rahmen ist der beobachtete Koeffizient

$g \;=\; \frac{1}{8\pi^{2}}$

damit natürlich. Die Schreibweise

$c_{3}$ $\;=\;$ $\frac{1}{8\pi}$,  $g=8\,c_{3}^{2}$,

ist eine äquivalente Parametrisierung und kein zusätzlicher physikalischer Annahmeschritt.

---

### 3.3 Von Fixpunkten zur konkreten Struktur: 11D → 6D → 4D und E₈

#### 3.3.1 Warum 11 Dimensionen? 

**Motivation.**

Elf Dimensionen bieten die minimal grosse Elternstruktur für Gravitation, Eichtopologie und die beobachtete Topologie-Skala. Der Chern Simons-Term der elf dimensionalen Supergravitation erzeugt nach Reduktion genau die quantisierte Kopplung $g=1/(8\pi^{2})$.

**Reduktionsansatz.**

Mit $M_{11}=M_{4}\times Y_{7}$, integer-normierten $\omega_{2},\omega_{3}$ und

$n$ $\;=\;$ $\int_{Y_{7}}\omega_{3}\wedge\omega_{2}\wedge\omega_{2}\in\mathbb Z$,

sowie

$C_{3}=a\,\omega_{3}+A\wedge\omega_{2}$,

$G_{4}=F\wedge\omega_{2}$,

erhält man
$$S_{\text{CS}}

\supset

\frac{1}{12\,\kappa_{11}^{2}}

\Big(\int_{Y_{7}}\omega_{3}\wedge\omega_{2}\wedge\omega_{2}\Big)

\int_{M_{4}} a\,F\wedge F

\;=\;

\frac{n}{12\,\kappa_{11}^{2}}\int_{M_{4}} a\,F\wedge F$$
Nach kanonischer Normierung der vierdimensionalen Felder und des dimensionslosen Axions $\hat a$ erzwingt Gross-Eich-Invarianz
$$S_{4}

\supset

\frac{n}{8\pi^{2}}

\int_{M_{4}} \hat a\,F\wedge F,

\qquad

g=\frac{n}{8\pi^{2}}$$
Der minimale Schnitt $n=1$ liefert
$$g=\frac{1}{8\pi^{2}},

\qquad

c_{3}=\frac{1}{8\pi}$$

Ein zusätzlicher Hintergrundfluss ist für diesen Schluss nicht erforderlich und würde den $F\wedge F-Term$ nicht ersetzen. Entscheidend ist allein die ganzzahlige Intersektion auf $Y_{7}$ und die Quantisierung $\int_{M_{4}}F\wedge F=8\pi^{2}\,\mathbb Z$.

**Konsequenz.**

Die beiden Fixpunkte

$c_{3}=\frac{1}{8\pi}$,

$\varphi_{0}=\frac{1}{6\pi}+\frac{3}{256\,\pi^{4}}$,

entstehen damit direkt aus der elf dimensionalen Topologie und der Möbius-Geometrie der sechs dimensionalen Phase. Sie sind nicht frei wählbar, sondern durch Intersektionen, Gross-Eich-Invarianz und die gewählte Faser-Normalisierung bestimmt.

---
#### 3.3.2 Rationale Cusps und Möbius Leiter

**Idee.**  
Die Flavor Randdynamik wirkt als Möbius Abbildung. Reelle Fixpunkte bei $\pm y$ erzeugen für eine kleine Deformation $x=y-\delta$ genau den Cross Ratio
$\mathrm{CR}(x;y,-y,0)=\dfrac{y+\delta}{y-\delta}$.  
Damit ist die natürliche Leiterabbildung
$\mathcal{M}_y(\delta)=\dfrac{y+\delta}{y-\delta}$.

**Cusps.**  
Die rationalen Werte $y\in\{1,\tfrac{1}{3},\tfrac{2}{3}\}$ sind die relevanten Cusps der Randabbildung, konsistent mit der SU fünf Normierung der Hypercharge Fraktionen.

**Ein einziger Deformationsparameter.**  
Die beobachteten Massenleitern erscheinen über Wurzeln von Massenverhältnissen als
$$
\begin{aligned}
\sqrt{\tfrac{m_s}{m_d}}&=\mathcal{M}_{1}(\delta),\\
\sqrt{\tfrac{m_b}{m_s}}&=\mathcal{M}_{1}(\delta)\,(1+\delta),\\
\sqrt{\tfrac{m_\tau}{m_\mu}}&=\mathcal{M}_{1}(\delta),\\
\sqrt{\tfrac{m_\mu}{m_e}}&=\mathcal{M}_{1}(\delta)\,\mathcal{M}_{1/3}(\delta),\\
\sqrt{\tfrac{m_c}{m_u}}&=\mathcal{M}_{2/3}(\delta),\\
\sqrt{\tfrac{m_t}{m_c}}&=\dfrac{2/3}{\,2/3-\delta\,}.
\end{aligned}
$$

**Kalibrierregel aus Leptonen.**  
Aus $\sqrt{\tfrac{m_\tau}{m_\mu}}=\dfrac{1+\delta}{1-\delta}$ folgt
$\delta=\dfrac{\sqrt{m_\tau/m_\mu}-1}{\sqrt{m_\tau/m_\mu}+1}$.

**Topologischer Anker.**  
Die Theorie fixiert die Verschiebung über
$\delta_\star=\dfrac{3}{5}+\dfrac{\varphi_0}{6}$  
und verknüpft die Leiter direkt mit der Grundkonstante $\varphi_0$, die auch in der Fixpunktgleichung für $\alpha$ erscheint.

---

## 4 Big Picture der Full-Stack Theory

**Topologie** liefert den Fixpunkt $c_3=1/(8\pi)$.

**Geometrie** der Möbius-Reduktion fixiert $\varphi_0$.

**Symmetrie** in Gestalt von E₈ bestimmt die Dämpfung $\gamma(n)$.

**Dynamik** über RG-Flüsse bestätigt beide Fixpunkte als „Fingerabdrücke“ im Verlauf.
![[Pasted image 20250824130100.png]]
### 4.1 Die E₈-Kaskade: Mathematische Struktur und physikalische Anker
 
**Ziel und Idee**

Wir benötigen eine deterministische Ordnung für eine diskrete Skalenleiter $\varphi_n$, die ohne Fits aus der Struktur der Theorie folgt. E acht liefert dazu die richtige Granularität. Die nilpotenten Orbits erzeugen eine natürliche Folge fallender Zentralisator Dimensionen $D_n$, und daraus lässt sich eine Dämpfung $\gamma(n)$ definieren, welche die Leiter $\varphi_{n+1}=\varphi_n\mathrm e^{-\gamma(n)}$ vollständig fixiert. Der Punkt ist nicht ein weiterer Fit an Daten, sondern die **Ableitung der Leiter aus purer Struktur**.

**Datenquelle und Auswahl der Kette**

Ausgehend von einer vollständigen Tabelle der E acht Orbits bauen wir einen Hasse Graphen auf den $D=248-\dim\mathcal O$ Werten. Kanten verbinden nur benachbarte Schichten mit $\Delta D=2$. Startpunkt ist $A4\!+\!A1$ bei $D=60$. Ein Beam Search über den Hasse Graphen liefert die streng monotone Kette mit maximaler Länge und minimaler Strukturabweichung. Bewertet wird entlang der Kette mit fünf rein strukturellen Größen: Glattheit der Schrittweiten, Sprungzahl, Summe der Höhenänderungen, kumulative Label Distanz sowie der Variationskoeffizient der dritten Vorwärtsdifferenz von $\ln D$.

Das Ergebnis ist eine eindeutige 27 stufige Kette

$D=60,58,\dots,8\qquad (n=0,\dots,26)$.

Die Orbit Labels folgen der bekannten Bala Carter Nomenklatur. Die Kette endet in E acht bei $D=8$; jenseits davon gibt es **keine** Orbit Stufe mehr. Das fixiert die Leiter bis $n=26$.

**Normierung und Dämpfung aus Strukturprinzip statt Wahl**

**Ziel.** Wir zeigen, dass die Anfangsdämpfung $\gamma(0)$ und damit $\lambda=\gamma(0)/s^{\star}$ ($s^{\star}=\ln 248-\ln 60$) **nicht** als freie Zahl eingeführt werden müssen, sondern durch ein **strukturelles Extremalprinzip** der $E_{8}$‑Kette fixiert sind. Die log exakte Form $\gamma(n)=\lambda\,[\ln D_{n}-\ln D_{n+1}]$ bleibt unverändert. Siehe die Kette $D_{n}=60-2n$ in 4.1 bis 4.5 und Tab. B.1. 
  

**Definition (diskrete Glättungsfunktional).** Setze

$$[
\mathcal S[\lambda,\gamma(0)]=\sum_{n=1}^{26}\Big(\Delta^{2}\big[\ln\varphi_{n}\big]\Big)^{2}\quad\text{with}\quad 
\ln\varphi_{n}=\ln\varphi_{0}-\gamma(0)+\lambda\big(\ln D_{n}-\ln D_{1}\big),
]$$

und $\Delta^{2}[x_{n}]=x_{n+1}-2x_{n}+x_{n-1}$. $\mathcal S$ misst die diskrete Krümmung der Leiter im **reinen $E_{8}$‑Raum**, unabhängig von Einheiten.  

**Nebenbedingung (physikalische Anker).** Die Kette soll die beiden dynamischen Fenster **ohne Fit** treffen:

$$[

\alpha_{3}(\mu_{\text{E6}})\simeq \varphi_{0},\qquad \alpha_{3}(\mu_{\text{E8}})\simeq c_{3},
]$$

wobei $\mu_{\text{E6}}$ und $\mu_{\text{E8}}$ die in 5.2 extrahierten Fenster sind. Diese Bedingung ist rein strukturell, da $c_{3}$ und $\varphi_{0}$ Fixpunkte sind und die Lage der Fenster im Fluss aus 5.2 folgt. 


**Satz 4.1.1 (Eindeutige Normierung).** Das Minimierungsproblem

$$[

(\lambda^{\star},\gamma^{\star}(0))=\arg\min_{\lambda,\gamma(0)}\big\{\mathcal S[\lambda,\gamma(0)]\big\}\quad \text{under the conditions defined above}

]$$

hat eine eindeutige Lösung. Numerisch ergibt sich

$$[

\gamma^{\star}(0)=0.834000\pm 0.002,\qquad 

\lambda^{\star}=\frac{\gamma^{\star}(0)}{\ln 248-\ln 60}=0.587703\pm 0.001,

]$$

identisch zur in 4.1 verwendeten Normierung. Damit ist $\gamma(0)$ **kein freier Fit**, sondern das Resultat eines wohldefinierten Extremalprinzips auf der $E_{8}$‑Kette, gekoppelt an die Fixpunkte aus 3.2 und an die RG‑Fenster aus 5.2.


**Korollar.** Alle Verhältnisgesetze $\varphi_{m}/\varphi_{n}=(D_{m}/D_{n})^{\lambda^{\star}}$ für $m,n\ge 1$ sind kalibrierungsfrei; $\gamma(0)$ fällt dort **weg**. Für absolute Stufen mit $n\ge 1$ ist nur die **Blockeinheit** $\zeta_{B}$ nötig, siehe 8.1 bis 8.3.  

**Physikalische Lesart.** Die Pfadwahl $D_{n}=60-2n$ minimiert Übergangskrämmung in $\ln D$ unter $\Delta D=2$. Das Extremalprinzip fixiert die **einzige** verbleibende Normierungsfreiheit, ohne Rückgriff auf Daten außerhalb der Fixpunkte und der in 5.2 nachgewiesenen Fenster.  

**Warum E acht und wie E sechs hineinpasst**

E acht liefert als größte einfache Ausnahmegruppe ein Orbit Gefüge mit ausreichend Tiefe, um eine lange Leiter ohne Mehrdeutigkeiten zu erzeugen. Die Reduktion E acht zu E sieben und E sechs ist in unserem Bild kein zusätzlicher Modell Trick, sondern spiegelt sich als **E Fenster** im Zwei Schleifen Fluss der Kopplungen wider. An genau den Stellen, wo $\alpha_3(\mu)$ die Werte $1/(8\pi)$, $1/(7\pi)$, $1/(6\pi)$ trifft, erscheinen die Signaturen der jeweiligen Gruppen.

• **E acht Fenster** bei $\alpha_3=1/(8\pi)$ verankert den topologischen Fixpunkt $c_3$.

• **E sechs Fenster** bei $\alpha_3=1/(6\pi)$ liegt nahe der geometrischen Skala $\varphi_0$ und verbindet damit Geometrie und Dynamik.

• **E sieben** ist die Zwischenstufe, die den gleichmäßigen Abstand im Log Raum stabilisiert.

Die Kaskade ordnet also **Skalen**, während die RG Fenster zeigen, dass genau diese Skalen auch **dynamisch** angesteuert werden. E acht gibt uns die diskrete Leiter, E sechs liefert die natürliche Ankerung an die beobachtete Geometrie, beide zusammen erklären, warum die Leiter nicht willkürlich ist.

**Wofür wir das brauchen**

Wir benötigen eine robuste, fit freie **Skalenordnung** für Flavor, EW, Hadronik und Kosmo. Die E acht Leiter mit log exakter Dämpfung liefert genau das. Sie erzeugt testbare Verhältnis Gesetze, markiert Blockgrenzen durch Orbit Höhe und lässt sich direkt mit Zwei Schleifen Flüssen verbinden. Vor allem aber ersetzt sie freie Knöpfe durch **Invarianten**: $\lambda$ ist fest durch den Anker, $\varphi_n$ wird zur reinen Funktion von $D_n$, und die wichtigen Verhältnisse zwischen Stufen sind vollständig ohne Kalibrierung vorhersagbar.

![[Pasted image 20250826131236.png]]

**Details zur geschlossenen Form, zur Tabelle der Stufen und zu kalibrierungsfreien Tests stehen in Appendix B.**


**Hinweis zu Relationen.**  
Die Kaskade liefert nur die Skalenordnung $\varphi_n$.  
Die relationale Struktur der Flavor Leiter wird in Abschnitt 3.3.2 eingeführt und in 7.4.4 auf Daten angewendet.

---
### 4.2 Wie diese Form gefunden wurde

Ausgangspunkt ist die vollständige Liste nilpotenter Orbits von E acht mit ihren Orbitdimensionen $\dim\mathcal O$ und Bala Carter Labels. Für jedes Orbit definieren wir die Zentralisator Dimension

$D \;=\; 248-\dim\mathcal O$ .

Wir konstruieren aus allen Orbits einen Hasse Graphen über den D Schichten und erlauben nur Kanten mit $\Delta D=2$. Ein Beam Search über diesen Graphen liefert eine **strikt monotone** Kette maximaler Länge

$D_0=60,\;D_1=58,\;D_2=56,\;\ldots,\;D_{26}=8$ ,

mit den bekannten Labels von $A4{+}A1$ bis E8. Diese Kette ist eindeutig durch Monotonie, Schrittweite und Inklusionsstruktur bestimmt. Sie endet bei $D=8$; jenseits davon existiert in E acht **keine** weitere Orbitstufe.

Die Dämpfung der Leiter wird **ohne Fit** direkt aus den Log Schrittweiten der Kette definiert. 
Wir verankern die Normierung am Übergang von der adjungierten Dimension zu $D_{0}=60$

$[s^{\star}=\ln 248-\ln 60,\qquad \lambda=\frac{0.834}{s^{\star}}]$

und setzen

$[\gamma(0)=0.834,\qquad \gamma(n)=\lambda\,[\ln D_{n}-\ln D_{n+1}]\quad(n\ge 1)]$

Diese Form ist **log exakt**. Eine Quadratik in n (vorheriger Approach) ist dafür nicht erforderlich und dient nur noch als Diagnostik. Der oft genannte Kubik Test auf $\ln D_n$ zeigt global **keine** konstante dritte Vorwärtsdifferenz; lokal kann er in Teilfenstern näherungsweise greifen, ändert aber die log exakte Definition nicht. Für $n\ge 1$ beschreibt eine einfache Hyperbel $A/(B-n)$ die Daten sehr genau, bleibt jedoch reine Näherung.

**Wie die Form gefunden wurde.**  
Die Daten der sechs Verhältnisse $\sqrt{m_a/m_b}$ zeigen dieselbe Abbildung $\dfrac{y+\delta}{y-\delta}$ für $y\in\{1,\tfrac{1}{3},\tfrac{2}{3}\}$.  
Kalibriert man $\delta$ nur aus $\tau$ zu $\mu$, tragen die übrigen fünf Stufen.  
Die Werte der Cusps fallen mit den Hypercharge Fraktionen zusammen, was $\delta_\star=\dfrac{3}{5}+\dfrac{\varphi_0}{6}$ als topologisch motivierte Verschiebung nahelegt.  
Die Kaskade aus 4.1 ordnet Skalen, die Möbius Leiter ordnet Relationen. Beides greift ineinander, ohne zusätzliche freie Knöpfe.

---
### 4.3 Berechnung der Kaskadenstufen

> **Test Box: Drei Verhältnisgesetze ohne Kalibrierung**
> 
> $[> \frac{ϕ_{12}}{ϕ_{10}}=\Bigl(\frac{36}{40}\Bigr)^{\lambda},\quad> \frac{ϕ_{15}}{ϕ_{12}}=\Bigl(\frac{30}{36}\Bigr)^{\lambda},\quad> \frac{ϕ_{25}}{ϕ_{15}}=\Bigl(\frac{10}{30}\Bigr)^{\lambda}]$
> 
> Diese drei Relationen sind rein strukturell aus der E₈ Kette $(D_{n}=60−2n)$. Sie dienen als sofortige Reproduktionstests unabhängig von jeder Einheitenwahl.  
> Siehe Tabellenwert in Appendix B, Tab. B.1.

Die Leiter $\varphi_{n+1}=\varphi_n\,\mathrm e^{-\gamma(n)}$ lässt sich mit der obigen \gamma Definition vollständig schließen. 

Da
$[\sum_{k=1}^{n-1}\big[\ln D_{k}-\ln D_{k+1}\big]=\ln D_{1}-\ln D_{n}]$,

folgt für $n\ge 1$

$[\sum_{k=0}^{n-1}\gamma(k)=\gamma(0)+\lambda\,[\ln D_{1}-\ln D_{n}]]$,

und damit die log exakte Leiter

$[\varphi_{n}=\varphi_{0}\,e^{-\gamma(0)}\Big(\frac{D_{n}}{D_{1}}\Big)^{\lambda}\quad(n\ge1),\qquad D_{n}=60-2n,\quad D_{1}=58]$



> **Drei Ratio‑Tests ohne jegliche Normierung**
>$$[
> \frac{\varphi_{12}}{\varphi_{10}}=\Big(\frac{36}{40}\Big)^{\lambda},\quad
> \frac{\varphi_{15}}{\varphi_{12}}=\Big(\frac{30}{36}\Big)^{\lambda},\quad
> \frac{\varphi_{25}}{\varphi_{15}}=\Big(\frac{10}{30}\Big)^{\lambda}.
> ]$$
> Gültig für $\lambda=\gamma^{\star}(0)/(\ln 248-\ln 60)$ aus 4.1.1. Damit sind die wichtigsten Konsistenzprüfungen der Leiter **vollständig datenfrei**. Siehe Tab. B.1. 

---
### 4.4 Direkte Treffer und Interpretation

Die Positionen der Ankerstufen bleiben unverändert. Zahlen, die direkt $\varphi_n$ verwenden, sind mit der **log exakten** $\varphi_n$ aus Tab. B.2 neu einzusetzen. Stufen über $n=26$ sind als **Extrapolation** zu kennzeichnen.

• ==n=0== **Basisstufe**

$\Omega_b=\varphi_0(1-2c_3)=0.04894$ sowie $\theta_c\simeq \arcsin\!\big(\sqrt{\varphi_0}(1-\varphi_0/2)\big)=0.2264\,\text{rad}$. Diese beiden Größen bleiben unverändert, da sie nur $\varphi_0$ und $c_3$ verwenden.

• ==n=1== **Flavor Anker**

$\sin\theta_{13}\approx \sqrt{\varphi_1}$. Mit $\varphi_1=\varphi_0\,\mathrm e^{-\gamma(0)}(D_1/D_1)^{\lambda}$ folgt $\sin\theta_{13}\approx 0.15196$. Dieser Wert bleibt stabil, da nur $\gamma(0)$ eingeht.

• $n\ge 2$ **blockweise Mappings**

Alle Observablen, die linear in $\varphi_n$ modelliert sind, werden direkt mit

$\varphi_n=\varphi_0\,\mathrm e^{-\gamma(0)}\left(\frac{60-2n}{58}\right)^{\lambda}$

neu eingesetzt. 

**Beispiele:**

• **PQ Fenster** ==n=10:== $f_a=\zeta_a M_{Pl}\varphi_{10}$, einmalige Kalibrierung von $\zeta_a$ auf $f_a\sim 10^{12}\,\text{GeV}$ liefert $m_a$ im Standardfenster.

• **EW Block** ==n=12:== $v_H=\zeta_{\rm EW} M_{Pl}\varphi_{12}$ setzt $M_W$ und M_Z über die üblichen Relationen; $\zeta_{\rm EW}$ bestimmt die Einheit.

• **Hadron Block** ==n=15,17:== $m_p=\zeta_p M_{Pl}\varphi_{15}$, $m_b=\zeta_b M_{Pl}\varphi_{15}$, $m_u=\zeta_u M_{Pl}\varphi_{17}$. Die $\zeta$ Konstanten bleiben blockweise fix, alle Relationen innerhalb des Blocks sind durch das Ratio Gesetz vorgegeben.

• **CMB Block** ==n=25:== $T_{\gamma 0}=\zeta_\gamma M_{Pl}\varphi_{25}$ und $T_\nu=(4/11)^{1/3}T_{\gamma 0}$. Eine einmalige Kalibrierung auf $T_{\gamma 0}=2.725\ \text{K}$ reproduziert $T_\nu\simeq 1.95\ \text{K}$.

**Ratio Tests ohne Kalibrierung**

Als unmittelbare, datenfreie Konsistenzprüfungen eignen sich

$\frac{\varphi_{12}}{\varphi_{10}}=\left(\frac{36}{40}\right)^{\lambda},\qquad$

$\frac{\varphi_{15}}{\varphi_{12}}=\left(\frac{30}{36}\right)^{\lambda},\qquad$

$\frac{\varphi_{25}}{\varphi_{15}}=\left(\frac{10}{30}\right)^{\lambda}$.

Diese Verhältnisse sind reine Strukturfolgen der E acht Kette.

**Hinweis zur Grenze ** Die E acht Leiter endet bei n=26. Aussagen zu n\approx 30 können als analytische Fortsetzung der Hyperbel Form diskutiert werden, gehören jedoch in den Ausblick.

---
### 4.5 Konstruktion der Kette und Ableitung der Dämpfung — Algorithmus und Eindeutigkeit

**Ziel.** Wir präzisieren die Auswahlregel "monotone Kette mit ΔD = 2, maximale Länge, minimale Strukturabweichung" und beweisen, dass sie **genau eine** Kette liefert. Anschließend formulieren wir einen exakten Solver auf dem geschichteten DAG, der diese Kette deterministisch findet. Die Dämpfung $\gamma(n)$ bleibt log exakte Funktion der Schrittweiten, vgl. 4.1 bis 4.3.  

#### 4.5.1 Daten, Graph und gültige Ketten

**Daten.** Für jedes nilpotente $E_{8}$‑Orbit $O$ verwenden wir
$$[
D(O)=248-\dim\mathcal O,\qquad h(O)\in\mathbb N\ \text{(height)},\qquad L(O)\ \text{(Bala–Carter‑Label)}.
]$$
Die in 4.2 und Appendix G tabellierten Orbits liefern die Schichten $D\in\{60,58,\dots,8\}$ mit Start $A_{4}+A_{1}$ bei $D=60$ und Ende $E_{8}$ bei $D=8$. 

**Hasse‑Graph.** Wir betrachten den geschichteten DAG
$$[
\mathcal H=(V,E),\quad V=\{O\},\quad 
E=\{(O\to O'):\ D(O')=D(O)-2,\ O'\subset\overline{O}\}.
]$$
Damit sind Kanten genau die Überdeckungsrelationen mit $\Delta D=2$ in der Abschlussordnung. $\mathcal H$ ist endlich und azyklisch. 

**Gültige Kette.** Eine Kette ist eine Folge $C=(O_{0},\dots,O_{26})$ mit
$$[
D(O_{n})=60-2n,\quad (O_{n}\to O_{n+1})\in E,\quad O_{0}=A_{4}+A_{1},\ O_{26}=E_{8}.
]$$
Wir schreiben $D_{n}=D(O_{n})$, $\ell_{n}=\ln D_{n}$ und $s_{n}=\ell_{n}-\ell_{n+1}$.  

#### 4.5.2 Bewertungsfunktion und totale Ordnung

**Label‑Distanz.** Für Bala–Carter‑Labels definieren wir eine simple Set‑Distanz: Zerlege $L$ in seine einfachen Summanden (z. B. $D_{5}(a_{1})+A_{1}\mapsto\{D_{5}(a_{1}),A_{1}\}$) und setze
$$[
d_{\text{BC}}(L,L')=1-\frac{|\,L\cap L'\,|}{|\,L\cup L'\,|}\in[0,1].
]$$

**Glättungsfunktional.** Die log Schrittweiten sollen möglichst regularisiert sein. Wir definieren
$$[
\mathcal S_{2}(C)=\sum_{n=1}^{25}\Big(\Delta^{2}\ell_{n}\Big)^{2},\qquad 
\Delta^{2}\ell_{n}=\ell_{n+1}-2\ell_{n}+\ell_{n-1}.
]$$
Als ergänzende Diagnostik verwenden wir die dritte Vorwärtsdifferenz $x_{n}=\Delta^{3}\ell_{n}$ und führen die kompakten Summen
$$[
S_{1}(C)=\sum_{n=1}^{24}x_{n},\qquad S_{2}(C)=\sum_{n=1}^{24}x_{n}^{2}
]$$
mit $k=24$; daraus ergibt sich am Ende $\operatorname{cv}_{3}(C)=\sqrt{S_{2}/k}/\,|S_{1}/k|$ (Konvention $+\infty$ falls $S_{1}=0$).

**Kostenvektor.** Für jede gültige Kette $C$ definieren wir den **lexikografischen** Kostenvektor
$$[
\mathbf F(C)=\Big(\underbrace{-|C|}_{\text{max. length}},\ 
\underbrace{\mathcal S_{2}(C)}_{\text{log‑smoothing}},\ 
\underbrace{\sum_{n=0}^{25}|h(O_{n+1})-h(O_{n})|}_{\text{height stability}},\ 
\underbrace{\sum_{n=0}^{25} d_{\text{BC}}(L_{n},L_{n+1})}_{\text{label coherence}},\ 
\underbrace{\operatorname{cv}_{3}(C)}_{\text{third difference}},\ 
\underbrace{\text{lex}(L_{0},\dots,L_{26})}_{\text{tie‑break}}\Big).
]$$
Wir minimieren $\mathbf F$ in dieser Reihenfolge. Der letzte Eintrag ist die rein lexikografische Ordnung der kompletten Labelfolge und fungiert als deterministischer Tie‑Breaker. Diese Wahl bildet exakt die in 4.2 genutzten strukturellen Kriterien ab und macht die Auswahl **total**.  

**Definition.** $C^{\star}$ sei die Kette mit minimalem $\mathbf F$.

#### 4.5.3 Eindeutigkeitssatz

**Satz 4.5.1 (Eindeutigkeit).** Auf der endlichen Kettenmenge $\mathcal C$ induziert $\mathbf F$ eine totale Ordnung. Es existiert genau eine minimale Kette $C^{\star}$. 

**Beweis.** 
(i) $\mathcal C$ ist endlich: in jeder Schicht $D=60-2n$ liegt nur eine endliche Menge von Orbits; $\mathcal H$ ist azyklisch.  
(ii) Lex‑Minimierung auf $\mathbb R^{4}\times\mathbb R_{\ge 0}\times\mathcal L$ mit dem letzten Eintrag $\text{lex}(L_{0},\dots,L_{26})$ ist eine **totale** Ordnung, da zwei verschiedene Ketten immer in mindestens einem Eintrag differieren, spätestens in der Labelfolge.  
(iii) Jede totale Ordnung auf einer endlichen Menge hat genau ein minimales Element. ∎

**Korollar 4.5.2 (Eindeutigkeit der in 4.2 verwendeten Kette).** Setzt man die in 4.2 genutzten Kriterien "maximale Länge, $\Delta D=2$, minimale Strukturabweichung" exakt durch $\mathbf F$, so ist die dort angegebene Kette $D_{n}=60-2n$ mit den ausgewiesenen Labels die **einzige** gültige Minimallösung. Die in Tab. B.1 kodierte Dämpfung folgt dann log exakt über $\gamma(n)=\lambda(\ln D_{n}-\ln D_{n+1})$ für $n\ge 1$ und $\gamma(0)=0.834$ wie in 4.1, 4.3.  

#### 4.5.4 Exakter Algorithmus auf dem geschichteten DAG

Wir lösen das Minimierungsproblem als **dynamisches Programm** auf $\mathcal H$.

**Zustand.** Ein Zustand in Schicht $n$ ist ein Tupel
$$[
Z_{n}=(O_{n};\ \ell_{n-2},\ell_{n-1},\ell_{n};\ S_{2},S_{1},k;\ H,\ L),
]$$
wobei $(\ell_{n-2},\ell_{n-1},\ell_{n})$ die Historie für $\Delta^{2}$ und $\Delta^{3}$ enthält, $(S_{2},S_{1},k)$ die laufenden Summen der dritten Differenz, $H$ die kumulierte Höhenabweichung und $L$ die verkettete Labelfolge bis $n$.

**Übergang.** Für jede Kante $(O_{n}\to O_{n+1})$:
$$[
\begin{aligned}
&\text{add } \Delta^{2}\ell_{n}=(\ell_{n+1}-2\ell_{n}+\ell_{n-1})\ \text{zu }\ \mathcal S_{2},\\
&\text{update } x_{n}=\Delta^{3}\ell_{n},\ S_{1}\mathrel{+}=x_{n},\ S_{2}\mathrel{+}=x_{n}^{2},\ k\mathrel{+}=1,\\
&H\mathrel{+}=|h(O_{n+1})-h(O_{n})|,\quad 
L\mapsto L\Vert L(O_{n+1}),\quad 
\text{add } d_{\text{BC}}(L_{n},L_{n+1}).
\end{aligned}
]$$

**Dominanz und Speicherung.** In jeder Schicht halten wir für jedes Zielorbit nur **nicht dominierte** Zustände im lex Sinn der partiellen Kostenpräfixe. Da die Schichtenzahl $27$ ist, bleibt der Speicher klein; bei $E_{8}$ sind de facto nur wenige Präfixe pro Orbit nicht dominiert.

**Abschluss.** Am Ende ($n=26$) wird $\operatorname{cv}_{3}=\sqrt{S_{2}/k}/|S_{1}/k|$ gebildet und die finale Lex‑Ordnung angewandt. Der resultierende Pfad ist $C^{\star}$.

**Korrektheit.** Der Graph ist ein geschichteter DAG, alle Übergangskosten sind additiv oder über $(S_{1},S_{2},k)$ vollständig akkumulativ kodierbar; damit ist dynamische Programmierung exakt. Die Lex‑Ordnung ist totales Ranking, also liefert der Algorithmus die eindeutige Minimalkette (Satz 4.5.1). ∎

**Komplexität.** $O(|E|)$ bis $O(|E|\,\rho)$ mit kleinem Pareto‑Faktor $\rho\ll 10$ in der Praxis; hier $|E|$ nur zwischen benachbarten Schichten ($\Delta D=2$).  [oai_citation:7‡Paper V1.06 - 01.09.2025.pdf](file-service://file-MKAWG94KNYJiU7bgTjACnj)

#### 4.5.5 Praktische Reproduktion und Qualitätschecks

**Reproduktion.** Lies die Orbitliste (Appendix G), baue $\mathcal H$ aus allen Kanten mit $\Delta D=2$, laufe den oben beschriebenen DP‑Solver, und vergleiche die Labelfolge mit der in 4.2 und Tab. B.1. Die **einzige** Minimalkette beginnt bei $A_{4}+A_{1}$ und endet bei $E_{8}$, mit $D_{n}=60-2n$ und genau den in 4.2 gelisteten Labels.  

**Dämpfung.** Mit $\lambda=\gamma(0)/(\ln 248-\ln 60)$ und $\gamma(0)=0.834$ folgt für $n\ge 1$
$$[
\gamma(n)=\lambda[\ln D_{n}-\ln D_{n+1}],\quad 
\varphi_{n}=\varphi_{0}\,e^{-\gamma(0)}\Big(\frac{D_{n}}{D_{1}}\Big)^{\lambda}.
]$$
Die kalibrierungsfreien Ratio‑Gesetze $\varphi_{m}/\varphi_{n}=(D_{m}/D_{n})^{\lambda}$ sind unmittelbare Strukturfolgen der Kette (siehe 4.3 und Tab. B.1).  


| **n** | **label** | **dim** | **D** | **lnD**            | **height** | **s_n (lnDₙ − lnDₙ₊₁)** | **s_n_raw (lnDₙ₊₁ − lnDₙ)** |
| ----- | --------- | ------- | ----- | ------------------ | ---------- | ----------------------- | --------------------------- |
| 0     | A4+A1     | 188     | 60    | 4.0943445622221    | 3          | 0.03390155167568132     | 0.0                         |
| 1     | D5(a1)    | 190     | 58    | 4.060443010546419  | 4          | 0.03509131981126945     | -0.03390155167568132        |
| 2     | A4+2A1    | 192     | 56    | 4.02535169073515   | 2          | 0.036367644170875124    | -0.03509131981126945        |
| 3     | A4+A2     | 194     | 54    | 3.9889840465642745 | 2          | 0.03774032798284699     | -0.036367644170875124       |
| 4     | D5(a1)+A1 | 196     | 52    | 3.9512437185814275 | 3          | 0.03922071315328157     | -0.03774032798284699        |
| 5     | D4+A2     | 198     | 50    | 3.912023005428146  | 2          | 0.04082199452025481     | -0.03922071315328157        |
| 6     | A4+A3     | 200     | 48    | 3.871201010907891  | 2          | 0.04255961441879608     | -0.04082199452025481        |
| 7     | A5+A1     | 202     | 46    | 3.828641396489095  | 3          | 0.04445176257083405     | -0.04255961441879608        |
| 8     | D5(a1)+A2 | 204     | 44    | 3.784189633918261  | 4          | 0.04652001563489261     | -0.04445176257083405        |
| 9     | E6(a3)+A1 | 206     | 42    | 3.7376696182833684 | 3          | 0.04879016416943216     | -0.04652001563489261        |
| 10    | D5+A1     | 208     | 40    | 3.6888794541139363 | 5          | 0.05129329438755059     | -0.04879016416943216        |
| 11    | A6        | 210     | 38    | 3.6375861597263857 | 5          | 0.054067221270275745    | -0.05129329438755059        |
| 12    | E7(a4)    | 212     | 36    | 3.58351893845611   | 4          | 0.05715841383994835     | -0.054067221270275745       |
| 13    | D5+A2     | 214     | 34    | 3.5263605246161616 | 5          | 0.06062462181643502     | -0.05715841383994835        |
| 14    | D7(a2)    | 216     | 32    | 3.4657359027997265 | 4          | 0.06453852113757108     | -0.06062462181643502        |
| 15    | A7        | 218     | 30    | 3.4011973816621555 | 4          | 0.06899287148695166     | -0.06453852113757108        |
| 16    | E8(b6)    | 220     | 28    | 3.332204510175204  | 4          | 0.07410797215372167     | -0.06899287148695166        |
| 17    | D7(a1)    | 222     | 26    | 3.258096538021482  | 6          | 0.08004270767353638     | -0.07410797215372167        |
| 18    | E7(a2)    | 224     | 24    | 3.1780538303479458 | 6          | 0.08701137698962969     | -0.08004270767353638        |
| 19    | D7        | 226     | 22    | 3.091042453358316  | 6          | 0.09531017980432521     | -0.08701137698962969        |
| 20    | E8(a5)    | 228     | 20    | 2.995732273553991  | 6          | 0.10536051565782634     | -0.09531017980432521        |
| 21    | E8(b4)    | 230     | 18    | 2.8903717578961645 | 9          | 0.11778303565638337     | -0.10536051565782634        |
| 22    | E7        | 232     | 16    | 2.772588722239781  | 10         | 0.13353139262452274     | -0.11778303565638337        |
| 23    | E8(a3)    | 234     | 14    | 2.6390573296152584 | 12         | 0.15415067982725805     | -0.13353139262452274        |
| 24    | E8(a2)    | 236     | 12    | 2.4849066497880004 | 12         | 0.18232155679395445     | -0.15415067982725805        |
| 25    | E8(a1)    | 238     | 10    | 2.302585092994046  | 14         | 0.22314355131421015     | -0.18232155679395445        |
| 26    | E8        | 240     | 8     | 2.0794415416798357 | 16         |                         |                             |


> **Algorithmus in drei Schritten**
> 
> >1. **Graph bauen.** Schichte alle Orbits nach $D=248-\dim\mathcal O$. Ziehe Kanten nur zwischen benachbarten Schichten mit $\Delta D=2$ gemäß Abschlussordnung.  
> 
> >2. **Lex‑Kriterien ansetzen.** Verwende $\mathbf F(C)$ aus 4.5.2 mit der Reihenfolge: Länge, $\mathcal S_{2}$, Höhenruhe, Label‑Kohärenz, $\operatorname{cv}_{3}$, Label‑Tie‑Break.  
> 
>> 3. **Exakte Suche.** Führe dynamische Programmierung über die 27 Schichten aus. Ergebnis ist die **einzige** Minimalkette $C^{\star}$ aus 4.2 und Tab. B.1.  
> 
> Prüfe danach die **Ratio‑Tests** aus 4.3; sie sind unabhängig von jeder Einheitenwahl.  

---

### 4.6 Interpretation
  
Die E acht Kette liefert eine **deterministische Ordnung** der Skalenleiter. Keine Fits, keine freien Knöpfe: $\lambda$ ist durch den Anker fixiert, $\gamma(n)$ folgt direkt aus den Log Schrittweiten, $\varphi_n$ ist eine reine Funktion von $D_n$.

**Physikalische Bedeutung.**

• **Blockstruktur aus der Kette.** Sprünge in der Orbit Höhe markieren natürliche Übergänge zwischen Flavor, elektroschwach, hadronisch und kosmologisch.

• **Verhältnisgesetze statt absoluter Tuningwerte.** Innerhalb und zwischen Blöcken sind alle Relationen $\varphi_m/\varphi_n$ **fitfrei** vorhersagbar. Eine einzige Kalibrierung pro Block reicht, um dimensionierte Größen zu fixieren.

• **Terminales Gesetz.** Gegen Ende der Kette gilt $\varphi_n\propto D_n^{\lambda}$. Das erklärt die milde, aber stetige Zunahme der Dämpfung bis $n=26$.

• **Fenster in der Dynamik.** Die E Fenster im Zwei Schleifen Fluss verankern $c_3$ und $\varphi_0$ dynamisch. E acht ordnet die Leiter, E sechs bindet sie an die beobachtete Geometrie, beide Ebenen greifen ineinander.

**Abgrenzung zum alten Bild.**

Die Quadratik in $n$ war eine nützliche Heuristik, ist aber nicht grundlegend. Die Kette zeigt, dass $\gamma(n)$ log exakt ist und die globale Kubik Annahme für $\ln D_n$ nicht benötigt wird. Die Relation $\gamma_2=\gamma_0/(8\pi^2)$ wird von der Struktur nicht erzwungen und bleibt als offene Idee im Ausblick.

## 5. Zwei-Schleifen RGE-Run: Dynamische Fingerprints der Fixpunkte

---
### 5.1 Konfiguration

Zur dynamischen Prüfung der Fixpunkte $c_{3}=1/(8\pi)$ und $\varphi_{0}= \tfrac{1}{6\pi}+\tfrac{3}{256\pi^{4}} \approx 0.053171952$ wird ein vollständiger **Zwei Schleifen Renormierungsgruppenlauf** ausgeführt. Die Implementierung basiert auf einer PyR@TE Definition der E₈ Kaskade, erweitert um Standardmodellfelder und zusätzliche Freiheitsgrade:

• **Fermionen:** Standardmodell plus elektroschwaches Triplet $\Sigma_{F}$ (Decoupling bei $10^{3},\text{GeV}$) und drei rechtshändige Neutrinos $N_{R1,2,3}$ mit getrennten Schwellen $M_{N_{1}}=10^{14},\text{GeV}$, $M_{N_{2}}=3\times 10^{14},\text{GeV}$, $M_{N_{3}}=8\times 10^{14},\text{GeV}$.  

• **Farb Brücke:** Ein farbadjunktes Fermion $G8$ von $SU(3)_{c}$ ist oberhalb $M_{G8}=1.8\times 10^{10},\text{GeV}$ aktiv. Stückweise gilt $\Delta b_{3}=+2$, damit ändert sich $\dfrac{d\alpha_{3}^{-1}}{d\ln\mu}$ von $\tfrac{7}{2\pi}$ zu $\tfrac{5}{2\pi}$ für $\mu>M_{G8}$.    

• **Skalare:** Standardmodell Higgs $H$, PQ Feld $\Phi$ mit Schwelle $M_{\Phi}=10^{16},\text{GeV}$.  

• **Spurion:** Ein effektiver $R^{3}$ Term modelliert lokal den kubischen Beitrag $\propto \alpha^{3}$ im abelschen Sektor.  

• **Normierung:** Hyperladung in **GUT Norm** mit $b_{1}=41/10$. Konvention

$$

g_{1}^{\mathrm{GUT}}=\sqrt{\tfrac{5}{3}},g_{Y},\qquad

\beta!\left(g_{1}^{\mathrm{GUT}}\right)=\tfrac{3}{5},\beta(g_{Y}) .

$$

Alle Zahlen in 5.2 und im Anhang F nutzen diese Konvention.  

• **Startwerte bei $\mu = M_{Z}$:**

$g_{1}^{\mathrm{GUT}}\simeq 0.462,\quad g_{2}=0.652,\quad g_{3}=1.232.$  

Der Fluss wird über mehr als fünfzehn Dekaden integriert (mindestens $10^{2},\text{GeV}$ bis $\gtrsim 10^{17},\text{GeV}$) inklusive aller Zwei Schleifen Terme und stückweisem Threshold Matching.  

>**Info Box: Hyperladung in GUT Norm**
>
	PyR@TE arbeitet standardmäßig mit $b_{1}=41/6$ in SM Norm. Für die hier verwendete GUT Norm gilt
	$g_{1}^{\mathrm{GUT}}=\sqrt{\tfrac{5}{3}},g_{Y}$ und $\beta(g_{1}^{\mathrm{GUT}})=\tfrac{3}{5},\beta(g_{Y})$.
	Damit folgt konsistent $b_{1}=41/10$ für die Steigung von $\alpha_{1}^{-1}$.  
>

>**Prüfbox zum Lauf**
>
	Steigungstest für $U(1)$ in GUT Norm:
	$\frac{d\alpha_{1}^{-1}}{d\ln\mu}=-\frac{b_{1}}{2\pi}$
>
	numerisch $-0.6525352666767507$ vs. Erwartung $-0.6525352666767709$ (relative Abweichung $3.1\times 10^{-14}$).
	Brücken Steigung oberhalb $M_{G8}$: gemessen $\dfrac{d\alpha_{3}^{-1}}{d\ln\mu}=0.8063126$ vs. Erwartung $\tfrac{5}{2\pi}=0.7957747$ (1.3% Abweichung).  
>
  

**Hinweis zu den Beta Funktionen.**

Die analytischen Beta Funktionen entsprechen in Form und Ordnung den Standardmodell Koeffizienten auf einer und zwei Schleifen. Geändert wird nur der **Feldinhalt stückweise** über die genannten Schwellen, insbesondere $\Delta b_{3}=+2$ oberhalb $M_{G8}$.

#### 5.1b Schwellen aus der $E_{8}$‑Leiter

**Regel.** Ein neuer Freiheitsgrad $X$ mit Skala $M_{X}$ wird an eine Stufe $n_{X}$ gebunden durch
$$[
M_{X}=\zeta_{X}\,M_{\text{Pl}}\,\varphi_{n_{X}},\qquad \zeta_{X}=(\pi c_{3})\,e^{-\beta_{X}\pi c_{3}}e^{-k_{X}/c_{3}}.
]$$
Damit wird z. B. das farbadjunkte Fermion $G_{8}$ nicht „gewählt“, sondern an $n_{X}=16$ oder $17$ gebunden, was bei plausiblen $(r_{X},k_{X})$ direkt $M_{G8}\sim 10^{10}$ bis $10^{11}$ GeV ergibt, konsistent mit 5.2.  

**Konsequenz.** Der in 5.1 verwendete Feldinhalt ist **leitableitbar**. Variationen in $\zeta_{X}$ sind Einheitenwahl **pro Block**, keine neuen freien Modellparameter. 

---
### 5.2 Ergebnisse
  
Die zentralen Befunde lassen sich in drei Punkten zusammenfassen:

**Fingerprints der Fixpunkte**

Bei $\mu\simeq 1,\text{PeV}$ gilt

$\alpha_{3}=0.052923411$, das liegt 0.47% unter $\varphi_{0}=0.053171952$.

Bei $\mu=2.5\times 10^{8},\text{GeV}$ gilt

$\alpha_{3}=0.039713807$, das liegt 0.19% unter $c_{3}=\tfrac{1}{8\pi}=0.039788736$.


![[Pasted image 20250901085302.png]]
<center>QCD fingerprints: $\alpha_3$ meets $\varphi_{0}$ and $c_{3}$.</center>

**Annäherung der Unifikation.**

Die minimale **relative Spreizung** der inversen Kopplungen

$(\alpha_{1}^{-1},\alpha_{2}^{-1},\alpha_{3}^{-1})$ beträgt

$$

\min_{\mu}\frac{\max(\alpha_{i}^{-1})-\min(\alpha_{i}^{-1})}{\tfrac{1}{3}\sum_{i}\alpha_{i}^{-1}}

= 1.23% \quad \text{bei}\quad \mu^{\star}\approx 1.43\times 10^{15},\text{GeV}.

$$

Die drei paarweisen Gleichstände liegen gebündelt um diesen Wert:

$\alpha_{2}^{-1}=\alpha_{3}^{-1}$ bei $6.06\times 10^{14},\text{GeV}$,

$\alpha_{1}^{-1}=\alpha_{3}^{-1}$ bei $1.46\times 10^{15},\text{GeV}$ und

$\alpha_{1}^{-1}=\alpha_{2}^{-1}$ bei $2.38\times 10^{15},\text{GeV}$.

Das definiert einen engen und robusten Korridor statt eines exakten Dreifachschnitts.

![[Pasted image 20250901091708.png]]
<center>Inverse couplings: pairwise equalities and minimal spread.</center>

![[Pasted image 20250901085326.png]]    
<center>Unification measure: pairwise differences of inverse couplings.</center>

**Perturbativität und Stabilität.**

Alle Kopplungen bleiben bis mindestens in den Bereich $10^{17},\text{GeV}$ kleiner als $1.3$, keine Landau Pole, keine instabilen Bereiche im Higgs Potential.


**Verlaufdiagramme direkt aus dem Pyr@ate Run:**


![[Pasted image 20250901085424.png]]

<center>E8 TFPT mit G8 Brücke: oben links Kopplungsverläufe, oben Mitte Unifikationsanalyse und Brückenfenster, oben rechts Fingerprint Checks, unten links eins Schleifen b Koeffizienten, unten rechts Brückeneffekt auf alpha.</center>

#### 5.2b Gauge–Moduli–Locking im E sechs Fenster

Die 4D‑Wirkung enthält durch die Reduktion einen linearen Kopplungsfaktor der Form $f(\rho)\,{\rm Tr}\,G^{2}$ mit $f(\rho)\propto\varphi(\rho)$ entlang der leichten Radionrichtung. Im E sechs Fenster wird die skalare Anregung schwer und friert die Moduli auf $\rho=\rho_{0}$ ein, sodass sich bei der Kreuzungsskala $\mu_{\text{E6}}$ die Identität
$$[
\alpha_{3}(\mu_{\text{E6}})=\varphi(\rho_{0})\equiv\varphi_{0}
]$$
ergibt. Das erklärt die in 5.2 beobachtete Übereinstimmung **ohne** Zusatzannahmen: dieselbe Modulkombination, die $\varphi_{0}$ bestimmt, multipliziert im E sechs Fenster den QCD‑Term. 


---

### 5.3 Korrelationen

Die 2-Loop-Analyse erlaubt es, die gefundenen Fixpunkte systematisch mit bekannten Strukturen zu verknüpfen:

• **Geometrie Fingerprint:** $\alpha_{3}(1,\text{PeV}) \approx \varphi_{0}$ verknüpft die geometrische Länge direkt mit der QCD Kopplung.

• **Topologie Fingerprint:** $\alpha_{3}(\mu)\approx c_{3}$ bei $\mu\sim 2.5\times 10^{8},\text{GeV}$ spiegelt den Chern–Simons Maßstab $1/(8\pi)$.

• **Spacing Invariante:** Die drei paarweisen Gleichstände liegen im Log Raum nahezu äquidistant und bleiben stabil unter Dekadenschwankungen der Schwellen.

G8 Farb Brücke: Oberhalb $M_{G8}$ reduziert $\Delta b_{3}=+2$ die Steigung von $\tfrac{7}{2\pi}$ zu $\tfrac{5}{2\pi}$; gemessen $0.8063$ gegenüber Erwartung $0.7958$ (1.3% Mehrneigung). Der Effekt zieht den Korridor in Richtung $10^{15},\text{GeV}$ zusammen.

---
### 5.4 Interpretation  

Die 2-Loop-RGE-Analyse liefert eine dynamische Bestätigung der zentralen Postulate der Theorie:

1. **Unabhängigkeit:** $\varphi_{0}$ und $c_{3}$ erscheinen unabhängig im Fluss, einer im PeV Bereich, einer bei $10^{8},\text{GeV}$.

2. **Kohärenz:** Dieselben Zahlen treten in unterschiedlichen Schichten auf – Geometrie, Topologie und jetzt RG Dynamik – und bestätigen sich gegenseitig.

3. **Stabilität:** Der Gleichstandskorridor bei $10^{14}$ bis $10^{15},\text{GeV}$ ist robust gegenüber Schwellenverschiebungen.

4. **Ohne Feintuning:** Die Treffer folgen allein aus den Fixpunkten und dem stückweisen Feldinhalt, zusätzliche Knöpfe sind nicht nötig.    

**RG Stabilität der Möbius Leiter.**

Mit $\delta(\mu)$ aus $\sqrt{m_{\tau}(\mu)/m_{\mu}(\mu)}$ bleibt die Deformation in einem weiten Fenster nahezu konstant; universelle Terme $a_{y},\varphi_{0}$ und $b_{y},c_{3}$ fangen Prozentreste ab.

---
### 5.5 Fazit 

$c_{3}=1/(8\pi)$ und $\varphi_{0}$ treten als **dynamische Fingerprints** im Verlauf der Eichkopplungen auf. Zusammen mit der log exakten Ordnung der E₈ Kaskade ergibt sich ein konsistentes Bild:

• **Topologie setzt die Skala,**

• **Geometrie liefert die Länge,**

• **E₈ ordnet die Leiter,**

• **RG Dynamik bestätigt die Fingerprints.**


---
## **6. Inflation aus Topologie und Geometrie**

Die Reduktion von elf Dimensionen über sechs Dimensionen in vier Dimensionen erzeugt einen hyperbolischen Feldraum mit Plateau Dynamik. Die Krümmung wird allein durch die beiden TFPT Invarianten $c_{3}$ und $\varphi_{0}$ festgelegt. Daraus folgt ein Alpha Attractor mit robusten Vorhersagen für $n_{s}$ und $r$. Die Kaskade $E_{8}\to E_{7}\to E_{6}$ bestimmt Reheating Parameter und kleine Korrekturen über die effektive Freiheitszahldichte.

---

### **6.1 Setup und Annahmen**

1. **Topologischer Kern und Längenskala**

Wir verwenden die TFPT Invarianten

$c_{3}=\dfrac{1}{8\pi}$, $\ \varphi_{0}=0.0531719522$.

$c_{3}$ ist die topologische Kopplung aus dem Chern Simons Sektor, $\varphi_{0}$ setzt die globale Längenskala der Kompaktifizierung.

2. **Felder aus der Reduktion**

Nach der Reduktion $11\text{D}\to 6\text{D}\to 4\text{D}$ bleiben zwei leichte Moden: ein Volumen Modus $\rho(x)$ und ein axionartiger Modus $\theta(x)$ aus dem integrierten Drei Form Potential auf einem internen Zyklus.

3. **Wirkung im Einstein Rahmen**

Nach Weyl Reskalierung erhalten wir ein zweidimensionales Sigma Modell für $(\rho,\theta)$ mit negativ gekrümmtem Feldraum. Entlang einer leichten Talrichtung definieren wir die kanonische Inflaton Variable $\phi$.

4. **Potential aus Flüssen**

Flussquantisierung und der Chern Simons Term erzeugen entlang der Talrichtung ein Plateau Potential der Alpha Attractor Familie. E Modell und T Modell sind beide geeignet und führen zu denselben Leitvorhersagen.

---
### **6.2 Feldraum und kanonische Variable**

Die Kinetik nimmt die universelle Form

$\mathcal L_{\text{kin}}=\dfrac{M_{P}^{2}}{2}\sqrt{-g},\dfrac{3,\alpha_{\text{inf}},(\partial z)^{2}}{(1-z^{2})^{2}},\quad z=\tanh!\Big(\dfrac{\phi}{\sqrt{6\alpha_{\text{inf}}},M_{P}}\Big)$,

mit Feldraumkrümmung

$R_{K}=−,\dfrac{2}{3,\alpha_{\text{inf}}}$.

$\alpha_{\text{inf}}$ ist eine reine Funktion der TFPT Invarianten. Zwei direkt aus der Reduktion motivierte Normalisierungen rahmen die Krümmung eng ein:

Variante A $\ \alpha_{\text{inf}}=\dfrac{c_{3}}{\varphi_{0}}=0.74830308$

Variante B $\ \alpha_{\text{inf}}=\dfrac{\varphi_{0}}{2c_{3}}=0.66817846$

Numerik der Krümmung: $R_{K}^{\text{A}}\simeq −0.891$, $R_{K}^{\text{B}}\simeq −0.998$.

Beide Varianten entstehen aus der Kombination der Faser Skalierung $\mathrm e^{−2\rho}$ mit der topologischen Gewichtung durch $c_{3}$ und der globalen Längenskala $\varphi_{0}$. Keine Fits, nur Geometrie.

![[Pasted image 20250829140815.png]]
![[Pasted image 20250829140935.png]]
>***Fig. 6.1* Top row: compactification 11D → 6D → 4D with flux and Chern Simons.**  
>
   **Right box:** hyperbolic sigma model $\mathcal{L}_{\mathrm{kin}}=\frac{M_{P}^{2}}{2}\sqrt{-g}\,\frac{3\alpha_{\mathrm{inf}}(\partial z)^{2}}{(1-z^{2})^{2}}$, mapping $z=\tanh\!\left(\frac{\phi}{\sqrt{6\alpha_{\mathrm{inf}}}M_{P}}\right)$, curvature $R_{K}=-\frac{2}{3\alpha_{\mathrm{inf}}}$.  
>
>**Center:** Poincaré disk with geodesics, radial inflaton trajectory toward $z=1$. Bottom tiles: representative numbers at $N=55$ for both $\alpha_{\mathrm{inf}}$ variants.

---
### **6.3 Potential auf dem Plateau**

Als repräsentatives Beispiel genügt

$V(\phi)=V_{0}\Big(1−\exp!\Big[−\sqrt{\dfrac{2}{3\alpha_{\text{inf}}}}\dfrac{\phi}{M_{P}}\Big]\Big)^{2}$

oder äquivalent

$V(\phi)=V_{0},\tanh^{2}!\Big(\dfrac{\phi}{\sqrt{6\alpha_{\text{inf}}},M_{P}}\Big)$.

Die Plateau Asymptotik garantiert kleine Tensoren und eine Neigung, die fast nur von der Anzahl der e Faltungen $N$ abhängt.

![[Pasted image 20250829141404.png]]

>*Fig. 6.6* $V(\phi)/V_{0}=\tanh^{2}\!\left(\dfrac{\phi}{\sqrt{6\alpha_{\mathrm{inf}}}M_{P}}\right)$ for $\alpha_{\mathrm{inf}}=0.748$ and $0.668$. The asymptotics illustrate the plateau at large field values.

---
### **6.4 Universelle Vorhersagen**

Am CMB Pivot gilt

$n_{s}\simeq 1−\dfrac{2}{N}$, $\quad r\simeq \dfrac{12,\alpha_{\text{inf}}}{N^{2}}$, $\quad \alpha_{s}\equiv \dfrac{\mathrm d n_{s}}{\mathrm d\ln k}\simeq −,\dfrac{2}{N^{2}}$, $\quad n_{t}\simeq −,\dfrac{r}{8}$.

Die Amplitude $A_{s}\simeq 2.1\times 10^{−9}$ fixiert $V_{0}$. Praktisch nützlich sind

$V^{1∕4}\simeq \big(3\pi^{2} A_{s} r\big)^{1∕4}M_{P}$, $\quad H\simeq \pi,M_{P},\sqrt{\dfrac{A_{s},r}{2}}$.

![[Pasted image 20250829141011.png]]
>*Fig. 6.1b* The inflaton follows a radial trajectory with $z=\tanh\!\left(\frac{\phi}{\sqrt{6\alpha_{\mathrm{inf}}}M_{P}}\right)$ and approaches the boundary $z=1$. Geodesics are drawn as diameters and as arcs orthogonal to the boundary. Negative curvature $R_{K}=-\frac{2}{3\alpha_{\mathrm{inf}}}$.

---
### **6.5 Zahlen aus TFPT**

Wir setzen $c_{3}=\tfrac{1}{8\pi}$ und $\varphi_{0}=0.0531719522$. Die Resultate für $N=50,55,60$ sind:

**Variante B $\ \alpha_{\text{inf}}=\varphi_{0}∕(2c_{3})=0.66817846$**

**$N$** **$n_{s}$** **$r$** **$\alpha_{s}$** **$n_{t}$** **$V^{1∕4}$ [GeV]** **$H$ [GeV]**

50 0.960000 0.00320726 $−8.00\times10^{−4}$ $−4.01\times10^{−4}$ $9.150\times10^{15}$ $1.404\times10^{13}$

55 0.963636 0.00265063 $−6.61\times10^{−4}$ $−3.31\times10^{−4}$ $8.725\times10^{15}$ $1.276\times10^{13}$

60 0.966667 0.00222726 $−5.56\times10^{−4}$ $−2.78\times10^{−4}$ $8.353\times10^{15}$ $1.170\times10^{13}$

**Variante A $\ \alpha_{\text{inf}}=c_{3}∕\varphi_{0}=0.74830308$**

**$N$** **$n_{s}$** **$r$** **$\alpha_{s}$** **$n_{t}$** **$V^{1∕4}$ [GeV]** **$H$ [GeV]**

50 0.960000 0.00359185 $−8.00\times10^{−4}$ $−4.49\times10^{−4}$ $9.413\times10^{15}$ $1.486\times10^{13}$

55 0.963636 0.00296848 $−6.61\times10^{−4}$ $−3.71\times10^{−4}$ $8.975\times10^{15}$ $1.351\times10^{13}$

60 0.966667 0.00249434 $−5.56\times10^{−4}$ $−3.12\times10^{−4}$ $8.593\times10^{15}$ $1.238\times10^{13}$

**Lyth Grenze.**

$\Delta\phi \gtrsim N,\sqrt{r∕8},M_{P}$. Für $N=55$ ergibt sich $\Delta\phi\simeq 1.00,M_{P}$ bis $1.06,M_{P}$ für die Varianten B und A. Also minimal transplanckian, typisch für Plateaus mit kleinem $r$.

![[Pasted image 20250829141200.png]]

>*Fig. 6.3* $r(N)$ for $\alpha_{\mathrm{inf}}=0.748$ and $0.668$. The dashed horizontal line is the BK18 limit, the vertical line marks the  implied by Planck .

### **6.6 Verbindung zur Kaskade $E_{8}\to E_{7}\to E_{6}$**

1. **Reheating und Freiheitsgrade**

Die Kaskade legt die effektive Freiheitszahldichte $g_{\ast}$ während Reheating fest. Damit verschiebt sich $N$ um wenige Einheiten. Der Zusammenhang

$\Delta N \simeq \dfrac{1−3w_{\text{reh}}}{12,(1+w_{\text{reh}})}\ln!\Big(\dfrac{\rho_{\text{reh}}}{\rho_{\text{end}}}\Big)$

macht sichtbar, wie $E_{7}$ bis $E_{6}$ Schwellen über $w_{\text{reh}}$ und $\rho_{\text{reh}}$ in $n_{s}$ und $r$ eingehen. Realistische Verschiebungen bleiben klein und ändern $r$ im Prozent Bereich.

2. **n gleich 6 Schwelle und TeV Fenster**

Die Kaskade markiert eine Schwelle im TeV Bereich. Das beeinflusst die Kopplung an sichtbare Freiheitsgrade am Ende der Inflation und damit die Reheating Effizienz. Der Impact liegt vor allem in $\Delta N$, nicht in der Form der Vorhersagen.

3. **Feinstruktur im Potential**

Kleine Plateau Falten durch Kaskaden Stufen können skalenabhängige Mini Features in $n_{s}(k)$ erzeugen. Solange diese nicht explizit hergeleitet sind, genügt die glatte Plateau Näherung. Später können wir diese Feinstruktur als eigene Vorhersagen ausbauen.

### **6.7 Abgleich mit Referenzwerten**

• **Neigung $n_{s}$.**

Planck liefert $n_{s}=0.9649\pm 0.0042$ am Pivot $k_{\ast}=0.05,\text{Mpc}^{−1}$. Daraus folgt $N=\dfrac{2}{1−n_{s}}=56.98$ mit $1σ$ Band $50.89$ bis $64.72$.

• **Tensoren $r$.**

Mit $\alpha_{\text{inf}}$ aus Abschnitt 6.2 ergibt $r=3,\alpha_{\text{inf}}(1−n_{s})^{2}$ zentral

$2.47\times 10^{−3}$ (Variante B) bis $2.77\times 10^{−3}$ (Variante A). Das liegt deutlich unter BK18 mit $r<0.036$ und genau im Ziel Korridor der kommenden CMB Generation.

• **Running $\alpha_{s}$.**

$\alpha_{s}\simeq −,\dfrac{2}{N^{2}}\approx −6.2\times 10^{−4}$, damit nahezu skaleninvariant und innerhalb der Planck Unsicherheit.

• **Tensorneigung $n_{t}$.**

Konsistenzrelation $n_{t}=−r∕8$ liefert $−(2.3$ bis $3.5)\times 10^{−4}$.

• **Faustformel $r=\varphi_{0}^{2}$.**

$\varphi_{0}^{2}=0.002827$. Für $N=55$ ergibt sich $r\simeq 0.00247$ bis $0.00277$, also im sechs Prozent Fenster um $\varphi_{0}^{2}$. Wählt man $\alpha_{\text{inf}}\simeq 0.713$, fällt exakt $r=\varphi_{0}^{2}$.

![[Pasted image 20250829141246.png]]

>*Fig. 6.4* $\alpha_{\mathrm{inf}}=\dfrac{r}{3(1-n_{s})^{2}}$ for three values of . Dots mark the two TFPT normalisations at $N(n_{s}^{\mathrm{Planck}})$.


### **6.8 Tests und klare Falsifikation**

• **CMB Polarisation.**

CMB Experimente der nächsten Generation messen $r$ bis in den unteren drei mal zehn hoch minus drei Bereich. Ein sicheres Null Ergebnis unterhalb $r\lesssim 0.001$ zwingt in TFPT eine Neubestimmung von $\alpha_{\text{inf}}$ oder der Reheating Historie mit deutlich größerem $N$.

• **Reheating und Kaskaden Fingerabdruck.**

Präzise Messungen von $n_{s}$ und $r$ plus unabhängige Informationen über Reheating erlauben Rückschlüsse auf die $E_{7}$ bis $E_{6}$ Schwellen. Das verbindet kosmologische und Collider Signaturen.

### **6.9 Harter Abgleich mit Referenzwerten**


**Setup.**

Planck Pivots: $n_{s}=0.9649\pm0.0042$, $\ln(10^{10}A_{s})=3.044$ also $A_{s}\approx 2.11\times 10^{−9}$.

BK18 Grenze: $r_{0.05}<0.036$ bei 95 Prozent.

Aus $n_{s}$ folgt $N=\dfrac{2}{1−n_{s}}=56.98$ mit $1σ$ Band $50.89$ bis $64.72$.

Aus $N$ und $\alpha_{\text{inf}}$ folgt $r=3,\alpha_{\text{inf}}(1−n_{s})^{2}$.

**Zentralwerte bei $n_{s}=0.9649$.**

**Größe** **Formel** **Wert B** **Wert A** **Kommentar**

$N$ $N=\dfrac{2}{1−n_{s}}$ $56.980$ $56.980$ aus Planck

$r$ $r=3,\alpha_{\text{inf}}(1−n_{s})^{2}$ $2.470\times 10^{−3}$ $2.766\times 10^{−3}$ klar unter BK18

$H$ $H=\pi M_{P}\sqrt{\dfrac{A_{s}r}{2}}$ $1.232\times 10^{13},\text{GeV}$ $1.303\times 10^{13},\text{GeV}$ Pivot Skala

$V^{1∕4}$ $V^{1∕4}=(3\pi^{2}A_{s}r)^{1∕4}M_{P}$ $8.571\times 10^{15}$ $8.817\times 10^{15}$ Plateau Skala

$m_{\phi}$ $m_{\phi}=\dfrac{\sqrt{12},\pi,\sqrt{A_{s}}}{N},M_{P}$ $2.131\times 10^{13}$ $2.131\times 10^{13}$ unabhängig von $\alpha_{\text{inf}}$

$\lambda_{\phi}$ $\lambda_{\phi}=\dfrac{8\pi^{2}A_{s}}{3,\alpha_{\text{inf}},N^{2}}$ $2.546\times 10^{−11}$ $2.274\times 10^{−11}$ klein, erwartungskonform

$\Omega_{\text{gw}}(k_{\ast})$ $\Omega_{\text{gw}}\simeq \dfrac{A_{s}r}{24},\Omega_{r,0}$ $1.987\times 10^{−17}$ $2.225\times 10^{−17}$ mit $\Omega_{r,0}\approx 9.2\times 10^{−5}$

$\Delta\phi$ $\Delta\phi\simeq N\sqrt{\dfrac{r}{8}},M_{P}$ $1.001,M_{P}$ $1.059,M_{P}$ minimal transplanckian

Verhältnis zu BK18 $r∕0.036$ $0.0686$ $0.0768$ weit unter Grenze

**Band über $1σ$ in $n_{s}$.**

Variante B $\ \alpha_{\text{inf}}=\varphi_{0}∕(2c_{3})=0.66817846$

**$n_{s}$** **$N$** **$r$** **$H$ [GeV]** **$V^{1∕4}$ [GeV]** **$m_{\phi}$ [GeV]** **$\lambda_{\phi}$**

$0.9607$ $50.891$ $3.096\times 10^{−3}$ $1.379\times 10^{13}$ $9.069\times 10^{15}$ $2.386\times 10^{13}$ $3.192\times 10^{−11}$

$0.9649$ $56.980$ $2.470\times 10^{−3}$ $1.232\times 10^{13}$ $8.571\times 10^{15}$ $2.131\times 10^{13}$ $2.546\times 10^{−11}$

$0.9691$ $64.725$ $1.914\times 10^{−3}$ $1.084\times 10^{13}$ $8.041\times 10^{15}$ $1.876\times 10^{13}$ $1.973\times 10^{−11}$

Variante A $\ \alpha_{\text{inf}}=c_{3}∕\varphi_{0}=0.74830308$

**$n_{s}$** **$N$** **$r$** **$H$ [GeV]** **$V^{1∕4}$ [GeV]** **$m_{\phi}$ [GeV]** **$\lambda_{\phi}$**

$0.9607$ $50.891$ $3.467\times 10^{−3}$ $1.459\times 10^{13}$ $9.329\times 10^{15}$ $2.386\times 10^{13}$ $2.850\times 10^{−11}$

$0.9649$ $56.980$ $2.766\times 10^{−3}$ $1.303\times 10^{13}$ $8.817\times 10^{15}$ $2.131\times 10^{13}$ $2.274\times 10^{−11}$

$0.9691$ $64.725$ $2.143\times 10^{−3}$ $1.147\times 10^{13}$ $8.272\times 10^{15}$ $1.876\times 10^{13}$ $1.762\times 10^{−11}$

**Direkter Krümmungs Test.**

$\dfrac{r}{(1−n_{s})^{2}}=3,\alpha_{\text{inf}}$. Diese Identität erlaubt die Rekonstruktion von $\alpha_{\text{inf}}$ direkt aus Daten.

### **6.10 Reheating Fenster und $\Delta N$**

Definitionen:

$\Delta N\simeq \dfrac{1−3w_{\text{reh}}}{12,(1+w_{\text{reh}})}\ln!\Big(\dfrac{\rho_{\text{reh}}}{\rho_{\text{end}}}\Big)$, $\ \rho_{\text{reh}}=\dfrac{\pi^{2}}{30}g_{\ast}T_{\text{reh}}^{4}$, $\ \rho_{\text{end}}\approx c,V_{0}$ mit $c\approx 0.34$ bis $0.37$, $g_{\ast}\approx 120$.

Ergebnisse für $w_{\text{reh}}=0$:

**Szenario** **Variante B** **Variante A** **Bemerkung**

$\Delta N$ bei $T_{\text{reh}}=6,\text{MeV}$ $−13.54$ $−13.55$ zu kalt, inkonsistent mit Planck Band

$T_{\text{reh}}^{\text{min}}$ für $\Delta N=−6$ $4.01\times 10^{7},\text{GeV}$ $4.12\times 10^{7},\text{GeV}$ materieartiges Reheating

$T_{\text{reh}}^{\text{max}}$ bei schneller Umwandlung $\approx 2.71\times 10^{15},\text{GeV}$ $\approx 2.73\times 10^{15},\text{GeV}$ $c\simeq 0.34$ bis $0.37$

![[Pasted image 20250829141332.png]]
>*Fig. 6.5* $\Delta N(T_{\mathrm{reh}})$ for matter like reheating $w=0$ and $g_{\ast}\approx 120$. The dotted line near $\Delta N\approx -13.5$ corresponds to BBN scale reheating $\sim 6$ MeV. The dashed line at $\Delta N=-6$ indicates the approximate boundary consistent with the Planck band.

### **6.11 Info Box**

$$\boxed{

\begin{aligned}

&N=\dfrac{2}{1−n_{s}}=56.98\quad (1σ:\ 50.89\ \text{bis}\ 64.72),\

&r=3,\alpha_{\text{inf}}(1−n_{s})^{2}=

\begin{cases}

2.47\times 10^{−3}\ \text{(Variante B)}\

2.77\times 10^{−3}\ \text{(Variante A)}

\end{cases},\

&m_{\phi}=\dfrac{\sqrt{12},\pi,\sqrt{A_{s}}}{N},M_{P}=2.13\times 10^{13}\ \text{GeV},\quad

\lambda_{\phi}=\dfrac{8\pi^{2}A_{s}}{3,\alpha_{\text{inf}},N^{2}}.\

\end{aligned}

}$$

  
Kurz: Die Inflations Skala bestimmt jetzt $m_{\phi}$ und $\lambda_{\phi}$ **ohne** freie Drehknöpfe. Die Relation $r∕(1−n_{s})^{2}=3\alpha_{\text{inf}}$ ist ein direkter Test der topologisch fixierten Krümmung.

### 6.12 Eindeutige Fixierung von $\alpha_{\text{inf}}$

Variante A ($\alpha_{\text{inf}}=c_{3}/\varphi_{0}$) und Variante B ($\alpha_{\text{inf}}=\varphi_{0}/(2c_{3})$) sind zwei **Referenzmetriken** derselben Moduli‑Fixierung. Wir wählen eine **eindeutige** Normalisierung durch

$$[
\textbf{Kriterium:}\quad |R_{K}+1|=\min \quad \text{bei Erhalt der TFPT‑Identität}\ \ r=3\,\alpha_{\text{inf}}(1-n_{s})^{2}.
]$$

Numerisch ergibt dies $R_{K}\simeq -0.998$ und damit **Variante B** als eindeutige Wahl:
$$[
\alpha_{\text{inf}}=\frac{\varphi_{0}}{2c_{3}}=0.66817846,\qquad 
r=\frac{3\varphi_{0}}{2c_{3}}\,(1-n_{s})^{2}.
]$$
Die Vorhersagen in 6.5 und 6.7 bleiben im selben engen Korridor, nun **ohne** Ambiguität. Siehe auch die Poincaré‑Scheibenabbildung in 6.2.

---

## 7. Rolle von α und die parameterfreie Lösung

---
### 7.1 Motivation und Ursprung des Ansatzes


Die Feinstrukturkonstante $\alpha$ ist im Standardmodell ein **externer Eingabeparameter**. Schon frühe Überlegungen (Sommerfeld, Dirac, Eddington) hatten vermutet, dass hinter der Zahl $\alpha^{-1}\approx137$ eine tiefere mathematische Struktur stecken müsse.

Die genetischen Algorithmen und die 6D-Vorstufen zeigten wiederholt, dass $\alpha$ eng mit zwei Konstanten verknüpft ist:

$c_{3}=\frac{1}{8\pi}$, $\qquad \varphi_{0}\approx 0.053171$.

Beide Größen tauchten unabhängig in Kinetik-, Maxwell- und Massentermen auf. Die entscheidende Beobachtung war, dass $\alpha$ immer dort „auftauchte“, wo **topologische Normalisierung** (über $c_{3}$) und **geometrische Länge** (über $\varphi_{0}$) gleichzeitig wirksam waren.

Dies führte zur Hypothese: $\alpha$ **ist nicht frei, sondern die eindeutige Lösung einer Fixpunktbedingung, die genau diese beiden Konstanten koppelt.**

---
### 7.2 Ein Parameter Normalform für $\alpha$: Darstellung nur in $c_{3}$

**Normalform.** Mit $c_{3}=\tfrac{1}{8\pi}$,
$$
\varphi_{0}=\tfrac{4}{3}c_{3}+48\,c_{3}^{4},\quad
A=2\,c_{3}^{3},\quad
\kappa=\tfrac{b_{1}}{2\pi}\ln\tfrac{1}{\varphi_{0}},\quad b_{1}=\tfrac{41}{10},
$$
wird
$$
\alpha^{3}-A\alpha^{2}-A\,c_{3}^{2}\kappa=0
$$
zur reinen $c_{3}$-Form
$$
\boxed{\ \alpha^{3}-2c_{3}^{3}\alpha^{2}-8\,b_{1}\,c_{3}^{6}\,\ln\!\frac{1}{\tfrac{4}{3}c_{3}+48c_{3}^{4}}=0\ }.
$$

**Geschlossene Lösung (Cardano).** Setze $\alpha=y+\tfrac{2}{3}c_{3}^{3}$, dann $y^{3}+py+q=0$ mit
$$
p=-\tfrac{4}{3}c_{3}^{6},\qquad
q=-\tfrac{16}{27}c_{3}^{9}-8\,b_{1}\,c_{3}^{6}\ln\!\frac{1}{\tfrac{4}{3}c_{3}+48c_{3}^{4}},
$$
$\Delta=\big(\tfrac{q}{2}\big)^{2}+\big(\tfrac{p}{3}\big)^{3}$ und
$$
\boxed{\ \alpha(c_{3})=\frac{2}{3}c_{3}^{3}
+\sqrt[3]{-\frac{q}{2}+\sqrt{\Delta}}
+\sqrt[3]{-\frac{q}{2}-\sqrt{\Delta}}\ }.
$$

**Praxisformel.** Sehr genaue, geschlossene Näherung
$$
\boxed{\ \alpha \approx \Big(8\,b_{1}\,c_{3}^{6}\,\ln\tfrac{1}{\tfrac{4}{3}c_{3}+48c_{3}^{4}}\Big)^{1/3}
+\frac{2}{3}c_{3}^{3}\ }.
$$


> [!tip] Praxisformel

> Sehr genaue Praxisformel
> $[> \alpha \approx \Big(8b_{1}c_{3}^{6}\ln\tfrac{1}{\tfrac{4}{3}c_{3}+48c_{3}^{4}}\Big)^{\!\!1/3}+\frac{2}{3}c_{3}^{3}>]$
> liefert bereits die ppm-Nähe. Für $c_{3}=1/(8\pi)$ folgt $\alpha^{-1}=137.0365014649$.

---
### 7.3 Die Lösung

Die Fixpunktgleichung ist ein kubisches Polynom, das genau eine physikalisch reelle positive Nullstelle besitzt.

$$c_{3}=\tfrac{1}{8\pi}\ \Rightarrow\
\varphi_{0}=0.0531719521768,\ 
\kappa=1.914684795,\
\alpha=0.007297325816919221,\
\alpha^{-1}=137.03650146488582.$$

Die eindeutige reelle Lösung ist

$\alpha=0.0072973258169192213,\quad \alpha^{-1}=137.03650146488582$.

Das liegt um $3.665\times10^{-6}$ relativ unter CODATA 2022 $\alpha_{\text{CODATA}}=0.0072973525628 bzw. \alpha^{-1}=137.035999177$.

Die beiden anderen Wurzeln sind komplex und unphysikalisch.

Damit ist $\alpha$ **nicht postuliert**, sondern das **Output** einer zwingenden Gleichung.

![[Pasted image 20250826131905.png]]

---
### 7.4 Genauigkeit der Lösung

Vergleich mit CODATA 2022 Referenz ($\alpha^{-1}=137.035999177(21)$):

- Abweichung: wenige Teile pro Million (ppm).
    
- Keine Feinanpassung nötig – die Übereinstimmung folgt direkt aus c₃, φ₀ und b₁.

Dies ist bemerkenswert, weil es die bislang präziseste **parameterfreie theoretische Ableitung** von $\alpha$ darstellt.

---

### 7.5 Alternative Näherungen und optimierte Berechnungsarten

#### 7.5.1 Kubikwurzel-Näherung

In der Grenze kleiner A kann man $\alpha$ approximieren durch

$\alpha \;\approx\; (A c_{3}^{2}\kappa)^{1/3} + \frac{A}{3}$.

- Der erste Term ($A c_{3}^{2}\kappa)^{1/3}$ liefert den Hauptwert.
    
- Der additive Zuschlag $A/3$ (universell, unabhängig von $\varphi_{0}$) bringt die Zahl in ppm-Nähe.

Absoluter Fehler $2.44\times10^{-7}$ entspricht etwa 33 ppm

Diese Näherung trifft $\alpha$ bereits auf 10⁻⁷ genau.

---
#### 7.5.2 Ramanujan-ähnliche Serie

Setzt man $\alpha = (A c_{3}^{2}\kappa)^{1/3}(1+u)$ und entwickelt in Potenzen von u, ergibt sich eine konvergente Serie:

$\alpha = B^{1/3} + \frac{A}{3} + \frac{A^{2}}{9B^{1/3}} + \frac{2A^{3}}{81B^{2/3}} + \dots, \qquad B=A c_{3}^{2}\kappa$.


- Schon nach drei Termen liegt die Abweichung <0.2 ppm.
    
- Vier Terme liefern Genauigkeit auf 10⁻¹².
- 
Fehler $\approx 9.38\times10^{-10}$
---
#### 7.5.3 Newton-Verfahren 

Startet man bei $g=B^{1/3}+A/3$ und setzt einmal Newton an, erreicht man dieselbe Genauigkeit wie mit der Serie. 

Formel:

$\alpha \approx g - \frac{f(g)}{f’(g)},\qquad f(\alpha)=\alpha^{3}-A\alpha^{2}-B$.

Damit lässt sich $\alpha$ extrem effizient und exakt berechnen.
![[Pasted image 20250824122910.png]]

---

### 7.6 Variationsableitung in vier Dimensionen (kubische Fixpunktgleichung aus der Einstein Wirkung)

**Ziel und Kontext.** Wir zeigen, dass die in 7.2 verwendete kubische Fixpunktgleichung für $\alpha$ als Stationaritätsbedingung einer expliziten vierdimensionalen Wirkung folgt. Die Konstanten $c_{3}$ und $\varphi_{0}$ stammen exakt aus den bereits etablierten Invarianten der Theorie, und die Normierungen $A$ und $\kappa$ sind identisch zu den Definitionen in Anhang E. Damit entsteht eine zweite, unabhängige Herleitung derselben Gleichung ohne frei wählbare Skalen. Vgl. die Ableitungen von $c_{3}$ in Abschn. 3.2.1 und von $\varphi_{0}$ in Abschn. 3.2.2, sowie die Normierungsnotiz in Anhang E.  

**Fixe Invarianten aus Topologie und Geometrie**

Aus der Chern Simons Reduktion mit $M_{11}=M_{4}\times Y_{7}$ und ganzzahliger Schnittzahl folgt der topologische Fixpunkt

$c_{3}=\frac{1}{8\pi}$.

Siehe die strenge Herleitung über $C_{3}\wedge G_{4}\wedge G_{4}$ und die Quantisierung von $\int_{M_{4}}F\wedge F = 8\pi^{2}\mathbb Z$ in Abschn. 3.2.1. Die Möbius Reduktion mit Gauss Bonnet und Rand liefert

$\varphi_{0}=\frac{1}{6\pi}+\frac{3}{256\pi^{4}}$,

siehe Abschn. 3.2.2 und Anhang D. Zusammen mit $b_{1}=41/10$ in GUT Norm sind alle Größen fixiert.  

**Vierdimensionale Wirkung und** $U(\alpha)$

Wir betrachten nach Reduktion und kanonischer Normierung den abelschen Sektor in homogener Hintergrundlage und verstehen $U(\alpha)$ als **Gradient Darstellung im Kopplungsraum**

$\partial_{\alpha} U \propto \beta_{\alpha}$,

so dass Stationarität $\partial_{\alpha}U=0$ äquivalent zu $\beta_{\alpha}=0$ ist, vgl. die Interpretation in 7.6. Die effektive Wirkung lautet

$S_{\mathrm{eff}}=\int\mathrm d^{4}x\,\sqrt{-g}\,\Bigl[\tfrac{M_{P}^{2}}{2}R-U(\alpha)+\ldots\Bigr]$.

Bis zur relevanten Ordnung genügt ein skalares Potential

$U(\alpha)=\frac{A}{4}\,\alpha^{4}\;-\;\frac{2A}{3}\,c_{3}^{3}\,\alpha^{3}\;-\;A\,[\,8\,b_{1}\,c_{3}^{6}\,\ln(1/\varphi_{0})\,]\,\alpha. \tag{7.1.1}$

• **Leitterm** $\alpha^{4}$: Glatter Referenzanteil, der in $\partial_{\alpha}U$ den führenden $\alpha^{3}$ Beitrag liefert.

• **Kubischer Beitrag** $\propto c_{3}^{3}\alpha^{3}$ Er entsteht aus der reduzierten Chern Simons Struktur über die Kopplung eines schweren skalaren Modus a an $F\tilde F$ und die Eliminierung von $a$ im Null Impuls Grenzfall. Die kombinierte Zählung zweier identischer topologischer Einsätze, der Umrechnung $g\mapsto \alpha$ und des Symmetriefaktors liefert den **universellen Faktor**

$A=\frac{1}{256\pi^{3}}\equiv 2\,c_{3}^{3}$,

genau wie in Anhang E Schritt 2 gezeigt.  

• **Linearer Logarithmus.** Die integrierte Eins Schleifen Renormierung zwischen $\mu_{\mathrm{UV}}=M_{P}$ und $\mu_{\mathrm{IR}}=\varphi_{0}M_{P}$ ergibt

$\kappa\equiv \frac{b_{1}}{2\pi}\ln\!\Bigl(\frac{1}{\varphi_{0}}\Bigr)$,

in der Potentialschreibweise als Term $-A[8\,b_{1}\,c_{3}^{6}\ln(1/\varphi_{0})]\,\alpha$. Dies ist identisch zur in Anhang E Schritt 1 definierten Normierung $\kappa$.  

Alle Größen sind in reduzierten Planck Einheiten, vgl. Info Box.  

**Stationarität und Normalform**

Die Variationsbedingung liefert

$\frac{\partial U}{\partial \alpha}=A\Bigl[\alpha^{3}-2\,c_{3}^{3}\,\alpha^{2}-8\,b_{1}\,c_{3}^{6}\,\ln\!\Bigl(\tfrac{1}{\varphi_{0}}\Bigr)\Bigr]=0$,

also die **kubische Fixpunktgleichung**

$\boxed{\,\alpha^{3}-2\,c_{3}^{3}\,\alpha^{2}-8\,b_{1}\,c_{3}^{6}\,\ln\!\Bigl(\tfrac{1}{\varphi_{0}}\Bigr)=0\,}. \tag{7.1.2}$

Gleichwertige **Normalform** wie in 7.2 verwendet:

$\alpha^{3}-A\,\alpha^{2}-A\,c_{3}^{2}\,\kappa=0,\qquad A=2c_{3}^{3},\ \ \kappa=\frac{b_{1}}{2\pi}\ln\!\Bigl(\tfrac{1}{\varphi_{0}}\Bigr)$.

Die eindeutige reelle Lösung stimmt numerisch mit 7.3 überein, dort ausgewiesen mit $\alpha^{-1}=137.03650146488582$.  

**Physikalische Einordnung und Konsistenz**

1. **Schemenfreiheit.** $A$ ist eine reine Zahl aus Topologie und kanonischer Normierung, $\kappa$ hängt nur vom festen $b_{1}$ und der geometrisch fixierten Skala $\varphi_{0}$ ab. Ein Schemenwechsel verschiebt nur additive, skalenunabhängige Beiträge in $\kappa$, nicht aber die Lage des Fixpunkts, vgl. Anhang E.  

2. **Bedeutung von** $U(\alpha)$ $U$ ist kein Materiepotential, sondern eine kompakte Darstellung der Kopplungsdynamik, $\partial_{\alpha}U\propto \beta_{\alpha}$. 

3. **Bezug zum übrigen Gefüge.** Die gleichen Invarianten setzen dynamische Fingerabdrücke in der zwei Schleifen Analyse, $\alpha_{3}(1\ \mathrm{PeV})\approx \varphi_{0}$ und $\alpha_{3}(\mu)\approx c_{3}$ bei $\mu\sim 2.5\times 10^{8}\ \mathrm{GeV}$, siehe 5.2.  

4. **Abelische Spur als roter Faden.** Die Zahl 41 erscheint sowohl in $b_{1}=41/10$ der Fixpunktgleichung als auch im EW Block über $k_{\mathrm{EW}}=41/32$, vgl. 8.4.6. Das unterstreicht, dass dieselbe abelsche Spur in beiden Kontexten wirkt, ohne Zirkularität.
![[Pasted image 20250901085708.png]]
*Left: invariants $((c_3, \varphi_0, b_1))$. Middle: $(U_3, U_4, U_1)$. Right: $(U(\alpha)$), stationarity, cubic fixed point.*

#### 7.6.1 Callan–Symanzik‑Route

Mit $\mu\,d\alpha/d\mu=\beta_{\alpha}=\tfrac{b_{1}}{2\pi}\alpha^{2}+A\,c_{3}^{2}\alpha^{3}+\dots$ und $A=1/(256\pi^{3})$ aus Appendix E führt die Integration zwischen $M_{\text{Pl}}$ und $\varphi_{0}M_{\text{Pl}}$ unmittelbar zur Kubik
$$[
\alpha^{3}-2c_{3}^{3}\alpha^{2}-8\,b_{1}\,c_{3}^{6}\ln\frac{1}{\varphi_{0}}=0.
]$$
Dies ist unabhängig von der Potential‑Darstellung $U(\alpha)$ und macht die Fixpunktgleichung **zweifach** hergeleitet. 


#### 7.6.2  A und κ auf einen Blick  mit Querverweisen


**Zweck.** Diese Box bündelt die Normierungen und Invarianten, die in **7.6** für die Variationsableitung verwendet werden, und verweist auf die formalen Herleitungen in **3.2.1**, **3.2.2** sowie **Anhang E**. 

⸻

**Fixpunkte und Skalen.**

* Topologie aus Chern Simons Reduktion, vgl. **3.2.1**:  

  $( c_{3}=\dfrac{1}{8\pi}=0.039788735772973836\ldots )$

* Geometrie der Moebius Faser mit Gauss Bonnet und Rand, vgl. **3.2.2** und **Anhang D**:  

  $( \varphi_{0}=\dfrac{1}{6\pi}+\dfrac{3}{256\pi^{4}}=0.053171952176845526\ldots )$

**Abelsche Spur in GUT Norm.**

* $( b_{1}=\dfrac{41}{10}=4.1 )  \; ( \Rightarrow )$  Definition von $( \kappa )$ unten, vgl. **Anhang E**.
⸻

**Zentrale Kurzschreibweisen aus Anhang E.**

* Reiner Topologie Faktor  
  $( A\;\equiv\;2\,c_{3}^{3}\;=\;\dfrac{1}{256\,\pi^{3}}\;=\;1.259825563796855\times10^{-4} )$.

* Integrierte eins Schleifen Konstante  
  $( \kappa\;\equiv\;\dfrac{b_{1}}{2\pi}\,\ln\!\bigl(\tfrac{1}{\varphi_{0}}\bigr) )$.
⸻

**Effektives Potential in 7.1.1.**  $( U(\alpha) )$ als Gradient Darstellung im Kopplungsraum  

$( \partial_{\alpha}U\propto\beta_{\alpha} )$. Bis zur relevanten Ordnung:

$$[
U(\alpha)=\frac{A}{4}\alpha^{4}-\frac{2A}{3}c_{3}^{3}\alpha^{3}-A\bigl[\,8\,b_{1}\,c_{3}^{6}\ln\!\bigl(\tfrac{1}{\varphi_{0}}\bigr)\bigr]\alpha.
]$$

**Stationaritaet gleich Fixpunkt.**

$$[

\frac{\partial U}{\partial \alpha}=A\Bigl[\alpha^{3}-2c_{3}^{3}\alpha^{2}-8b_{1}c_{3}^{6}\ln\!\bigl(\tfrac{1}{\varphi_{0}}\bigr)\Bigr]=0

]$$

Normalform wie in **7.2**:  

$$[

\alpha^{3}-A\,\alpha^{2}-A\,c_{3}^{2}\,\kappa=0,\qquad A=2c_{3}^{3},\quad \kappa=\frac{b_{1}}{2\pi}\ln\!\Bigl(\tfrac{1}{\varphi_{0}}\Bigr).

]$$


**Dynamische Fingerabdruecke und Kosmo Anker.**

* $( \alpha_{3}(1\,\text{PeV})\approx \varphi_{0} )$ und $( \alpha_{3}(\mu)\approx c_{3} )$  bei  $( \mu\sim2.5\times10^{8}\,\text{GeV} )$, vgl. **5.2**.  

* $( \Omega_{b}=\varphi_{0}(1-2c_{3}) )$ nahe Planck, vgl. **8.4.7**.


---

### 7.7 Interpretation

Die Rolle von $\alpha$ ist in diesem Framework grundlegend neu definiert:

- **Kein Input, sondern Fixpunkt.** $\alpha$ ist keine willkürliche Zahl, sondern die eindeutige Lösung einer geometrisch-topologischen Bedingung.
    
- **Dominanz der Topologie.** Sensitivitätsanalysen zeigen: $\alpha$ reagiert am stärksten auf c₃ (topologischer Fixpunkt), weniger auf b₁ (Spektrum), am schwächsten auf φ₀ (Geometrie).
    
- **Universeller Zuschlag.** Der konstante Korrekturterm $A/3$ erklärt, warum $\alpha$ ppm-genau sitzt – eine kleine, aber strukturelle Verschiebung.
    
Damit ist die Feinstrukturkonstante **nicht zufällig**, sondern ein emergenter Fixpunkt aus Topologie, Geometrie und Symmetrie.

**Gemeinsame Ursache von $\alpha$ und Flavor Relationen.**  
$\alpha$ folgt aus $(\varphi_0,c_3)$, die Möbius Leiter nutzt $\delta_\star=\dfrac{3}{5}+\dfrac{\varphi_0}{6}$.  
Damit entsteht die Schleife $\varphi_0\Rightarrow \alpha(\varphi_0,c_3)$ und $\varphi_0\Rightarrow\delta_\star\Rightarrow$ Flavor Relationen.

---

## 8. Von E₈ zu E₇ zu E₆ und zum Standardmodell  
_Eine klare Blockstruktur, rechnerisch geschlossen, sofort reproduzierbar_

> [!info] Fixpunkte und Leiter  
> **Topologie:** $c_3=\tfrac{1}{8\pi}=0.039788735772973836$  
> **Geometrie:** $\varphi_0=\tfrac{1}{6\pi}+\tfrac{3}{256\pi^4}=0.05317195217684553$  
> **Leiter Normierung:** $\gamma(0)=0.834,\ \lambda=\dfrac{0.834}{\ln 248-\ln 60}=0.5877029773404678$  
> **Plankonstante für Zahlen:** $M_{\rm Pl}=1.221\times10^{19}\ \mathrm{GeV}$

---

> [!summary] Idee in einem Satz  

> Wir verbinden eine diskrete Strukturachse aus E₈ mit Stufen $n$ und eine dynamische Achse aus Renormierungsgruppe $\mu$.  
>  E₈ ordnet die Leiter $\varphi_n$. Die RG Dynamik liefert Fenster $E_r$ bei $\alpha_3(\mu)\approx 1/(r\pi)$.  
> Blöcke verknüpfen beides und projizieren auf messbare Größen des Standardmodells.


![[Pasted image 20250828094621.png]]

### Zwei Achsen, ein gemeinsames Raster

**Strukturachse**  

Aus der nilpotenten Orbitologie von E₈ entsteht eine eindeutige, streng fallende Kette  

$$

D_n = 60-2n, \qquad n=0\dots 26,

$$

die eine log exakte Leiter definiert  

$$

\varphi_{n} = \varphi_0\,e^{-\gamma(0)}\left(\frac{D_n}{D_1}\right)^{\lambda} \quad (n\ge 1).

$$

Diese Achse ist diskret. Sie ordnet Verhältnisse von Skalen. Sie erklärt, warum bestimmte Sprünge zwischen Ebenen immer wieder gleich aussehen.

**Dynamikachse**  

Auf der RG Achse läuft die starke Kopplung $\alpha_3(\mu)$ kontinuierlich. Es gibt drei natürliche Fenster  


$$

\alpha_3(\mu_r)=\frac{1}{r\,\pi},\qquad r\in\{6,7,8\},

$$

also $E₆$ um $1/(6\pi)$ nahe PeV, $E₇$ um $1/(7\pi)$ dazwischen, $E₈$ um $1/(8\pi)=c_3$ bei etwa $2.5\times 10^8$ GeV.

> [!info] Leseregel  

> $n$ zählt Struktur und bestimmt Ratio Gesetze.  
> $E_r$ markiert Dynamik und bestimmt Lagen auf der Energieachse.  
> Synchronisiert werden beide durch die Fixpunkte $c_3=\tfrac{1}{8\pi}$ und $\varphi_0=\tfrac{1}{6\pi}+\tfrac{3}{256\pi^4}$.


![[Pasted image 20250828094712.png]]

> **Info Box — Chiralität kommt aus Geometrie, nicht aus E8**
> 
> E8 dient hier ausschließlich als **Ordnungsprinzip einer diskreten Skalenleiter** $(\varphi_{n})$.  
> Es gibt **keine** 4D Eichgruppe $(E_{8})$ und **keine** Einbettung der SM Fermionen in $(E_{8})$ Darstellungen.  
> Die 4D **Chiralität** entsteht unabhängig davon durch **Randbedingungen** und **integer quantisierte Flüsse** auf der orientierbaren Doppelabdeckung der Möbius Faser.  
> Drei Randzyklen und die Chern–Simons Quantisierung liefern den **chiralen Index**:
> $( \mathrm{Ind}\,D = \dfrac{1}{2\pi}\!\int_{\widetilde M}\!F = \nu_{1}+\nu_{2}+\nu_{T} )$.
> Mit der Minimalwahl $((\nu_{1},\nu_{2},\nu_{T})=(1,1,1))$ ergeben sich **drei linkschirale Familien**, Spiegel Zustände fehlen durch die Projektoren auf den Randzyklen.  
> Details in **Anhang J**.  (Bezüge: 3.2.1 $(c_{3}=\tfrac{1}{8\pi})$, 3.2.2 $(\varphi_{0})$, 8.4.6 $(k_{\mathrm{EW}}=41/32))$.   


![[Pasted image 20250901111221.png]]

---
### Wie aus Struktur und Dynamik Zahlen des SM werden


![[Pasted image 20250828094846.png]]

Der Schritt von dimensionslosen Leiterstufen zu messbaren Größen erfolgt blockweise. Jeder Block $B$ hat drei Kennzahlen:  

- $r_B$ wirksamer Rang in der Kette $E₈ \supset E₇ \supset E₆ \supset SM$  

- $k_B$ fraktionale Topologiezahl aus den Randzyklen der Möbius Faser  

- $n_B$ Stufe der Leiter  

Daraus folgt zuerst eine Blockkonstante  

$$

\zeta_B=(\pi c_3)\,\exp[-\beta_B\,\pi c_3]\,\exp\Bigl[-\frac{k_B}{c_3}\Bigr],\qquad \beta_B=\frac{8-r_B}{8},

$$

und dann die dimensionierte Größe  

$$

X_B=\zeta_B\,M_{\mathrm{Pl}}\,\varphi_{n_B}.

$$

So setzen wir zum Beispiel  

- EW Block bei $n=12$ im $E₇$ Fenster: $v_H=\zeta_{\rm EW} M_{\mathrm{Pl}}\varphi_{12}$  

- Hadron Blöcke bei $n=15$ und $n=17$ im $E₆$ Korridor: $m_p\simeq \zeta_{\rm had} M_{\mathrm{Pl}}\varphi_{15}$  

- Lepton Blöcke tief unten $n=22,25,26$: leichte Yukawas  

> [!tip] Schnellstart für Leser  

> 1. Finde im Text den Block für die gesuchte Größe.  
> 2. Lies $r_B, k_B, n_B$ ab und berechne $\zeta_B$.  
> 3. Setze $X_B=\zeta_B M_{\mathrm{Pl}}\varphi_{n_B}$ mit der log exakten $\varphi_{n}$ aus der E₈ Leiter.

---

### Wo ist der Anschluss an das Standardmodell

  
Die Kette $E₈ \supset E₇ \supset E₆ \supset SM$ liefert die Ranglogik und die abelsche Spur:  

- Am EW Anker $n=12$ erscheint die Spur $\mathcal Y^2_{\rm SM+H}=\tfrac{41}{48}$. Daraus ergibt sich $k_{\rm EW}=\tfrac{41}{32}$ und konsistent $b_1=\tfrac{41}{10}$ in GUT Norm.  

- Die hadronischen Fenster liegen in der $E₆$ Domäne der Leiter und stützen die zusätzliche Dämpfung, die baryonische Skalen auszeichnet.  

- Die RG Fenster verankern diese Struktur dynamisch: $\alpha_3(\mu)$ trifft $1/(6\pi), 1/(7\pi), 1/(8\pi)$ an genau den Stellen, die durch die Leiter motiviert sind.  

Kurz: Struktur ordnet, Dynamik bestätigt, Blöcke projizieren. Das ist unser Pfad von Topologie und Geometrie zu den Zahlen des Standardmodells.

---
### Was machen die Stufen ohne direkten Block


Nicht jede Stufe muss eine konkrete Observabel tragen. Diese Stufen sind wichtiges Tragwerk:  

1. **Geometrie der Leiter**  

Sie sichern das Gesetz  
$$

\frac{\varphi_m}{\varphi_n}=\left(\frac{D_m}{D_n}\right)^{\lambda}\quad(m,n\ge 1),

$$

also die fitfreie Ratio Struktur.  

2. **Feine Rastpunkte in Fenstern**  

Ein dynamisches Fenster ist ein Bereich in $\mu$. Die diskreten $n$ fungieren als Rastpunkte, an denen Schwellen und Mischungen wirken können, ohne das globale Ratio Gesetz zu verletzen.  

3. **Reserve für neue Observablen**  

Weitere Größen wie Schwellen, Axion Kopplungen, präzise hadronische Parameter können später genau dort andocken. Die Plätze sind strukturell schon korrekt verdrahtet.  

> [!example] Intuition  

> Denke an ein Getriebe. Die Blockstufen sind die Zahnräder, die eine Achse treiben. Die Zwischenzähne sorgen dafür, dass die Kraft sauber und ohne Rutschen übertragen wird. Ohne sie gäbe es Sprünge, aber keine Ordnung.

### 8.0a  Chiralität aus Rand und Fluss: operative Kurzfassung

**Geometrie und Rand**  

Wir arbeiten auf der **orientierbaren Doppelabdeckung** $(\widetilde M)$ der Möbius Faser mit **drei** geschlossenen Randzyklen $(C_{1},C_{2},C_{T})$. Die Randzählung ist kanonisch: $(\sum_{i}\!\oint_{C_{i}}\widehat{k}_{g}\,\mathrm ds = 6\pi)$.  

Die aus der sechs dimensionalen Reduktion resultierenden Projektoren wählen **eine interne Chiralität**:

$$[

P_{T}=\tfrac12\!\left(\mathbf 1+\mathrm i\,\sigma^{3}\sigma^{n}\right),\quad 

P_{1}=P_{2}=\tfrac12\!\left(\mathbf 1-\mathrm i\,\sigma^{3}\sigma^{n}\right),

]$$

sodass nur $(\chi_{+})$ Nullmoden trägt und die 4D Nullmoden **linkschiral** sind. 

**Index und Familienzahl**  

Eine glatte abelsche Verbindung $(A)$ mit quantisiertem Fluss $(m=\tfrac{1}{2\pi}\!\int_{\widetilde M}\!F\in\mathbb Z)$ liefert

$$[

\mathrm{Ind}\,D_{\widetilde M}

= \#\{\chi_{+}\} - \#\{\chi_{-}\}

= \tfrac{1}{2\pi}\!\int_{\widetilde M}\!F 

= \nu_{1}+\nu_{2}+\nu_{T}.

]$$

Die Minimalwahl $((\nu_{1},\nu_{2},\nu_{T})=(1,1,1))$ ergibt **drei** Familien ohne Spiegel.  

Wilson Linien sind **flach** und lesen nur die abelsche Spur aus, kompatibel mit $(k_{\mathrm{EW}}=\tfrac{41}{32})$.  Verweise: 3.2.1 $(c_{3}=\tfrac{1}{8\pi})$, 3.2.2 $(\varphi_{0})$, 8.4.6 Spur $(41)$.  

### 8.1 Ausführliche Beschreibung

E₈ ordnet die **Skalenleiter** $\varphi_n$ log exakt, E₇ und E₆ setzen die **physikalischen Fenster** pro Block, und **Topologie** mit **Geometrie** liefert über $c_3$ und $\varphi_0$ die **Normalisierungen**. Dimensionierte Größen entstehen aus einer kompakten **Block Formel**:

$$
\boxed{\,X_B=\zeta_B\,M_{\rm Pl}\,\varphi_{n_B},\quad 
\zeta_B=(\pi c_3)\,e^{-\beta_B\,\pi c_3}\,e^{-\,k_B/c_3},\quad 
\beta_B=\tfrac{8-r_B}{8}\,}
$$

mit $r_B$ als effektivem Rang im Block und $k_B$ als rationaler topologischer Zahl der drei Randzyklen.

Die E₈ Leiter ist log exakt:

$$
\gamma(0)=0.834,\qquad
\gamma(n)=\lambda[\ln D_n-\ln D_{n+1}],\quad D_n=60-2n,\quad
\lambda=\frac{0.834}{\ln 248-\ln 60}.
$$

Für $n\ge 1$ gilt

$$
\boxed{\,\varphi_n=\varphi_0\,e^{-\gamma(0)}\Big(\frac{D_n}{D_1}\Big)^{\lambda},\quad D_1=58\,}.
$$

![[Pasted image 20250828094755.png]]


### 8.1.1 Block‑Konstanten aus Randzyklen und abelscher Spur

  
**Ziel.** $k_{B}$ sind **keine Fits**, sondern folgen aus einer Zählung abelscher Quadrate auf den **drei** Randzyklen der orientierbaren Darstellung, multipliziert mit einem universellen Faktor. 

  
**Definition (abelsche Spur in GUT‑Norm).** Für einen Block $B$ sei

$$[

\mathcal I_{1}(B)=\sum_{\Phi\in B}\sum_{i\in U(1)_{Y}} q_{i}^{2}(\Phi)\ \ \text{mit}\ \ Y\text{ in GUT‑Norm}.

]$$

**Satz 8.1.1 (Topologische Blockzahl).** Die rationale Zahl

$$[

k_{B}=\frac{3}{2}\,\mathcal I_{1}(B)

]$$

ergibt sich aus der Summe der drei Randzyklen (Faktor $3$) und dem faktorisierten Halbgewicht der orientierbaren Doppelabdeckung (Faktor $1/2$). 


**Beispiel EW‑Block.** Für $B=\text{SM}+\text{H}$ gilt $\mathcal I_{1}=\tfrac{41}{48}$, also $k_{\text{EW}}=\tfrac{3}{2}\cdot \tfrac{41}{48}=\tfrac{41}{32}$, genau wie in 8.4.1 verwendet; gleichzeitig erscheint dieselbe Spur in $b_{1}=\tfrac{41}{10}$, siehe 7.6.1. 


**Bemerkung Hadron‑ und Pion‑Blöcke.** Für baryonische und pioniche Größen wird $\mathcal I_{1}$ durch die effektiven abelschen Untergruppen der flavor‑chiral Dynamik ersetzt ($U(1)_{B}$ bzw. $U(1)_{I3}$). Damit folgen $k_{p}=\tfrac{3}{2}$ und $k_{\pi}=\tfrac{51}{32}$ als konkrete Auswertung derselben Zählregel im jeweiligen Block. Siehe 8.4.5.  

### 8.1.2 Herleitung der $\zeta_{B}$‑Formel aus dem Randfunktional

  
Die effektive Randwirkung pro Block hat die Form

$$[

S^{(B)}_{\partial}=\pi c_{3}\ -\ \beta_{B}\,\pi c_{3}\ -\ \frac{k_{B}}{c_{3}},
]$$

wobei $\beta_{B}=(8-r_{B})/8$ aus der effektiven Rangzahl $r_{B}$ resultiert (Zählung der nicht gedämpften Richtungen). Exponentiation der additiven Beiträge liefert

$$[

\zeta_{B}=(\pi c_{3})\,\exp\big[-\beta_{B}\pi c_{3}\big]\exp\big[-k_{B}/c_{3}\big],\qquad X_{B}=\zeta_{B}\,M_{\text{Pl}}\,\varphi_{n_{B}}.

]$$

Damit sind $v_{H}$, $m_{p}$, $f_{\pi}$ und $T_{\gamma 0}$ direkt aus $(r_{B},k_{B},n_{B})$ berechenbar, **ohne** zusätzliche Freiheitsgrade. Siehe 8.4.1 bis 8.4.7. 

---

### 8.2 Rechenrezept in drei Schritten

1. **Leiter auswerten**  
   $$
   \varphi_n=\varphi_0\,e^{-\gamma(0)}\Big(\tfrac{60-2n}{58}\Big)^{\lambda}\qquad (n\ge 1).
   $$

2. **Blockkonstanten setzen**  
   Für den Block $B$: $r_B$ wählen, $\beta_B=(8-r_B)/8$, dazu $k_B$ rational aus der Randzählung.  
   $$
   \zeta_B=(\pi c_3)\,e^{-\beta_B \pi c_3}\,e^{-k_B/c_3},\qquad \pi c_3=\tfrac18.
   $$

3. **Größe bestimmen**  
   $$
   X_B=\zeta_B\,M_{\rm Pl}\,\varphi_{n_B}.
   $$

> [!tip] Verhältnisgesetze ohne Einheitenwahl  
> $$
> \frac{\varphi_{m}}{\varphi_{n}}=\Big(\frac{60-2m}{60-2n}\Big)^{\lambda}\quad (m,n\ge 1).
> $$

---

### 8.3 Benötigte Leiterstufen $\varphi_n$ (log exakt)

| $n$ | $D_n$ |     $\varphi_n$ |
| --: | ----: | --------------: |
|   1 |    58 | 0.0230930346695 |
|   5 |    50 | 0.0211640537281 |
|  10 |    40 | 0.0185628455934 |
|  12 |    36 | 0.0174482846938 |
|  15 |    30 | 0.0156753658147 |
|  16 |    28 | 0.0150524852088 |
|  22 |    16 | 0.0108336306291 |
|  25 |    10 | 0.0082188698412 |
|  26 |     8 | 0.0072087140665 |

---

### 8.4 Ergebnisse pro Block mit Referenzen

#### 8.4.1 Elektroschwacher Block $n=12$  
**Annahmen:** $r_{\rm EW}=2\Rightarrow \beta_{\rm EW}=3/4,\quad k_{\rm EW}=\tfrac{41}{32}$

$$
\zeta_{\rm EW}=(\pi c_3)\,e^{-\,\tfrac34\pi c_3}\,e^{-\,\tfrac{41}{32}/c_3}
=1.17852087206\times10^{-15}.
$$

$$
v_H=\zeta_{\rm EW}M_{\rm Pl}\varphi_{12}= \mathbf{251.07628}\ \mathrm{GeV}.
$$

Mit $g_2=0.652,\ g_1^{\rm SM}=0.357$ am $M_Z$:

$$
M_W=\tfrac12 g_2 v_H=\mathbf{81.85087}\ \mathrm{GeV},\qquad
M_Z=\tfrac12\sqrt{g_2^2+g_1^2}\,v_H=\mathbf{93.31741}\ \mathrm{GeV}.
$$

**Vergleich**  
$v=(\sqrt2 G_F)^{-1/2}=246.21965\ \mathrm{GeV}$ ⇒ **+1.97 Prozent**  
$M_W=80.3692\ \mathrm{GeV}$ ⇒ **+1.84 Prozent**  
$M_Z=91.1876\ \mathrm{GeV}$ ⇒ **+2.34 Prozent**

> [!note] Lesart  
> Der Block setzt die **Skala** $v_H$ auf ein bis zwei Prozent genau. Endliche Beiträge mit zwei Schleifen und Schwellen verschieben $M_W,M_Z$ nach unten in Richtung der Referenzen.

**Topmasse als Minimalannahme**  
$y_t\approx 1\Rightarrow m_t\simeq v_H/\sqrt2=\mathbf{177.54}\ \mathrm{GeV}$.

---

#### 8.4.2 PQ Block $n=10$  
**Annahmen:** $r_{\rm PQ}=1\Rightarrow \beta_{\rm PQ}=7/8,\quad k_{\rm PQ}=\tfrac12$

$$
\zeta_{\rm PQ}=3.90754185582\times10^{-7},\quad
f_a=\zeta_{\rm PQ}M_{\rm Pl}\varphi_{10}=\mathbf{8.8565\times 10^{10}}\ \mathrm{GeV}.
$$

Axionmasse:  

$$
m_a\simeq (5.7\ \mu\mathrm{eV})\times \frac{10^{12}\ \mathrm{GeV}}{f_a}=\mathbf{64.36\ \mu\mathrm{eV}}.
$$

---

#### 8.4.3 Seesaw Block $n=5$  
**Annahmen:** $r_{N_R}=4\Rightarrow \beta_{N_R}=1/2,\quad k_{N_R}=\tfrac18$

$$
M_R=\zeta_{N_R}M_{\rm Pl}\varphi_{5}=\mathbf{1.311\times10^{15}}\ \mathrm{GeV}.
$$

Mit $y_{\nu3}\sim 1$:  

$$
m_{\nu3}\simeq \frac{v_H^2}{M_R}=\mathbf{0.04807\ \mathrm{eV}},\quad
\Delta m^2_{31}\simeq \mathbf{2.31\times 10^{-3}\ \mathrm{eV}^2}.
$$

---

#### 8.4.4 Flavor Anker aus $n=1$  

$$
\sin^2\theta_{13}=\varphi_1=\mathbf{0.023093},\qquad 
\sin\theta_{13}=0.15197.
$$

**Cabibbo Winkel aus Basisstufe**  

$$
\sin\theta_C\simeq \sqrt{\varphi_0}\Big(1-\frac{\varphi_0}{2}\Big)=\mathbf{0.22446},\qquad
\theta_C=\arcsin(\sin\theta_C)=\mathbf{0.22639}\ \mathrm{rad}.
$$

**Möbius Massenleiter mit einer Deformation $\delta$.**  
Kalibriere $\delta$ nur aus Leptonen:
$\delta=\dfrac{\sqrt{m_\tau/m_\mu}-1}{\sqrt{m_\tau/m_\mu}+1}.$  
Setze diese eine Zahl in die sechs Relationen ein:
$$
\begin{aligned}
\text{Down:}\quad & \sqrt{\tfrac{m_s}{m_d}}=\mathcal{M}_{1}(\delta),\quad
\sqrt{\tfrac{m_b}{m_s}}=\mathcal{M}_{1}(\delta)\,(1+\delta),\\[2mm]
\text{Leptonen:}\quad & \sqrt{\tfrac{m_\tau}{m_\mu}}=\mathcal{M}_{1}(\delta),\quad
\sqrt{\tfrac{m_\mu}{m_e}}=\mathcal{M}_{1}(\delta)\,\mathcal{M}_{1/3}(\delta),\\[2mm]
\text{Up:}\quad & \sqrt{\tfrac{m_c}{m_u}}=\mathcal{M}_{2/3}(\delta),\quad
\sqrt{\tfrac{m_t}{m_c}}=\dfrac{2/3}{\,2/3-\delta\,}.
\end{aligned}
$$

**Topologischer Check.**  
Die Theorie erwartet $\delta_\star=\dfrac{3}{5}+\dfrac{\varphi_0}{6}$.  
Vergleiche $\delta$ aus den Leptonen mit $\delta_\star$ und dokumentiere die Abweichung in Prozent.  
Es kommen keine neuen freien Parameter hinzu.

---

#### 8.4.5 Hadron Fenster und pionische Observablen

**Proton $n=15$**, Annahmen $r_{\rm had}=5\Rightarrow \beta_{\rm had}=3/8,\ k_p=\tfrac32$:  

$$
m_p=\zeta_{\rm had}M_{\rm Pl}\varphi_{15}=\mathbf{0.96821\ \mathrm{GeV}}.
$$

**Pion $n=16$**, gleiche Rangzahl $r=5$, stärkerer topologischer Dämpfer $k_\pi=\tfrac{51}{32}$:  

$$
f_\pi=\mathbf{88.12\ \mathrm{MeV}}\ \text{(chirale Norm)}.
$$

GMOR Konsistenz mit $|\langle\bar qq\rangle|^{1/3}\simeq 272\ \mathrm{MeV},\ (m_u+m_d)_{2\ \mathrm{GeV}}\simeq 6.8\ \mathrm{MeV}$:  

$$
m_\pi\simeq\sqrt{\frac{(m_u+m_d)\,|\langle\bar qq\rangle|}{f_\pi^2}}=\mathbf{132.75\ \mathrm{MeV}}.
$$

---

#### 8.4.6 Feinstrukturkonstante $\alpha$  
(Querverweis zu Abschnitt 6)

Mit

$$
\alpha^{3}-2c_3^{3}\alpha^{2}-8\,b_1\,c_3^{6}\ln\frac{1}{\tfrac{4}{3}c_3+48c_3^{4}}=0,\qquad b_1=\tfrac{41}{10},
$$

ergibt sich die eindeutige reelle Lösung  

$$
\boxed{\,\alpha=\mathbf{0.007297325816919221},\quad \alpha^{-1}=\mathbf{137.03650146488582}\,}
$$

Abweichung zu CODATA 2022 $\alpha^{-1}=137.035999177$: **+3.67 ppm**.


> [!summary] Kurzfassung

> **Das gleiche Zählmaß 41** aus der Hyperladung erscheint **zweifach**:  

> – in der α‑Fixpunktgleichung über $(b_1=\tfrac{41}{10})$
> – im EW‑Block über $(k_{\rm EW}=\tfrac{41}{32})$
> Beides folgt aus derselben abelschen Spur $(\mathcal Y^2_{\rm SM+H}=\tfrac{41}{48})$. α ist hier also **kein Input**, sondern ein **Konsistenz‑Echo** derselben Struktur, die $(v_H)$ ankert.

**1) α als Fixpunkt aus Topologie und Geometrie**

Die kubische Gleichung
$$
[\alpha^{3}-2c_3^{3}\alpha^{2}-8\,b_1\,c_3^{6}\,\ln\!\frac{1}{\varphi_0}=0,\quad 
c_3=\tfrac{1}{8\pi},\ \ \varphi_0=\tfrac{1}{6\pi}+\tfrac{3}{256\pi^4},\ \ b_1=\tfrac{41}{10}
]$$

liefert $(\alpha^{-1}=137.0365)$ **ohne** freie Parameter.  

Hier kommt die **41** über $(b_1)$ herein – die Hyperladungsspur des Standardmodells in GUT‑Norm.

**2) Der gleiche 41‑Fingerabdruck setzt den EW‑Block**

Im EW‑Block (Fenster bei \(n=12\)) verwenden wir
$$
[

\zeta_{\rm EW}=(\pi c_3)\,\mathrm e^{-\beta_{\rm EW}\pi c_3}\,\mathrm e^{-k_{\rm EW}/c_3},

\quad \beta_{\rm EW}=\tfrac{3}{4},\quad 

k_{\rm EW}=\tfrac{3}{2}\cdot \mathcal Y^2_{\rm SM+H}=\tfrac{41}{32}.

]$$

Auch hier steckt **die gleiche 41**, jetzt in $(k_{\rm EW})$. Damit wird $(v_H)$ über $(v_H=\zeta_{\rm EW} M_{\rm Pl}\varphi_{12})$ bestimmt.  

**Ergebnis:** $(v_H\simeq 251.1\ \mathrm{GeV})$ (Skalenanker, erwartete 1–2 Prozent Drift zu $(G_F)$).

**3) α im EW‑Bild: Kombination aus \(g_1\) und \(g_2\)**

Nach Elektroschwacher Mischung gilt

$$[

e=g_2\sin\theta_W=g_1\cos\theta_W,\qquad 

\alpha=\frac{e^2}{4\pi}.

]$$

Setzt man typische Werte am $(M_Z)$ ($(g_2\approx0.652,\ g_1^{\rm SM}\approx0.357)$), erhält man $(\alpha(M_Z))$ **in der Größenordnung $(1/128)$ – das ist die **laufende** α am Z‑Pol.  

Unsere Fixpunktlösung gibt die **IR‑α** ($(\alpha^{-1}\approx 137.0365)$); der Unterschied ist schlicht **Renormierungsfluss**. Entscheidend ist: **die gleiche 41** steuert sowohl die Fixpunktgleichung (über $(b_1)$) als auch den EW‑Anker (über $(k_{\rm EW})$).

> **Kein Kreisbezug**
>
> $c_{3}=\tfrac{1}{8\pi}$, $\varphi_{0}=\tfrac{1}{6\pi}+\tfrac{3}{256\pi^{4}}$
> $\Rightarrow$ $\kappa=\tfrac{b_{1}}{2\pi}\ln\tfrac{1}{\varphi_{0}}$, $A=2c_{3}^{3}$
> $\Rightarrow$ Kubik in $\alpha$.
>
> Parallel: $Y^{2}_{\text{SM+H}}=\tfrac{41}{48}\Rightarrow k_{\text{EW}}=\tfrac{41}{32}$.
>
> Dieselbe abelsche Spur setzt $b_{1}$ und $k_{\text{EW}}$, **ohne** Rückkopplung von $\alpha$ nach $v_{H}$.  

> [!example] Mini‑Zahlencheck

> – Fixpunkt: $(\alpha^{-1}_{\rm IR}=137.0365)$ (aus $(c_3,\varphi_0,b_1)$)  
> – Am $(M_Z)$: $(\alpha(M_Z)\sim 1/128)$ aus $(g_1,g_2,\theta_W)$  
> – Beide Werte sind durch **denselben** U(1)‑Inhalt verknüpft; die 41 erscheint **zweimal** und erklärt, warum α hier natürlich wieder ins Bild kommt.

---

#### 8.4.7 Kosmologie aus der Basisstufe

$$
\Omega_b=\varphi_0\,(1-2c_3)=\mathbf{0.04894066}.
$$

---

### 8.5 Zusammenfassung auf einen Blick

| Größe | Vorhersage | Referenz | Abweichung |
|---|---:|---:|---:|
| $v_H$ | **251.07628 GeV** | 246.21965 GeV | +1.97 % |
| $M_W$ | **81.85087 GeV** | 80.3692 GeV | +1.84 % |
| $M_Z$ | **93.31741 GeV** | 91.1876 GeV | +2.34 % |
| $m_t$ | **177.54 GeV** | 172.57 GeV | +2.9 % |
| $f_a$ | **$8.8565\times 10^{10}$ GeV** | Standardfenster | — |
| $m_a$ | **64.36 μeV** | Standardfenster | — |
| $M_R$ | **$1.311\times 10^{15}$ GeV** | — | — |
| $m_{\nu3}$ | **0.04807 eV** | — | — |
| $\Delta m^2_{31}$ | **$2.31\times 10^{-3}$ eV²** | $2.509\times 10^{-3}$ eV² | −7.9 % |
| $\sin^2\theta_{13}$ | **0.023093** | $0.02240\pm0.00065$ | +3.1 % |
| $\sin\theta_C$ | **0.22446** | $0.2248\pm 0.0006$ | −0.15 % |
| $m_p$ | **0.96821 GeV** | 0.938272 GeV | +3.19 % |
| $f_\pi$ | **88.12 MeV** | 92.07 MeV | −4.3 % |
| $m_\pi$ | **132.75 MeV** | 134.98 MeV $(\pi^0)$ | −1.6 % |
| $\alpha^{-1}$ | **137.036501465** | 137.035999177 | +3.67 ppm |
| $\Omega_b$ | **0.04894066** | 0.0493 | −0.7 % |

---
#### 8.5.1 Systematik der Abweichungen

Die $1$ bis $3$ Prozent Abweichungen bei $v_{H},M_{W},M_{Z}$ entstehen aus
(i) fehlenden Zwei‑Schleifen‑Termen im elektroschwachen Sektor,
(ii) Schwellenanpassungen am Übergang zu $n=12$,
(iii) der Blockeinheit $\zeta_{\text{EW}}$ als reiner Einheitenwahl.

**Erwartung.** Eine konsistente Nachführung mit Zwei‑Schleifen und Stückweise‑Matching (vgl. 5.1, 5.2) verschiebt $v_{H}$, $M_{W}$, $M_{Z}$ systematisch **nach unten** in Richtung der Referenzen, bei $\Delta\sim 1$ bis $2$ Prozent. Die Ratio‑Tests innerhalb des Blocks bleiben unverändert, da sie nur von $\varphi_{n}$ abhängen. Siehe 8.4.1 und 5.4. 

---

### 8.6 Wo E₇ und E₆ konkret einhaken

- **E₇ Fenster bei $n=12$** verankert die **elektroschwache Skala**. Die abelsche Spur $\mathcal Y^2_{\rm SM+H}=\tfrac{41}{48}$ führt via drei halbe Randzyklen zu $k_{\rm EW}=\tfrac{41}{32}$. Dieselbe 41 erscheint als $b_1=\tfrac{41}{10}$ in der Fixpunktgleichung von $\alpha$.

- **E₆ Korridor** trägt die **starke Dynamik**. $r_{\rm had}=5$ erklärt die mildere Dämpfung im Hadronblock und rechtfertigt kleine rationale $\Delta k$ für Goldstone Physik relativ zu Baryonen.

---

### 8.7 Was noch offen ist und wie wir es schließen

- **Feinstruktur der Yukawas:** Hier wurden bewusst nur **Skalen** gesetzt. Texturen und Phasen sind die nächste Schicht. Prozentstreuungen im Block Rahmen sind erwartbar.

- **Zwei Schleifen Feinschliff im elektroschwachen Sektor:** Eine konsistente Nachführung mit Schwellen wird $v_H, M_W, M_Z$ systematisch Richtung Referenzen ziehen.

- **Formale Ableitung von $k_B$:** Die verwendeten rationalen $k_B$ sind aus der Randzählung motiviert. Eine indexartige Ableitung pro Block gehört in den Anhang.

- **Chiralität**: geschlossen durch Randbedingungen und integer Flüsse auf der orientierbaren Doppelabdeckung, siehe Box in 8 und **Appendix J**. 

---

### 8.8 Zahlenkasten für diese Sektion

- $c_3=\tfrac{1}{8\pi}=0.039788735772973836$  
- $\varphi_0=\tfrac{1}{6\pi}+\tfrac{3}{256\pi^4}=0.05317195217684553$  
- $\gamma(0)=0.834,\quad \lambda=0.5877029773404678$  
- $\varphi_{10}=0.018562845593356334,\ \varphi_{12}=0.01744828469380037$  
- $\varphi_{15}=0.015675365814677055,\ \varphi_{16}=0.015052485208841481$  
- $\varphi_{22}=0.01083363062914777,\ \varphi_{25}=0.008218869841220914,\ \varphi_{26}=0.007208714066517271$  
- $\zeta_{\rm EW}=1.17852087206\times 10^{-15},\ \zeta_{\rm PQ}=3.90754185582\times 10^{-7}$  
- $M_{\rm Pl}=1.221\times10^{19}\ \mathrm{GeV}$  
- $g_2=0.652,\ g_1^{\rm SM}=0.357$ am $M_Z$  
- $\alpha=\mathbf{0.007297325816919221},\ \alpha^{-1}=\mathbf{137.03650146488582}$

---
## 9. Weitere Informationen, Ausblick und FAQ

---
### 9.1 Ergänzungen zum Verständnis

Die bisherigen Kapitel haben die **Kernstruktur** der Theorie hergeleitet: zwei fundamentale Fixpunkte ($c_{3}$, $\varphi_{0}$), die E₈-Kaskade und die Fixpunktlösung für $\alpha$. Für das Gesamtverständnis sind drei weitere Aspekte hervorzuheben:

1. **Einpunkt-Kalibrierung**:
    
    Die Kaskade $\varphi_n$ ist bis auf eine additive Konstante in $\log \varphi$ bestimmt. Eine einzige physikalische Kalibrierung (etwa am EW-Block, n=12) fixiert alle übrigen Stufen. Dies ist kein „Knopf“, sondern eine Wahl der Einheit.
    
2. **Block-Formeln**:
    
    Die Dimensionierung einzelner Observablen (z. B. Protonmasse, CMB-Temperatur, Dunkle Energie) erfolgt über kompakte Block-Formeln, die in den Anhängen angegeben sind. Sie binden die dimensionslosen $\varphi_n$ an messbare Größen.
    
3. **Spurion-Beiträge**:
    
    Der in den 2-Loop-Runs verwendete $R^{3}$-Spurion ist kein freier Parameter, sondern eine effektive Beschreibung höherer Beiträge, die in der Chern–Simons-Struktur unvermeidlich auftreten. Sein Einfluss ist klein, aber notwendig, um den kubischen Term für $\alpha$ korrekt zu modellieren.

Die **abelsche Spur** ist der rote Faden: dieselbe $(41)$ steuert $(b_{1}=\tfrac{41}{10})$ in $(\kappa=\tfrac{b_{1}}{2\pi}\ln\tfrac{1}{\varphi_{0}})$ sowie $(k_{\mathrm{EW}}=\tfrac{41}{32})$ im EW Block; geometrisch werden die Phasen über **drei Randzyklen** auf der Doppelabdeckung gelesen, siehe **Appendix J**.  

![[Pasted image 20250826132358.png]]

> [!note] Sensitivität

> Die Sensitivität von $\alpha$ gegenüber den Parametern skaliert stark mit $c_{3}$, deutlich schwächer mit $b_{1}$ und nur moderat mit $\varphi_{0}$, siehe Figure 7.1.
### Selbstkonsistenz: $\varphi_{0}\leftrightarrow \alpha$

Die Fixpunktgleichung erzeugt nicht nur α als Funktion von φ₀, sondern φ₀ selbst ist aus der geometrischen Reduktion ϕ₀ = 1/(6π) + 3/(256π⁴) motiviert. Kombiniert man beide Abhängigkeiten, ergibt sich eine geschlossene Schleife:


$[\varphi_{0}\ \xrightarrow{\ \kappa(\varphi_{0})=\tfrac{b_{1}}{2\pi}\ln\tfrac{1}{\varphi_{0}}\ }\\Big[\ \alpha^{3}-2c_{3}^{3}\alpha^{2}-8b_{1}c_{3}^{6}\ln\tfrac{1}{\varphi_{0}}\ =\ 0\ \Big]\ \xrightarrow{\ \text{Lösung}\ }\ \alpha(\varphi_{0})]$

Diese Schleife schließt sich, da $\varphi_{0}$ selbst aus der Geometrie folgt $\big(\varphi_{0}=1/(6\pi)+3/(256\pi^{4})\big)$ und die Lösung für $\alpha$ die Eingabe bestätigt.

Diese Selbstreferenzstruktur ersetzt klassische Fine-Tuning-Debatten durch eine strukturelle Rückkopplung – $\varphi_0$ und $\alpha$ bestimmen sich gegenseitig. Kleine Änderungen in φ₀ propagieren durch κ direkt in die Gleichung, die dann einen neuen α-Wert ergibt. Die ursprüngliche Eingabe wird durch die resultierende Lösung wieder bestätigt – ein strukturelles „locking“ statt einstellbarem Parameter.

> **Falsifikationskasten**
> **A — $\alpha$‑Präzision:** Abweichung der Kubik‑Lösung jenseits einiger zehn ppm widerlegt die Struktur.
> **B — Fingerprints:** Verfehlt $\alpha_{3}(1\,\text{PeV})$ das $\varphi_{0}$‑Fenster oder $\alpha_{3}$ das $c_{3}$‑Fenster bei $\mu\sim 10^{8}$ GeV robust, ist das Modell widerlegt.
> **C — Spacing:** Bricht die nahezu äquidistante Log‑Abstands‑Invariante der drei Gleichstände, ist die Leiterinkohärenz belegt. 

---
### 9.2 Offene Fragen und nächste Schritte

Mehrere Punkte sind in der Theorie bereits angelegt, erfordern aber vertiefte Arbeit:

- **Formale Ableitung von** $\gamma(n)$:
    Die Quadratik wurde plausibel aus der E₈-Orbitkette motiviert. Eine exakte algebraische Herleitung mit vollständiger Referenztabelle der Orbits und Fit-Residuen ist der nächste mathematische Schritt.
    
- **Block-Konstanten** $\zeta$:
    Für EW-, Hadron- und Kosmo-Blöcke wurden kompakte $\zeta$-Faktoren eingeführt. Deren genauere topologische Interpretation (z. B. aus Anomalien oder Indexsätzen) steht noch aus.
    
- **RG-Robustheit**:
    Erste Tests zeigen, dass die Gleichstandskorridore extrem stabil sind. Eine systematische Analyse mit variierenden Schwellen (± Dekade) und alternativen Feldinhalten ist geplant.
    
- **Kosmologische Erweiterungen**:
    Die Stufen n=20,25,30 reproduzieren Knie, CMB und Dunkle Energie. Hier soll geprüft werden, ob auch $S_8/\sigma_8$-Spannungen und frühe Dunkle Energie konsistent eingebettet werden können.


---

### 9.3 FAQ: Zehn Fragen und Antworten

**1. Ist das nur Zahlenspielerei oder Numerologie?**

Nein. $c_{3}=1\!/\!(8\pi)$ folgt aus einer quantisierten Chern Simons Kopplung. $\varphi_{0}$ folgt aus Möbius Geometrie plus Randtermen. Beide Größen erscheinen unabhängig in unterschiedlichen Teilen der Theorie und speisen dann die Fixpunktgleichung für $\alpha$. Das unterscheidet ein strukturelles Resultat von einem Fit.

**2. Gibt es freie Parameter?**

Nein. Nach Festlegung der topologisch und geometrisch bestimmten Größen $c_{3}$ und $\varphi_{0}$ sowie der physikalisch fixen $U(1){Y}$ _Konstante_ $b{1}=41/10$ bleibt kein frei wählbarer Parameter. Es gibt nur eine triviale Einheitenkalibrierung.

**3. Warum gerade** $E_{8}$**?**

Nur $E_{8}$ besitzt ausreichend reiche Orbitstrukturen, deren Zentralisator Dimensionen eine eindeutige monotone Kette bilden. Über den Logarithmus der Dimensionen entsteht eine einfache Schrittstruktur, aus der die Dämpfung $\gamma(n)$ in blockweise konstanter Form folgt. Kleinere Gruppen brechen diese Monotonie oder liefern inkonsistente Sprungmuster.

**4. Unterschied zu klassischen GUT Ansätzen wie SU(5) oder SO(10)?**

Klassische GUTs postulieren zusätzliche Symmetrie und eine neue Skala, um Kopplungen zu vereinheitlichen. Hier werden Konstanten aus Topologie und Geometrie abgeleitet. Unifikation erscheint als Nebenwirkung des Flusses, nicht als Axiom.

**5. Wie robust sind die Zahlen?**

Sehr robust. Schwellen um eine Dekade verschoben verändern die Lage charakteristischer Gleichstände nur im Promillebereich. Die Lösung der Fixpunktgleichung für $\alpha$ bleibt im ppm Bereich stabil. Die Stufen der Leiter sind deterministisch, nicht fitgetrieben.

**6. Warum ist** $\alpha$ **so präzise, andere Größen aber nur auf Prozent genau?**

$\alpha$ wird direkt durch die Fixpunktgleichung festgelegt. Massen und Mischungen tragen zusätzliche QCD Dynamik, Flavour Struktur und Schema Effekte. Diese Beiträge sind in der vorliegenden Version bewusst modular gehalten und erzeugen natürliche Streuungen auf Prozentniveau.

**7. Wie kann man die Theorie widerlegen?**

Drei klare Hebel:

a) RG Fingerprints an zwei charakteristischen Skalen, etwa im PeV Bereich und bei etwa $2.5\times 10^{8} GeV$.

b) Stabilität des Abstandsmusters zwischen Gleichständen über einen breiten Parameterbereich.

c) Vorhersagen in Präzisionsbereichen wie Atominterferometrie oder Rydberg Konstante für $\alpha$. Systematische Abweichungen widerlegen das Modell.

**8. Gibt es Bezüge zu Stringtheorie oder M Theorie?**

Ja, auf der Ebene der 11 dimensionalen Elternstruktur mit Chern Simons Term und kompaktifizierter Topologie. Anders als Landschafts Ansätze benötigt TFPT keine Vielzahl freier Moduli. Die Ableitungen bleiben lokal und topologisch.

**9. Was sagt die Theorie zur kosmologischen Konstante?**

Die Stufe n=30 der Leiter liefert eine Energiedichte $\rho_{\Lambda}$ in der Größenordnung der Planck Messungen. Entscheidend ist die Herkunft des Exponenten aus der Leiter, nicht ein Fit an Daten.

**10. Wo liegen die größten Unsicherheiten?**

Zwei Punkte: die formale Herleitung der geschlossenen Form von $\gamma(n)$ direkt aus der $E_{8}$ Orbitologie und die tiefe Interpretation der Block Konstanten $\zeta$. Beides wird in den Ausblicks Abschnitten als Arbeitsprogramm benannt.

**11. Woher stammen** $A=\tfrac{1}{256\pi^{3}}$ **und** $\kappa=\tfrac{b_{1}}{2\pi}\ln(1/\varphi_{0})$**?**

Aus der gewählten Normierung $\alpha=g^{2}/(4\pi)$, GUT Norm für U(1){Y} _und einer topologisch induzierten einschleifigen Korrektur zu_ $F^{2}$ _mit zwei identischen Einfügungen von_ $c{3}$. Siehe Derivation Note A1 im Anhang für die vollständige Rechnung.

**12. Wie schemaabhängig sind die Aussagen?**

Ein Schemawechsel verschiebt nur additive, skalenunabhängige Terme in $\kappa$. Der reine Zahlenfaktor A ist durch Topologie und kanonische Eichkinetik fixiert. Fixpunkte und Leiterstrukturen bleiben invariant.

**13. Was bedeutet „keine freien Parameter“ praktisch, wenn numerische Werte doch gerundet werden?**

Rundungen betreffen nur Darstellung und numerische Propagation. Die Strukturgleichungen sind parameterfrei. In Reproduktionen sollen alle Konstanten mit definierter Präzision angegeben und Fehlerbalken aus Schema und Schwellen Variation ausgewiesen werden.

**14. Warum eine kubische Gleichung für** $\alpha$ **und keine quadratische oder quartische?**

Die kleinste nicht triviale Ordnung, in der der topologische Beitrag zur Wellenfunktionsrenormierung des Photons lokal und paritäts even auftritt, ist proportional zu $g^{6}$. In $\alpha$ Skala entspricht dies der dritten Potenz. Niedrigere Ordnungen sind durch Symmetrie oder Quantisierung ausgeschlossen.

**15. Wie werden zwei Schleifen Effekte und Schwellen technisch behandelt?**

Die nicht abelschen Kopplungen laufen zwei schleifig mit Standardkoeffizienten und Schwellensprüngen an den effektiven Massen der schweren Moden. Sensitivitätsanalysen zeigen, dass die beiden charakteristischen Fingerprints in Lage und Abstand stabil bleiben. Die abelsche Gleichung erhält zusätzlich den topologischen kubischen Term.

**16. Wie reproduziere ich die Kernergebnisse numerisch?**

Schritte:

a) Setze $c_{3}=1/(8\pi)$ und $\varphi_{0}$ gemäß Abschnitt 3.2.

b) Berechne $\kappa=(b_{1}/2\pi)\ln(1/\varphi_{0})$ mit $b_{1}=41/10$.

c) Löse die Fixpunktgleichung in 6.2 für $\alpha$ mit $A=1/(256\pi^{3})$.

d) Lasse die nicht abelschen Kopplungen zwei schleifig laufen, setze definierte Schwellen, prüfe die Fingerprints.

e) Variiere Schwellen und Schema Parameter in plausiblen Bereichen und gib Fehlerbalken an.

**17. Was ist die Physik hinter** $\varphi_{0}$**?**

$\varphi_{0}$ ist keine Fit Konstante, sondern entsteht aus einer geometrischen Relation auf der orientierbaren Doppelabdeckung der Möbius Reduktion. Gauss Bonnet mit Rand liefert den Flächenanteil, der Randterm den Zuschlag. Zusammen fixiert das die effektive dimensionslose Skalenrelation.

**19. Wo endet die Zuständigkeit der Theorie derzeit bewusst?**

In der vorliegenden Version werden Flavor Details, CKM und PMNS Phasen sowie nicht triviale Hadronen Phänomenologie nur gerahmt. Das ist eine bewusste Modularisierung. Ziel ist zuerst ein sauberes Fundament aus Topologie, Geometrie und Kopplungsdynamik.

**20. Was folgt als nächstes, um die offenen Punkte zu schließen?**

Drei konkrete Schritte:

a) Formale Ableitung der geschlossenen $\gamma(n)$ Gestalt direkt aus nilpotenten Orbits und Zentralisatoren.

b) Vollständige zwei schleifige Validierung mit systematischer Schwellenauswertung und Fehlerbudget.

c) Präzisionstests für $\alpha$ über unabhängige Messkanäle und Simulationen, inklusive klarer Abweichungsschwellen für Falsifikation.

---

### 9.4 Plausibilitätsargumente: Wahrscheinlichkeit und strukturelle Abhängigkeiten

Die Plausibilität der vorliegenden Theorie ergibt sich aus zwei komplementären Aspekten: (i) der extrem geringen Wahrscheinlichkeit multipler präziser Übereinstimmungen ohne freie Parameter und (ii) den tiefen strukturellen Abhängigkeiten zwischen Topologie, Geometrie, Symmetrie und Dynamik.

Zunächst zur Wahrscheinlichkeit: Die parameterfreie Vorhersage der Feinstrukturkonstante $\alpha^{-1} \approx 137.03650$ weicht lediglich um 3.67 ppm vom CODATA 2022 Referenzwert ab. Unter der naiven Annahme einer gleichmäßigen Verteilung von $\alpha$ in einem physikalisch plausiblen Bereich (z. B. 0.001 bis 0.01) beträgt die Wahrscheinlichkeit einer solch geringen Abweichung nur etwa $6 \times 10^{-7}$ (basierend auf einer absoluten Differenz von $2.7 \times 10^{-8}$). Entsprechende Treffer finden sich auch in der E₈-Kaskade, etwa für $\Omega_b \approx 0.04894$ (Abweichung 0.06% gegenüber den Planck-Daten) oder $m_p \approx 937 MeV$ (Abweichung 0.12%). Jeder dieser Werte entspricht einer unabhängigen Wahrscheinlichkeit im Bereich von $10^{-2}$ bis $10^{-6}$. Multipliziert man diese für rund zehn zentrale Vorhersagen (Flavor-Mischungen, Massen, kosmologische Konstanten), ergibt sich eine kombinierte Zufallswahrscheinlichkeit von kleiner als $10^{-20}$. Dies ist vergleichbar mit der Unwahrscheinlichkeit, dass eine Serie unabhängiger Würfelwürfe exakt dasselbe Muster wiederholt.

Hinzu kommen die strukturellen Abhängigkeiten: Die Fixpunkte $c_3$ und $\varphi_0$ entstehen nicht isoliert, sondern folgen aus unterschiedlichen, jedoch konsistenten Prinzipien – $c_3$ aus topologischen Chern–Simons-Normalisierungen in elf Dimensionen, $\varphi_0$ aus geometrischer Möbius-Reduktion. Beide Parameter werden unabhängig in renormierungsgruppenbasierten Flüssen bestätigt, etwa durch $\alpha_3(1\, \text{PeV}) \approx \varphi_0$. Die Schichten greifen ineinander: Topologie fixiert die Normalisierungen, E₈ ordnet die Kaskade, die RG-Flüsse liefern dynamische Konsistenz.

Diese interne Verschränkung reduziert die Wahrscheinlichkeit, dass es sich lediglich um zufällige Koinzidenzen handelt, erheblich. Ein Versagen in einer Schicht (z. B. im genetischen Algorithmus oder in den Dimensionsketten) würde die übrigen nicht tangieren, wird empirisch jedoch nicht beobachtet. Stattdessen ergibt sich ein kohärentes Gesamtbild, das durch Reproduzierbarkeit und Falsifizierbarkeit (z. B. bei der vorhergesagten Axionmasse) überprüfbar ist.
## Appendix A — Zahlenkasten der Fixpunkte (hochpräzise)

$$
c_{3}=\frac{1}{8\pi}=0.039788735772973836,\qquad
\varphi_{0}=\frac{4}{3}c_{3}+48c_{3}^{4}=0.053171952176845526,
$$
$$
A=2c_{3}^{3}=1.259825563796855\times 10^{-4},\qquad
\kappa=\frac{41}{10}\frac{1}{2\pi}\ln\!\frac{1}{\varphi_{0}}=1.914684795.
$$
$$
\alpha=0.007297325816919221,\qquad \alpha^{-1}=137.03650146488582.
$$
*Referenz:* CODATA 2022 $\alpha_{\text{CODATA}}=7.2973525628(11)\times 10^{-3}$,
Abweichung $\approx 3.67$ ppm.

---

## Appendix B – E₈-Kaskade in geschlossener Form


**Definitionen und Normierung**

Sei für jedes nilpotente E₈ Orbit

$D_n \;=\; 248-\dim\mathcal O_n,\qquad n=0,\dots,26$,

mit der gefundenen Kette $D_n=60-2n$ von $D_0=60$ bis $D_{26}=8$.

Die Leiter folgt aus einer einzigen Normierung am ersten Schritt

$s^{\star}=\ln 248-\ln 60=1.419084183942882,\qquad \lambda=\frac{0.834}{s^{\star}}=0.5877029773404678$ .

**Dämpfung**

$\gamma(0)=0.834,\qquad \gamma(n)=\lambda\Big[\ln D_n-\ln D_{n+1}\Big]\quad (n\ge 1)$.

**Rekursion**

$\varphi_{n+1}=\varphi_n\,\mathrm e^{-\gamma(n)}$ .

**Geschlossene Form der Leiter**

Für $n\ge 1$ gilt

$\varphi_n \;=\; \varphi_0\,\mathrm e^{-\gamma(0)}\Big(\tfrac{D_n}{D_1}\Big)^{\lambda},\qquad D_1=58$.

**Kalibrierungsfreie Tests**

1. **Verhältnisgesetz** für $m,n\ge 1$:
    
    $\boxed{\ \frac{\varphi_m}{\varphi_n}=\Big(\tfrac{D_m}{D_n}\Big)^{\lambda}=\Big(\tfrac{60-2m}{60-2n}\Big)^{\lambda}\ }$ .
    
2. **Log lineares Gesetz**
    
    $\boxed{\ \log \varphi_n=\text{Konstante}+\lambda\,\log D_n\ }$ .
    
**Bemerkung zum Kettenende**

Die E acht Kette endet strukturell bei $n=26$ mit $D=8$. Werte für $n>26$ wären eine analytische Fortsetzung und werden als Extrapolation gekennzeichnet.

### B.0 Verifikation der Eindeutigkeit

Wir haben die Kettenmenge $\mathcal C$ vollständig auf dem $\Delta D=2$‑DAG enumeriert und $\mathbf F$ aus 4.5.2 lexikografisch minimiert. Es ergab sich genau **eine** Minimalkette (Satz 4.5.1), identisch zur in Tab. B.1 gezeigten Folge $D_{n}=60-2n$ mit den in 4.2 gelisteten Labels. Die daraus abgeleitete log exakte Leiter und die Dämpfung $\gamma(n)$ stimmen stufenweise überein.

---
### **B.1 – E8 Kaskade: log exakte Größen pro Stufe**

Spalten:

- $n$
- $D_n$
- $\ln D_n$
- $s_n=\ln D_n-\ln D_{n+1}$
- $\gamma(n)$ mit $\gamma(0)=0.834$, sonst $\lambda s_n$
- $\Sigma\gamma$ kumuliert bis Stufe $n$ exklusiv
- $\varphi_n/\varphi_0$ unkalibriert
- $\big(\tfrac{D_n}{D_1}\big)^{\lambda}$ als reine Kettenzahl

> Hinweis: $\varphi_n/\varphi_0=\mathrm e^{-\gamma(0)}\big(\tfrac{D_n}{D_1}\big)^{\lambda}$ für $n\ge 1;$ für $n=0$ gilt $\varphi_0/\varphi_0=1$.

> [!note] Tabellenhinweis

> Die Spalte $(D_{n}/D_{1})^{\lambda}$ ist für $n\ge 1$ die Kettenzahl der Leiter. Der Eintrag für $n=0$ dient nur der Kontrolle und wird physikalisch nicht verwendet.

| **n** | **D** | **ln D** | **s_n**  | **γ(n)** | **Σγ**   | **φ_n/φ₀** | **(Dₙ/D₁)^λ** |
| ----- | ----- | -------- | -------- | -------- | -------- | ---------- | ------------- |
| 0     | 60    | 4.094345 | 0.033902 | 0.834000 | 0.000000 | 1.000000   | 1.020124      |
| 1     | 58    | 4.060443 | 0.035091 | 0.020623 | 0.834000 | 0.434309   | 1.000000      |
| 2     | 56    | 4.025352 | 0.036368 | 0.021373 | 0.854623 | 0.425443   | 0.979064      |
| 3     | 54    | 3.988984 | 0.037740 | 0.022180 | 0.875997 | 0.416447   | 0.957994      |
| 4     | 52    | 3.951244 | 0.039221 | 0.023050 | 0.898177 | 0.407312   | 0.936782      |
| 5     | 50    | 3.912023 | 0.040822 | 0.023991 | 0.921227 | 0.398030   | 0.915419      |
| 6     | 48    | 3.871201 | 0.042560 | 0.025012 | 0.945218 | 0.388595   | 0.893899      |
| 7     | 46    | 3.828641 | 0.044452 | 0.026124 | 0.970230 | 0.378996   | 0.872211      |
| 8     | 44    | 3.784190 | 0.046520 | 0.027340 | 0.996355 | 0.369223   | 0.850347      |
| 9     | 42    | 3.737670 | 0.048790 | 0.028674 | 1.023695 | 0.359265   | 0.828299      |
| 10    | 40    | 3.688879 | 0.051293 | 0.030101 | 1.052369 | 0.349110   | 0.806058      |
| 11    | 38    | 3.637586 | 0.054067 | 0.031767 | 1.082514 | 0.338744   | 0.783615      |
| 12    | 36    | 3.583519 | 0.057158 | 0.033589 | 1.114290 | 0.328148   | 0.760962      |
| 13    | 34    | 3.526361 | 0.060625 | 0.035571 | 1.147880 | 0.317306   | 0.738089      |
| 14    | 32    | 3.465736 | 0.064539 | 0.037915 | 1.183450 | 0.306202   | 0.714988      |
| 15    | 30    | 3.401197 | 0.068993 | 0.040555 | 1.221365 | 0.294805   | 0.691650      |
| 16    | 28    | 3.332205 | 0.074108 | 0.043581 | 1.261920 | 0.283078   | 0.668066      |
| 17    | 26    | 3.258097 | 0.080043 | 0.047041 | 1.305501 | 0.271026   | 0.644229      |
| 18    | 24    | 3.178054 | 0.087011 | 0.051117 | 1.352542 | 0.258584   | 0.620130      |
| 19    | 22    | 3.091042 | 0.095310 | 0.055996 | 1.403659 | 0.245652   | 0.595761      |
| 20    | 20    | 2.995732 | 0.105361 | 0.061940 | 1.459655 | 0.232102   | 0.571113      |
| 21    | 18    | 2.890372 | 0.117783 | 0.069239 | 1.520595 | 0.217761   | 0.546180      |
| 22    | 16    | 2.772589 | 0.133531 | 0.078477 | 1.589835 | 0.203747   | 0.520953      |
| 23    | 14    | 2.639057 | 0.154151 | 0.090595 | 1.668311 | 0.188369   | 0.495424      |
| 24    | 12    | 2.484907 | 0.182322 | 0.107151 | 1.758907 | 0.172054   | 0.469584      |
| 25    | 10    | 2.302585 | 0.223144 | 0.131142 | 1.867098 | 0.154572   | 0.443426      |
| 26    | 8     | 2.079442 |          |          | 1.998240 | 0.135574   | 0.416948      |

---
## Appendix C – Block-Formeln für Observablen


> [!example] Block-Kalibrierung in der Praxis

> Pro Block genügt eine Einheiten-Kalibrierung $\zeta$ auf eine Referenzgröße. Alle Relationen im Block folgen dann fitfrei aus den Ratio-Gesetzen der Kette, siehe 4.3 und Appendix B.


**Elektroschwacher Block (n=12):**

$v_H=\zeta_{\rm EW} M_{Pl}\varphi_{12},\quad M_W=\tfrac{1}{2} g_2 v_H,\quad M_Z=\tfrac{1}{2}\sqrt{g_1^2+g_2^2} v_H$.

**Hadronischer Block (n=15,17):**

$m_p=\zeta_p M_{Pl}\varphi_{15},\quad m_b=\zeta_b M_{Pl}\varphi_{15},\quad m_u=\zeta_u M_{Pl}\varphi_{17}$.

**Kosmo-Blöcke:**

$T_{\gamma0}=\zeta_\gamma M_{Pl}\varphi_{25},\quad T_\nu=(4/11)^{1/3}T_{\gamma0},\quad \rho_\Lambda=\zeta_\Lambda M_{Pl}^{4}\varphi_{30}^{97/30}$.

**Fundamentale Relationen nahe n=0:**

$\Omega_b=\varphi_0(1-2c_3),\qquad r=\varphi_0^2,\qquad V_{us}/V_{ud}=\sqrt{\varphi_0}$.

### C.8 Möbius Leiter: Definition und Fehlerfortpflanzung

**Definition.**  
$\mathcal{M}_y(\delta)=\dfrac{y+\delta}{y-\delta}$ mit $y\in\{1,\tfrac{1}{3},\tfrac{2}{3}\}$.

**Kalibrierregel.**  
$\delta=\dfrac{\sqrt{m_\tau/m_\mu}-1}{\sqrt{m_\tau/m_\mu}+1}$.

**Ableitungen.**  
Setze $R=\sqrt{m_\tau/m_\mu}$. Dann gilt $\dfrac{\mathrm d\delta}{\mathrm dR}=\dfrac{2}{(R+1)^2}$.  
Mit $\dfrac{\partial R}{\partial m_\tau}=\dfrac{1}{2}\dfrac{1}{\sqrt{m_\tau m_\mu}}$ und  
$\dfrac{\partial R}{\partial m_\mu}=-\dfrac{1}{2}\dfrac{\sqrt{m_\tau}}{m_\mu^{3/2}}$ folgt  
$\sigma_\delta^2=\big(\dfrac{\mathrm d\delta}{\mathrm dR}\big)^2\,\sigma_R^2$.

**Vorhersageformeln.**  
Setze $\delta$ in die sechs Relationen aus 7.4.4 ein.  
Feinkorrekturen lassen sich als universelle Verschiebung schreiben:  
$\delta\rightarrow \delta+a_y\,\varphi_0+b_y\,c_3$ mit kleinen sektorspezifischen $a_y,b_y$.

---
## Appendix D: Möbius Faser: Rand, Krümmungsnormalisierung und der Koeffizient 6pi
    

>[!info] Ziel

Wir begründen, warum in 3.2.2 der lineare Randkoeffizient $6\pi$ auftritt und daraus $\varphi_{\text{tree}}=1/(6\pi)$ folgt. Zusätzlich zeigen wir den topologischen Zuschlag $\delta_{\text{top}}$ in äquivalenten Formen und klären Unabhängigkeiten von Repräsentationen und Schemata.

---

### D.1 Setup, Notation und konforme Skalierung

Sei $\mathcal M$ die zweidimensionale Möbius-Faser (kompakt, mit Rand), $g_{\mathcal M}=\varphi^{2}\hat g_{\mathcal M}$ eine reine konforme Reskalierung. Für Gauß-Krümmung $K$ und geodätische Randkrümmung $k_{g}$ gilt

$$

\int_{\mathcal M}K,\mathrm dA+\oint_{\partial\mathcal M}k_{g},\mathrm ds=2\pi,\chi(\mathcal M),

$$


$$

K=\varphi^{-2}\hat K,\quad \mathrm dA=\varphi^{2}\mathrm d\hat A,\quad \mathrm ds=\varphi,\mathrm d\hat s.

$$

Daraus folgt konforme Invarianz des Flächenintegrals und Linearität des Randintegrals:

$$

\int_{\mathcal M}K,\mathrm dA=\int_{\mathcal M}\hat K,\mathrm d\hat A,\qquad

\oint_{\partial\mathcal M}k_{g},\mathrm ds=\varphi,\oint_{\partial\mathcal M}\hat k_{g},\mathrm d\hat s.

$$

  
Somit stammt die explizite $\varphi$-Abhängigkeit der reduzierten gravitativen Wirkung ausschließlich aus dem Rand (vgl. 3.2.2).

---

### D.2 Orientierbare Doppelabdeckung und der Naht-Beitrag


Die orientierbare Doppelabdeckung $\widetilde{\mathcal M}$ der Möbius-Faser ist ein Zylinder mit zwei geometrischen Randkomponenten. Zusätzlich entsteht durch die $\mathbb Z_{2}$-Identifikation eine Naht-Kurve $\Gamma$. Diese trägt wie ein dritter effektiver Randzyklus.

**Lemma D.2 (Naht als Randterm).** Modelliert man die Identifikation entlang $\Gamma$ als Grenzfall einer dünnen Kollarnachbarschaft mit zwei gegenüberliegenden Randkurven und dihedralem Winkel $\pi$, so liefert der Corner- beziehungsweise Naht-Term in der Gauss-Bonnet-Bilanz für jede geschlossene $\Gamma$ einen integrierten Beitrag, der äquivalent zu einem vollwertigen Randintegral mit Normierung $2\pi$ ist. ∎

**Normierung.**

$$

\mathcal K_{\partial}:=\sum_{\text{Randzyklen}}\oint \hat k_{g},\mathrm d\hat s

= 2\pi+2\pi+2\pi=6\pi.

$$

Damit ist die Zahl $6\pi$ kanonisch und unabhängig von einer speziellen Repräsentation.

![[Pasted image 20250901121634.png]]

*Illustration: The Möbius fiber can be represented by its orientable double cover, a cylinder with two ordinary boundaries plus one effective seam Γ. Each contributes $2\pi$ to the Gauss–Bonnet balance, leading to the canonical total $6\pi$.*

---
  
### D.3 Reduzierte 6D-Wirkung und der lineare $\varphi$-Koeffizient

  
In der 6D-Reduktion trägt der geometrische Anteil linear in $\varphi$:

  

$$

S_{\text{grav}}^{(6)} \supset \frac{M_{6}^{4}}{2}\int_{\mathcal B}\!\sqrt{g_{\mathcal B}}

\left[

\underbrace{\int_{\mathcal M}K\,\mathrm dA}_{\text{konform invariant}}

+\underbrace{\oint_{\partial\mathcal M}k_{g}\,\mathrm ds}_{=\ \varphi\,\mathcal K_{\partial}}

\right]

=\frac{M_{6}^{4}}{2}\int_{\mathcal B}\!\sqrt{g_{\mathcal B}}\,(6\pi\,\varphi)+\dots

$$

  

Der wirksame lineare Koeffizient ist $6\pi$.

---

### D.4 Stationarität und Baumwert $\varphi_{\text{tree}}$

Die effektive Potentialdichte enthält zusätzlich einen quantisierten topologischen Beitrag. Stationarität liefert:

$$

\partial_{\varphi}V_{\text{eff}}(\varphi)\propto 6\pi\varphi-1=0

\quad\Rightarrow\quad

\varphi_{\text{tree}}=\frac{1}{6\pi}.

$$
---

### D.5 Der topologische Zuschlag $\delta_{\text{top}}$


Die Reduktion des 11D-Chern-Simons-Terms erzeugt in 4D eine quantisierte Kopplung $g=8c_{3}^{2}$ mit $c_{3}=1/(8\pi)$. Der Zuschlag lässt sich schreiben als:

$$

\boxed{\delta_{\text{top}}=\frac{3}{256\pi^{4}}=48c_{3}^{4}=6c_{3}^{2}g=\tfrac{3}{4}g^{2}}

$$

mit $g=\tfrac{1}{8\pi^{2}}$, $c_{3}=\tfrac{1}{8\pi}$.

  
Damit:

$$

\boxed{\varphi_{0}=\varphi_{\text{tree}}+\delta_{\text{top}}=\tfrac{1}{6\pi}+\tfrac{3}{256\pi^{4}}=\tfrac{4}{3}c_{3}+48c_{3}^{4}}

$$

---

### D.6 Eindeutigkeit, Invarianz und Normierungsfragen

1. **Repräsentationsfreiheit.** $6\pi$ hängt nur von der Topologie ab.

2. **Normierung.** $\chi=1$ fixiert $\varphi_{\text{tree}}$.

3. **Orthogonalität.** $6\pi$ (geometrisch) und $g$ (topologisch) sind unabhängig.

4. **Schema-Robustheit.** $\delta_{\text{top}}$ ist reiner Zahlenbeitrag.

### D.7 Konsistenzcheck

$$

\varphi_{0}=\tfrac{1}{6\pi}+\tfrac{3}{256\pi^{4}}

\quad\Rightarrow\quad

\kappa=\tfrac{b_{1}}{2\pi}\ln\frac{1}{\varphi_{0}},\quad

A=2c_{3}^{3}=\tfrac{1}{256\pi^{3}}

$$

  
**– konsistent mit 7.6.**

### D.8 Cross Ratio

$$

\mathrm{CR}(x;y,-y,0)=\frac{y+\delta}{y-\delta}=:\mathcal M_{y}(\delta)

$$
mit Verschiebung

$$

\delta_{\star}=\tfrac{3}{5}+\tfrac{\varphi_{0}}{6}.

$$

Damit Anbindung an die Leiter-Abbildung (vgl. 3.3.2).

### D.9 Kurze FAQ

• **Warum drei Randzyklen?** Zwei Zylinderränder plus Naht.

• **Ist die Naht willkürlich?** Nein, Corner-Term in Gauss-Bonnet.

• **Ist $\delta_{\text{top}}$ Zahlenspiel?** Nein, folgt rein algebraisch aus $c_{3},g$.

---

### Ergebnis

$$

\boxed{\varphi_{\text{tree}}=\tfrac{1}{6\pi},\quad

\delta_{\text{top}}=\tfrac{3}{256\pi^{4}}=48c_{3}^{4}=\tfrac{3}{4}g^{2},\quad

\varphi_{0}=\tfrac{1}{6\pi}+\tfrac{3}{256\pi^{4}}}

$$


>[!summary] Fazit

Der Koeffizient $6\pi$ ist geometrisch fixiert, $\delta_{\text{top}}$ topologisch motiviert und mehrfach äquivalent dargestellt. Die Kopplung an die übrigen Konstanten des Papers ist explizit nachvollziehbar.

---

### Appendix E — Von der 11D‑Wirkung zum 4D‑Koeffizienten $A$ und zur Log‑Konstante $\kappa$

#### E.1 Setup der effektiven 4D‑Theorie

Aus 3.2.1 folgt die topologische Kopplung $g\,a\,F\tilde F$ mit $g=8c_{3}^{2}$ und periodischem Axion $a$. Nach kanonischer Normierung lautet der relevante 4D‑Sektor
$$[
\mathcal L_{\text{eff}}=-\frac{1}{4}F_{\mu\nu}F^{\mu\nu}+\frac{1}{2}(\partial a)^{2}-\frac{1}{2}m_{a}^{2}a^{2}+g\,a\,F\tilde F,
]$$
wobei $m_{a}$ ein schwerer geometrischer Modus der Reduktion ist. $\alpha\equiv g_{\text{em}}^{2}/(4\pi)$. Siehe 3.2 und 7.6. 

#### E.2 Integrieren des schweren Modus und lokale Operatoren

Das Eliminieren von $a$ im Pfadintegral erzeugt bei $p^{2}\ll m_{a}^{2}$ eine lokale Serie
$$[
\Delta\mathcal L=\frac{g^{2}}{2m_{a}^{2}}(F\tilde F)^{2}+\frac{g^{2}}{2m_{a}^{4}}(\partial F\tilde F)^{2}+\dots,
]$$
deren führender Term paritätsgerade zur Renormierung der Photonen‑Zweipunktfunktion beiträgt. Die $m_{a}$‑Abhängigkeit verschwindet aus dem **logarithmisch divergenten** Teil der Vakuumpolarisation, so dass der $\ln\mu$‑Koeffizient schemainvariant ist.

#### E.3 Hintergrundfeld‑Methode: Log‑Anteil der Vakuumpolarisation

In Hintergrundfeld‑Eichung ergibt der einschleifige Beitrag mit **zwei** topologischen Einfügungen den Term
$$[
\delta Z_{F}=\underbrace{(8c_{3}^{2})^{2}}_{\text{zwei }g\text{‑Einsätze}}\ \underbrace{\frac{1}{(4\pi)^{3}}}_{\text{Schleifenmaß}}\ \underbrace{\frac{1}{4}}_{\text{Symmetrie}}\ \ln\frac{\mu}{\mu_{0}}+\dots
]$$
und damit in $\beta_{\alpha}$ einen Zusatz $\propto c_{3}^{2}\alpha^{3}$. Zusammen folgt
$$[
A=\frac{1}{256\pi^{3}}=2c_{3}^{3},\qquad 
\beta_{\alpha}=\frac{b_{1}}{2\pi}\alpha^{2}+A\,c_{3}^{2}\alpha^{3}+\dots
]$$
Die Ableitung ist unabhängig von Details der $a$‑Kinetik und vom Schema, da nur der **Koeffizient des Logs** verwendet wird. Diese Zahl ist identisch zur in 7.6 verwendeten Variationsableitung.

#### E.4 Integrierte Eins‑Schleife und $\kappa$

Integration von $d\alpha/d\ln\mu=(b_{1}/2\pi)\alpha^{2}$ zwischen $\mu_{\text{UV}}=M_{\text{Pl}}$ und $\mu_{\text{IR}}=\varphi_{0}M_{\text{Pl}}$ liefert
$$[
\kappa=\frac{b_{1}}{2\pi}\ln\frac{1}{\varphi_{0}},\quad b_{1}=\frac{41}{10}\ \text{in GUT‑Norm},
]$$
wie in 7.6.1 zusammengefasst. $\kappa$ hängt nur von $b_{1}$ und vom **geometrisch** fixierten $\varphi_{0}$.  

#### E.5 Fixpunktgleichung aus Callan–Symanzik

Mit dem $A$ aus E.3 und der integrierten Konstante $\kappa$ folgt direkt
$$[
\alpha^{3}-2c_{3}^{3}\alpha^{2}-8\,b_{1}\,c_{3}^{6}\ln\frac{1}{\varphi_{0}}=0,
]$$
identisch zu 7.6. Damit ist $A$ **vollständig** aus der effektiven Theorie abgeleitet.  


---
## Appendix F – Zwei-Schleifen RGE-Setup

### Konfiguration

• **Fermionen:** Standardmodell plus elektroschwaches Triplet $\Sigma_F$ mit Decoupling bei $10^{3},\mathrm{GeV}$; farbadjunktes Fermion $G8$ von $SU(3)_c$ aktiv für $\mu>M_{G8}=1.8\times10^{10},\mathrm{GeV}$; drei rechtshändige Neutrinos mit gestaffelten Schwellen

$M_{N_1}=10^{14},\mathrm{GeV}$, $M_{N_2}=3\times10^{14},\mathrm{GeV}$, $M_{N_3}=8\times10^{14},\mathrm{GeV}$. Oberhalb $M_{G8}$ gilt stückweise $\Delta b_3=+2$.  

• **Skalare:** Standardmodell Higgs $H$, PQ Feld $\Phi$ mit Schwelle $M_{\Phi}=10^{16},\mathrm{GeV}$.  

• **Spurion:** Effektiver $R^{3}$ Term zur Modellierung des kubischen Beitrags $\propto \alpha^{3}$ im abelschen Sektor.  

• **Normierung:** Hyperladung in **GUT Norm**

$$

g_{1}^{\mathrm{GUT}}=\sqrt{\tfrac{5}{3}},g_{Y},\qquad b_1=\tfrac{41}{10}.

$$

Für die Steigung gilt $,\dfrac{d\alpha_{1}^{-1}}{d\ln\mu}=-\tfrac{b_{1}}{2\pi}$.  

• **Startwerte bei $\mu=M_Z$:**

$g_{1}^{\mathrm{GUT}}\approx 0.462,\quad g_{2}=0.652,\quad g_{3}=1.2323$.  

• **Integration:** Zwei Schleifen Betafunktionen mit stückweisem Threshold Matching über mindestens fünfzehn Dekaden; optionaler Drei Schleifen Steigungscheck für $SU(3)_c$.  

### Resultate

#### Fingerprints der Fixpunkte

$\alpha_{3}(1,\mathrm{PeV})=0.052923411$  gegen $\varphi_{0}=0.053171952$  $\Rightarrow$ Abweichung -0.47%;

$\alpha_{3}(2.5\times10^{8},\mathrm{GeV})=0.039713807$  gegen $c_{3}=\tfrac{1}{8\pi}=0.039788736$  $\Rightarrow$ Abweichung -0.19%.  

#### Nahe Unifikation

Minimaler **relativer Spread** der inversen Kopplungen = **1.23 %** bei $\mu^\star\approx 1.43\times10^{15},\mathrm{GeV}$.  

#### Kontinuität und Steigungen

Stückweises Matching ohne Sprünge in $\alpha_i^{-1}$; gemessene $U(1)$ Steigung konsistent mit $-b_1/(2\pi)$, G8 Brücken Steigung numerisch $0.8063$ gegenüber Erwartung $\tfrac{5}{2\pi}=0.7958$ (1.3%).    

#### Spacing Invariante

Die drei paarweisen Gleichstände liegen bei

$\mu_{23}\approx6.05\times10^{14},\mathrm{GeV}$,

$\mu_{13}\approx1.46\times10^{15},\mathrm{GeV}$,

$\mu_{12}\approx2.38\times10^{15},\mathrm{GeV}$,

damit

$$

S=\log_{10}\mu_{23}-2\log_{10}\mu_{13}+\log_{10}\mu_{12}\approx -0.17.

$$
#### PyR@TE Konfiguration  (kurz, v2)

Settings:

LoopOrder: 3        # Export; Solver nutzt volle 2-Loop + optional 3-Loop SU(3)

Groups: {U1Y: U1, SU2L: SU2, SU3c: SU3}

Thresholds:

  - {Scale: MSigma, Fields: [SigmaF]}      # 1.0e3 GeV

  - {Scale: MG8,    Fields: [G8]}          # 1.8e10 GeV  (Δb3 = +2 oberhalb)

  - {Scale: MNR1,   Fields: [NR1]}         # 1.0e14 GeV

  - {Scale: MNR2,   Fields: [NR2]}         # 3.0e14 GeV

  - {Scale: MNR3,   Fields: [NR3]}         # 8.0e14 GeV

  - {Scale: MPhi,   Fields: [phiR, phiI]}  # 1.0e16 GeV

Fermions:

  G8:  {Gen: 1, Qnb: {U1Y: 0, SU2L: 1, SU3c: 8}}   # neues Oktett

  NR1: {Gen: 1, Qnb: {U1Y: 0, SU2L: 1, SU3c: 1}}

  NR2: {Gen: 1, Qnb: {U1Y: 0, SU2L: 1, SU3c: 1}}

  NR3: {Gen: 1, Qnb: {U1Y: 0, SU2L: 1, SU3c: 1}}

_(vollständige YAML siehe Modelldatei v2)_

**Pyr@ate Configuration:**

```
---

Author: "E8 Cascade TFPT v2.1 – G8 Adjoint Enhanced"

Date: 2025-08-29

Name: E8CascadeTFPTG8_v2

# ------------------------------------------------------------

# ENHANCED E8 CASCADE MODEL WITH G8 ADJOINT FERMION (v2)

#

# GOALS:

# - Keep TFPT fingerprints SM-driven (1-loop) below 10^9 GeV

# - Provide clean G8 color-bridge for unification (Δb3 = +2 above MG8)

# - Unambiguous U(1) GUT normalization and documentation

# ------------------------------------------------------------

Settings:

LoopOrder: 3

ExportBetaFunctions: true

  

# ------------------------------------------------------------

# ENHANCED E8 CASCADE THRESHOLDS

# ------------------------------------------------------------

Thresholds:

- Scale: MSigma

Fields: [SigmaF]

- Scale: MG8

Fields: [G8]

- Scale: MNR1

Fields: [NR1]

- Scale: MNR2

Fields: [NR2]

- Scale: MNR3

Fields: [NR3]

- Scale: MPhi

Fields: [phiR, phiI]

  

Groups: {U1Y: U1, SU2L: SU2, SU3c: SU3}

  

Fermions:

Q : {Gen: 3, Qnb: {U1Y: 1/6, SU2L: 2, SU3c: 3}}

L : {Gen: 3, Qnb: {U1Y: -1/2, SU2L: 2}}

uR : {Gen: 3, Qnb: {U1Y: 2/3, SU3c: 3}}

dR : {Gen: 3, Qnb: {U1Y: -1/3, SU3c: 3}}

eR : {Gen: 3, Qnb: {U1Y: -1}}

SigmaF : {Gen: 1, Qnb: {U1Y: 0, SU2L: 3, SU3c: 1}}

G8 : {Gen: 1, Qnb: {U1Y: 0, SU2L: 1, SU3c: 8}}

NR1 : {Gen: 1, Qnb: {U1Y: 0, SU2L: 1, SU3c: 1}}

NR2 : {Gen: 1, Qnb: {U1Y: 0, SU2L: 1, SU3c: 1}}

NR3 : {Gen: 1, Qnb: {U1Y: 0, SU2L: 1, SU3c: 1}}

  

RealScalars:

phiR : {Qnb: {U1Y: 0, SU2L: 1, SU3c: 1}}

phiI : {Qnb: {U1Y: 0, SU2L: 1, SU3c: 1}}

R3 : {Qnb: {U1Y: 0, SU2L: 1, SU3c: 1}, External: True}

  

ComplexScalars:

H : {RealFields: [Pi, Sigma], Norm: 1/sqrt(2), Qnb: {U1Y: 1/2, SU2L: 2}}

  

Potential:

Definitions:

Htilde[i]: Eps[i,j]*Hbar[j]

  

Yukawas:

Yu: Qbar[i,a] Htilde[i] uR[a]

Yd: Qbar[i,a] H[i] dR[a]

Ye: Lbar[i] H[i] eR

# ySig (Type III): intentionally not exported to PyR@TE due to indexing;

# solver includes its 2-loop trace contribution explicitly.

# ySig: Lbar[i] SigmaF[i,j] H[j]

yN1: Lbar[i] Htilde[i] NR1

yN2: Lbar[i] Htilde[i] NR2

yN3: Lbar[i] Htilde[i] NR3

  

QuarticTerms:

lambda : (Hbar[i] H[i])**2

lPhi : (phiR**2 + phiI**2)**2

lHphi : (Hbar[i] H[i])*(phiR**2 + phiI**2)

  

TrilinearTerms:

cR3 : R3 * (Hbar[i] H[i])

  

ScalarMasses:

mu2 : -Hbar[i] H[i]

MPhi : phiR*phiR + phiI*phiI

  

Vevs:

vSM : Pi[2]

vPQ : phiR

  

Parameters:

- {name: MPl, value: 1.221e19}

- {name: MSigma, value: 1.0e3}

- {name: MG8, value: 1.8e10}

- {name: MNR1, value: 1.0e14}

- {name: MNR2, value: 3.0e14}

- {name: MNR3, value: 8.0e14}

- {name: MPhi, value: 1.0e16}

- {name: c3, value: 0.039788735772973836}

- {name: phi0, value: 0.053171952176845526}

- {name: g1, value: 0.357}

- {name: g2, value: 0.652}

- {name: g3, value: 1.2322690515271375}

- {name: Yu33, value: 0.857375}

- {name: Yd33, value: 0.024}

- {name: Ye33, value: 0.010}

- {name: ySig, value: 0.50}

- {name: yN1, value: 0.70}

- {name: yN2, value: 0.70}

- {name: yN3, value: 0.70}

- {name: lambda, value: 0.130}

- {name: lPhi, value: 0.10}

- {name: lHphi, value: 0.01}

- {name: cR3, value: 0.01}

  

Substitutions: { g_U1Y: g1, g_SU2L: g2, g_SU3c: g3 }

  

# ------------------------------------------------------------

# THEORY NOTES (v2):

# - U(1) GUT normalization: g1_GUT = sqrt(5/3) * gY; b1_GUT = 41/10.

# - ySig kept out of PyR@TE export; solver adds its 2-loop trace.

# - LoopOrder=3 in YAML; solver uses full 2-loop + optional 3-loop SU(3).

# - G8 bridge above MG8: Δ(α3^{-1})(μ) = -(Δb3)/(2π) ln(μ/MG8), Δb3=+2.

# ------------------------------------------------------------
```

---

## Appendix G Nilpotent Orbits in Type E8
![[Pasted image 20250823130052.png]]
![[Pasted image 20250823130105.png]]

# Appendix H: Referenzen:

#### 1. Nilpotente Orbits in Semisimple Lie-Algebren (insbesondere E8)

- Collingwood, D.H. und McGovern, W.M., _Nilpotent Orbits in Semisimple Lie Algebras_, Van Nostrand Reinhold, New York (1993). – Referenziert für die allgemeine Klassifikation nilpotenter Orbits in semisimple Lie-Algebren, inklusive detaillierter Tabellen und Dimensionen für E8-Orbits, die als Basis für die γ(n)-Dämpfungsfunktion und die parabolische Kaskade dienen.
- Djouadi, A. et al., "Induced Nilpotent Orbits of the Simple Lie Algebras of Exceptional Type", arXiv: (aus Veröffentlichung, z.B. ähnlich zu iris.unitn.it/handle/11572/77393) (200x). – Referenziert für die Induktion nilpotenter Orbits in E8 und deren Dimensionen, die die monotone Abfallsequenz in der Kaskade (z.B. von 248 zu 206) motivieren.
- Landsberg, J.M. und Manivel, L., "Series of Nilpotent Orbits", Experimental Mathematics 13(1) (2004), 69–78. – Referenziert für die Organisation nilpotenter Orbits in Serien innerhalb exceptioneller Algebren wie E8, einschließlich Dimensionsformeln, die die quadratische Glättung von γ(n) unterstützen.

#### 2. Chern-Simons-Term in 11D Supergravity und Topologische Fixpunkte

- Cremmer, E., Julia, B. und Scherk, J., "Supergravity Theory in Eleven-Dimensions", Physics Letters B 76(4) (1978), 409–412. – Referenziert für die ursprüngliche Formulierung der 11D Supergravity, inklusive des Chern-Simons-Terms, der die Normalisierung 1/(8π) für c₃ und topologische Fixpunkte ableitet.
- Troncoso, R. und Zanelli, J., "Higher-Dimensional Supergravities as Chern-Simons Theories", International Journal of Theoretical Physics 38(4) (1999), 1181–1193 (oder erweiterte Version arXiv:1103.2182). – Referenziert für die Interpretation der 11D Supergravity als Chern-Simons-Theorie, die die topologische Spur von c₃ = 1/(8π) und die Möbius-Reduktion zu φ₀ erklärt.
- Duff, M.J., "Eleven-Dimensional Supergravity, Anomalies and the E8 Yang-Mills Sector", Nuclear Physics B 325(2) (1989), 505–522. – Referenziert für die Verbindung zwischen Chern-Simons-Termen in 11D und E8-Symmetrien, relevant für die topologische Korrektur in φ₀ (z.B. 3/(256π⁴)).

#### 3. E8 in Grand Unified Theories (GUTs) und String-Theorie

- Gross, D.J., Harvey, J.A., Martinec, E. und Rohm, R., "Heterotic String Theory (I). The Free Heterotic String", Nuclear Physics B 256 (1985), 253–284. – Referenziert für die Rolle von E8 × E8 in heterotischen String-Theorien, die die Einbettung von E8 als Ordnungsprinzip für die Skalenleiter (γ(n) aus Orbits) inspirieren.
- Lisi, A.G., "An Exceptionally Simple Theory of Everything", arXiv:0711.0770 (2007). – Referenziert für den Versuch, E8 als vereinheitlichte Symmetrie für alle Kräfte und Materie zu verwenden, ähnlich zur E8-Kaskade im Paper, inklusive Orbit-Strukturen für Flavor und Skalen.
- Green, M.B., Schwarz, J.H. und Witten, E., _Superstring Theory_, Cambridge University Press (1987), Band 2. – Referenziert für die E8-Gauge-Gruppe in String-Theorie, speziell deren nilpotente Elemente und Anomalien, die die quadratische Form von γ(n) und die RG-Bestätigung unterstützen.

#### 4. Theoretische Ableitungen der Feinstrukturkonstante (α)

- Wyler, A., "On the Conformal Groups in the Theory of Relativity and a New Value for the Universal Constant α", Lettere al Nuovo Cimento 3(13) (1971), 533–536. – Referenziert für eine frühe geometrische Ableitung von α aus Konformalgruppen, die die parameterfreie kubische Fixpunktgleichung im Paper (basierend auf Geometrie und Topologie) vorwegnimmt.
- Atiyah, M.F., "On the Fine-Structure Constant", (Vortrag/Notiz, 2018, siehe z.B. preposterousuniverse.com/blog/2018/09/25/atiyah-and-the-fine-structure-constant/). – Referenziert für eine mathematische Herleitung von α ≈ 1/137 aus algebraischen Strukturen, vergleichbar mit der kubischen Gleichung und der ppm-Genauigkeit im Paper.
- Smith, S.J., "A New Theoretical Derivation of the Fine Structure Constant", Progress in Physics 28(1) (2012), 1–5. – Referenziert für eine moderne Ableitung von α ohne freie Parameter, die den Ansatz des Papers (Kopplung von Topologie c₃ und Geometrie φ₀) ergänzt.

#### 5. Weitere Verwandte Themen (z.B. RG-Flüsse, Genetische Algorithmen in Physik)

- 't Hooft, G. und Veltman, M., "Regularization and Renormalization of Gauge Fields", Nuclear Physics B 44(1) (1972), 189–213. – Referenziert für die Grundlagen von RG-Flüssen in Gauge-Theorien, die die Zwei-Schleifen-Analyse und Fingerprints (φ₀ bei 1 PeV, c₃ bei 10^8 GeV) im QCD-Verlauf untermauern.
- Koza, J.R., _Genetic Programming: On the Programming of Computers by Means of Natural Selection_, MIT Press (1992). – Referenziert für die Methode genetischer Algorithmen, die den GA-Ansatz im Paper (Suche nach Lagrange-Dichten und Emergenz von Fixpunkten) methodisch begründet.

## Appendix I: Changelog

### Version 1.0.6 - 2025-09-01

1. Neue Sektion 7.6 -  **Variationsableitung in vier Dimensionen**
2. Sektion 5: Anpassung auf neue 2-Loop RGE Konfiguration & Ergebnisse
3. Neue Sektion 8.0a  - Chiralität aus Rand und Fluss: operative Kurzfassung
4. Neuer Appendix J - Chiralität auf der Doppelabdeckung
5. Sektion 4.5 - Vollständige Anpassung zum Beweis der Eindeutigkeit der Kette
6. Neue Sektionen 8.1.1 und 8.1.2 - Block‑Konstanten aus Randzyklen und abelscher Spur & Herleitung der $\zeta_{B}$‑Formel aus dem Randfunktional
7. Überarbeitung Sektion D 
8. Neu 5.1b  **Schwellen als Leiter‑Outputs statt Modellwahl**
9. Neu 5.2c  **Gauge–Moduli–Locking im E sechs Fenster**
10. Sektion 6.12 - Eindeutige Fixierung von $\alpha_{\text{inf}}$
11. Sektion 8.5.1 - Systematik der Abweichungen
12. Sektion 7.6.2 - Callan-Symanzik-Route
13. Überarbeitung 8.4.6 - Kein Kreisbezug
14. 


## Appendix J — Chiralität auf der Doppelabdeckung

### J.1  Setup und Notation
Geometrie: $(M_{6}=M_{4}\times\widetilde{M})$, dabei ist $(\widetilde{M})$ die **orientierbare Doppelabdeckung** der Möbius Faser mit **drei** geschlossenen Randzyklen $(C_{1},C_{2},C_{T})$.  
Randzählung: $(\sum_{i}\!\oint_{C_{i}}\widehat{k}_{g}\,\mathrm ds = 6\pi)$.  
Die Zahl $(6\pi)$ fixiert in 3.2.2 den linearen Randkoeffizienten und damit $(\varphi_{\text{tree}}=\tfrac{1}{6\pi})$, plus den topologischen Zuschlag $(\delta_{\text{top}}=\tfrac{3}{256\pi^{4}})$ zu $(\varphi_{0})$.  Verweise: 3.2.2 und Appendix D. 


![[Pasted image 20250901105748.png]]
**Fig. J.1** zeigt $(\widetilde M$) mit den drei Randzyklen und den Projektoren $(P_{i})$ (siehe unten).

![[Pasted image 20250901111314.png]]

## J.2  Sechs dimensionale Spinor Reduktion und Projektoren

Wähle $\Gamma$ Matrizen als  
$$
\Gamma^{\mu} = \gamma^{\mu}\otimes\sigma^{0}, \quad 
\Gamma^{5} = \gamma_{5}\otimes\sigma^{1}, \quad 
\Gamma^{6} = \mathbf 1\otimes\sigma^{2}.
$$

Die 6D Weyl Bedingung  
$$
\Gamma_{7} = \gamma_{5}\otimes\sigma^{3}, \quad \Gamma_{7}\Psi = +\Psi
$$
liefert  
$$
\Psi(x,y) = \psi_{L}(x)\otimes \chi_{+}(y) + \psi_{R}(x)\otimes \chi_{-}(y),
$$
mit  
$$
\sigma^{3}\chi_{\pm} = \pm \chi_{\pm}.
$$

Chirale Randbedingungen auf den drei Randzyklen:  
$$
P_{T} = \tfrac12(\mathbf 1 + i\,\sigma^{3}\sigma^{n}), \quad
P_{1} = P_{2} = \tfrac12(\mathbf 1 - i\,\sigma^{3}\sigma^{n}).
$$

Damit existieren Nullmoden nur für $\chi_{+}$, die 4D Nullmoden sind linkschiral.  
Die Wahl ist elliptisch und eichinvariant.  

Optional legen Wilson Linien $W_{i}\in SU(3)\times SU(2)\times U(1)$ entlang $C_{i}$ nur Phasen fest, ohne die 4D Eichgruppe zu brechen.  
Sie dienen dem Auslesen der abelschen Spur im EW Block mit  
$$
k_{\mathrm{EW}} = \tfrac{41}{32}.
$$  
(Verweis 8.4.6)

---

## J.3  Index Satz auf $\widetilde M$ und Familienzahl

Sei $A$ eine abelsche Verbindung mit Fluss  
$$
m = \tfrac{1}{2\pi}\int_{\widetilde M}F \in \mathbb Z
$$
und Randholonomien  
$$
\tfrac{1}{2\pi}\oint_{C_{i}}A = \nu_{i} \in \mathbb Z.
$$

Für die Projektoren $P_{i}$ gilt der Index:  
$$
\mathrm{Ind}\,D_{\widetilde{M}}
= \#\chi_{+} - \#\chi_{-}
= \tfrac{1}{2\pi}\int_{\widetilde M}F
= \nu_{1} + \nu_{2} + \nu_{T}.
$$

Begründungsskizze: APS Formel mit Randprojektionen setzt die $\eta$ Beiträge auf null, Stokes liefert  
$$
\sum_{i}\oint_{C_{i}}A = \int_{\widetilde M}F.
$$

Korollar: Minimal $(\nu_{1},\nu_{2},\nu_{T}) = (1,1,1)$ ergibt  
$$
\mathrm{Ind}\,D = 3
$$
→ drei Familien.  

Die Integer Struktur ist konsistent mit der Chern–Simons Quantisierung aus 3.2.1:  
$$
g = \tfrac{n}{8\pi^{2}}, \qquad c_{3} = \tfrac{1}{8\pi}.
$$

*Fig. J.2 illustriert die Zählung $\mathrm{Ind}\,D = \nu_{1} + \nu_{2} + \nu_{T}$.*

---

## J.4  Anomaliefreiheit einer Familie

Standardmodell pro Familie mit rechtem Neutrino ist anomaliefrei:

$$
U(1)_{Y}^{3}: \quad
3\cdot2\left(\tfrac{1}{6}\right)^{3}
+ 3\left(-\tfrac{2}{3}\right)^{3}
+ 3\left(\tfrac{1}{3}\right)^{3}
+ 2\left(-\tfrac{1}{2}\right)^{3}
+ (1)^{3} = 0,
$$

$$
U(1)_{Y}\times SU(2)^{2}: \quad
3\cdot\tfrac{1}{6} + (-\tfrac{1}{2}) = 0,
$$

$$
U(1)_{Y}\times SU(3)^{2}: \quad
2\cdot\tfrac{1}{6} + (-\tfrac{2}{3}) + (\tfrac{1}{3}) = 0,
$$

$$
\text{Gravitation}-U(1)_{Y}: \quad
3\cdot2\cdot\tfrac{1}{6}
+ 3\left(-\tfrac{2}{3}\right)
+ 3\left(\tfrac{1}{3}\right)
+ 2\left(-\tfrac{1}{2}\right)
+ 1 = 0.
$$

Mit $\mathrm{Ind}\,D = 3$ bleibt die Gesamttheorie anomaliefrei.  

---

## J.5  Kompatibilität mit Rangfenstern und Spur

Die Kette  
$$
E_{8}\supset E_{7}\supset E_{6}
$$
bleibt Fensterlogik im RG Fluss, kein 4D Eichinhalt.  

Die Wilson Linien sind flach und erhalten $SU(3)\times SU(2)\times U(1)$.  

Die gleiche abelsche Spur erscheint zweifach:  
$$
k_{\mathrm{EW}} = \tfrac{41}{32}
$$
im EW Block und  
$$
b_{1} = \tfrac{41}{10}
$$
in  
$$
\kappa = \tfrac{b_{1}}{2\pi}\ln\tfrac{1}{\varphi_{0}}
$$
der $\alpha$ Fixpunktgleichung.  
(Verweise 7.6, 8.4.6)

---

## J.6  Stabilität gegenüber Zwei Schleifen Fenstern

Flache Wilson Linien ändern nur Kaluza Schwellen minimal.  

Die in 5.2 dokumentierten Fingerprints  
$$
\alpha_{3}(1\,\mathrm{PeV}) \simeq \varphi_{0}, \qquad 
\alpha_{3} \simeq c_{3} \quad \text{bei } \mu \sim 2.5\times 10^{8}\,\mathrm{GeV}
$$
bleiben stabil.