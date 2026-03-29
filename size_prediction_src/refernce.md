A Taguchi-Fuzzy Multi-Response Optimization Approach to the Synthesis of Ni–Ti (Nitinol) Nanoparticles via Plasma-Assisted Micro-EDM
Plasma-Assisted Micro-EDM for Controlled Nitinol Nanoparticle Synthesis: A Hybrid Taguchi–Fuzzy Multi-Objective Optimization Approach

Abstract
This study examines the controlled fabrication of Nitinol nanoparticles using a specially developed plasma-assisted micro-electrical discharge machining (micro-EDM) system. A Taguchi L9 experimental design was employed to evaluate the influence of key parameters voltage, peak current, duty cycle, and discharge frequency on nanoparticle formation. Particle morphology and size were analyzed using SEM, while zeta potential measurements were carried out to assess dispersion stability.The results showed noticeable variations in particle size and stability under different operating conditions, highlighting the difficulty of achieving both smaller size and improved suspension behavior simultaneously. To address this, a hybrid optimization strategy combining the Taguchi method with fuzzy logic was applied. The fuzzy inference approach integrated multiple response variables into a single Fuzzy Performance Index (FPI) to identify optimal conditions.The optimum parameter combination was found to be 25 V, 3 A, 30% duty cycle, and 5 kHz. Validation experiments conducted at these settings produced nanoparticles with an average size of 227.8 nm and improved uniformity. Energy dispersive X-ray spectroscopy (EDS) confirmed the formation of Ti–Ni nanoparticles with only minor surface oxidation. Overall, the findings demonstrate that micro-EDM is a viable non-conventional technique for synthesizing Nitinol nanoparticles and provide a structured hybrid optimization framework for advanced nanomaterial production.
Keywords:
Micro-EDM, Nitinol nanoparticles, Multi-objective optimization (fuzzy logic), Nanofluids.

Introduction
Smart nanomaterials are rapidly moving from research to industrial use, with the market growing from USD 12.87 billion in 2025 to a projected USD 29.97 billion by 2032, demon-starting a CAGR of 12.83% driven by rising demand for adaptive, high-performance materials [1] Smart materials have emerged as a game-changing class of engineering materials due to their inherent capacity to respond reliably to environmental stimuli such as temperature, stress, electric fields, and magnetic fields. Nickel-titanium (NiTi) shape memory alloy, also known as Nitinol, stands out because to its distinct shape memory effect (SME), superplastic behavior, high damping capacity, outstanding corrosion resistance, and exceptional biocompatibility [2], [3], [4]. Nitinol’s multifunctional properties have enabled it to be widely used in a variety of industries, including biomedical implants, micro-actuators, aerospace components, robotics, and adaptive structural systems [5], [6], [7].
In recent years, materials research has increasingly transitioned from bulk smart alloys to nano-structured and nanoparticle-based smart materials, motivated by the understanding that nanoscale materials frequently display markedly different physical, thermal, and functional properties relative to their bulk equivalents [8], [9] Nitinol nanoparticles exhibit size-dependent phase transformation behavior, altered martensitic transformation temperatures, increased surface reactivity, and adjustable functional responses, rendering them suitable for applications in smart coatings, nano-actuation, targeted drug delivery, and responsive sensors [10], [11], [12] However, the regulated production of compositionally stable and functionally reliable Nitinol nanoparticles remains a substantial difficulty.
NiTi-based nanostructures have been investigated using conventional synthesis techniques as sol-gel processing, chemical reduction, mechanical alloying, laser ablation, and physical vapor deposition [13], [14], [15], [16] Although these technologies are capable of producing nanoscale particles, they typically have drawbacks such as oxidation, chemical contamination, compositional divergence because nickel or titanium preferentially evaporates, high processing costs, and limited scalability [17], [18], [19]. Furthermore, chemical synthesis methods typically involve risky precursors and lengthy post-processing steps, which limits their applicability for ecologically friendly and biomedical-grade applications [20] These limitations highlight the need for alternative, precision-based, non-chemical methods of creating smart nitinol nanoparticles.
Electrical discharge machining (EDM) is a non-traditional manufacturing technology widely recognized for its capability to machine electrically conductive, hard-to-machine materials irrespective of their mechanical strength or hardness [21], [22], [23] In EDM, material removal happens through a series of regulated electrical discharges between a tool electrode and a conductive workpiece submerged in a dielectric liquid. The intense concentrated thermal energy created during each discharge results in fast melting and partial vaporization of the workpiece material, followed by quenching due to the surrounding dielectric fluid [24] This unusual thermophysical environment naturally favors the development of small trash particles, many of which fall within the micro- and nano-scale range.
Improvements in micro-electrical discharge machining (micro-EDM) have increased the precision and control of the discharge process, allowing for feature sizes in the micrometer range and improving plasma channel stability [25], [26], [27]. Micro-EDM is particularly ideal for processing advanced functional materials such as shape memory alloys, metallic glasses, and intermetallic compounds, which are otherwise difficult to machine using conventional techniques [28] Micro-EDM has significant advantages for Nitinol because it reduces cutting forces, mechanical deformation, and allows for precise control over heat input [29].
Recent research indicates that the debris produced during EDM should not be considered as waste but rather a potential reservoir of useful micro- and nano-particles. In every discharge pulse, molten bits from the workpiece get blasted out and freeze almost instantly in the dielectric fluid, cooling at rates over 106 k/s. This forms round or nearly round particles with highly refined internal structures. The extreme heat shifts and rapid cooling in micro-EDM create ideal conditions for nanoparticles featuring unstable phases, lots of defects, and altered surface properties[30], [31], [32]
This rapid solidification is critical for smart alloys such as Nitinol, as it may preserve essential martensitic or austenitic phases even at very small scales [24]. Nitinol nanoparticle production in micro-EDM is significantly influenced by discharge current, pulse-on time, pulse-off time, open-circuit voltage, dielectric fluid, electrode type, and polarity [33], [34], [35]. They have a direct impact on particle size spread, morphology, composition, and phase retention by controlling discharge energy, plasma behavior, and cooling rate.
Prior research has primarily focused on traditional metals such as copper, aluminum, silver, and tungsten, however the synthesis of EDM nanoparticles has gained popularity [36], [37], [38], [39], [40], and [41]. Nanoparticles of smart alloys, particularly nitinol, are still not well understood. Furthermore, rather than focusing on the characteristics and performance of the nanoparticles themselves, research tends to focus on machining outcomes including removal rates, tool wear, and surface finish [4, 43]. A significant gap in the field is represented by this oversight.
Micro-EDM on Nitinol is difficult due to its heat sensitivity and phase shifts. Excess discharge energy increases the likelihood of selective nickel evaporation, titanium oxidation, or brittle intermetallic, all of which reduce the shape memory and super elastic properties of both the workpiece surface and the resultant nanoparticles. Tight control and optimization of process parameters are thus required to ensure optimal composition and functionality in these EDM-produced Nitinol nanoparticles [24], [44], and [45].
Full characterisation is necessary to understand how structure, characteristics, and processes influence nanoparticle generation in micro-EDM. SEM and TEM provide particle shapes and sizes, whereas XRD reveals phase composition and crystal structures [46], [47], [48]. EDS manages elemental breakdown, while DSC studies nanoscale martensitic changes and heat reactions [49]. Together, these methods validate if Nitinol retains its smart properties after micro-EDM.

Nitinol nanoparticles produced using micro-EDM have significant potential for the development of multifunctional smart systems. Because of their high surface-to-volume ratio and inherent shape memory behavior, they are ideal for reactive coatings, adaptive composite materials, and biomedical devices requiring controlled drug release or localized actuation [50], [51]. Furthermore, by limiting chemical waste, allowing material recovery from machining debris, and merging fabrication and nanomaterial synthesis, micro-EDM-based nanoparticle formation fits with sustainable manufacturing principles [52]. 

The use of micro-EDM to create nitinol nanoparticles holds great promise for multipurpose smart applications. Adaptive composite structures, reactive surface coatings, and biomedical systems needing controlled drug delivery or localized actuation can all benefit from their large surface-to-volume ratio and intrinsic shape memory [24], [53].
The controlled synthesis of nitinol nanoparticles is examined in this work in relation to the parameters of the micro-EDM process. After conducting a thorough morphological characterization and dispersion stability analysis, a combined Taguchi–fuzzy logic approach was used to produce particles with a smaller size and better stability. The study illustrates the possibility of creating functional Nitinol nanoparticles through unconventional manufacturing methods by directly linking discharge parameters and nanoparticle characteristics.
Experimentation
Nitinol nanoparticles were prepared using a custom-built, laboratory-scale plasma-assisted micro/nano-EDM system. The setup provided accurate control over electrode motion, discharge consistency, and processing parameters. A specially fabricated X–Y–Z positioning stage, operated by micro-stepper motors and managed through an Arduino controller, ensured precise and repeatable alignment of the tool electrode.
A closed-loop control mechanism was implemented to maintain a stable and appropriate spark gap between the 0.8 mm diameter copper electrode and the machining region. The electrical circuit permitted independent regulation of gap voltage, peak current, duty cycle, and discharge frequency. The system was also coupled with an ultrasonic unit, and deionized water served as the dielectric medium. This dielectric fluid facilitated spark formation, enabled rapid quenching, and assisted in dispersing the synthesized Nitinol nanoparticles, resulting in the formation of stable nanofluids.

Figure 1: Schematic diagram of the experimental setup
Design of Experiments and Synthesis Procedure
The Design of Experiments (DOE) approach was applied to methodically assess the way process parameters affected the synthesis of nitinol nanoparticles. Based on initial testing and previous investigation, they determined four main parameters: duty cycle (%), frequency (kHz), peak current (A), and gap voltage (V). Each experiment was performed at three levels to find non-linear patterns, as illustrated in Table 1. 
Table 1: Factors and levels considered for the experimental study
Factors
Leve 1
Level 2
Level 3
Gap Voltage (V)
20
25
30
Current (A)
1
2
3
Duty Cycle (%)
30
40
50
Frequency (kHz)
3
4
5


A Taguchi L9 orthogonal array (Table 2) was employed to organize the experimental matrix, allowing for a comprehensive analysis of the parameter space through nine distinct runs. This method ensures a balanced comparison of factor levels while minimizing experimental cost, effort and redundancy. The synthesis was performed by submerging Nitinol workpiece and Nitinol electrode in deionized water, coupled with ultrasonicator operating at the frequency of 40kHz (PS-06A). During synthesizing the pulse discharge waveform were monitored and captured by ROGOL oscilloscope ( DHO924). Following each run, the resulting colloidal samples were collected in borosilicate vials as shown in Figure 1 and monitored for initial physical stability and sedimentation behavior over a period of twenty-four hours.

Table 2: Taguchi L9 experimental layout for nanoparticle synthesis
Sl. No.
Gap Voltage (V)
Current (A)
Duty Cycle (%)
Frequency kHz
1
20
1
30
3
2
20
2
40
4
3
20
3
50
5
4
25
1
40
5
5
25
2
50
3
6
25
3
30
4
7
30
1
50
4
8
30
2
30
5
9
30
3
40
3

Characterization Techniques:
The synthesized Nitinol nanofluids produced from each experimental runs were characterized to evaluate particle size and dispersion stability. The particle size distribution and morphology were analyzed using scanning electron microcopy (SEM). SEM analysis study was performed using a GEMINI 300 scanning electron microscope (Carls Zeiss, Germany). Prior to imaging the nanofluids samples were ultrasonicate for the period of 2 hours using PreCiSonic (Model – USB 1.5 L/H/DTC) for the well dispersion of particles and were prepared by placing small droplet of nanofluid onto a carbon-coated copper grid and allowed to dry.
The dispersion stability of the synthesized nanofluids was evaluated through zeta potential test. Zeta potential test was performed using a Litesizer 500 (Anton Paar, Austria). The measurements were conducted at room temperature and the obtained results were used as an indicator of the electrostatic stability of the synthesized nanoparticles suspensions.
Results and discussion
SEM-Based Particle Size Analysis
 The shape and particle size properties of the produced nitinol nanoparticles derived from the nine experimental combinations were examined using scanning electron microscopy (SEM). To guarantee repeatability, the synthesis experiment was carried out twice for each run, and samples taken from both trails were subjected to SEM imaging.


From the SEM images particle size dimensions were extracted using image analysis software. All the particle sizes were compiled and subsequent processed with statistical analysis. The synthesized particle size distributions were evaluated and frequency distribution analysis was conducted to identify the number of particles falling within specific size intervals.
A particle size distribution histogram with log-normal fitting for selected experimental runs are shown in Figure 4, which is frequently used to characterize nanoparticles size distributions produced under rapid nucleation and quenching circumstances, was employed to match the experimental particle size data.  Important statistical characteristics such as minimum particle size, maximum particle size, mean particle size and typica; percentile diameters like D10, D20, D50, D80, D90 could be extracted from the fitted distributions. 


Table 3: Summary of Statistical Metrics
Experiment
Average Size (nm)
Minimum (nm)
Median
(nm)
Maximum
(nm)
D10
D20
D50
D80
D90
RUN 1
381.49
115
351
1010
164
206
351
802
861
RUN 2
388.5178
91
430
1071
165
211
430
982
1050
RUN 3
317.9062
65
365
752
130
193
365
548
604
RUN 4
240.4507
54
270.875
722
83
129
270.875
478
544.5
RUN 5
295.4024
30
279
868
159
193
279
478
709
RUN 6
225.217
110
237
484
114
130
237
415
436
RUN 7
361.9433
192
330
893
224
269
330
505
655
RUN 8
258.8177
33
228
1130
66
142
228
560
783
RUN 9
426.6988
202
408
1159
210.5
238
408
830
973


Particle size distribution varies significantly between experimental runs, according to the statistical parameters provided in Table 3. Run 6 exhibited comparatively smaller particle size, characterized by a lower average particle size and a lower D50 value in comparison to previous runs. Run 6 also had much lower D80 and D90 values, pointing to a narrower size spread. In contrast, other trials produced bigger particles and wider distributions, likely from clumping or bigger molten blobs during discharges. This shows micro-EDM settings play a big role in nanoparticle formation and size patterns.
Zeta potential Analysis and Dispersion Stability
The zeta potential of the Nitinol nanofluids was measured in all nine experimental runs to determine their stability. For each run, two independent samples were prepared, and the average values and standard deviations are displayed in Figure 5.
Under various operating conditions, the obtained zeta potential values ranged from near-neutral to rather high magnitudes. Run 6 demonstrated the highest average zeta potential (−16.35 mV) of the nine experiments, demonstrating the nanofluid's electrostatic stability. Run 5 had zeta potential values close to zero, indicating a higher tendency for particle agglomeration and decreased dispersion stability. The other tests showed intermediate stability levels, demonstrating that surface charge generation is extremely sensitive to the treatment settings applied.
The differences in zeta potential can be explained by changes in discharge energy, plasma temperature, and cooling rate, which change the size of the particles, the structure of the surface, and the amount of oxidation on the surface. Changes in particle size distribution and surface charge density from one run to the next cause changes in measured values because zeta potential depends heavily on surface chemistry and effective surface area [54]. Thus, the observed variations result from discrepancies in synthesis conditions rather than from experimental error.














Representative zeta potential distribution plots for selected experimental runs (run6, run5 and run9) are presented in Figure 6 to illustrate stable, unstable and intermediate dispersion behavior, respectively. A summary of average zeta potential values for all nine runs, along with error is presented in Figure 5. 
The drastic variation observed in the seta potential values across the experimental runs, in combination with various particle size characteristics obtained from SEM analysis, highlights the multi-response nature of the nanoparticles synthesis process. This necessitates the adoption of a multi-criteria optimization approach, which is addressed in the subsequent section using a fuzzy logic based framework.
Multi-Objective Optimization Using Fuzzy Logic
Rationale
The traditional Taguchi technique is best suited for single-response optimization. However, the current study includes two opposing performance characteristics: particle size determined by SEM analysis and dispersion stability determined by zeta potential research. Direct Taguchi optimization is ineffective due to smaller particle size and higher magnitude of zeta potential. To overcome this constraint, a fuzzy logic-based multi-objective optimization methodology was combined with the Taguchi method. Fuzzy logic aggregates many performance responses into a single scalar index known as the Fuzzy Performance Index (FPI), which can then be utilized for Taguchi optimization. 
Fuzzy Inference System Configuration
The Matlab Fuzzy Logic Toolbox was utilized to carry out fuzzy logic operations. A mamdani-type fuzzy inference system (FIS) was created, which included fuzzification, rule assessment, and defuzzification procedures. Triangular membership functions were used for all variables due of their simplicity and computational efficiency. Based on the experimental results, we developed two input variables and one output variable: 

Input 1: SEM Average Particle Size (mm)
Smaller particle sizes are preferred, this input was classified into three linguistics categories: small, medium, and large.
Input 2: Zeta potential (absolute value, mV)
Higher absolute zeta potential indicates better dispersion stability, this input was categorized as low, medium and high.
Output: Fuzzy Performance Index (FPI)
The FPI reflects total synthesis performance and was standardized between 0 and 1.
 The SEM average particle size was defined between 220 nm and 440 nm, and the zeta potential range was selected between 0 mV and 17 mV based on experimental extremes. The output FPI was divided into nine linguistics levels to ensure proper resolution in performance rating. Figure 7 depicts the membership functions used for the input and output variables.
For the combination of three particle size levels and three zeta potential levels (3 * 3), nine inference rules were developed. Smaller particle size and greater dispersion stability were given priority in the design of the rule base, which was done in compliance with the experimental goals. The rule base's rationale can be explained as follows:

Best case: Small particle size and high zeta potential yield a high FPI.
Worst case: Large particle size and low zeta potential result in a low FPI.























The complete rule set used in the fuzzy logic inference system is presented in Table 4. 
             Table 4: Fuzzy Inference Rules for Nitinol Nanoparticle Optimization
Rule No.
IF SEM Size is...
AND Zeta Potential is...
THEN FPI is...
Logic Justification
1
Small
High
Huge
Best Case (High Quality)
2
Small
Medium
Very Large
Excellent stability, good size
3
Small
Low
Large
Good size, poor stability
4
Medium
High
Medium Large
Average size, excellent stability
5
Medium
Medium
Medium
Neutral / Average Performance
6
Medium
Low
Small Medium
Average size, poor stability
7
Large
High
Small
Poor size, excellent stability
8
Large
Medium
Very Small
Poor size, average stability
9
Large
Low
Tiny
Worst Case (Low Quality)


Defuzzification was performed using the centroid (center of gravity) method, which converts the aggregated fuzzy output into a single crisp FPI value ranging from 0 to 1. This method provides a balanced and physically meaningful representation of the overall performance.
Fuzzy Performance Index Evaluation
Table 5 summarized the estimated Fuzzy Performance Index values from all nine testing runs. These FPI values represent the combination effect of particle size and dispersion stability, and were then used as a single repose parameter in Taguchi optimization.
Table 5: Fuzzy Performance Index
Exp No:
Voltage (V)
Current (A)
Duty Cycle (%)
Frequency kHz
SEM Images (Average Size) (nm)
Zeta Potential mV
FPI
RUN1
20
1
30
3
381.48995
4.15
0.0853
RUN2
20
2
40
4
388.51781
5.05
0.124
RUN3
20
3
50
5
317.90622
1.5
0.374
RUN4
25
1
40
5
240.45073
2.55
0.746
RUN5
25
2
50
3
295.40238
0.3
0.374
RUN6
25
3
30
4
225.217
16.35
0.959
RUN7
30
1
50
4
361.94327
10.4
0.514
RUN8
30
2
30
5
258.81772
9.4
0.874
RUN9
30
3
40
3
426.69883
10.5
0.146


Among the experimental runs, run 6 exhibited the highest FPI, showing greater overall performance due to the combination of lower particle size and much higher dispersion stability. Lower FPI values, on the other hand, are indicative of undesirable particle size and zeta potential combinations.
Taguchi Analysis using Fuzzy performance Index
The Fuzzy Performance Index (FPI) values derived from the fuzzy logic analysis were used as the single response variable for Taguchi optimization. Since higher FPI equates to improved overall nanoparticle quality, combining smaller particle size and enhanced dispersion stability, the “larger-the-better” signal-to-noise (S/N) ratio criterion was employed for the analysis.

Figure 8: Main effects plot for the S/N ratio
 The impact of each process parameter on the multiple objective performance represented by the FPI is depicted in Figure 7 main effects plot for the S/N ratios. The level that corresponds to the maximum S/N ratio for each control factor represents the ideal setting since it shows enhanced synthesis performance and resilience. The following process parameters were found to be ideal for the manufacture of nitinol nanoparticles based on the S/N ratio analysis:
Voltage: 25 V (Level 2)
Current: 3 A (Level 3)
Duty Cycle: 30% (Level 1)
Frequency: 5 kHz (Level 3)
The hybrid fuzzy logic-Taguchi approach effectively identified an optimal parameter combination of 25 V, 3A, 30% duty cycle, and 5 kHz which resulted in the highest Fuzzy Performance Index. This parameter set represents the best compromise between minimizing particle size and maximizing dispersion stability among the investigated experimental conditions.
Validation of Optimized Process parameters
SEM and EDS Analysis
To validate the effectiveness of the optimized micro-EDM process parameters identified through the fuzzy logic-Taguchi approach, SEM analysis was performed on nanoparticles synthesized under the optimal condition of 25 V, 3 A, 30 % duty cycle and 5 kHz frequency. Representative SEM micrographs and the corresponding particle size distribution are shown in Figure 9.

(a)

(b)

(c)

Figure 9: (a) and (b) SEM images and (c) statistical analysis for optimized parameter
The SEM images reveal the formation of predominantly nanoscale particles with the reduced agglomeration compared to the several non-optimized experimental runs. Quantitative particle size was carried out measuring 125 individual particles, and the resulting data were statistically analyzed using a log-normal distribution model. The optimized condition yielded an average particle size of approximately 220.93 nm, with median (D50) values of 161 nm, indicating a clear shift toward finer particle formation.  
Table 6: Statistical Metrics for Optimized Parameter
Experiment
Average Size (nm)
Minimum (nm)
Median
(nm)
Maximum
(nm)
D10
D20
D50
D80
D90
Optimized Parameter
220.93
47
161
956
90
97
161
306
435


The particle size distribution indicates that many nanoparticles are in the lower size range, with D10 and D20 values of 90 and 97 nm, respectively. Furthermore, D80 and D90 values of 306 nm and 435 nm show a reduced size spread compared to non-optimized synthesis circumstances. The distribution of plasma-assisted nanoparticles is highly skewed towards smaller particle sizes, despite the presence of a tiny population of bigger particles (as evidenced by the greatest measured size).
The fitted log-normal distribution parameters demonstrate the efficacy of the modified discharge conditions in promoting successful nanoparticle nucleation, rapid quenching, and regulated material degradation. The optimization method is validated by the reduction in average particle size and enhancement of distribution uniformity, indicating that the micro-EDM process can produce nitinol nanoparticles with superior size attributes under optimal operating conditions. 

Energy Dispersive X-ray Spectroscopy (EDS) analysis was performed to verify the elemental composition of nanoparticles synthesized under the optimized micro-EDM parameters. The corresponding EDS spectrum and elemental composition is presented in Figure 10.

Figure 10: EDS Spectrum and composition table
Titanium (Ti) and nickel (Ni), the two main alloying elements, are visible in the spectrum, indicating that the electrode material was successfully converted into nanoparticle form during the micro-EDM process. The optimized nanoparticles had 51.31 weight percent Ti and 27.01 weight percent Ni, suggesting that the basic makeup of nitinol is preserved throughout plasma-assisted synthesis. 
There was also a 14.70 weight percent O oxygen peak found. Rapid quenching in the deionized water dielectric medium and surface oxidation of nanoparticles during high-temperature discharge are responsible for the oxygen's presence. Such oxide formation is expected given the high surface-to-volume ratios of nanoparticles and is common in plasma-based nanoparticle synthesis approaches.

The spectrum showed a slight copper signal (6.98 weight percent Cu). Cu may be present in the micro-EDM system due to secondary transfer processes, such as material migration from fixtures, electrical connections, or auxiliary conductive components that are exposed to the discharge environment. The comparatively low concentration of copper indicates that it represents trace inclusion throughout the synthesis process rather than dominating the composition.

The feasibility of using micro-EDM as an unconventional method to produce nitinol nanoparticles is validated by the EDS analysis, which shows that the synthesized nanoparticles are primarily composed of Ti and Ni, with minor elemental contributions coming from the high-energy plasma environment.
Zeta Potential
Zeta potential analysis was performed on the nanofluid synthesized under the optimized process parameter (25V, 3A, 30% duty cycle and 5 kHz). Figure 11 depicts the zeta potential distribution. The optimized sample’s zeta potential of -19.6 mV shows an improvement in magnitude in comparison to the values obtained under the non-optimized experimental conditions and this comparison can be attributed to the higher negative surface charge indicated that is improved electrostatic repulsion, indicating improved initial dispersion stability. 
However, when the nanofluid were stored in the glass vail and time-dependent observation showed that the electrostatic stabilization magnitude was insufficient to sustain long-term suspension stability. Partial coagulation and sedimentation were observed after about an hour when sample were kept in still position without ultrasonication, suggesting a slow decrease in effective interparticle repulsive interactions. Van der Waals attraction forces may finally outweigh electrostatic repulsion in metallic nanoparticles made in pure deionized water without the use of surfactants.

Figure 11: Zeta Potential for optimized process parameter.
Conclusion
This study utilized a specially designed plasma-assisted micro-EDM system to fabricate Nitinol nanoparticles. A Taguchi L9 experimental design was used to systematically test how the voltage, peak current, duty cycle, and discharge frequency of the process affect the shape and behaviour of the particles.
SEM observations demonstrated distinct differences in the distribution of particle sizes under various operating conditions, with discharge energy having significant effects on particle size. Important differences in dispersion stability between the produced nanofluids were also revealed by zeta potential measurements, underscoring the multi-response character of the optimization problem.
A combined Taguchi and fuzzy logic optimization strategy was used to simultaneously reduce particle size and enhance dispersion stability. Effective process condition selection was made possible by the fuzzy inference system's integration of several response parameters into a single Fuzzy Performance Index (FPI). The highest FPI value was obtained by identifying the ideal parameter set, which was 25 V, 3 A, 30% duty cycle, and 5 kHz.
Confirmation tests conducted under these ideal conditions revealed better particle properties, with a median particle diameter (D50) of 161 nm and an average particle size of roughly 227.8 nm. According to EDS results, nickel and titanium were present, confirming the formation of Nitinol nanoparticles with low oxide content and trace amounts of elemental elements that are usually connected to plasma-assisted processing.
The effectiveness of the fuzzy logic-Taguchi optimization in improving nanoparticles dispersion characteristics was demonstrated by the synthesize nanofluid under optimized conditions enhanced zeta potential of -19.6 mV, which was greater in magnitude than the non-optimized experimental runs. To achieve long-term colloidal stability, however, further stabilizer insertion may be needed, as seen by the agglomeration and sedimentation seen after about an hour.   
These results indicate that micro-EDM can be a useful alternative technique for the regulated synthesis of nitinol nanoparticles. Fuzzy logic combined with the Taguchi method provides a dependable approach to multi-objective optimization in the synthesis of nanoparticles. In order to increase the potential applications of the developed nanofluids, future research may concentrate on evaluating their functional performance and long-term stability.

References
[1]	“Nanomaterials Market Size, Competitors & Forecast to 2032.” Accessed: Feb. 02, 2026. [Online]. Available: https://www.researchandmarkets.com/report/nanomaterials
[2]	“K. Otsuka, C. M. Wayman, Shape memory materials, Cambridge university press,”.
[3]	D. C. . Lagoudas, Shape memory alloys : modeling and engineering applications. Springer, 2008.
[4]	S. Miyazaki and K. Otsuka, “Development of Shape Memory Alloys,” ISIJ International, vol. 29, no. 5, pp. 353–377, 1989, doi: 10.2355/ISIJINTERNATIONAL.29.353.
[5]	N. B. Morgan, “Medical shape memory alloy applications—the market and its products,” Materials Science and Engineering: A, vol. 378, no. 1–2, pp. 16–23, Jul. 2004, doi: 10.1016/J.MSEA.2003.10.326.
[6]	J. A. Shaw and S. Kyriakides, “Thermomechanical aspects of NiTi,” J. Mech. Phys. Solids, vol. 43, no. 8, pp. 1243–1281, Aug. 1995, doi: 10.1016/0022-5096(95)00024-D.
[7]	T. Duerig, A. Pelton, and D. Stöckel, “An overview of nitinol medical applications,” Materials Science and Engineering: A, vol. 273–275, pp. 149–160, Dec. 1999, doi: 10.1016/S0921-5093(99)00294-4.
[8]	H. Gleiter, “Nanostructured materials: basic concepts and microstructure,” Acta Mater., vol. 48, no. 1, pp. 1–29, Jan. 2000, doi: 10.1016/S1359-6454(99)00285-2.
[9]	“Nanostructured Materials: Processing, Properties, and Applications: Second Edition | Request PDF.” Accessed: Feb. 02, 2026. [Online]. Available: https://www.researchgate.net/publication/298017216_Nanostructured_Materials_Processing_Properties_and_Applications_Second_Edition
[10]	J. Frenzel, E. P. George, A. Dlouhy, C. Somsen, M. F. X. Wagner, and G. Eggeler, “Influence of Ni on martensitic phase transformations in NiTi shape memory alloys,” Acta Mater., vol. 58, no. 9, pp. 3444–3458, May 2010, doi: 10.1016/J.ACTAMAT.2010.02.019.
[11]	T. Waitz, T. Antretter, F. D. Fischer, and H. P. Karnthaler, “Size effects on martensitic phase transformations in nanocrystalline NiTi shape memory alloys,” Materials Science and Technology, vol. 24, no. 8, pp. 934–940, 2008, doi: 10.1179/174328408X302620;WGROUP:STRING:PUBLICATION.
[12]	Y. C. Cao, “Nanomaterials for biomedical applications,” Nanomedicine, vol. 3, no. 4, pp. 467–469, Aug. 2008, doi: 10.2217/17435889.3.4.467;CTYPE:STRING:JOURNAL.
[13]	Y. Uslugil, A. Feyzullah, M. Ekrem, and A. Avcı, “Mechanical alloying of NiTi alloys with Ag and Cu nanoparticles: synthesis and characterization,” Phys. Scr., vol. 100, no. 7, Jul. 2025, doi: 10.1088/1402-4896/ADE018.
[14]	A. Balachandran, S. P. Sreenilayam, K. Madanan, S. Thomas, and D. Brabazon, “Nanoparticle production via laser ablation synthesis in solution method and printed electronic application - A brief review,” Results in Engineering, vol. 16, p. 100646, Dec. 2022, doi: 10.1016/J.RINENG.2022.100646.
[15]	T. Mousavi, F. Karimzadeh, and M. H. Abbasi, “Synthesis and characterization of nanocrystalline NiTi intermetallic by mechanical alloying,” Materials Science and Engineering: A, vol. 487, no. 1–2, pp. 46–51, Jul. 2008, doi: 10.1016/J.MSEA.2007.09.051.
[16]	V. Amendola and M. Meneghetti, “Laser ablation synthesis in solution and size manipulation of noble metal nanoparticles,” Physical Chemistry Chemical Physics, vol. 11, no. 20, pp. 3805–3821, 2009, doi: 10.1039/B900654K.
[17]	K. H. Ho and S. T. Newman, “State of the art electrical discharge machining (EDM),” Int. J. Mach. Tools Manuf., vol. 43, no. 13, pp. 1287–1300, Oct. 2003, doi: 10.1016/S0890-6955(03)00162-7.
[18]	M. Kunieda, B. Lauwers, K. P. Rajurkar, and B. M. Schumacher, “Advancing EDM through Fundamental Insight into the Process,” CIRP Annals, vol. 54, no. 2, pp. 64–87, Jan. 2005, doi: 10.1016/S0007-8506(07)60020-1.
[19]	N. Mohd Abbas, D. G. Solomon, and M. Fuad Bahari, “A review on current research trends in electrical discharge machining (EDM),” Int. J. Mach. Tools Manuf., vol. 47, no. 7–8, pp. 1214–1228, Jun. 2007, doi: 10.1016/J.IJMACHTOOLS.2006.08.026.
[20]	T. Masuzawa, “State of the Art of Micromachining,” CIRP Annals, vol. 49, no. 2, pp. 473–488, Jan. 2000, doi: 10.1016/S0007-8506(07)63451-9.
[21]	H. C. Lin, K. M. Lin, and I. S. Cheng, “Electro-discharge machining characteristics of TiNi shape memory alloys,” J. Mater. Sci., vol. 36, no. 2, pp. 399–404, Jan. 2001, doi: 10.1023/A:1004824511973/METRICS.
[22]	S. F. Hsieh, M. H. Lin, S. L. Chen, S. F. Ou, T. S. Huang, and X. Q. Zhou, “Surface modification and machining of TiNi/TiNb-based alloys by electrical discharge machining,” The International Journal of Advanced Manufacturing Technology 2016 86:5, vol. 86, no. 5, pp. 1475–1485, Jan. 2016, doi: 10.1007/S00170-015-8257-X.
[23]	Y. Lü, X. Gao, and Y. Hou, “Material removal and microstructure evolution of nitinol shape memory alloy by magnetic field assisted-wire electrical discharge machining,” Mater. Today Commun., vol. 48, p. 113681, Sep. 2025, doi: 10.1016/J.MTCOMM.2025.113681.
[24]	J. W. Mwangi, L. Weisheit, V. D. Bui, M. Y. Zanjani, and A. Schubert, “Influence of Micro-EDM on the Phase Transformation Behaviour of Medical-Grade Nitinol,” Shape Memory and Superelasticity 2018 4:4, vol. 4, no. 4, pp. 450–461, Nov. 2018, doi: 10.1007/S40830-018-00195-1.
[25]	M. P. S. K. Kiran and S. S. Joshi, “Modeling of Surface Roughness and the Role of Debris in Micro-EDM,” J. Manuf. Sci. Eng., vol. 129, no. 2, pp. 265–273, Apr. 2007, doi: 10.1115/1.2540683.
[26]	P. Kumar, P. K. Singh, D. Kumar, V. Prakash, M. Hussain, and A. K. Das, “A novel application of micro-EDM process for the generation of nickel nanoparticles with different shapes,” Materials and Manufacturing Processes, vol. 32, no. 5, pp. 564–572, Apr. 2017, doi: 10.1080/10426914.2016.1244832.
[27]	J. Wang, J. A. Sánchez, Z. Wang, B. Izquierdo, and I. Ayesta, “Observations on debris composition and size distribution in WEDM,” Procedia CIRP, vol. 95, pp. 331–336, 2020, doi: 10.1016/J.PROCIR.2020.02.323.
[28]	M. Gostimirovic, M. Radovanovic, M. Madic, D. Rodic, and N. Kulundzic, “Inverse electro-thermal analysis of the material removal mechanism in electrical discharge machining,” International Journal of Advanced Manufacturing Technology, vol. 97, no. 5–8, pp. 1861–1871, Jul. 2018, doi: 10.1007/S00170-018-2074-Y.
[29]	D. E. Nicholson et al., “Standardization of Shape Memory Alloys from Material to Actuator,” Shape Memory and Superelasticity, vol. 9, no. 2, pp. 353–363, Jun. 2023, doi: 10.1007/S40830-023-00431-3.
[30]	J. W. Murray, J. Sun, D. V. Patil, T. A. Wood, and A. T. Clare, “Physical and electrical characteristics of EDM debris,” J. Mater. Process. Technol., vol. 229, no. 9–11, pp. 54–60, Mar. 2016, doi: 10.1016/j.jmatprotec.2015.09.019.
[31]	A. K. Khanra, L. C. Pathak, and M. M. Godkhindi, “Microanalysis of debris formed during electrical discharge machining (EDM),” Journal of Materials Science 2007 42:3, vol. 42, no. 3, pp. 872–877, Feb. 2007, doi: 10.1007/s10853-006-0020-0.
[32]	P. Kumar, P. K. Singh, D. Kumar, V. Prakash, M. Hussain, and A. K. Das, “A novel application of micro-EDM process for the generation of nickel nanoparticles with different shapes,” Materials and Manufacturing Processes, vol. 32, no. 5, pp. 564–572, Apr. 2017, doi: 10.1080/10426914.2016.1244832.
[33]	J. Vora et al., “Machining parameter optimization and experimental investigations of nano-graphene mixed electrical discharge machining of nitinol shape memory alloy,” Journal of Materials Research and Technology, vol. 19, pp. 653–668, Jul. 2022, doi: 10.1016/j.jmrt.2022.05.076.
[34]	R. Chaudhari et al., “Experimental Investigations and Effect of Nano-Powder-Mixed EDM Variables on Performance Measures of Nitinol SMA,” Materials 2022, Vol. 15, vol. 15, no. 20, Oct. 2022, doi: 10.3390/ma15207392.
[35]	R. R. Kolhapure and D. S. Badkar, “Study of Machining Characteristics of Non-conventional Methods,” pp. 288–305, May 2023, doi: 10.2991/978-94-6463-136-4_26.
[36]	R. K. Sahu, S. H. Somashekhar, and P. V. Manivannan, “Investigation on Copper Nanofluid Obtained through Micro Electrical Discharge Machining for Dispersion Stability and Thermal Conductivity,” Procedia Eng., vol. 64, no. 2, pp. 946–955, Jan. 2013, doi: 10.1016/j.proeng.2013.09.171.
[37]	K. H. Tseng, H. L. Lee, C. Y. Liao, K. C. Chen, and H. S. Lin, “Rapid and Efficient Synthesis of Silver Nanofluid Using Electrical Discharge Machining,” J. Nanomater., vol. 2013, no. 1, p. 174939, Jan. 2013, doi: 10.1155/2013/174939.
[38]	C. Y. Chang, K. H. Tseng, J. T. Chang, M. Y. Chung, and Z. Y. Lin, “A Study of Nano-Tungsten Colloid Preparing by the Electrical Spark Discharge Method,” Micromachines (Basel)., vol. 13, no. 11, p. 2009, Nov. 2022, doi: 10.3390/mi13112009.
[39]	R. K. Sahu, S. S. Hiremath, P. V. Manivannan, and M. Singaperumal, “An Innovative Approach for Generation of Aluminium Nanoparticles Using Micro Electrical Discharge Machining,” Procedia Materials Science, vol. 5, pp. 1205–1213, Jan. 2014, doi: 10.1016/j.mspro.2014.07.424.
[40]	K. H. Tseng, J. C. Huang, K. C. Chen, and C. Y. Liao, “Preparation of metal nano-fluid via electrical discharge machining,” Proceedings of the 2011 6th IEEE Conference on Industrial Electronics and Applications, ICIEA 2011, pp. 1085–1089, 2011, doi: 10.1109/ICIEA.2011.5975748.
[41]	K. H. Tseng, M. Y. Chung, and C. Y. Chang, “Parameters for Fabricating Nano-Au Colloids through the Electric Spark Discharge Method with Micro-Electrical Discharge Machining,” Nanomaterials, vol. 7, no. 6, p. 133, Jun. 2017, doi: 10.3390/nano7060133.
[42]	B. Nahak and A. Gupta, “A review on optimization of machining performances and recent developments in electro discharge machining,” Manuf. Rev. (Les Ulis)., vol. 6, p. 2, 2019, doi: 10.1051/mfreview/2018015.
[43]	R. Çakıroğlu and M. Günay, “Comprehensive analysis of material removal rate, tool wear and surface roughness in electrical discharge turning of L2 tool steel,” Journal of Materials Research and Technology, vol. 9, no. 4, pp. 7305–7317, Jul. 2020, doi: 10.1016/j.jmrt.2020.04.060.
[44]	M. S. Mahdieh, “Recast layer and heat-affected zone structure of ultra-fined grained low-carbon steel machined by electrical discharge machining,” Proc. Inst. Mech. Eng. B J. Eng. Manuf., vol. 234, no. 5, pp. 933–944, Apr. 2020, doi: 10.1177/0954405419889202.
[45]	X. M. Wang, X. Q. Zeng, G. S. Wu, S. S. Yao, and L. B. Li, “Surface oxidation behavior of MgNd alloys,” Appl. Surf. Sci., vol. 253, no. 22, pp. 9017–9023, Sep. 2007, doi: 10.1016/j.apsusc.2007.05.023.
[46]	M. Honarvar, B. Konh, T. K. Podder, A. P. Dicker, Y. Yu, and P. Hutapea, “X-ray Diffraction Investigations of Shape Memory NiTi Wire,” J. Mater. Eng. Perform., vol. 24, no. 8, pp. 3038–3048, Aug. 2015, doi: 10.1007/s11665-015-1574-2.
[47]	J. I. Goldstein, D. E. Newbury, J. R. Michael, N. W. M. Ritchie, J. H. J. Scott, and D. C. Joy, “Scanning electron microscopy and x-ray microanalysis,” Scanning Electron Microscopy and X-ray Microanalysis, pp. 1–550, Jan. 2017, doi: 10.1007/978-1-4939-6676-9.
[48]	D. B. Williams and C. B. Carter, “The Transmission Electron Microscope,” Transmission Electron Microscopy, pp. 3–22, 2009, doi: 10.1007/978-0-387-76501-3_1.
[49]	P. Gill, T. T. Moghadam, and B. Ranjbar, “Differential Scanning Calorimetry Techniques: Applications in Biology and Nanoscience,” J. Biomol. Tech., vol. 21, no. 4, p. 167, Dec. 2010, Accessed: Feb. 03, 2026. [Online]. Available: https://pmc.ncbi.nlm.nih.gov/articles/PMC2977967/
[50]	K. Otsuka and C. Wayman, Shape memory materials. 1999. Accessed: Feb. 03, 2026. [Online]. Available: https://books.google.co.in/books?hl=en&lr=&id=DvItE9XUlN8C&oi=fnd&pg=PR11&dq=K.+Otsuka,+C.+M.+Wayman,+Shape+Memory+Materials,+Cambridge+University+Press,+1998.&ots=xmgXDdgf-P&sig=fpynzyOMF60P4i1on4R8OQHbR_8
[51]	T. Duerig, A. Pelton, and D. Stöckel, “An overview of nitinol medical applications,” Materials Science and Engineering: A, vol. 273–275, no. 1, pp. 149–160, Dec. 1999, doi: 10.1016/s0921-5093(99)00294-4.
[52]	P. Kumar, P. K. Singh, D. Kumar, V. Prakash, M. Hussain, and A. K. Das, “A novel application of micro-EDM process for the generation of nickel nanoparticles with different shapes,” Materials and Manufacturing Processes, vol. 32, no. 5, pp. 564–572, Apr. 2017, doi: 10.1080/10426914.2016.1244832.
[53]	B. Nahak and A. Gupta, “A review on optimization of machining performances and recent developments in electro discharge machining,” Manuf. Rev. (Les Ulis)., vol. 6, p. 2, 2019, doi: 10.1051/mfreview/2018015.
[54]	F. Yu et al., “Dispersion stability of thermal nanofluids,” Oct. 01, 2017, Elsevier B.V. doi: 10.1016/j.pnsc.2017.08.010.
 

