# MSS references

[Handbook lecture notes](https://wiley.fossen.biz/) list the edition for each
chapter.

The books and papers below cover the models, methods, and examples used in
the Marine Systems Simulator (MSS). See the
[complete publications list](https://www.fossen.biz/publications/) for
additional work by Thor I. Fossen.

## Foundations of the marine craft model

- T. I. Fossen (1991). *Nonlinear Modelling and Control of Underwater Vehicles*.
  Dr. Ing. thesis, Department of Engineering Cybernetics, Norwegian Institute
  of Technology (NTH), Trondheim. See the
  [historical and mathematical overview of Fossen’s marine craft model](https://fossen.biz/html/marineCraftModel.html).
- T. I. Fossen (1994). [*Guidance and Control of Ocean Vehicles*](https://www.wiley.com/en-us/guidance-and-control-of-ocean-vehicles-p-9780471941132).
  Wiley.
- T. I. Fossen (2011). [*Handbook of Marine Craft Hydrodynamics and Motion Control*](https://doi.org/10.1002/9781119994138).
  1st ed., Wiley.
- T. I. Fossen (2021). [*Handbook of Marine Craft Hydrodynamics and Motion Control*](https://www.wiley.com/en-us/handbook-of-marine-craft-hydrodynamics-and-motion-control-2nd-edition-p-9781119575030).
  2nd ed., Wiley.
- T. I. Fossen (2027). [*Handbook of Marine Craft Hydrodynamics and Motion Control*](https://wiley.fossen.biz/).
  3rd ed., Wiley, in preparation.
- S. I. Sagatun and T. I. Fossen (1991). [Lagrangian formulation of underwater vehicles’ dynamics](https://doi.org/10.1109/ICSMC.1991.169823).
  *IEEE International Conference on Systems, Man, and Cybernetics*,
  Charlottesville, VA, 1029–1034.
  See [`m2c`](../LIBRARY/modeling/m2c.m) for the Coriolis–centripetal matrix
  computed from the rigid-body inertia or added-mass matrix.

## Marine craft modeling and hydrodynamics

- T. I. Fossen (2025). [Maneuvering coefficient estimation from frequency-dependent added mass and damping: A power-based approach](https://doi.org/10.1016/j.oceaneng.2025.122494).
  *Ocean Engineering*, 341, 122494.
  See [`computeManeuveringModel`](../HYDRO/computeManeuveringModel.m) and
  [`exManeuveringModel`](../mssExamples/exManeuveringModel.m).
- T. I. Fossen (2005). [A nonlinear unified state-space model for ship maneuvering and control in a seaway](https://doi.org/10.1142/S0218127405013691).
  *International Journal of Bifurcation and Chaos*, 15(9), 2717–2746.
  [Author manuscript](https://www.fossen.biz/publications/2005%20Fossen%20IJBC.pdf).
- T. I. Fossen and O.-E. Fjellstad (1995). [Nonlinear modelling of marine vehicles in 6 degrees of freedom](https://doi.org/10.1080/13873959508837004).
  *Journal of Mathematical Modelling of Systems*, 1(1), 17–28.
  [Author PDF](https://www.fossen.biz/publications/1995%20Fossen%20and%20Fjellstad%20JMMS.pdf).
- T. Perez and T. I. Fossen (2009). [A MATLAB tool for parametric identification of radiation-force models of ships and offshore structures](https://doi.org/10.4173/mic.2009.1.1).
  *Modeling, Identification and Control*, 30(1), 1–15.
  See the [FDI toolbox and its supporting papers](../FDI_TOOLBOX/FDI-Toolbox.md).
- J. J. Jensen, A. E. Mansour, and A. S. Olsen (2004). [Estimation of ship motions using closed-form expressions](https://doi.org/10.1016/S0029-8018(03)00108-2).
  *Ocean Engineering*, 31(1), 61–85.
  See [`waveresponse345`](../LIBRARY/environment/waveresponse345.m).
- T. I. Fossen, S. I. Sagatun, and A. J. Sørensen (1996). [Identification of dynamically positioned ships](https://doi.org/10.1016/0967-0661%2896%2900014-7).
  *Control Engineering Practice*, 4(3), 369–376.
  See the [supply vessel model](../CRAFT/SHIP/models/supply.m).
- K. Nomoto, T. Taguchi, K. Honda, and S. Hirano (1957). [On the steering qualities of ships](https://doi.org/10.3233/ISP-1957-43504).
  *International Shipbuilding Progress*, 4(35), 354–370.
  See [`nomoto.m`](../GNC/nomoto.m) for the first- and second-order steering models.
- D. Clarke, P. Gedling, and G. Hine (1983). [The application of manoeuvring criteria in hull design using linear theory](https://cir.nii.ac.jp/crid/1571698600852988416).
  *Transactions of the Royal Institution of Naval Architects*, 125, 45–68.
  See [`clarke83`](../CRAFT/SHIP/models/clarke83.m).
- T. Prestero (2001). [Verification of a six-degree of freedom simulation model for the REMUS autonomous underwater vehicle](https://dspace.mit.edu/entities/publication/727deeca-8cf5-441a-9287-707a0d957fa6).
  Master’s thesis, MIT–WHOI Joint Program.
  See the [REMUS 100 model](../CRAFT/AUV/models/remus100.m).

## Guidance

- T. I. Fossen, M. Breivik, and R. Skjetne (2003). [Line-of-sight path following of underactuated marine craft](https://doi.org/10.1016/S1474-6670%2817%2937809-6).
  *IFAC Proceedings Volumes*, 36(21), 211–216.
- M. Breivik and T. I. Fossen (2005). [Principles of guidance-based path following in 2D and 3D](https://doi.org/10.1109/CDC.2005.1582226).
  *IEEE Conference on Decision and Control and European Control Conference*,
  627–634.
  [Author manuscript](https://fossen.biz/publications/2005%20Breivik%20and%20Fossen%20CDC.pdf).
- T. I. Fossen and K. Y. Pettersen (2014). [On uniform semiglobal exponential stability (USGES) of proportional line-of-sight guidance laws](https://doi.org/10.1016/j.automatica.2014.10.018).
  *Automatica*, 50(11), 2912–2917.
  [Author manuscript](https://fossen.biz/publications/2014%20Fossen%20and%20Pettersen%20Automatica.pdf).
- A. M. Lekkas and T. I. Fossen (2014). [Integral LOS path following for curved paths based on a monotone cubic Hermite spline parametrization](https://doi.org/10.1109/TCST.2014.2306774).
  *IEEE Transactions on Control Systems Technology*, 22(6), 2287–2301.
  [Author manuscript](https://fossen.biz/publications/2014%20Lekkas%20and%20Fossen%20IEEE%20TCST.pdf).
- T. I. Fossen, K. Y. Pettersen, and R. Galeazzi (2015). [Line-of-sight path following for Dubins paths with adaptive sideslip compensation of drift forces](https://doi.org/10.1109/TCST.2014.2338354).
  *IEEE Transactions on Control Systems Technology*, 23(2), 820–827.
  [Author manuscript](https://backend.orbit.dtu.dk/ws/files/129469382/FossenPettersenGaleazzi2014.pdf).
- T. I. Fossen (2023). [An adaptive line-of-sight (ALOS) guidance law for path following of aircraft and marine craft](https://doi.org/10.1109/TCST.2023.3259819).
  *IEEE Transactions on Control Systems Technology*, 31(6), 2887–2894.
  See [`ALOSpsi`](../GNC/ALOSpsi.m) and [`ALOSpsiHermite`](../GNC/ALOSpsiHermite.m).
- T. I. Fossen and A. P. Aguiar (2024). [A uniform semiglobal exponential stable adaptive line-of-sight (ALOS) guidance law for 3-D path following](https://doi.org/10.1016/j.automatica.2024.111556).
  *Automatica*, 163, 111556.
  See [`ALOS3D`](../GNC/ALOS3D.m).
- E. M. Coates and T. I. Fossen (2026). [A spherical amplitude–phase formulation for 3-D adaptive line-of-sight (ALOS) guidance with USGES stability guarantees](https://doi.org/10.1016/j.automatica.2026.113149).
  *Automatica*, 192, 113149.
  See [`SIMremus100`](../CRAFT/AUV/SIMremus100.m).
- E. Børhaug, A. Pavlov, and K. Y. Pettersen (2008). [Integral LOS control for path following of underactuated marine surface vessels in the presence of constant ocean currents](https://doi.org/10.1109/CDC.2008.4739352).
  *IEEE Conference on Decision and Control*, 4984–4991.
  See [`ILOSpsi`](../GNC/ILOSpsi.m).
- H. Xu, T. I. Fossen, and C. Guedes Soares (2020). [Uniformly semiglobally exponential stability of vector field guidance law and autopilot for path-following](https://doi.org/10.1016/j.ejcon.2019.09.007).
  *European Journal of Control*, 53, 88–97.
  [Author manuscript](https://www.fossen.biz/publications/2021%20Xu%20et%20al%20EJC.pdf).
- T. I. Fossen (2022). [Line-of-sight path-following control utilizing an extended Kalman filter for estimation of speed and course over ground from GNSS positions](https://doi.org/10.1007/s00773-022-00872-y).
  *Journal of Marine Science and Technology*, 27, 806–813.
  See [`SIMmariner`](../CRAFT/SHIP/SIMmariner.m).

## Navigation and estimation

- T. I. Fossen (2023). [Feedback error-state Kalman filter with time-delay compensation for hydroacoustic-aided inertial navigation of underwater vehicles](https://doi.org/10.1016/j.conengprac.2023.105603).
  *Control Engineering Practice*, 138, 105603.
- S. Fossen and T. I. Fossen (2021). [Five-state extended Kalman filter for estimation of speed over ground (SOG), course over ground (COG) and course rate of unmanned surface vehicles (USVs): Experimental results](https://doi.org/10.3390/s21237910).
  *Sensors*, 21(23), 7910.
  See [`EKF_5states`](../GNC/EKF_5states.m).
- R. Mahony, T. Hamel, and J.-M. Pflimlin (2008). [Nonlinear complementary filters on the special orthogonal group](https://doi.org/10.1109/TAC.2008.923738).
  *IEEE Transactions on Automatic Control*, 53(5), 1203–1218.
  See the quaternion attitude observer [`quatObserver`](../INS/functions/quatObserver.m)
  and its simulation [`SIMquatObserver`](../INS/SIMquatObserver.m).
- H. F. Grip, T. I. Fossen, T. A. Johansen, and A. Saberi (2013). [Nonlinear observer for GNSS-aided inertial navigation with quaternion-based attitude estimation](https://doi.org/10.1109/ACC.2013.6579849).
  *American Control Conference*, 272–279.
  See [`quatObserver`](../INS/functions/quatObserver.m).
- J. L. Crassidis, F. L. Markley, and Y. Cheng (2007). [Survey of nonlinear attitude estimation methods](https://doi.org/10.2514/1.22452).
  *Journal of Guidance, Control, and Dynamics*, 30(1), 12–28.
  See [`quatMEKF`](../INS/functions/quatMEKF.m).
- M. D. Shuster and S. D. Oh (1981). [Three-axis attitude determination from vector observations](https://doi.org/10.2514/3.19717).
  *Journal of Guidance and Control*, 4(1), 70–77.
  See [`quest`](../LIBRARY/kinematics/quest.m).

## Control

### Passivity, wave filtering, and dynamic positioning

- T. I. Fossen and J. P. Strand (1999). [Passive nonlinear observer design for ships using Lyapunov methods: Full-scale experiments with a supply vessel](https://doi.org/10.1016/S0005-1098%2898%2900121-6).
  *Automatica*, 35(1), 3–16.
  [Author manuscript](https://fossen.biz/publications/1999%20Fossen%20and%20Strand%20Automatica.pdf).
  See [`exPassiveWavefilterDP`](../mssExamples/exPassiveWavefilterDP.m).
- T. I. Fossen (2000). [Nonlinear passive control and observer design for ships](https://doi.org/10.4173/mic.2000.3.1).
  *Modeling, Identification and Control*, 21(3), 129–184.
  Covers the passivity framework for ship control, observers, and dynamic positioning.
- T. I. Fossen and J. P. Strand (2001). [Nonlinear passive weather optimal positioning control (WOPC) system for ships and rigs: Experimental results](https://doi.org/10.1016/S0005-1098%2801%2900006-1).
  *Automatica*, 37(5), 701–715.
  Original passivity-based WOPC design.
- T. I. Fossen (2026). [Nonlinear PID control for weather-optimal positioning of fully actuated and underactuated marine craft](https://doi.org/10.2139/ssrn.7350329).
  SSRN preprint. See the nonlinear PID example [`exWOPC`](../mssExamples/exWOPC.m).

### Control allocation

- T. A. Johansen and T. I. Fossen (2013). [Control allocation — A survey](https://doi.org/10.1016/j.automatica.2013.01.035).
  *Automatica*, 49(5), 1087–1103.
  See [`allocPseudoinverse`](../GNC/allocPseudoinverse.m) and [`SIMosv`](../CRAFT/SHIP/SIMosv.m).
