# College_Research


## I. Gravity Spectral Decomposition


### Intro

Gravity modelling has high ambiguity, even though it has been controlled with complementary data. Gravity anomaly measured at the surface is a superposition of many element from shallow to deep anomalies with different density and geometry. Suppose that we can separate gravity data from each depth of every element causing the superposition on surface anomaly, objective data from every depth can be obtained. And therefore, modelling ambiguity can also be mitigated. Decomposition of this anomaly signal based on its frequency component can be applied to separate deep anomaly with lower wavenumber, from the shallow anomaly with higher wavenumber. Hence, profiles from different depth can modeled for more objective interpretation.

Previous research also applied Spectral Decomposition to model structure of subduction in Eastern area of Java Island. Using 1D operator filter, satellite gravity data are decomposed into anomalies at different depths. The result is 2D profile of subduction structure with near surface depth. The model is validated with its geometry fits earthquake hypocenter distribution at crustal interface. The resolution of model can also be increased into 3D model with wider range area. Spectral analysis from resulting objective model can be done more easily in 3D. This research is done to analyze model obtained by Spectral Decomposition of gravity data. 


### Data and Method

Area of the research is in Java Island. The geological system is impacted by subduction of Indo-Australia Oceanic plate from South to Eurasia Continental plate in the North. This research used satellite gravity data from WGM2012. It is a 2° by 2° grid world Bouguer Gravity Anomaly data in geographical coordinate. Reprojection is done to UTM coordinate in 49S zone to extract distance parameter. The resolution of projected map is 3.7 km by 3.7 km. Then South-North 600 km long 1D profiles are then digitized in 2 km by 2 km resolution, hence 301 profiles are digitized.

These profiles are then transformed into frequency domain using Fast Fourier Transform algorithm. Earthquake hypocenter locations downloaded from USGS are also used as a validator for obtained model. From previous research, we also used filter parameters. The research used one-dimensional Lowpass Interpolated Finite Impulse Response (IFIR) with passband magnitude of 0 dB, stopband magnitude of -30 dB, with passband edge frequency from the i-th wavenumber on i-th iteration, and stopband edge frequency from the i+1 on i-th iteration. Cutoff frequency of filter would determine the depth of spectra. Hence, IFIR with narrowband characteristic is used as filter operator in this research. Filtering of two-dimensional surface anomaly map is done repetitively, as every iteration goes, the cutoff frequency increases.

Resulting spectra are saved into ASCII format for modelling. Modelling is done by interpolating spectra at its respective depth into two-dimensional Bouguer Anomaly maps. The maps are then analyzed with earthquake hypocenter locations as its validator. Then, spectra are gridded into three-dimensional cube of Bouguer Anomaly model. The cube is also analyzed using earthquake hypocenter locations.


### Result and Discussion

Result from filtering is Bouguer Gravity Anomaly Spectrum, which frequency components above cutoff frequency are attenuated. Shown in figure 2 and 3, spectrum with lower cutoff frequency has relatively smoother trend than that of higher cutoff frequency. This is caused by the limited amount of frequency component passed in filtering iteration in lower cutoff frequency.

From surface Bouguer Gravity Anomaly map, high anomaly indicates higher density. It is caused by relatively higher density of Oceanic Indo-Australia plate in Southern area. Lower anomaly indicates lower density, caused by relatively lower density of Continental Eurasia plate in the North. The transition area is indicated by light blue to greenish color in the middle section of map. This area’s Bouguer Gravity value is controlled by both Oceanic and Continental plate. As it is the area of subduction zone, where both plates converge, assimilate and creating area of high stress. The high stress causes the Bouguer Gravity Anomaly value at the area to be relatively higher.

Obtained three-dimensional model is a fitting model. The geometry resembles subduction zone. As the density varies by the stress at the location, the anomaly is differentiated into high, moderate, and low value. As shown by figure 4, 5 and 6, objective modelling of subsurface structure using Bouguer Gravity Anomaly data using Spectrum Decomposition is proven to be viable. Model would give structural geology insight of the area, which varies with density.


### Conclusions 

With wide ranged satellite gravity data and proper filter operator, subsurface structure can be modelled. Spectral analysis shows that Bouguer Gravity Anomaly, which varies with density, can give insight of structural geology. Higher anomaly is caused either by the density of rock subsurface material, and high stress which act on the area. Lower cutoff frequency from filter operator would have relatively poor feature than that of higher cutoff frequency. Limitation of this method is that to reach deeper anomaly, the length of digitized sample needs to be increased. By so, longer wavelength from low wavenumber frequency component can be obtained. Hence more detailed profile can be achieved by increasing the amount of total frequency component.


### References
Ashadi, A.L., Harmoko, U., Yuliyanto, G. and Kaka, S.I., 2015, Bulletin of the Seismological Society of America, Vol.105, 1711–1720.
Biswas, A. 2015, Geoscience Frontiers, Vol.6, 875–893.
Handyarso, A. and Kadir, W.G.A., 2017, Journal of Engineering and Technological Sciences, Vol.49, 423–437.
Mandal, A. and Niyogi, S., 2018, Journal of Applied Geophysics, Vol.159, 218–227.
Pradana, F.H., 2017, Thesis, Institut Teknologi Sepuluh Nopember.
Sandwell, D.T., Müller, R.D., Smith, W.H., Garcia, E. and Francis, R., 2014, Science, Vol.346, 65–67.
Sandwell, D.T. dan Smith, W.H., 2014, Journal of Geodesy, Vol.88, 765–771.
Setiawan, M.R. and Setiawan, A., 2017, Jurnal Fisika Indonesia, Vol.19.


